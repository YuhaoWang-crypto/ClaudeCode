"""One-command evidence pipeline for ANY target (Problem B).

Mines all three sources for a UniProt id, runs the molecular + fragment +
protein-layer enrichment analyses against decoys, and writes data, figures, and a
markdown report.

Examples
--------
    python analysis/run_target.py --uniprot P43220 --name GLP1R
    python analysis/run_target.py --uniprot P07306 --name ASGR1 \
        --drugclip-decoy P00533 --ppi-decoys P04406,P00533

Defaults: DrugCLIP decoy = EGFR (P00533); PPI decoys = GAPDH (P04406) + EGFR.
Everything is free (network + CPU). Outputs under data/targets/<NAME>/ and
results/figures/<NAME>_*.png ; report at results/reports/<NAME>.md .
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from lipidlib import targetpipe as tp

BLUE, ORANGE, GREY, GREEN = "#0072B2", "#E69F00", "#8A8A8A", "#009E73"


def _save(df, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--uniprot", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--min-pchembl", type=float, default=6.0)
    ap.add_argument("--drugclip-decoy", default="P00533", help="UniProt of DrugCLIP decoy target")
    ap.add_argument("--ppi-decoys", default="P04406,P00533",
                    help="comma-separated UniProt ids for PPI decoys")
    args = ap.parse_args()

    tdir = ROOT / "data" / "targets" / args.name
    fdir = ROOT / "results" / "figures"
    rdir = ROOT / "results" / "reports"
    for d in (tdir, fdir, rdir):
        d.mkdir(parents=True, exist_ok=True)
    report: list[str] = [f"# {args.name} ({args.uniprot}) — evidence report\n"]

    # --- mine ------------------------------------------------------------
    print(f"[1/5] ChEMBL …", flush=True)
    chembl_id = tp.resolve_chembl_target(args.uniprot)
    chembl = tp.mine_chembl(chembl_id, args.min_pchembl) if chembl_id else pd.DataFrame()
    if len(chembl):
        _save(chembl, tdir / f"{args.name}_chembl_ligands.csv")
    n_sm = int((chembl.get("class") == "small_molecule").sum()) if len(chembl) else 0
    report.append(f"- **ChEMBL** ({chembl_id}): {len(chembl)} actives "
                  f"({n_sm} small-molecule).")
    print(f"      {chembl_id}: {len(chembl)} actives")

    print(f"[2/5] DrugCLIP …", flush=True)
    drug = tp.mine_drugclip(args.uniprot)
    if len(drug):
        _save(drug, tdir / f"{args.name}_drugclip_predicted.csv")
    report.append(f"- **DrugCLIP**: {len(drug)} predicted hits "
                  f"(MW median {drug['mw'].median():.0f})." if len(drug)
                  else "- **DrugCLIP**: no coverage for this UniProt.")
    print(f"      {len(drug)} predicted ligands")

    print(f"[3/5] humanPPI …", flush=True)
    ppi = tp.mine_humanppi(args.uniprot)
    if len(ppi):
        _save(ppi, tdir / f"{args.name}_humanppi_partners.csv")
    report.append(f"- **humanPPI**: {len(ppi)} predicted partners "
                  f"({int(ppi['is_cell_surface'].sum())} cell-surface)." if len(ppi)
                  else "- **humanPPI**: no partners found.")
    print(f"      {len(ppi)} partners")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "figure.dpi": 130})

    # --- molecular + fragment enrichment ---------------------------------
    print(f"[4/5] molecular/fragment enrichment …", flush=True)
    if len(drug) and len(chembl):
        actives_sm = chembl[chembl["class"] == "small_molecule"]["smiles"].tolist()
        decoy_drug = tp.mine_drugclip(args.drugclip_decoy)
        decoy_sm = decoy_drug["smiles"].tolist() if len(decoy_drug) else []
        me = tp.molecular_enrichment(drug["smiles"], actives_sm, decoy_sm)
        fe = tp.fragment_enrichment(drug["smiles"], decoy_sm, actives_sm)
        report.append(
            f"\n## Molecular cross-comparison (DrugCLIP vs measured)\n"
            f"- whole-molecule NN-Tanimoto: target median {me['fg_median']:.3f} vs "
            f"decoy {me['decoy_median']:.3f}; MWU p={me['mwu_p']:.2g}; "
            f"EF@0.35={me['enrichment_factor']:.2f}×\n"
            f"- fragment-in-active containment: {fe['fg_hit']}/{fe['fg_tot']} vs "
            f"decoy {fe['decoy_hit']}/{fe['decoy_tot']}; Fisher p={fe['fisher_p']:.2g}, "
            f"OR={fe['odds']:.2f}\n"
            f"- MW mismatch: DrugCLIP {drug['mw'].median():.0f} vs actives "
            f"{chembl[chembl['class']=='small_molecule']['mw'].median():.0f} Da")

        fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.4))
        bins = np.linspace(0, 1, 26)
        ax[0].hist(me["decoy_sim"], bins=bins, density=True, color=GREY, alpha=.75,
                   label=f"decoy (n={me['n_decoy']})")
        ax[0].hist(me["fg_sim"], bins=bins, density=True, color=BLUE, alpha=.75,
                   label=f"{args.name} (n={me['n_fg']})")
        ax[0].set_xlabel("Nearest Tanimoto to a measured active")
        ax[0].set_ylabel("Density")
        ax[0].set_title(f"A  Whole-molecule enrichment\nEF@0.35={me['enrichment_factor']:.1f}×, "
                        f"p={me['mwu_p']:.1e}")
        ax[0].legend(frameon=False, fontsize=8.5)
        mb = np.linspace(100, 950, 30)
        ax[1].hist(chembl[chembl["class"] == "small_molecule"]["mw"].dropna(), bins=mb,
                   density=True, color=ORANGE, alpha=.75, label="measured actives")
        ax[1].hist(drug["mw"].dropna(), bins=mb, density=True, color=BLUE, alpha=.75,
                   label="DrugCLIP hits")
        ax[1].set_xlabel("Molecular weight (Da)"); ax[1].set_ylabel("Density")
        ax[1].set_title("B  Chemical-space (size) check")
        ax[1].legend(frameon=False, fontsize=8.5)
        fig.tight_layout(); fig.savefig(fdir / f"{args.name}_molecular.png", bbox_inches="tight")
        plt.close(fig)
    else:
        report.append("\n## Molecular cross-comparison\n- skipped (missing ChEMBL or DrugCLIP data).")

    # --- protein-layer enrichment ----------------------------------------
    print(f"[5/5] protein-layer enrichment …", flush=True)
    if len(ppi):
        decoys = {}
        for up in [d.strip() for d in args.ppi_decoys.split(",") if d.strip()]:
            dp = tp.mine_humanppi(up)
            if len(dp):
                decoys[up] = dp
        pe = tp.ppi_surface_enrichment(ppi, decoys)
        vs = "; ".join(f"vs {k}: {v['frac']:.1%} (p={v['p']:.1g})" for k, v in pe["vs"].items())
        report.append(
            f"\n## Protein-layer enrichment (humanPPI)\n"
            f"- {args.name} cell-surface partners: {pe['target_surface']}/{pe['target_total']} "
            f"= {pe['target_frac']:.1%}\n- {vs}"
            + (f"; vs genome ~25%: p={pe['vs_genome']['p']:.1g}" if "vs_genome" in pe else ""))
        surf = ppi[ppi["is_membrane"]].sort_values("AF", ascending=False)
        _save(surf, tdir / f"{args.name}_surface_partners.csv")

        fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
        labels = [args.name] + list(pe["vs"].keys()) + ["genome"]
        fracs = [pe["target_frac"]] + [v["frac"] for v in pe["vs"].values()] + \
                [pe.get("vs_genome", {}).get("frac", tp.GENOME_MEMBRANE_FRAC)]
        ax[0].bar(labels, fracs, color=[BLUE] + [GREY] * len(pe["vs"]) + ["#CCC"], width=.66)
        ax[0].set_ylabel("Cell-surface partner fraction")
        ax[0].set_title(f"A  {args.name} surface enrichment")
        ax[0].tick_params(axis="x", labelsize=8)
        top = surf.head(10).iloc[::-1]
        yp = np.arange(len(top))
        bc = [GREEN if str(c).lower().startswith("high") else BLUE for c in top["confident"]]
        ax[1].barh(yp, top["AF"], color=bc, height=.7)
        ax[1].set_yticks(yp); ax[1].set_yticklabels(top["Gene Name"], fontsize=8)
        ax[1].set_xlabel("AF_Score"); ax[1].set_title("B  Top surface partners")
        fig.tight_layout(); fig.savefig(fdir / f"{args.name}_humanppi.png", bbox_inches="tight")
        plt.close(fig)
    else:
        report.append("\n## Protein-layer enrichment\n- skipped (no humanPPI data).")

    (rdir / f"{args.name}.md").write_text("\n".join(report) + "\n")
    print(f"\nreport -> {rdir / f'{args.name}.md'}")
    print(f"figures -> {fdir}/{args.name}_molecular.png , {args.name}_humanppi.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
