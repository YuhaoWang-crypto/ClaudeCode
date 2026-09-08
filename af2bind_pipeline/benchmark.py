"""Run AF2BIND across a panel of drug targets and report calibrated behaviour.

    python -m af2bind_pipeline.benchmark --out results/panel

One Modal session, one AlphaFold2 pass per target, ground truth taken from each
structure's own crystallographic ligand. The point is not a leaderboard: it is
to show how the score behaves across target classes that differ in size, in how
buried the site is, and in how much non-drug HETATM clutter the structure holds,
so that a number on a new target can be put in context.

Every entry names its ligand explicitly. Automatic ligand selection is fine on a
clean soluble structure, but on a glycoprotein or a membrane protein it will
happily count glycans, cofactors and crystallisation lipids as the drug site.

Large targets need a bigger card than the A10G default:

    AF2BIND_GPU=L40S python -m af2bind_pipeline.benchmark --out results/panel
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

#: (pdb, chain, ligand codes, target class, description)
PANEL = [
    ("6w70", "A", ["GG2"], "designed binder",
     "ABLE, a de novo apixaban binder (the paper's own example)"),
    ("1stp", "A", ["BTN"], "small soluble",
     "streptavidin + biotin"),
    ("3hs4", "A", ["AZM"], "metalloenzyme",
     "carbonic anhydrase II + acetazolamide"),
    ("1hxw", "A", ["RIT"], "interface site",
     "HIV-1 protease + ritonavir, bound at the homodimer interface"),
    ("1m17", "A", ["AQ4"], "kinase",
     "EGFR kinase domain + erlotinib"),
    ("1iep", "A", ["STI"], "kinase",
     "ABL1 kinase + imatinib, DFG-out"),
    ("3ert", "A", ["OHT"], "nuclear receptor",
     "estrogen receptor alpha LBD + 4-hydroxytamoxifen"),
    ("2rh1", "A", ["CAU"], "GPCR",
     "beta2-adrenergic receptor (T4L fusion) + carazolol"),
    ("4eiy", "A", ["ZMA"], "GPCR",
     "A2A adenosine receptor (BRIL fusion) + ZM241385"),
    ("3ln1", "A", ["CEL"], "large glycoprotein",
     "COX-2 + celecoxib"),
]


def run_panel(out_root: str | Path, seeds=(0,), targets=None, top_n: int = 15):
    from .modal_app import app, pair_features
    from .run import run_from_features
    from . import structure

    out_root = Path(out_root)
    out_root.mkdir(parents=True, exist_ok=True)
    panel = [p for p in PANEL if targets is None or p[0] in targets]

    rows = []
    with app.run():
        for pdb, chain, codes, klass, desc in panel:
            d = out_root / f"{pdb}_{chain}"
            d.mkdir(parents=True, exist_ok=True)
            row = {"pdb": pdb, "chain": chain, "ligand": codes,
                   "class": klass, "description": desc}
            try:
                path = structure.fetch_structure(pdb, d)
                result = pair_features.remote(
                    pdb_text=Path(path).read_text(),
                    chain=chain,
                    mask_sidechains=True,
                )
                report = run_from_features(
                    result, path, d, seeds=seeds, top_n=top_n, ligand_codes=codes
                )
            except Exception as exc:  # noqa: BLE001
                row["error"] = f"{type(exc).__name__}: {exc}"
                print(f"[af2bind] {pdb}/{chain} FAILED: {row['error']}", flush=True)
                rows.append(row)
                continue

            v = report.get("ligand_validation", {})
            row.update(
                n_residues=report["n_residues"],
                n_contacts=v.get("n_contact_residues"),
                positive_rate=v.get("positive_rate"),
                roc_auc=v.get("roc_auc"),
                average_precision=v.get("average_precision"),
                precision_at_10=v.get("top_n", {}).get("10", {}).get("precision_at_10"),
                precision_at_15=v.get("top_n", {}).get("15", {}).get("precision_at_15"),
                enrichment_at_10=v.get("top_n", {}).get("10", {}).get(
                    "enrichment_over_random"),
                best_rank_of_a_true_site=v.get("best_rank_of_a_true_site"),
                max_p_bind=report["p_bind"]["max"],
                mean_p_bind=report["p_bind"]["mean"],
                n_above_0_5=report["p_bind"]["n_above_0.5"],
                n_pockets=len(report["pockets"]),
            )
            rows.append(row)
            print(
                f"[af2bind] {pdb}/{chain} L={row['n_residues']} "
                f"AUC={row['roc_auc']} AP={row['average_precision']} "
                f"P@10={row['precision_at_10']} max={row['max_p_bind']}",
                flush=True,
            )

    (out_root / "panel.json").write_text(json.dumps(rows, indent=2))
    (out_root / "panel.md").write_text(as_markdown(rows))
    return rows


def as_markdown(rows: list[dict]) -> str:
    head = (
        "| target | class | L | positives | ROC-AUC | AP | P@10 | "
        "enrich@10 | max p | mean p |\n"
        "|---|---|---|---|---|---|---|---|---|---|\n"
    )
    body = []
    for r in rows:
        if "error" in r:
            body.append(
                f"| {r['pdb']} {r['chain']} | {r['class']} | — | — | "
                f"FAILED: {r['error']} | | | | | |"
            )
            continue
        if r.get("positive_rate") is None:
            # the run succeeded but no ligand matched, so there is no self-check
            body.append(
                f"| {r['pdb']} {r['chain']} ({'+'.join(r['ligand'])}) | {r['class']} | "
                f"{r['n_residues']} | no ligand matched | — | — | — | — | "
                f"{r['max_p_bind']} | {r['mean_p_bind']} |"
            )
            continue
        body.append(
            f"| {r['pdb']} {r['chain']} ({'+'.join(r['ligand'])}) | {r['class']} | "
            f"{r['n_residues']} | {r['n_contacts']} ({100 * r['positive_rate']:.1f}%) | "
            f"{r['roc_auc']} | {r['average_precision']} | {r['precision_at_10']} | "
            f"{r['enrichment_at_10']}x | {r['max_p_bind']} | {r['mean_p_bind']} |"
        )
    return head + "\n".join(body) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="af2bind-benchmark")
    ap.add_argument("--out", default="results/panel")
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--targets", default=None,
                    help="comma-separated PDB ids to restrict the panel to")
    args = ap.parse_args(argv)
    rows = run_panel(
        args.out,
        seeds=tuple(int(s) for s in args.seeds.split(",") if s.strip()),
        targets=(
            {t.strip().lower() for t in args.targets.split(",")}
            if args.targets else None
        ),
    )
    print("\n" + as_markdown(rows))
    print(f"[af2bind] wrote {args.out}/panel.json and panel.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
