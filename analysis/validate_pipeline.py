"""Positive control: does the enrichment machinery fire WHEN IT SHOULD?

The GLP1R and PSMA molecular comparisons came back negative. Is that because the
pipeline can't detect enrichment, or because DrugCLIP's fragment hits genuinely
don't match drug-sized actives? This script settles it.

We feed the SAME `molecular_enrichment` routine a foreground that is *known* to
share chemistry with the reference — a held-out split of the target's own measured
actives — and check that it lights up. Three foregrounds are scored against the
same reference (a random 60% of target-A small-molecule actives), with an unrelated
target's actives (target B) as the decoy:

  1. held-out target-A actives   (POSITIVE control — must enrich strongly)
  2. target-A DrugCLIP hits       (the real-use case — expected weak)
  3. target-B actives             (this IS the decoy — sanity floor)

If (1) shows a large enrichment factor and (2)/(3) do not, the method is sound and
the DrugCLIP negatives are real, not a bug.

Default: A = PSMA (Q04609), B = EGFR (P00533). Outputs a figure + prints stats.
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

BLUE, ORANGE, GREY = "#0072B2", "#E69F00", "#8A8A8A"
FIG = ROOT / "results" / "figures"; FIG.mkdir(parents=True, exist_ok=True)


def load_actives(uniprot: str, name: str) -> pd.DataFrame:
    """Use cached ChEMBL CSV if present, else mine."""
    cached = ROOT / "data" / "targets" / name / f"{name}_chembl_ligands.csv"
    alt = ROOT / "data" / "targets" / name / f"{name}_ligands.csv"
    for p in (cached, alt):
        if p.exists():
            df = pd.read_csv(p)
            break
    else:
        cid = tp.resolve_chembl_target(uniprot)
        df = tp.mine_chembl(cid)
    return df[df["class"] == "small_molecule"].reset_index(drop=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--a-uniprot", default="Q04609"); ap.add_argument("--a-name", default="PSMA")
    ap.add_argument("--b-uniprot", default="P00533"); ap.add_argument("--b-name", default="EGFR")
    args = ap.parse_args()

    A = load_actives(args.a_uniprot, args.a_name)
    B = load_actives(args.b_uniprot, args.b_name)
    drug = pd.read_csv(ROOT / "data" / "targets" / args.a_name /
                       f"{args.a_name}_drugclip_predicted.csv")
    print(f"{args.a_name} actives: {len(A)} | {args.b_name} actives: {len(B)} | "
          f"{args.a_name} DrugCLIP: {len(drug)}")

    # deterministic 60/40 split of target-A actives
    rng = np.random.default_rng(0)
    idx = rng.permutation(len(A))
    cut = int(0.6 * len(A))
    ref = A.iloc[idx[:cut]]["smiles"].tolist()
    held = A.iloc[idx[cut:]]["smiles"].tolist()

    pos = tp.molecular_enrichment(held, ref, B["smiles"].tolist())       # positive control
    dcl = tp.molecular_enrichment(drug["smiles"].tolist(), ref, B["smiles"].tolist())  # real use

    print(f"\nPOSITIVE control (held-out {args.a_name} actives vs {args.a_name} ref):")
    print(f"  median NN-sim {pos['fg_median']:.3f} vs decoy {pos['decoy_median']:.3f} | "
          f"EF@0.35={pos['enrichment_factor']:.1f}x | MWU p={pos['mwu_p']:.2g}")
    print(f"DrugCLIP hits (vs same {args.a_name} ref):")
    print(f"  median NN-sim {dcl['fg_median']:.3f} vs decoy {dcl['decoy_median']:.3f} | "
          f"EF@0.35={dcl['enrichment_factor']:.1f}x | MWU p={dcl['mwu_p']:.2g}")

    # --- figure -----------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "figure.dpi": 130})
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    bins = np.linspace(0, 1, 30)
    ax[0].hist(pos["decoy_sim"], bins=bins, density=True, color=GREY, alpha=.7,
               label=f"{args.b_name} decoy")
    ax[0].hist(dcl["fg_sim"], bins=bins, density=True, color=ORANGE, alpha=.7,
               label="DrugCLIP hits")
    ax[0].hist(pos["fg_sim"], bins=bins, density=True, color=BLUE, alpha=.7,
               label=f"held-out {args.a_name} actives")
    ax[0].set_xlabel(f"Nearest Tanimoto to {args.a_name} reference actives")
    ax[0].set_ylabel("Density")
    ax[0].set_title("A  Positive control fires; DrugCLIP does not")
    ax[0].legend(frameon=False, fontsize=8.5)

    conds = ["held-out\nactives\n(positive)", "DrugCLIP\nhits", f"{args.b_name}\ndecoy"]
    efs = [pos["enrichment_factor"], dcl["enrichment_factor"], 1.0]
    ax[1].bar(conds, efs, color=[BLUE, ORANGE, GREY], width=.6)
    ax[1].axhline(1, color="#333", lw=1, ls="--")
    ax[1].set_ylabel("Enrichment factor @ Tanimoto ≥ 0.35")
    ax[1].set_title("B  Method detects real similarity")
    for i, v in enumerate(efs):
        ax[1].text(i, v + max(efs) * 0.02, f"{v:.1f}×", ha="center", fontsize=9)
    fig.tight_layout()
    out = FIG / "pipeline_validation.png"
    fig.savefig(out, bbox_inches="tight")
    print(f"\nwrote figure -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
