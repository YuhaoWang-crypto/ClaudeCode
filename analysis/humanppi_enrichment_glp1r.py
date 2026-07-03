"""Protein-layer enrichment of GLP1R's predicted interactome (humanPPI).

The molecular streams (ChEMBL, DrugCLIP) don't corroborate each other for GLP1R.
The protein layer does carry signal: is GLP1R's predicted interactome enriched for
cell-surface / membrane proteins (relevant handles/co-receptors on GLP1R+ cells)?

Tested against two decoy proteins pulled from the same humanPPI database:
  - GAPDH (P04406)  — soluble/cytosolic
  - EGFR  (P00533)  — another membrane receptor (a stringent control)
and a rough genome baseline (~25% of the human proteome is membrane-associated).

Honest caveat surfaced here too: humanPPI's binary "confident" flag is very strict
(few pass), so we rank the actionable shortlist by the continuous AF_Score.

Outputs
-------
  data/targets/GLP1R/GLP1R_surface_partners.csv   ranked membrane partners
  results/figures/glp1r_humanppi_enrichment.png
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, binomtest

ROOT = Path(__file__).resolve().parents[1]
TARGETS = ROOT / "data" / "targets"
FIG = ROOT / "results" / "figures"; FIG.mkdir(parents=True, exist_ok=True)
BLUE, ORANGE, GREY, GREEN = "#0072B2", "#E69F00", "#8A8A8A", "#009E73"
GENOME_MEMBRANE_FRAC = 0.25  # ~human proteome membrane-associated fraction


def load(name: str) -> pd.DataFrame:
    df = pd.read_csv(TARGETS / name / f"{name}_humanppi_partners.csv")
    loc = df["Subcellular Locality"].astype(str)
    df["is_membrane"] = loc.str.contains("membrane", case=False)
    df["is_cell_surface"] = loc.str.contains("Cell membrane", case=False)
    df["AF"] = pd.to_numeric(df["AF_Score"], errors="coerce")
    return df


def frac(df, col):
    return int(df[col].sum()), len(df)


def main() -> None:
    glp = load("GLP1R")
    gapdh = load("GAPDH_decoy")
    egfr = load("EGFR_ppi_decoy")

    # --- enrichment: cell-surface fraction --------------------------------
    gm, gt = frac(glp, "is_cell_surface")
    print(f"GLP1R cell-surface partners: {gm}/{gt} = {gm/gt:.1%}")
    print(f"GLP1R any-membrane partners: {glp['is_membrane'].sum()}/{gt} "
          f"= {glp['is_membrane'].mean():.1%}\n")

    for name, dec in [("GAPDH (soluble)", gapdh), ("EGFR (membrane)", egfr)]:
        dm, dt = frac(dec, "is_cell_surface")
        odds, p = fisher_exact([[gm, gt - gm], [dm, dt - dm]], alternative="greater")
        print(f"vs {name:18s}: {dm}/{dt}={dm/dt:.1%}  Fisher OR={odds:.2f} p={p:.2g}")

    bt = binomtest(gm, gt, GENOME_MEMBRANE_FRAC, alternative="greater")
    print(f"vs genome baseline {GENOME_MEMBRANE_FRAC:.0%}: binomial p={bt.pvalue:.2g}")

    # --- actionable shortlist: top cell-surface partners by AF_Score -------
    surf = glp[glp["is_membrane"]].sort_values("AF", ascending=False)
    keep = ["Gene Name", "Interacting Partner", "AF_Score", "RF_Score",
            "confident", "Subcellular Locality", "Function"]
    surf[keep].to_csv(TARGETS / "GLP1R" / "GLP1R_surface_partners.csv", index=False)
    print(f"\ntop membrane partners of GLP1R (by AF_Score):")
    print(surf.head(8)[["Gene Name", "AF_Score", "confident"]].to_string(index=False))
    print(f"wrote shortlist -> {TARGETS/'GLP1R'/'GLP1R_surface_partners.csv'}")

    # --- figure -----------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "figure.dpi": 130})
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))

    # Panel A: cell-surface fraction, GLP1R vs decoys vs baseline
    labels = ["GLP1R", "EGFR\n(membrane)", "GAPDH\n(soluble)", "genome\nbaseline"]
    fracs = [gm / gt, *[frac(d, "is_cell_surface")[0] / frac(d, "is_cell_surface")[1]
                        for d in (egfr, gapdh)], GENOME_MEMBRANE_FRAC]
    cols = [BLUE, GREY, GREY, "#CCCCCC"]
    bars = ax[0].bar(labels, fracs, color=cols, width=0.66)
    counts = [f"{gm}/{gt}", f"{frac(egfr,'is_cell_surface')[0]}/{len(egfr)}",
              f"{frac(gapdh,'is_cell_surface')[0]}/{len(gapdh)}", "~1/4"]
    for b, c in zip(bars, counts):
        ax[0].text(b.get_x() + b.get_width() / 2, b.get_height() + 0.008, c,
                   ha="center", fontsize=8.5)
    ax[0].set_ylabel("Fraction of partners at the cell surface")
    ax[0].set_title("A  GLP1R interactome is surface-enriched")

    # Panel B: top membrane partners by AF_Score (actionable)
    top = surf.head(10).iloc[::-1]
    ypos = np.arange(len(top))
    barcol = [GREEN if str(c).lower().startswith("high") else BLUE
              for c in top["confident"]]
    ax[1].barh(ypos, top["AF"], color=barcol, height=0.7)
    ax[1].set_yticks(ypos); ax[1].set_yticklabels(top["Gene Name"], fontsize=8)
    ax[1].set_xlabel("AF_Score (predicted interface confidence)")
    ax[1].set_title("B  Top surface partners (green = 'high' conf.)")
    fig.tight_layout()
    out = FIG / "glp1r_humanppi_enrichment.png"
    fig.savefig(out, bbox_inches="tight")
    print(f"wrote figure -> {out}")


if __name__ == "__main__":
    main()
