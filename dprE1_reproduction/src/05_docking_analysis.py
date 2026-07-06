#!/usr/bin/env python3
"""
Part 2 analysis: compare reproduced Vina affinities to the paper's GOLD ChemPLP
scores (Table 3), by rank correlation, and render a comparison figure.

GOLD ChemPLP: higher = better.  Vina affinity: lower (more negative) = better.
So a faithful reproduction shows a NEGATIVE Pearson / Spearman between the two.
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import spearmanr, pearsonr

OUT = Path(__file__).resolve().parent.parent / "results"
vina = json.load(open(OUT / "docking_results.json"))

# paper Table 3 GOLD ChemPLP scores
GOLD = {
    "DprE1_WT": {"TCA1": 68.743, "GTD_9.1": 90.471, "GTD_9.2": 90.365,
                 "GTD_9.3": 85.549, "GTD_9.4": 76.870, "GTD_9.5": 83.990,
                 "GTD_9.6": 93.041, "GTD_9.7": 90.455, "GTD_9.8": 84.702,
                 "GTD_9.9": 90.842, "GTD_9.10": 92.770},
    "DprE1_Y314C": {"TCA1": 68.959, "GTD_9.1": 97.765, "GTD_9.2": 87.700,
                    "GTD_9.3": 82.559, "GTD_9.4": 86.001, "GTD_9.5": 75.381,
                    "GTD_9.6": 79.481, "GTD_9.7": 94.694, "GTD_9.8": 77.748,
                    "GTD_9.9": 83.118, "GTD_9.10": 78.145},
    "CYP2C9": {"TCA1": 81.405, "GTD_9.1": 74.172, "GTD_9.2": 68.795,
               "GTD_9.3": 88.096, "GTD_9.4": 66.266, "GTD_9.5": 91.293,
               "GTD_9.6": 92.333, "GTD_9.7": 87.593, "GTD_9.8": 96.425,
               "GTD_9.9": 84.957, "GTD_9.10": 98.204},
}

fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
md_table = ["| Compound | " + " | ".join(
    f"{t} Vina / GOLD" for t in vina) + " |",
    "|" + "---|" * (len(vina) + 1)]

corr_lines = []
compounds = ["TCA1"] + [f"GTD_9.{i}" for i in range(1, 11)]
for ax, target in zip(axes, vina):
    names = [c for c in compounds if c in vina[target] and c in GOLD[target]]
    v = np.array([vina[target][c] for c in names])
    g = np.array([GOLD[target][c] for c in names])
    rho, p = spearmanr(v, g)
    r, _ = pearsonr(v, g)
    corr_lines.append(f"{target}: Spearman rho = {rho:+.2f} (p={p:.3f}), "
                      f"Pearson r = {r:+.2f}  [negative expected]")
    ax.scatter(g, v, c=["#d62728" if n == "TCA1" else "#1f77b4" for n in names],
               s=60, zorder=3)
    for n, gg, vv in zip(names, g, v):
        ax.annotate(n.replace("GTD_", ""), (gg, vv), fontsize=7,
                    xytext=(3, 3), textcoords="offset points")
    ax.set_title(f"{target}\nSpearman rho = {rho:+.2f}")
    ax.set_xlabel("Paper GOLD ChemPLP (higher = better)")
    ax.set_ylabel("Reproduced Vina affinity\n(kcal/mol, lower = better)")
    ax.grid(alpha=0.3)

plt.suptitle("Reproduced Vina docking vs paper GOLD ChemPLP (Table 3)  —  "
             "TCA1 in red", y=1.02, fontsize=11)
plt.tight_layout()
plt.savefig(OUT / "docking_comparison.png", dpi=130, bbox_inches="tight")

# build a markdown table for RESULTS.md
lines = ["| Compound | WT Vina | WT GOLD | Y314C Vina | Y314C GOLD | "
         "CYP2C9 Vina | CYP2C9 GOLD |",
         "|---|---|---|---|---|---|---|"]
for c in compounds:
    row = [c]
    for t in ("DprE1_WT", "DprE1_Y314C", "CYP2C9"):
        row.append(f"{vina[t].get(c, '—')}")
        row.append(f"{GOLD[t].get(c, '—')}")
    lines.append("| " + " | ".join(str(x) for x in row) + " |")

print("\n".join(corr_lines))
print()
wt = vina["DprE1_WT"]
print(f"DprE1_WT: {sum(1 for k,v_ in wt.items() if k!='TCA1' and v_ < wt['TCA1'])}"
      f"/10 GTD candidates dock better than TCA1 (paper: 10/10)")

(OUT / "docking_table.md").write_text(
    "\n".join(lines) + "\n\n" + "\n".join(corr_lines) + "\n")
print(f"\nWrote {OUT/'docking_comparison.png'} and {OUT/'docking_table.md'}")


if __name__ == "__main__":
    pass
