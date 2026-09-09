#!/usr/bin/env python
"""Summary figure: leave-one-cell-line-out AUROC against a no-proteome control.

Every model bar comes from the paper's own Supplementary Table S5C, averaged over
its repeated runs and restricted to the 15 cell lines that also appear in the
label matrix. The reference line is the drug-mean control computed from
Supplementary Table S1 sheet C_Efficacy over those same 15 cell lines.
"""

from __future__ import annotations

import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Validated categorical slots 1 and 2 (blue, orange) on the light surface.
BLUE, ORANGE = "#2a78d6", "#eb6834"
SURFACE, INK, INK_2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"

LABEL = {
    "ppODE (Non-SWA)": "ProteinTalks (non-SWA)",
    "ppODE (SWA)": "ProteinTalks (SWA)",
    "Linear Regression": "Linear regression",
    "Random Forest": "Random forest",
    "Deepsynergy": "DeepSynergy",
    "genecompass": "GeneCompass",
    "geneformer": "Geneformer",
    "uce": "UCE",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default=os.path.join(HERE, "results", "paired_cellline_test.json"))
    ap.add_argument("--out", default=os.path.join(HERE, "figures", "cellline_auroc.png"))
    args = ap.parse_args()

    with open(args.results) as f:
        res = json.load(f)
    per = res["per_model"]
    control = list(per.values())[0]["control_auroc"]

    rows = sorted(
        ((LABEL.get(k, k), v["model_auroc"], v["wilcoxon_p"], k) for k, v in per.items()),
        key=lambda r: r[1],
    )
    names = [r[0] for r in rows]
    vals = [r[1] for r in rows]
    is_pt = ["ppODE" in r[3] for r in rows]

    fig, ax = plt.subplots(figsize=(9.8, 5.4), dpi=200)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    colors = [ORANGE if p else BLUE for p in is_pt]
    bars = ax.barh(names, vals, height=0.62, color=colors, zorder=3)

    ax.axvline(control, color=INK, lw=1.6, ls="--", zorder=4)
    ax.annotate(
        f"drug-mean control, {control:.3f}\n(no proteomics, no model)",
        xy=(control, len(names) - 0.42), xytext=(control - 0.014, len(names) - 0.28),
        ha="right", va="center", fontsize=8.5, color=INK,
    )

    # Value labels sit in a fixed column clear of the reference line, so no
    # label can collide with it regardless of bar length.
    label_x = 1.010
    for b, v, p, pt in zip(bars, vals, [r[2] for r in rows], is_pt):
        # Show the paired p-value against the control for every model, so the
        # reader can see which differences are real in both directions.
        note = f"  (p = {p:.3f})" if p < 0.001 else (
            f"  (p = {p:.3f})" if p < 0.05 else f"  (p = {p:.2f}, n.s.)")
        if p < 0.001:
            note = "  (p < 0.001)"
        ax.text(label_x, b.get_y() + b.get_height() / 2, f"{v:.3f}{note}",
                va="center", fontsize=9, color=INK_2)

    ax.set_xlim(0.40, 1.30)
    ax.set_xticks([0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    ax.set_xlabel("AUROC, leave-one-cell-line-out (15 cell lines)", fontsize=9.5, color=INK_2)
    ax.set_title(
        "A no-proteomics control beats every comparator the paper reports",
        fontsize=12, color=INK, pad=34, loc="left", weight="bold",
    )
    ax.xaxis.grid(True, color=GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(colors=INK_2, labelsize=9.5, length=0)

    handles = [
        plt.Rectangle((0, 0), 1, 1, color=ORANGE),
        plt.Rectangle((0, 0), 1, 1, color=BLUE),
    ]
    ax.legend(handles, ["ProteinTalks", "Comparators reported by the paper"],
              loc="lower left", bbox_to_anchor=(0.0, 1.005), ncol=2, frameon=False,
              fontsize=9, labelcolor=INK_2, handlelength=1.1, handleheight=1.1,
              columnspacing=1.6)

    fig.text(0.01, 0.015,
             "Model values: Supplementary Table S5C, averaged over repeated runs. "
             "Control: Table S1 sheet C_Efficacy. p from Wilcoxon signed-rank, paired by cell line.",
             fontsize=7.5, color=INK_2)
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out, facecolor=SURFACE)
    print(f"[saved] {args.out}")


if __name__ == "__main__":
    main()
