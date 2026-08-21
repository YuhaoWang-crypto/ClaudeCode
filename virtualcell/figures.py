"""Figures for the zero-shot cross-cell-line benchmark.

Two charts, each doing one job:

``zero_shot_metrics``   magnitude comparison across models, one panel per metric
``context_ablation``    change in score as source contexts are added

Colours are the validated categorical slots 1-4 (blue, orange, aqua, yellow)
assigned in fixed order to models, never cycled.  Slots 3 and 4 fall below 3:1
against the light surface, so every bar carries a direct value label.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from . import metrics as M

# validated categorical palette, light mode, fixed order
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
INK = "#0b0b0b"
INK_2 = "#52514e"
SURFACE = "#fcfcfb"
GRID = "#dcdbd6"

SHORT = {
    "control (delta=0)": "control (Δ=0)",
    "global mean [challenge baseline]": "global mean (baseline)",
    "naive transfer": "naive transfer",
    "ContextTransfer (ours)": "ContextTransfer",
}
PANELS = [
    ("discrimination_score_l1", "Perturbation discrimination", 0.5),
    ("overlap_at_100", "DE overlap @100", 0.0),
    ("de_direction_match", "DE direction match", 0.0),
    ("pearson_delta", "Pearson (effect)", 0.0),
]


def _style(ax) -> None:
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_2, labelsize=8, length=0)
    ax.yaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def zero_shot_metrics(payload: dict, out: Path) -> None:
    """One panel per metric; bars are cell lines, colour is the model."""
    results, lines = payload["results"], payload["lines"]
    models = list(SHORT)

    fig, axes = plt.subplots(1, len(PANELS), figsize=(14, 4.0), facecolor=SURFACE)
    width = 0.2
    x = np.arange(len(lines) + 1)          # cell lines, plus their mean

    for ax, (key, title, floor) in zip(axes, PANELS):
        _style(ax)
        for m_i, model in enumerate(models):
            vals = [results[c][model][key] for c in lines]
            vals.append(float(np.mean(vals)))
            bars = ax.bar(x + (m_i - 1.5) * width, vals, width * 0.88,
                          color=SERIES[m_i], label=SHORT[model],
                          linewidth=0, zorder=3)
            # slots 3-4 are under 3:1 on this surface, so label every bar
            for b, v in zip(bars, vals):
                ax.text(b.get_x() + b.get_width() / 2, v + 0.012, f"{v:.2f}",
                        ha="center", va="bottom", fontsize=5.6, color=INK_2,
                        rotation=90)
        if floor:
            ax.axhline(floor, color=INK_2, linewidth=1, linestyle=(0, (4, 3)),
                       zorder=2)
            ax.text(0.99, floor, "chance ", transform=ax.get_yaxis_transform(),
                    va="bottom", ha="right", fontsize=7, color=INK_2)
        ax.set_title(title, fontsize=10, color=INK, pad=8)
        ax.set_xticks(x)
        ax.set_xticklabels(lines + ["mean"], fontsize=8)
        top = max(max(results[c][m][key] for c in lines for m in models), floor)
        ax.set_ylim(0, top * 1.28)

    axes[0].set_ylabel("score (higher is better)", fontsize=8.5, color=INK_2)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False,
               fontsize=9, bbox_to_anchor=(0.5, -0.02), labelcolor=INK_2)
    fig.suptitle("Zero-shot knockdown prediction in a cell line with no perturbation "
                 "training data", fontsize=12.5, color=INK, y=1.0)
    fig.tight_layout(rect=(0, 0.06, 1, 0.94))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def context_ablation(ablation: dict, out: Path) -> None:
    """Score against the number of source cell lines, one line per held-out line."""
    fig, ax = plt.subplots(figsize=(6.6, 4.4), facecolor=SURFACE)
    _style(ax)

    names = list(ablation)
    sizes = sorted(int(k) for k in ablation[names[0]])
    per_line = {}
    for i, name in enumerate(names):
        ys = [float(np.mean([r["vcc_score"] for r in ablation[name][str(k)]]))
              for k in sizes]
        per_line[name] = ys
        ax.plot(sizes, ys, marker="o", markersize=7, linewidth=2,
                color=SERIES[i], label=name, zorder=3,
                markeredgecolor=SURFACE, markeredgewidth=2)

    mean = [float(np.mean([per_line[n][j] for n in names]))
            for j in range(len(sizes))]
    ax.plot(sizes, mean, linewidth=3, color=INK, zorder=4, linestyle=(0, (5, 2)),
            label="mean")

    # No direct labels here: the four lines converge to within a couple of
    # percent at the right edge, so end-of-line labels overprint each other.
    # The legend carries identity instead, which is what it is for.
    ax.set_xticks(sizes)
    ax.set_xlabel("source cell lines used for training", fontsize=9, color=INK_2)
    ax.set_ylabel("VCC score vs challenge baseline", fontsize=9, color=INK_2)
    ax.set_title("More contexts, not more cells:\nthe perturbation panel is the "
                 "same in every source line", fontsize=11.5, color=INK, pad=10)
    ax.set_xlim(sizes[0] - 0.15, sizes[-1] + 0.15)
    ax.legend(frameon=False, fontsize=8.5, loc="upper left",
              labelcolor=INK_2, handlelength=1.6)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def magnitude_frontier(payload: dict, out: Path) -> None:
    """Error against discrimination as predicted effects are scaled.

    The single-source configuration cannot clear the challenge baseline on
    discrimination and on mean absolute error at the same time, and the shape of
    that constraint is more useful than any one tuned point.  Left panel: both
    metrics against effect scale, with the baseline crossing marked.  Right
    panel: the trade-off itself, one curve per held-out line.
    """
    betas = payload["betas"]
    # plot the no-prior sweep; the transferability-prior arms are a table
    folds = {k: v["arm"]["none:0"] for k, v in payload["folds"].items()}
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.3), facecolor=SURFACE)
    for ax in axes:
        _style(ax)

    ax = axes[0]
    for i, name in enumerate(folds):
        ax.plot(betas, folds[name]["pds"], marker="o", markersize=5,
                linewidth=2, color=SERIES[i], label=name, zorder=3,
                markeredgecolor=SURFACE, markeredgewidth=1.5)
    ax.axhline(0.5, color=INK_2, linewidth=1, linestyle=(0, (4, 3)), zorder=2)
    ax.text(betas[-1], 0.503, "chance", fontsize=8, color=INK_2, ha="right")
    ax.set_xlabel("effect scale β", fontsize=9, color=INK_2)
    ax.set_ylabel("perturbation discrimination", fontsize=9, color=INK_2)
    ax.set_title("Discrimination is bought with magnitude", fontsize=11.5,
                 color=INK, pad=10)
    ax.legend(frameon=False, fontsize=8.5, loc="lower right", labelcolor=INK_2,
              handlelength=1.6)

    ax = axes[1]
    for i, name in enumerate(folds):
        d = folds[name]
        ax.plot(d["mae_vs_base"], d["pds"], marker="o", markersize=5,
                linewidth=2, color=SERIES[i], label=name, zorder=3,
                markeredgecolor=SURFACE, markeredgewidth=1.5)
    ax.axvline(0.0, color=INK, linewidth=1.2, zorder=2)
    ax.set_xlabel("MAE minus the challenge baseline  (left of the line is better)",
                  fontsize=9, color=INK_2)
    ax.set_ylabel("perturbation discrimination", fontsize=9, color=INK_2)
    ax.set_title("…and paid for in error", fontsize=11.5, color=INK, pad=10)

    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def main() -> None:
    root = Path("results")
    figs = Path("figures/virtualcell")
    if (root / "gwps_frontier.json").exists():
        magnitude_frontier(json.loads((root / "gwps_frontier.json").read_text()),
                           figs / "magnitude_frontier.png")
    if (root / "context.json").exists():
        zero_shot_metrics(json.loads((root / "context.json").read_text()),
                          figs / "zero_shot_context.png")
    if (root / "double.json").exists():
        zero_shot_metrics(json.loads((root / "double.json").read_text()),
                          figs / "zero_shot_double.png")
    if (root / "ablation.json").exists():
        context_ablation(json.loads((root / "ablation.json").read_text()),
                         figs / "context_ablation.png")


if __name__ == "__main__":
    main()
