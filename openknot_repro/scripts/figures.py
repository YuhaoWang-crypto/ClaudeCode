"""Figures for the OpenKnot benchmark reproduction.

  fig_success_by_method.png -- per-method success rate per round (the paper's
                               Fig. 1H/1J and Fig. 2C/2F bar charts)
  fig_ai_vs_human.png       -- targets solved by any AI method vs by Eterna
                               participants vs the starting sequences
  fig_score_validation.png  -- recomputed minus released OpenKnot score

Usage:  python scripts/figures.py   (after validate_score.py and success_rates.py)
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"

# Validated categorical slots (dataviz reference palette, light surface).
BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SOFT = "#52514e"
GRID = "#dedcd5"

METHOD_ORDER = [
    "Starting sequence",
    "Rosetta",
    "Rosetta-LoRes",
    "3DRNA",
    "MPNN-fixbb",
    "MPNN-RFdiff",
    "codesign-RFdiff",
    "gRNAde",
    "gRNAde-no3d",
    "Struct2SeQ",
    "Struct2SeQ-SHAPE",
    "Eterna",
]
ROUND_TITLES = {
    1: "Round 1 - 17 targets",
    2: "Round 2 - same 17 targets, RNet available",
    3: "Round 3 - 20 new targets, up to 100 nt",
    4: "Round 4 - 20 new targets, 117-240 nt",
}


def style(ax) -> None:
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_SOFT, labelsize=8, length=0)
    ax.xaxis.grid(True, color=GRID, linewidth=0.6)
    ax.yaxis.grid(False)
    ax.set_axisbelow(True)


def figure_success_by_method(per_method: pd.DataFrame) -> Path:
    data = per_method[per_method["score_column"] == "OKS_recomputed"]
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 9.5), facecolor=SURFACE)
    for ax, rnd in zip(axes.ravel(), (1, 2, 3, 4)):
        style(ax)
        block = data[data["round"] == rnd]
        methods = [m for m in METHOD_ORDER if m in set(block["method"])]
        y = np.arange(len(methods))
        height = 0.38
        for offset, flag, color, label in (
            (height / 2, False, ORANGE, "all designs"),
            (-height / 2, True, BLUE, "RNet F1 >= 0.8 only"),
        ):
            sub = block[block["rnet_filter"] == flag].set_index("method").reindex(methods)
            ax.barh(
                y + offset,
                sub["success_rate"].to_numpy(),
                height=height,
                color=color,
                label=label,
                xerr=sub["standard_error"].to_numpy(),
                error_kw={"ecolor": INK_SOFT, "elinewidth": 0.9, "capsize": 2},
            )
            for yi, (rate, n_ok, n_sub) in enumerate(
                zip(sub["success_rate"], sub["n_success"], sub["n_targets_submitted"])
            ):
                if np.isnan(rate):
                    continue
                se = sub["standard_error"].to_numpy()[yi]
                ax.text(
                    rate + (0 if np.isnan(se) else se) + 3.0,
                    y[yi] + offset,
                    f"{int(n_ok)}/{int(n_sub)}",
                    va="center",
                    fontsize=7,
                    color=INK_SOFT,
                )
        ax.set_yticks(y)
        ax.set_yticklabels(methods, fontsize=8, color=INK)
        ax.invert_yaxis()
        ax.set_xlim(0, 126)
        ax.set_xticks([0, 25, 50, 75, 100])
        ax.set_title(ROUND_TITLES[rnd], fontsize=10, color=INK, loc="left", pad=8)
        if rnd in (3, 4):
            ax.set_xlabel("% of submitted targets with a design scoring > 90", fontsize=8.5,
                          color=INK_SOFT)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles, labels, loc="upper right", frameon=False, fontsize=8.5,
        labelcolor=INK_SOFT, ncols=2, bbox_to_anchor=(0.99, 0.995),
    )
    fig.suptitle(
        "Pseudoknot design success by method, recomputed from released SHAPE data",
        fontsize=12.5, color=INK, x=0.01, ha="left", y=0.985,
    )
    fig.text(
        0.01, 0.005,
        "Success: at least one design for the target with OpenKnot score > 90 "
        "(designs failing the release signal-to-noise gate excluded). Error bars: binomial standard error.",
        fontsize=7.5, color=INK_SOFT,
    )
    fig.tight_layout(rect=(0, 0.02, 1, 0.965))
    path = FIGURES / "fig_success_by_method.png"
    fig.savefig(path, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    return path


def figure_ai_vs_human(per_group: pd.DataFrame) -> Path:
    data = per_group[
        (per_group["score_column"] == "OKS_recomputed") & per_group["rnet_filter"]
    ]
    groups = ["AI (any method)", "Eterna (human)", "Starting sequence"]
    colors = {groups[0]: BLUE, groups[1]: ORANGE, groups[2]: AQUA}
    rounds = [1, 2, 3, 4]
    fig, ax = plt.subplots(figsize=(8.5, 4.6), facecolor=SURFACE)
    style(ax)
    ax.xaxis.grid(False)
    ax.yaxis.grid(True, color=GRID, linewidth=0.6)
    width = 0.26
    x = np.arange(len(rounds))
    for k, group in enumerate(groups):
        sub = data[data["group"] == group].set_index("round").reindex(rounds)
        positions = x + (k - 1) * width
        ax.bar(positions, sub["success_rate"].to_numpy(), width=width * 0.92,
               color=colors[group], label=group)
        for xi, (rate, solved, n) in enumerate(
            zip(sub["success_rate"], sub["n_solved"], sub["n_targets"])
        ):
            if np.isnan(rate):
                continue
            ax.text(positions[xi], rate + 2, f"{int(solved)}/{int(n)}", ha="center",
                    fontsize=8, color=INK_SOFT)
    ax.set_xticks(x)
    ax.set_xticklabels([f"Round {r}" for r in rounds], fontsize=9, color=INK)
    ax.set_ylim(0, 112)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylabel("% of targets solved", fontsize=9, color=INK_SOFT)
    ax.legend(
        frameon=False, fontsize=8.5, labelcolor=INK_SOFT, ncols=3,
        loc="lower left", bbox_to_anchor=(0.0, 1.0),
    )
    ax.set_title(
        "Targets solved by at least one design (OpenKnot score > 90, RNet F1 >= 0.8)",
        fontsize=11, color=INK, loc="left", pad=28,
    )
    fig.tight_layout()
    path = FIGURES / "fig_ai_vs_human.png"
    fig.savefig(path, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    return path


def figure_score_validation(scores: pd.DataFrame) -> Path:
    delta = scores["delta"].dropna()
    exact = float((delta.abs() < 1e-6).mean())
    fig, ax = plt.subplots(figsize=(8.0, 4.2), facecolor=SURFACE)
    style(ax)
    bins = np.linspace(-3, 3, 121)
    ax.hist(delta.clip(-3, 3), bins=bins, color=BLUE)
    ax.set_yscale("log")
    ax.set_xlabel("recomputed minus released OpenKnot score", fontsize=9, color=INK_SOFT)
    ax.set_ylabel("designs (log scale)", fontsize=9, color=INK_SOFT)
    ax.set_title(
        f"Score reproduction across {len(delta):,} released designs "
        f"- {100 * exact:.1f}% bit-exact",
        fontsize=11, color=INK, loc="left", pad=10,
    )
    ax.text(
        0.99, 0.94,
        "non-zero residuals are single residues whose released\n"
        "reactivity was rounded onto a scoring threshold,\n"
        "plus the two puzzles rescored against a second target",
        transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color=INK_SOFT,
    )
    fig.tight_layout()
    path = FIGURES / "fig_score_validation.png"
    fig.savefig(path, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    return path


def main() -> int:
    FIGURES.mkdir(exist_ok=True)
    per_method = pd.read_csv(RESULTS / "success_rates_by_method.csv")
    per_group = pd.read_csv(RESULTS / "success_rates_by_group.csv")
    scores = pd.read_csv(RESULTS / "recomputed_scores.csv")
    for path in (
        figure_success_by_method(per_method),
        figure_ai_vs_human(per_group),
        figure_score_validation(scores),
    ):
        print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
