"""Figures for the OpenKnot benchmark reproduction.

  fig_success_by_method.png -- per-method success rate per round (the paper's
                               Fig. 1H/1J and Fig. 2C/2F bar charts)
  fig_ai_vs_human.png       -- targets solved by any AI method vs by Eterna
                               participants vs the starting sequences
  fig_score_validation.png  -- recomputed minus released OpenKnot score
  fig_rnet_validation.png   -- RNet reproduction, and simulated vs measured score

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




def figure_rnet_validation() -> Path | None:
    """RNet reproduction, and what the simulated score is worth on each molecule.

    The middle and right panels are the same designs scored two ways, so the
    difference between them is the choice of molecule: the isolated design, or
    the padded construct the experiment actually measured.
    """
    path = RESULTS / "rnet_predictions.csv"
    if not path.exists():
        return None
    data = pd.read_csv(path)
    padded_path = RESULTS / "rnet_padded_check.csv"
    padded = (
        pd.read_csv(padded_path).drop_duplicates(subset=["design_sequence"])
        if padded_path.exists()
        else None
    )
    n_panels = 3 if padded is not None else 2
    fig, axes = plt.subplots(1, n_panels, figsize=(5.4 * n_panels, 4.9), facecolor=SURFACE)

    ax = axes[0]
    style(ax)
    ax.yaxis.grid(True, color=GRID, linewidth=0.6)
    ax.scatter(data["released_RNet_F1"], data["rnet_F1"], s=16, color=BLUE,
               alpha=0.55, linewidths=0)
    ax.plot([0, 1], [0, 1], color=INK_SOFT, linewidth=1, linestyle=(0, (4, 3)))
    ax.axvline(0.8, color=GRID, linewidth=1)
    ax.axhline(0.8, color=GRID, linewidth=1)
    ax.set_xlabel("released RNet F1", fontsize=9, color=INK_SOFT)
    ax.set_ylabel("RNet F1 recomputed here", fontsize=9, color=INK_SOFT)
    agree = float(((data["rnet_F1"] >= 0.8) == (data["released_RNet_F1"] >= 0.8)).mean())
    ax.set_title(
        f"Secondary structure\nsame filter call on {100 * agree:.0f}% of {len(data)} designs",
        fontsize=10, color=INK, loc="left", pad=10,
    )

    def score_panel(ax, measured, simulated, title, note):
        style(ax)
        ax.yaxis.grid(True, color=GRID, linewidth=0.6)
        ax.scatter(measured, simulated, s=16, color=ORANGE, alpha=0.55, linewidths=0)
        limits = (30, 102)
        ax.plot(limits, limits, color=INK_SOFT, linewidth=1, linestyle=(0, (4, 3)))
        ax.axvline(90, color=GRID, linewidth=1)
        ax.axhline(90, color=GRID, linewidth=1)
        ax.set_xlim(*limits)
        ax.set_ylim(*limits)
        ax.set_xlabel("experimental OpenKnot score", fontsize=9, color=INK_SOFT)
        ax.set_ylabel("simulated from RNet reactivity", fontsize=9, color=INK_SOFT)
        rho = pd.DataFrame({"a": simulated, "b": measured}).corr(method="spearman").iloc[0, 1]
        ax.set_title(f"{title}\nSpearman {rho:.2f}", fontsize=10, color=INK, loc="left", pad=10)
        ax.text(0.03, 0.96, note, transform=ax.transAxes, va="top", fontsize=7.5,
                color=INK_SOFT)

    if padded is None:
        score_panel(axes[1], data["target_openknot_score"], data["simulated_openknot_score"],
                    "Scoring without an experiment",
                    "points above the diagonal are\ndesigns RNet scores too generously")
    else:
        score_panel(axes[1], padded["target_openknot_score"], padded["simulated_unpadded"],
                    "Predicting on the design alone",
                    f"{len(padded)} designs\nnot the molecule that was measured")
        score_panel(axes[2], padded["target_openknot_score"], padded["simulated_padded"],
                    "Predicting on the padded construct",
                    "the same designs, and the same\nmolecule the experiment saw")

    fig.suptitle(
        f"RibonanzaNet reproduced on released designs",
        fontsize=12.5, color=INK, x=0.01, ha="left", y=0.98,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out = FIGURES / "fig_rnet_validation.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    return out


def figure_grnade_designs() -> Path | None:
    """Per-target scores: designs made here, the paper's design, and its measurement.

    A dot plot rather than bars: every value sits between 80 and 100, and bars
    would have to start at zero and waste the range that carries the signal.
    """
    path = RESULTS / "design_comparison.csv"
    if not path.exists():
        return None
    table = pd.read_csv(path).sort_values(["round", "puzzle"]).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(12.0, 5.0), facecolor=SURFACE)
    style(ax)
    ax.yaxis.grid(True, color=GRID, linewidth=0.6)
    ax.xaxis.grid(False)
    x = np.arange(len(table))

    ours = table["our_gRNAde_best_simulated"].to_numpy()
    theirs = table["released_gRNAde_simulated"].to_numpy()
    measured = table["released_gRNAde_experimental"].to_numpy()
    ax.vlines(x, np.minimum(ours, theirs), np.maximum(ours, theirs),
              color=GRID, linewidth=1.4, zorder=1)
    ax.scatter(x, theirs, s=34, color=ORANGE, zorder=3, label="the paper's design, simulated")
    ax.scatter(x, ours, s=34, color=BLUE, zorder=3, label="best design made here, simulated")
    ax.scatter(x, measured, s=40, color=AQUA, marker="D", zorder=2,
               label="the paper's design, measured")
    ax.axhline(90, color=INK_SOFT, linewidth=1, linestyle=(0, (4, 3)))

    ax.set_xticks(x)
    ax.set_xticklabels(table["puzzle"], fontsize=7.5, color=INK, rotation=90)
    ax.set_xlim(-0.8, len(table) - 0.2)
    ax.set_ylim(78.5, 101)
    ax.set_ylabel("OpenKnot score", fontsize=9, color=INK_SOFT)
    ax.text(-0.6, 90.4, "success cutoff", fontsize=7.5, color=INK_SOFT,
            ha="left", va="bottom")
    ax.legend(frameon=False, fontsize=8.5, labelcolor=INK_SOFT, ncols=3,
              loc="lower left", bbox_to_anchor=(0.0, 1.0))
    ax.set_title(
        "gRNAde designs for the Round 3 and Round 4 targets",
        fontsize=11, color=INK, loc="left", pad=30,
    )
    fig.text(
        0.01, 0.005,
        "Blue vs orange is like-for-like, both scored in silico: blue from 8 samples per target "
        "here, orange from the paper's search of up to a million.\nOrange vs green is the price "
        "of scoring the isolated design; predicting on the padded construct closes most of it.",
        fontsize=7.5, color=INK_SOFT,
    )
    fig.tight_layout(rect=(0, 0.065, 1, 1))
    out = FIGURES / "fig_grnade_designs.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    return out


def main() -> int:
    FIGURES.mkdir(exist_ok=True)
    per_method = pd.read_csv(RESULTS / "success_rates_by_method.csv")
    per_group = pd.read_csv(RESULTS / "success_rates_by_group.csv")
    scores = pd.read_csv(RESULTS / "recomputed_scores.csv")
    for path in (
        figure_success_by_method(per_method),
        figure_ai_vs_human(per_group),
        figure_score_validation(scores),
        figure_rnet_validation(),
        figure_grnade_designs(),
    ):
        if path is not None:
            print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
