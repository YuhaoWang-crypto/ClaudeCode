"""Step 9a -- figures.

Four figures, each written as PNG and SVG.  Colours follow one rule set: an
ordinal blue ramp for tiers, a single-hue sequential blue ramp for magnitude,
and the fixed status palette for pass/fail marks.  Every highlighted group
carries a marker shape as well as a colour, so nothing is identified by colour
alone.
"""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

from . import config as C

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
INK_MUTED = "#8a8985"
GRID = "#e4e3df"

TIER_COLOUR = {"Tier 1": "#104281", "Tier 2": "#2a78d6", "Tier 3": "#86b6ef"}
SERIES_BLUE = "#2a78d6"
SERIES_ORANGE = "#eb6834"
STATUS_CRITICAL = "#d03b3b"
STATUS_GOOD = "#0ca30c"
SEQ_BLUE = LinearSegmentedColormap.from_list(
    "seq_blue", ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])


def _style(ax):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_2, labelsize=9, length=3, color=GRID)
    ax.title.set_color(INK)


def _save(fig, name: str):
    for ext in ("png", "svg"):
        path = C.FIGURES / f"{name}.{ext}"
        fig.savefig(path, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {name}.png / {name}.svg")


def fig_ranking(ranked: pd.DataFrame, n: int = 25):
    top = ranked.head(n).iloc[::-1]
    fig, ax = plt.subplots(figsize=(9, 8.5))
    colours = [TIER_COLOUR.get(t, "#86b6ef") for t in top.tier]
    ax.barh(top.gene, top.final_score, color=colours, height=0.68)
    for y, (score, safety) in enumerate(zip(top.final_score, top.safety_score)):
        ax.text(score + 0.004, y, f"{score:.3f}", va="center", fontsize=8.5, color=INK_2)
    ax.axvline(C.TIER1_THRESHOLD, color=INK_MUTED, lw=1, ls=(0, (4, 3)))
    ax.axvline(C.TIER2_THRESHOLD, color=INK_MUTED, lw=1, ls=(0, (4, 3)))
    ax.text(C.TIER1_THRESHOLD, n + 0.2, f" Tier 1 >= {C.TIER1_THRESHOLD}",
            fontsize=8, color=INK_MUTED)
    ax.text(C.TIER2_THRESHOLD, n + 0.2, f" Tier 2 >= {C.TIER2_THRESHOLD}",
            fontsize=8, color=INK_MUTED)
    ax.set_xlabel("composite score  (geometric mean of tumour quality, safety, consensus)",
                  color=INK_2, fontsize=10)
    ax.set_title(f"Top {n} antibody-accessible surface antigens in lung adenocarcinoma",
                 fontsize=13, pad=26, loc="left")
    ax.xaxis.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    handles = [Line2D([], [], marker="s", ls="", color=c, label=t, markersize=9)
               for t, c in TIER_COLOUR.items() if (ranked.tier == t).any()]
    ax.legend(handles=handles, frameon=False, fontsize=9, loc="lower right", labelcolor=INK_2)
    _style(ax)
    _save(fig, "fig1_candidate_ranking")


def fig_compartments(ranked: pd.DataFrame, n: int = 25):
    top = ranked.head(n)
    cols = ["epi_mean", "caf_mean", "immune_mean", "endothelial_mean"]
    labels = ["Epithelial\n(incl. malignant)", "CAF", "Immune", "Endothelial"]
    mat = np.log1p(top[cols].to_numpy(dtype=float))
    fig, ax = plt.subplots(figsize=(6.6, 9))
    im = ax.imshow(mat, cmap=SEQ_BLUE, aspect="auto")
    ax.set_xticks(range(len(labels)), labels, fontsize=9, color=INK_2)
    ax.set_yticks(range(len(top)), top.gene, fontsize=9, color=INK_2)
    vmax = mat.max() if mat.size else 1
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax.text(j, i, f"{mat[i, j]:.1f}", ha="center", va="center", fontsize=7.5,
                    color="#ffffff" if mat[i, j] > 0.55 * vmax else INK_2)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label("log1p(mean expression, counts per 10,000)", color=INK_2, fontsize=9)
    cbar.ax.tick_params(colors=INK_2, labelsize=8)
    cbar.outline.set_visible(False)
    ax.set_title("Compartment expression of the top candidates\n"
                 "cross-atlas consensus, counts per 10,000 per cell",
                 fontsize=12, pad=14, loc="left", color=INK)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    _save(fig, "fig2_compartment_heatmap")


def fig_therapeutic_index(scored: pd.DataFrame, ranked: pd.DataFrame, n_top: int = 12):
    top = ranked.head(n_top)
    pos = ranked[ranked.gene.isin(C.VALIDATION_POSITIVES)]
    neg = ranked[ranked.gene.isin(C.NEGATIVE_CONTROLS)]
    fig, ax = plt.subplots(figsize=(9.2, 7.2))
    ax.scatter(ranked.tumour_quality, ranked.safety_score, s=12, c="#d8d7d2",
               edgecolors="none", label=f"all ranked candidates (n={len(ranked)})")
    ax.scatter(top.tumour_quality, top.safety_score, s=78, c=SERIES_BLUE, marker="o",
               edgecolors=SURFACE, linewidths=1.4, label=f"top {n_top} nominations", zorder=3)
    ax.scatter(pos.tumour_quality, pos.safety_score, s=86, c=SERIES_ORANGE, marker="s",
               edgecolors=SURFACE, linewidths=1.4, label="clinically validated antigens", zorder=4)
    ax.scatter(neg.tumour_quality, neg.safety_score, s=110, c=STATUS_CRITICAL, marker="^",
               edgecolors=SURFACE, linewidths=1.4, label="negative controls", zorder=5)
    # Greedy label placement: try a few offsets and keep the first that does not
    # collide with a label already placed.
    placed: list[tuple[float, float, float, float]] = []
    fig.canvas.draw()
    x_span = ranked.tumour_quality.max() - ranked.tumour_quality.min()
    y_span = ranked.safety_score.max() - ranked.safety_score.min()
    offsets = [(7, 4), (7, -11), (-9, 6), (-9, -13), (7, 15), (7, -24),
               (-9, 17), (-9, -26), (7, 27), (7, -37), (-9, 30), (-9, -40)]
    labelled: set[str] = set()
    # A gene can belong to two groups (MSLN is both a nomination and a validated
    # antigen); label it once, under the group that says more about it.
    for df, colour in ((pos, SERIES_ORANGE), (neg, STATUS_CRITICAL), (top, INK_2)):
        for r in df.sort_values("tumour_quality").itertuples():
            if r.gene in labelled:
                continue
            labelled.add(r.gene)
            w = 0.013 * x_span * max(4, len(r.gene))
            h = 0.045 * y_span
            chosen = None
            for dx, dy in offsets:
                bx = r.tumour_quality + dx * x_span / 620
                by = r.safety_score + dy * y_span / 480
                if not any(abs(bx - px) < (w + pw) / 2 and abs(by - py) < (h + ph) / 2
                           for px, py, pw, ph in placed):
                    chosen = (dx, dy, bx, by)
                    break
            # never drop a label: fall back to the farthest offset
            dx, dy, bx, by = chosen or (7, -48, r.tumour_quality + 7 * x_span / 620,
                                        r.safety_score - 48 * y_span / 480)
            placed.append((bx, by, w, h))
            ax.annotate(r.gene, (r.tumour_quality, r.safety_score),
                        textcoords="offset points", xytext=(dx, dy),
                        fontsize=8, color=colour,
                        arrowprops=dict(arrowstyle="-", color=GRID, lw=0.7,
                                        shrinkA=0, shrinkB=3) if abs(dy) > 14 else None)
    ax.set_xlabel("tumour quality  (specificity, intensity, uniformity, accessibility, druggability)",
                  color=INK_2, fontsize=10)
    ax.set_ylabel("safety coefficient  (1 = silent in normal tissue)", color=INK_2, fontsize=10)
    ax.set_title("Therapeutic index: tumour profile against normal-tissue risk",
                 fontsize=13, pad=14, loc="left")
    ax.margins(x=0.06, y=0.10)   # room for labels below the lowest points
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=9, loc="lower left", labelcolor=INK_2)
    _style(ax)
    _save(fig, "fig3_therapeutic_index")


def fig_validation(ranked: pd.DataFrame, positives: pd.DataFrame, negatives: pd.DataFrame):
    n = len(ranked)
    rows = []
    for r in positives.itertuples():
        rows.append((r.gene, r.rank, "validated antigen"))
    for r in negatives.itertuples():
        rows.append((r.gene, r.rank, "negative control"))
    rows = [(g, rk, k) for g, rk, k in rows if pd.notna(rk)]
    rows.sort(key=lambda x: x[1])
    genes = [r[0] for r in rows]
    ranks = [r[1] for r in rows]
    kinds = [r[2] for r in rows]

    fig, ax = plt.subplots(figsize=(9, 6.4))
    y = np.arange(len(rows))
    for yi, (rk, kind) in enumerate(zip(ranks, kinds)):
        colour = SERIES_ORANGE if kind == "validated antigen" else STATUS_CRITICAL
        marker = "s" if kind == "validated antigen" else "^"
        ax.hlines(yi, 1, rk, color=GRID, lw=1.4, zorder=1)
        ax.scatter([rk], [yi], s=96, color=colour, marker=marker,
                   edgecolors=SURFACE, linewidths=1.3, zorder=3)
        ax.text(rk * 1.12, yi, f"{int(rk)}", va="center", fontsize=8.5, color=INK_2)
    for k, style in ((10, ":"), (20, "--"), (C.NEGATIVE_CONTROL_FAIL_RANK, "-.")):
        ax.axvline(k, color=INK_MUTED, lw=1, ls=style)
        ax.text(k, len(rows) - 1.35, f"rank {k}", fontsize=8, color=INK_MUTED,
                ha="center", va="bottom")
    ax.set_xscale("log")
    ax.set_yticks(y, genes, fontsize=9.5)
    ax.set_xlim(1, n * 1.6)

    # validated antigens that never reached the ranking are stated, not omitted
    unranked = [r.gene for r in positives.itertuples() if pd.isna(r.rank)]
    if unranked:
        ax.text(1.05, len(rows) - 0.45,
                "not scored (below the antigen display threshold): " + ", ".join(unranked),
                fontsize=8.5, color=INK_2, va="center")
    ax.set_xlabel(f"rank among {n} scored candidates (log scale)", color=INK_2, fontsize=10)
    ax.set_title("Where the pre-registered validation set lands",
                 fontsize=13, pad=14, loc="left")
    handles = [
        Line2D([], [], marker="s", ls="", color=SERIES_ORANGE, markersize=9,
               label="clinically validated LUAD antigen"),
        Line2D([], [], marker="^", ls="", color=STATUS_CRITICAL, markersize=10,
               label="negative control (should rank low)"),
    ]
    ax.legend(handles=handles, frameon=False, fontsize=9, loc="lower right", labelcolor=INK_2)
    ax.xaxis.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    _style(ax)
    _save(fig, "fig4_validation_ranks")


def run() -> list[str]:
    ranked = pd.read_csv(C.RESULTS / "s6_ranked_candidates.csv")
    positives = pd.read_csv(C.RESULTS / "s7_validation_positives.csv")
    negatives = pd.read_csv(C.RESULTS / "s7_negative_controls.csv")
    print("[S9a] figures")
    fig_ranking(ranked)
    fig_compartments(ranked)
    fig_therapeutic_index(ranked, ranked)
    fig_validation(ranked, positives, negatives)
    names = ["fig1_candidate_ranking", "fig2_compartment_heatmap",
             "fig3_therapeutic_index", "fig4_validation_ranks"]
    C.write_provenance("s9_figures", {"figures": names,
                                      "formats": ["png", "svg"],
                                      "directory": str(C.FIGURES)})
    return names


if __name__ == "__main__":
    run()
