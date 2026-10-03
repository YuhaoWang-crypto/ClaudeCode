"""Figures for the GSK3-beta campaign report.

Colour follows the validated reference palette (blue/orange categorical slots on
the light chart surface; the full generated pool is drawn as neutral context
rather than as a series, because the three groups are nested rather than
independent). Chrome is recessive, every multi-series panel carries a legend,
and the series that sits below 3:1 contrast is direct-labelled.
"""

from __future__ import annotations

import json
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch
from rdkit import Chem
from rdkit.Chem import Draw

# --- palette (light mode, validated) ---------------------------------------
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"
S1 = "#2a78d6"   # blue
S2 = "#eb6834"   # orange
S3 = "#1baf7a"   # aqua (direct-labelled wherever used)
CONTEXT = "#b8b6af"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "axes.edgecolor": BASELINE,
    "axes.labelcolor": INK_2,
    "axes.titlecolor": INK,
    "text.color": INK,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "xtick.labelcolor": INK_2,
    "ytick.labelcolor": INK_2,
    "grid.color": GRID,
    "grid.linewidth": 0.8,
    "axes.grid": True,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "figure.dpi": 200,
    "savefig.dpi": 200,
    "savefig.facecolor": SURFACE,
})


def _clean(ax, title=None, xlabel=None, ylabel=None):
    if title:
        ax.set_title(title, fontsize=11, fontweight="bold", loc="left", pad=10)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=9)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=9)
    ax.set_axisbelow(True)
    ax.tick_params(length=0, labelsize=8)
    return ax


# --------------------------------------------------------------------------
def fig_convergence(history, out):
    gens = [h["generation"] for h in history]
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    series = [
        ("Best in population", [h["best"] for h in history], S1),
        ("Population mean", [h["mean"] for h in history], S2),
        ("Best novel (Tanimoto < 0.4)", [h["best_novel"] for h in history], S3),
    ]
    for label, ys, c in series:
        ax.plot(gens, ys, color=c, linewidth=2.0, solid_capstyle="round",
                label=label)
    # End labels, nudged apart so near-equal final values do not overprint.
    ends = sorted(((s[1][-1], s[2]) for s in series), key=lambda e: e[0])
    min_gap = 0.016
    placed: list[float] = []
    for value, colour in ends:
        y = value if not placed else max(value, placed[-1] + min_gap)
        placed.append(y)
        ax.annotate(f"{value:.3f}", (gens[-1], y),
                    textcoords="offset points", xytext=(7, 0),
                    fontsize=8, color=colour, va="center", fontweight="bold")
    _clean(ax, "GA convergence", "Generation", "Objective score")
    ax.set_xlim(1, max(gens) + 3.5)
    ax.set_ylim(0.55, 0.95)
    ax.legend(fontsize=8, loc="lower right", labelcolor=INK_2)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def fig_activity_qed(all_scored, survivors, top10, out):
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    a = np.array([[r["activity"], r["qed"]] for r in all_scored
                  if not r.get("gated")])
    s = np.array([[r["activity"], r["qed"]] for r in survivors])
    t = np.array([[r["activity"], r["qed"]] for r in top10])
    ax.scatter(a[:, 0], a[:, 1], s=7, c=CONTEXT, alpha=0.45, linewidths=0,
               label=f"All generated (n={len(a)})", rasterized=True)
    ax.scatter(s[:, 0], s[:, 1], s=26, c=S1, alpha=0.85, linewidths=0.6,
               edgecolors=SURFACE, label=f"Passed cascade (n={len(s)})")
    ax.scatter(t[:, 0], t[:, 1], s=95, c=S2, linewidths=1.6,
               edgecolors=SURFACE, label=f"Selected (n={len(t)})", zorder=5)
    for r in top10:
        ax.annotate(r["design_id"].replace("GSK3B-DN-", ""),
                    (r["activity"], r["qed"]), fontsize=6.5, color=SURFACE,
                    ha="center", va="center", zorder=6, fontweight="bold")
    ax.axvline(0.5, color=MUTED, linewidth=0.9, linestyle=(0, (4, 3)), zorder=1)
    ax.axhline(0.6, color=MUTED, linewidth=0.9, linestyle=(0, (4, 3)), zorder=1)
    ax.annotate("activity gate 0.50", (0.5, 0.03), fontsize=7, color=MUTED,
                rotation=90, ha="right", va="bottom")
    ax.annotate("QED gate 0.60", (0.015, 0.6), fontsize=7, color=MUTED,
                ha="left", va="bottom")
    _clean(ax, "Surrogate activity vs drug-likeness",
           "GSK3B oracle activity (surrogate)", "QED")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(0, 1.0)
    ax.legend(fontsize=8, loc="upper left", labelcolor=INK_2)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def fig_properties(gen_props, ref_props, out):
    panels = [("mw", "Molecular weight (Da)", (100, 650)),
              ("logp", "cLogP", (-2, 8)),
              ("tpsa", "TPSA (A$^2$)", (0, 200)),
              ("sa", "SA_Score (1 easy - 10 hard)", (1, 7))]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.0))
    for ax, (key, label, rng) in zip(axes.ravel(), panels):
        g = np.array([p[key] for p in gen_props if key in p and p[key] is not None])
        r = np.array([p[key] for p in ref_props if key in p and p[key] is not None])
        bins = np.linspace(rng[0], rng[1], 36)
        ax.hist(g, bins=bins, color=S1, alpha=0.75, density=True,
                label="Generated (cascade survivors)", linewidth=0)
        ax.hist(r, bins=bins, histtype="step", color=S2, linewidth=2.0,
                density=True, label="Known GSK3B actives (ChEMBL)")
        _clean(ax, None, label, "density")
        ax.set_yticks([])
    # Figure-level legend, outside the panels, so it cannot overprint the data.
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=8, loc="upper right",
               bbox_to_anchor=(0.995, 0.975), ncol=2, labelcolor=INK_2)
    fig.suptitle("Property distributions: designs vs known actives",
                 fontsize=11, fontweight="bold", x=0.02, ha="left", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def fig_novelty(all_sims, surv_sims, threshold, out):
    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    bins = np.linspace(0, 1, 51)
    ax.hist(all_sims, bins=bins, color=CONTEXT, linewidth=0,
            label=f"All generated (n={len(all_sims)})")
    ax.hist(surv_sims, bins=bins, color=S1, linewidth=0,
            label=f"Passed cascade (n={len(surv_sims)})")
    ax.axvline(threshold, color=S2, linewidth=2.0)
    ax.annotate(f"novelty cutoff {threshold:.2f}", (threshold, ax.get_ylim()[1] * 0.92),
                fontsize=8, color=S2, ha="right", va="top",
                xytext=(-6, 0), textcoords="offset points", fontweight="bold")
    _clean(ax, "Similarity of generated molecules to the nearest known active",
           "Max Tanimoto to any ChEMBL GSK3B active (ECFP4)", "count")
    ax.legend(fontsize=8, loc="upper right", labelcolor=INK_2)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def fig_cascade(log, out):
    labels = [e["stage"] for e in log]
    vals = [e["out"] for e in log]
    start = log[0]["in"]
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    ymax = start * 1.02
    for i, (lab, v) in enumerate(zip(labels, vals)):
        # rounded data-end bar, anchored at x=0, 2px visual gap between rows
        ax.add_patch(FancyBboxPatch(
            (0, i - 0.34), max(v, start * 0.004), 0.68,
            boxstyle="round,pad=0,rounding_size=%f" % (start * 0.006),
            mutation_aspect=0.0012, linewidth=0, facecolor=S1))
        ax.annotate(f"{v:,}", (v, i), xytext=(8, 0), textcoords="offset points",
                    fontsize=8.5, color=INK_2, va="center", fontweight="bold")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels([f"{l}" for l in labels], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(0, ymax * 1.12)
    ax.set_ylim(len(labels) - 0.5, -0.5)
    _clean(ax, f"Filter cascade attrition (from {start:,} unique molecules)",
           "Molecules remaining", None)
    ax.grid(axis="y", visible=False)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def fig_structures(records, out, cols=5, label_key="design_id", sub=True,
                   title=None):
    mols, legends = [], []
    for r in records:
        m = Chem.MolFromSmiles(r["smiles"] if isinstance(r, dict) else r)
        if m is None:
            continue
        mols.append(m)
        if isinstance(r, dict) and sub:
            legends.append(f"{r.get(label_key,'')}  score {r['score']:.3f}\n"
                           f"act {r['activity']:.2f} | QED {r['qed']:.2f} | "
                           f"SA {r['sa']:.2f}")
        elif isinstance(r, dict):
            legends.append(str(r.get(label_key, "")))
        else:
            legends.append("")
    img = Draw.MolsToGridImage(mols, molsPerRow=cols,
                               subImgSize=(300, 260), legends=legends,
                               returnPNG=False)
    img.save(out)
    return out


def build_all(result_dir: str):
    R = pathlib.Path(result_dir)
    figs = R / "figures"
    figs.mkdir(parents=True, exist_ok=True)
    history = json.load(open(R / "ga_history.json"))
    all_scored = json.load(open(R / "ga_all_scored.json"))
    survivors = json.load(open(R / "cascade_survivors.json"))
    top10 = json.load(open(R / "top10_designs.json"))
    log = json.load(open(R / "cascade_log.json"))

    fig_convergence(history, figs / "fig1_convergence.png")
    fig_activity_qed(all_scored, survivors, top10, figs / "fig2_activity_qed.png")
    fig_cascade(log, figs / "fig3_cascade.png")
    return figs
