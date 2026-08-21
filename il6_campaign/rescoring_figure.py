"""Figure: what the three scoring upgrades do to the IL-6 shortlist.

Left    the two predictors disagree — Boltz-2 ipSAE_min vs Chai-1 ipSAE_min
Middle  self-consistency DockQ of every shortlisted design, with the gate
Right   rank under the old ipTM score vs rank under the gated ensemble

Usage:
    python rescoring_figure.py --rescored rescored.csv --out figures/rescoring.png
"""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SURFACE, INK, INK2, INK3, GRID = "#1a1a19", "#ffffff", "#c3c2b7", "#8f8e86", "#3a3a37"
S1, S2, S4 = "#3987e5", "#d95926", "#c98500"  # validated dark-mode categorical slots
GATE = 0.23


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK3, labelsize=9)
    ax.grid(color=GRID, lw=0.6, alpha=0.55)
    ax.set_axisbelow(True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rescored", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    rows = list(csv.DictReader(args.rescored.open()))
    for r in rows:
        for k in ("ipsae_boltz", "ipsae_chai", "dockq", "ensemble_score"):
            r[k] = float(r[k])
        for k in ("rank_old_iptm", "rank_new_ensemble"):
            r[k] = int(r[k])
        r["pass"] = r["dockq"] >= GATE
    runs = {"siteI_v2": (S1, "site I run"), "siteII_v2": (S2, "site II run")}

    fig = plt.figure(figsize=(15.5, 5.4), facecolor=SURFACE)
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1], wspace=0.26,
                          left=0.055, right=0.985, top=0.76, bottom=0.15)

    # -- 1: predictor agreement
    ax = fig.add_subplot(gs[0, 0]); style(ax)
    ax.plot([0, 0.95], [0, 0.95], color=INK3, lw=1, ls="--", alpha=0.6)
    for run, (c, lab) in runs.items():
        sub = [r for r in rows if r["run"] == run]
        ax.scatter([r["ipsae_boltz"] for r in sub], [r["ipsae_chai"] for r in sub],
                   s=70, color=c, edgecolor=SURFACE, linewidth=1.5, label=lab, zorder=3)
    n_zero = sum(1 for r in rows if r["ipsae_chai"] == 0)
    ax.annotate(f"{n_zero} designs Chai-1 scores at zero\nthat Boltz-2 scored 0.46–0.83",
                xy=(0.62, 0.02), xytext=(0.30, 0.30), color=INK2, fontsize=9,
                arrowprops=dict(arrowstyle="->", color=INK3, lw=1))
    ax.set_xlabel("ipSAE$_{min}$ — Boltz-2 (designed the binders)", color=INK2, fontsize=9.5)
    ax.set_ylabel("ipSAE$_{min}$ — Chai-1 (orthogonal judge)", color=INK2, fontsize=9.5)
    ax.set_xlim(0, 0.95); ax.set_ylim(-0.03, 0.95)
    leg = ax.legend(frameon=False, fontsize=9, loc="upper left", labelcolor=INK2)
    for t in leg.get_texts():
        t.set_color(INK2)
    ax.set_title("the designer over-scores its own work", color=INK, fontsize=11, pad=9, loc="left")

    # -- 2: DockQ gate
    ax = fig.add_subplot(gs[0, 1]); style(ax)
    for i, (run, (c, lab)) in enumerate(runs.items()):
        sub = sorted([r for r in rows if r["run"] == run], key=lambda r: r["dockq"])
        x = i + np.linspace(-0.22, 0.22, len(sub))
        ax.scatter(x, [r["dockq"] for r in sub], s=70, color=c, edgecolor=SURFACE,
                   linewidth=1.5, zorder=3, alpha=[1.0 if r["pass"] else 0.35 for r in sub])
        n_ok = sum(r["pass"] for r in sub)
        ax.text(i, 0.94, f"{n_ok}/{len(sub)} pass", color=INK2, fontsize=9.5, ha="center")
    ax.axhline(GATE, color=S4, lw=1.6, ls="--")
    ax.text(1.45, 0.70, "gate: DockQ 0.23\n(acceptable pose)", color=S4, fontsize=9, ha="right")
    ax.set_xticks([0, 1]); ax.set_xticklabels([v[1] for v in runs.values()], color=INK2, fontsize=9.5)
    ax.set_ylabel("DockQ, design model vs Chai-1 prediction", color=INK2, fontsize=9.5)
    ax.set_xlim(-0.5, 1.5); ax.set_ylim(0, 1.0)
    ax.set_title("did the independent model reproduce the designed pose?",
                 color=INK, fontsize=11, pad=9, loc="left")

    # -- 3: rank churn
    ax = fig.add_subplot(gs[0, 2]); style(ax)
    ax.plot([0, 19], [0, 19], color=INK3, lw=1, ls="--", alpha=0.6)
    for run, (c, lab) in runs.items():
        sub = [r for r in rows if r["run"] == run]
        ax.scatter([r["rank_old_iptm"] for r in sub if r["pass"]],
                   [r["rank_new_ensemble"] for r in sub if r["pass"]],
                   s=70, color=c, edgecolor=SURFACE, linewidth=1.5, zorder=3)
        ax.scatter([r["rank_old_iptm"] for r in sub if not r["pass"]],
                   [r["rank_new_ensemble"] for r in sub if not r["pass"]],
                   s=70, facecolor="none", edgecolor=c, linewidth=1.6, zorder=3)
    worst = max((r for r in rows if not r["pass"]), key=lambda r: -r["rank_old_iptm"])
    ax.annotate("was rank 1 by ipTM,\nfails the pose gate",
                xy=(worst["rank_old_iptm"], worst["rank_new_ensemble"]),
                xytext=(4.5, 4.0), color=INK2, fontsize=9,
                arrowprops=dict(arrowstyle="->", color=INK3, lw=1))
    ax.scatter([], [], s=70, color=INK3, edgecolor=SURFACE, label="passes gate")
    ax.scatter([], [], s=70, facecolor="none", edgecolor=INK3, label="fails gate")
    leg = ax.legend(frameon=False, fontsize=9, loc="lower right", labelcolor=INK2)
    for t in leg.get_texts():
        t.set_color(INK2)
    ax.set_xlabel("rank by design-time ipTM (old)", color=INK2, fontsize=9.5)
    ax.set_ylabel("rank by two-predictor ipSAE ensemble (new)", color=INK2, fontsize=9.5)
    ax.set_xlim(0, 19); ax.set_ylim(0, 19)
    ax.set_title("the shortlist reshuffles", color=INK, fontsize=11, pad=9, loc="left")

    fig.suptitle("Re-scoring the IL-6 shortlist: orthogonal judge + ipSAE$_{min}$ + DockQ gate",
                 color=INK, fontsize=15.5, x=0.055, ha="left", y=0.955)
    fig.text(0.055, 0.865,
             "32 shortlisted designs · Boltz-2 designed them, Chai-1 re-folded them from sequence on a GPU · "
             "still no experimental data: every number is a prediction",
             color=INK3, fontsize=10, ha="left")
    fig.savefig(args.out, dpi=190, facecolor=SURFACE)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
