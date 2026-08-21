"""Campaign summary panel: interface-confidence distributions, epitope targeting,
and the independent-confirmation vs off-target-control comparison.

Usage:
    python campaign_summary.py --designs designs.json --iface v2_iface.json --out figures/campaign_summary.png
"""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SURFACE = "#1a1a19"
INK = "#ffffff"
INK2 = "#c3c2b7"
INK3 = "#8f8e86"
GRID = "#3a3a37"
# validated categorical slots (dark mode) - see dataviz/references/palette.md
S1, S2, S3, S4 = "#3987e5", "#d95926", "#199e70", "#c98500"

RUNS = [
    ("siteI", "site I, epitope index bug", S1, ":"),
    ("siteI_v2", "site I, corrected", S1, "-"),
    ("siteII", "site II, epitope index bug", S2, ":"),
    ("siteII_v2", "site II, corrected", S2, "-"),
]

# independent confirmation co-folds: 5 Boltz-2 samples each, no design template,
# binder run without an MSA (recorded from the Boltz API responses)
CONFIRM = [
    ("IL6-S2-01 vs IL-6", [0.9500, 0.9553, 0.9477, 0.9511, 0.9497], 0.765, S2),
    ("IL6-S1-01 vs IL-6", [0.8953, 0.8911, 0.8968, 0.8701, 0.8853], 0.211, S1),
    ("IL6-S1-02 vs IL-6", [0.8539, 0.7688, 0.8519, 0.8912, 0.8623], 0.0027, S1),
    ("IL6-S1-01 vs CLEC12A\n(off-target control)", [0.7450, 0.8131, 0.7552, 0.7646, 0.8472], 0.00067, S4),
]


def style(ax):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK3, labelsize=9)
    ax.grid(color=GRID, lw=0.6, alpha=0.55)
    ax.set_axisbelow(True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--designs", type=Path, required=True)
    ap.add_argument("--iface", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    designs = json.loads(args.designs.read_text())
    iface = json.loads(args.iface.read_text())

    fig = plt.figure(figsize=(15.5, 5.4), facecolor=SURFACE)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1, 1.05], wspace=0.28,
                          left=0.055, right=0.985, top=0.76, bottom=0.16)

    # -- panel 1: cumulative distribution of design-time interface confidence
    ax = fig.add_subplot(gs[0, 0]); style(ax)
    for key, label, color, ls in RUNS:
        v = np.sort([r["iptm"] for r in designs[key]])
        y = np.arange(1, len(v) + 1) / len(v)
        ax.plot(v, 1 - y, color=color, lw=2.0, ls=ls, label=label)
    ax.axvline(0.85, color=INK3, lw=1.0, ls="--", alpha=0.7)
    ax.text(0.843, 0.045, "ipTM 0.85", color=INK3, fontsize=8.5, ha="right")
    ax.set_xlabel("design-time interface confidence (ipTM)", color=INK2, fontsize=9.5)
    ax.set_ylabel("fraction of designs above", color=INK2, fontsize=9.5)
    ax.set_xlim(0.1, 1.0); ax.set_ylim(0, 1)
    leg = ax.legend(frameon=False, fontsize=8.5, loc="lower left", labelcolor=INK2)
    for t in leg.get_texts():
        t.set_color(INK2)
    ax.set_title("120 designs per run, 4 runs", color=INK, fontsize=11, pad=9, loc="left")

    # -- panel 2: epitope targeting of the 12 top-ranked designs of each corrected run
    ax = fig.add_subplot(gs[0, 1]); style(ax)
    s1 = [r for r in iface if r["name"].startswith("siteI_v2")]
    s2 = [r for r in iface if r["name"].startswith("siteII_v2")]
    for grp, color, x0, label in ((s1, S1, 0, "site I run"), (s2, S2, 1, "site II run")):
        key = "s1" if x0 == 0 else "s2"
        vals = sorted(r[key] for r in grp)
        jitter = np.linspace(-0.17, 0.17, len(vals))
        ax.scatter(x0 + jitter, vals, s=64, color=color, edgecolor=SURFACE, linewidth=1.6, zorder=3)
        med = float(np.median(vals))
        ax.plot([x0 - 0.28, x0 + 0.28], [med, med], color=INK2, lw=2.0, zorder=4)
        ax.text(x0 + 0.33, med, f"median {med:.2f}", color=INK2, fontsize=9, va="center")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["site I run\n(vs IL-6Rα footprint)", "site II run\n(vs gp130 footprint)"],
                       color=INK2, fontsize=9)
    ax.set_ylabel("fraction of the natural footprint contacted", color=INK2, fontsize=9.5)
    ax.set_ylim(0, 1.05); ax.set_xlim(-0.5, 1.75)
    ax.set_title("top 12 designs per run, epitope recall", color=INK, fontsize=11, pad=9, loc="left")

    # -- panel 3: independent confirmation vs off-target control
    ax = fig.add_subplot(gs[0, 2]); style(ax)
    ys = np.arange(len(CONFIRM))[::-1]
    for y, (label, vals, bc, color) in zip(ys, CONFIRM):
        ax.scatter(vals, [y] * len(vals), s=52, color=color, edgecolor=SURFACE, linewidth=1.5, zorder=3)
        ax.plot([min(vals), max(vals)], [y, y], color=color, lw=2.0, alpha=0.5, zorder=2)
        ax.text(0.995, y + 0.28, f"binding score {bc:.3g}", color=INK2, fontsize=8.5, ha="right")
    ax.set_yticks(ys)
    ax.set_yticklabels([c[0] for c in CONFIRM], color=INK2, fontsize=9)
    ax.set_xlabel("ipTM of 5 independent co-folds", color=INK2, fontsize=9.5)
    ax.set_xlim(0.70, 1.0); ax.set_ylim(-0.6, len(CONFIRM) - 0.35)
    ax.set_title("independent confirmation (no design template)", color=INK, fontsize=11, pad=9, loc="left")

    fig.suptitle("De novo miniprotein binders against human IL-6 — in-silico campaign summary",
                 color=INK, fontsize=15.5, x=0.055, ha="left", y=0.955)
    fig.text(0.055, 0.865,
             "480 designs across two epitopes · Boltz-2 design + independent co-fold confirmation · "
             "no experimental validation: every number is a prediction",
             color=INK3, fontsize=10, ha="left")

    fig.savefig(args.out, dpi=190, facecolor=SURFACE)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
