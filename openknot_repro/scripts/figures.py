"""Figures: does the in-silico (RibonanzaNet) OpenKnot score predict the measured one?"""
import pickle, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8983"
BLUE, ORANGE = "#2a78d6", "#eb6834"
SUCCESS = 90.0

res = pickle.load(open("/home/user/work/rnet_preds_r13.pkl", "rb"))
df = pd.DataFrame([{k: v for k, v in r.items()
                    if k not in ("pred_react", "exp_react", "seq", "struct")} for r in res])
d = df[(df.sn_filter == 1) & df.oks_exp_published.notna()].copy()
d["exp"] = d.oks_exp_published
d["hit"] = d.exp > SUCCESS

plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
                     "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
                     "axes.spines.top": False, "axes.spines.right": False})

# ---- Figure A: in-silico vs measured, one panel per round ----
rounds = sorted(d["round"].unique())
fig, axes = plt.subplots(1, len(rounds), figsize=(4.1 * len(rounds), 3.9), sharey=True)
for ax, rnd in zip(np.atleast_1d(axes), rounds):
    g = d[d["round"] == rnd]
    hb = ax.hexbin(g.oks_insilico, g.exp, gridsize=38, cmap="Blues", mincnt=1,
                   linewidths=0, extent=(40, 100, 40, 100))
    ax.axhline(SUCCESS, color=ORANGE, lw=1.5, ls="--", zorder=3)
    ax.text(41.5, SUCCESS + 1.2, "experimental success cutoff (90)",
            color=ORANGE, fontsize=7.5, va="bottom")
    rho = spearmanr(g.oks_insilico, g.exp).correlation
    ax.set_title(f"Round {rnd}   n={len(g):,}   Spearman rho = {rho:.2f}",
                 fontsize=9.5, color=INK, pad=8)
    ax.set_xlabel("in-silico OpenKnot score\n(RibonanzaNet-predicted SHAPE)")
    ax.grid(True, color="#e6e5e0", lw=0.6)
    ax.set_axisbelow(True)
    ax.set_xlim(40, 100); ax.set_ylim(40, 100)
np.atleast_1d(axes)[0].set_ylabel("measured OpenKnot score\n(experimental SHAPE)")
cb = fig.colorbar(hb, ax=np.atleast_1d(axes).tolist(), pad=0.02, fraction=0.035)
cb.set_label("designs per bin", color=INK2, fontsize=8)
cb.outline.set_edgecolor(MUTED)
fig.suptitle("In-silico screening score vs. what the wet lab actually measured",
             fontsize=11, color=INK, x=0.02, ha="left", y=0.99)
fig.savefig("/home/user/work/fig_insilico_vs_measured.png", dpi=190, bbox_inches="tight")

# ---- Figure B: success rate by in-silico score decile ----
fig2, ax = plt.subplots(figsize=(6.4, 3.9))
for rnd, color in zip(rounds, (BLUE, ORANGE)):
    g = d[d["round"] == rnd].copy()
    g["dec"] = pd.qcut(g.oks_insilico, 10, labels=False, duplicates="drop")
    s = g.groupby("dec").agg(x=("oks_insilico", "mean"), y=("hit", "mean"), n=("hit", "size"))
    ax.plot(s.x, 100 * s.y, "-o", color=color, lw=2, ms=6,
            markeredgecolor=SURFACE, markeredgewidth=1.5, label=f"Round {rnd}", zorder=3)
    ax.axhline(100 * g.hit.mean(), color=color, lw=1, ls=":", zorder=1)
    ax.text(57.5, 100 * g.hit.mean() + 1, f"Round {rnd} base rate {100 * g.hit.mean():.0f}%",
            color=color, fontsize=7.5, va="bottom")
ax.set_xlim(56, 101)
ax.set_xlabel("in-silico OpenKnot score (decile mean)")
ax.set_ylabel("designs that succeeded\nin the wet lab (%)")
ax.set_title("Screening on the in-silico score does enrich for real successes",
             fontsize=10.5, color=INK, pad=8, loc="left")
ax.grid(True, color="#e6e5e0", lw=0.6); ax.set_axisbelow(True)
ax.legend(frameon=False, loc="upper left", fontsize=8.5, bbox_to_anchor=(0.0, 1.0))
fig2.savefig("/home/user/work/fig_enrichment.png", dpi=190, bbox_inches="tight")
print("wrote /home/user/work/fig_insilico_vs_measured.png and fig_enrichment.png")
