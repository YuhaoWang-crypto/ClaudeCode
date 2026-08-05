"""Shared plotting style and paths for the Particle Life experiments."""

from __future__ import annotations

import json
import os
import sys

import matplotlib
import matplotlib.font_manager

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.abspath(os.path.join(HERE, "..", "figures"))
RES = os.path.abspath(os.path.join(HERE, "..", "results"))
os.makedirs(FIG, exist_ok=True)
os.makedirs(RES, exist_ok=True)

BG = "#0d1117"
FG = "#e6edf3"
GRID = "#30363d"
ACC = ["#39d353", "#ff5ec4", "#58a6ff", "#f0883e", "#a371f7", "#e3b341", "#f85149"]

# WenQuanYi Zen Hei is present in this image and covers CJK; keep DejaVu first so
# the Latin/maths glyphs stay unchanged and only the CJK falls through to it.
_CJK = [f for f in ("WenQuanYi Zen Hei", "Noto Sans CJK SC", "Source Han Sans SC")
        if f in {x.name for x in matplotlib.font_manager.fontManager.ttflist}]

plt.rcParams.update({
    "font.sans-serif": ["DejaVu Sans", *_CJK],
    "axes.unicode_minus": False,
    "figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG,
    "text.color": FG, "axes.labelcolor": FG, "axes.edgecolor": GRID,
    "xtick.color": FG, "ytick.color": FG, "grid.color": GRID,
    "axes.titlecolor": FG, "font.size": 9, "axes.grid": True,
    "grid.alpha": 0.25, "legend.framealpha": 0.15, "figure.dpi": 130,
    "axes.prop_cycle": plt.cycler(color=ACC),
})


def save(fig, name):
    path = os.path.join(FIG, name)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {os.path.relpath(path, os.path.dirname(HERE))}")
    return path


def dump(obj, name):
    path = os.path.join(RES, name)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, default=float)
    print(f"  -> {os.path.relpath(path, os.path.dirname(HERE))}")
    return path


def species_colors(S):
    return [ACC[i % len(ACC)] for i in range(S)]


def scatter_snapshot(ax, sim, size=1.6, alpha=0.9):
    cols = species_colors(sim.p.S)
    for s in range(sim.p.S):
        m = sim.types == s
        ax.scatter(sim.pos[m, 0], sim.pos[m, 1], s=size, c=cols[s],
                   linewidths=0, alpha=alpha)
    ax.set_xlim(0, sim.p.L[0])
    ax.set_ylim(0, sim.p.L[1])
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
