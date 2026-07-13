"""
Figure 1 — Photonic band structure of the Wu-Hu crystal (real PWE calculation):
double Dirac cone at the critical radius, and the gap that opens on either side,
with the C6 band-inversion diagnosis printed for the trivial vs topological gap.
"""
from __future__ import annotations
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import geometry as geo
import symmetry as sym

FIGDIR = os.path.join(os.path.dirname(__file__), "..", "figures")


def gamma_gap(pc):
    f = pc.solve_k([0, 0], 6)
    # gap between the two doublets (bands 1-2 and 3-4, 0-indexed)
    return f[3] - f[2]


def run(nmax=9, save=True):
    cases = [
        ("Shrunken  R=0.30a\n(trivial)", geo.R_SHRUNK),
        ("Critical  R=a/3\n(double Dirac cone)", geo.R_CRITICAL),
        ("Expanded  R=0.36a\n(topological)", geo.R_EXPAND),
    ]
    path, ticks, keys = geo.gamma_M_K_path(n_per_seg=60)
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.6), sharey=True)
    diagnostics = {}
    for ax, (title, R) in zip(axes, cases):
        pc = geo.wu_hu_crystal(R, nmax=nmax)
        b = pc.bands(path, nbands=7)
        for k in range(b.shape[1]):
            ax.plot(np.arange(len(path)), b[:, k], color="#1b3b6f", lw=1.6)
        for t in ticks:
            ax.axvline(t - 1, color="0.85", lw=.8)
        ax.set_xticks([max(t - 1, 0) for t in ticks])
        ax.set_xticklabels(keys)
        ax.set_title(title, fontsize=10.5)
        ax.set_xlim(0, len(path) - 1)
        ax.set_ylim(0.30, 0.75)
        # annotate gap / Dirac at Gamma (index 0)
        g = gamma_gap(pc)
        if g < 0.01:
            ax.plot(0, pc.solve_k([0, 0], 3)[1], "o", color="#c0392b", ms=8, zorder=5)
        else:
            f = pc.solve_k([0, 0], 6)
            ax.axhspan(f[2], f[3], color="#e8a33d", alpha=.28)
            d = sym.diagnose(pc)
            diagnostics[title] = d
            ax.text(2, (f[2] + f[3]) / 2, f"gap={g:.3f}\n{d['phase']}",
                    fontsize=8.5, va="center",
                    bbox=dict(fc="white", ec="0.6", alpha=.8))
    axes[0].set_ylabel(r"Frequency  $\omega a / 2\pi c$")
    fig.suptitle("Wu-Hu photonic topological crystal — TM bands "
                 "(Si rods ε=11.7, honeycomb cluster) : Dirac cone → gap opening → band inversion",
                 fontsize=11)
    fig.tight_layout()
    if save:
        out = os.path.join(FIGDIR, "fig1_bands_dirac_gap.png")
        fig.savefig(out, dpi=130)
        print("saved", out)
    return diagnostics


if __name__ == "__main__":
    d = run()
    print("\nBand-inversion diagnosis (from real PWE modes, C6 symmetry indicator):")
    for k, v in d.items():
        tag = k.replace("\n", " ")
        print(f"  {tag:38s} lower|l|={v['l_lower']}  upper|l|={v['l_upper']} -> {v['phase']}")
