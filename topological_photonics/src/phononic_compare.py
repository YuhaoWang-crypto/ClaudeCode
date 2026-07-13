"""
Figure 8 (complete package) — two-system comparison: PHOTONIC vs PHONONIC.

The identical Wu-Hu hexagonal-cluster lattice is solved with two different wave
kernels (Maxwell TM for light, scalar acoustic for sound). Both show a double
Dirac cone at R=a/3 and a band-inverted gap for R≠a/3, driven by the same C6
p(l=±1)/d(l=±2) mechanism — but the phononic gap sits at a higher normalized
frequency and on higher bands. This is the deliverable enabling the two-system
comparison (photonic chip [2,3] vs GHz elastic waveguide [4]).
"""
from __future__ import annotations
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import geometry as geo
import acoustic as ac
import symmetry as sym

FIGDIR = os.path.join(os.path.dirname(__file__), "..", "figures")


def _phase_label(l_low, l_high):
    return "topological (d below p)" if np.mean(l_low) > np.mean(l_high) \
        else "trivial (p below d)"


def run(nmax=7, save=True):
    path, ticks, keys = geo.gamma_M_K_path(n_per_seg=50)
    xt = [max(t - 1, 0) for t in ticks]
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    # ---- Photonic row (bands 0-6; Dirac doublets are bands 1-2 / 3-4) ----
    for ax, (ttl, R) in zip(axes[0], [("critical R=a/3", geo.R_CRITICAL),
                                      ("expanded R=0.36a", geo.R_EXPAND)]):
        pc = geo.wu_hu_crystal(R, nmax=nmax)
        b = pc.bands(path, 7)
        for k in range(b.shape[1]):
            ax.plot(np.arange(len(path)), b[:, k], color="#1b3b6f", lw=1.4)
        ax.set_ylim(0.30, 0.70)
        _decorate(ax, xt, keys, ttl, "PHOTONIC (TM light)")
        if R != geo.R_CRITICAL:
            f = pc.solve_k([0, 0], 6)
            ax.axhspan(f[2], f[3], color="#e8a33d", alpha=.3)
            lab = _phase_label(sym.c6_angular_momentum(pc, [1, 2]),
                               sym.c6_angular_momentum(pc, [3, 4]))
            ax.text(2, (f[2] + f[3]) / 2, f"gap={f[3]-f[2]:.3f}\n{lab}",
                    fontsize=8, va="center", bbox=dict(fc="white", alpha=.8, ec="0.6"))
        else:
            ax.plot(0, pc.solve_k([0, 0], 3)[1], "o", color="#c0392b", ms=8)

    # ---- Phononic row (Dirac doublets are bands 3-4 / 5-6). Its topological
    #      side is R<a/3 (opposite to photonic) for this material contrast. ----
    for ax, (ttl, R) in zip(axes[1], [("critical R=a/3", geo.R_CRITICAL),
                                      ("shrunken R=0.30a", geo.R_SHRUNK)]):
        c = ac.AcousticCrystal(R, nmax=nmax)
        b = c.bands(path, 8)
        for k in range(b.shape[1]):
            ax.plot(np.arange(len(path)), b[:, k], color="#6a1b3b", lw=1.4)
        ax.set_ylim(0.60, 1.30)
        _decorate(ax, xt, keys, ttl, "PHONONIC (scalar sound)")
        if R != geo.R_CRITICAL:
            f = c.solve_k([0, 0], 7)
            ax.axhspan(f[4], f[5], color="#e8a33d", alpha=.3)
            lab = _phase_label(sym.c6_angular_momentum(c, [3, 4]),
                               sym.c6_angular_momentum(c, [5, 6]))
            ax.text(2, (f[4] + f[5]) / 2, f"gap={f[5]-f[4]:.3f}\n{lab}",
                    fontsize=8, va="center", bbox=dict(fc="white", alpha=.8, ec="0.6"))
        else:
            ff = c.solve_k([0, 0], 7)
            ax.plot(0, ff[4], "o", color="#c0392b", ms=8)

    fig.suptitle("Engine 5 — Two-system comparison: same Wu-Hu lattice, two wave kernels.\n"
                 "Both host a double Dirac cone at R=a/3 and a C6-inverted topological gap "
                 "(phononic gap on higher bands, higher f).", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    if save:
        out = os.path.join(FIGDIR, "fig8_photonic_vs_phononic.png")
        fig.savefig(out, dpi=130)
        print("saved", out)


def _decorate(ax, xt, keys, ttl, sysname):
    for t in xt:
        ax.axvline(t, color="0.85", lw=.7)
    ax.set_xticks(xt); ax.set_xticklabels(keys)
    ax.set_xlim(0, xt[-1])
    ax.set_ylabel(r"$\omega a / 2\pi c$")
    ax.set_title(f"{sysname} — {ttl}", fontsize=10)


if __name__ == "__main__":
    run()
