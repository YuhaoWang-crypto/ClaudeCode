"""
Figures 12 & 13 — ref-[4]-class CASE STUDY: a full-vector elastic valley-Hall
topological waveguide of air holes in silicon (GHz-range on a real Si platform).

fig12 : elastic band structure (P-SV) — Dirac cone at K for the inversion-
        symmetric honeycomb, and the gap that opens at K when the two holes
        differ (broken inversion) -> valley-Hall phase.
fig13 : Berry curvature over the Brillouin zone (mass-metric FHS on the lower
        Dirac band), concentrated with opposite signs at the K and K' valleys;
        the valley Chern (C_v ~ ±½ in magnitude but NOT strictly quantized —
        it depends on the integration patch; total Chern ≈ 0). This is the
        topological invariant appropriate to the broken-inversion / T-symmetric
        class — NOT a global Chern number.
"""
from __future__ import annotations
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import elastic_pwe as ep

FIGDIR = os.path.join(os.path.dirname(__file__), "..", "figures")
LOWER_BAND = 1  # lower band of the K Dirac pair (band 0 is the acoustic branch)


def bands_figure(save=True):
    path, ticks, keys = ep.KMGamma_path(n=55)
    xt = [max(t - 1, 0) for t in ticks]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), sharey=True)
    for ax, (tag, (rA, rB)) in zip(
        axes, [("symmetric holes\nDirac cone at K", (0.28, 0.28)),
               ("broken inversion (r$_A$≠r$_B$)\ngap at K → valley-Hall", (0.38, 0.16))]):
        c = ep.ElasticCrystal(rA, rB, nmax=9)
        b = c.bands(path, 6)
        for k in range(b.shape[1]):
            ax.plot(np.arange(len(path)), b[:, k], color="#6a1b3b", lw=1.5)
        for t in xt:
            ax.axvline(t, color="0.85", lw=.7)
        ax.set_xticks(xt); ax.set_xticklabels(keys)
        ax.set_xlim(0, xt[-1]); ax.set_ylim(0, 0.55)
        ax.set_title(tag, fontsize=10.5)
        Kidx = ticks[1]
        fK = c.solve_k(path[Kidx], 3)
        if rA == rB:
            ax.plot(Kidx, fK[1], "o", color="#c0392b", ms=8, zorder=5)
        else:
            ax.axhspan(fK[1], fK[2], color="#e8a33d", alpha=.3)
            ax.text(2, (fK[1] + fK[2]) / 2, f"gap={fK[2]-fK[1]:.3f}", fontsize=8.5,
                    va="center", bbox=dict(fc="white", alpha=.8, ec="0.6"))
    axes[0].set_ylabel(r"Frequency  $f=\omega a/2\pi c_t$")
    fig.suptitle("Case study (ref-[4] class) — full-vector ELASTIC valley-Hall crystal: "
                 "air holes in silicon, GHz Lamb/P-SV waves", fontsize=11)
    fig.tight_layout()
    if save:
        out = os.path.join(FIGDIR, "fig12_elastic_bands.png")
        fig.savefig(out, dpi=130)
        print("saved", out)


def berry_patch(c, center, band, half=1.3, Nk=22):
    """
    Single-band (mass-metric) FHS Berry curvature on a small cartesian patch
    centred at `center`, AWAY from the Γ acoustic degeneracy. Returns the grid,
    the curvature, and the integrated flux / 2π (the valley Chern of this valley).
    """
    ax = np.linspace(-half, half, Nk)
    KX, KY = np.meshgrid(center[0] + ax, center[1] + ax, indexing="ij")
    vecs = np.empty((Nk, Nk), dtype=object)
    R = None
    for i in range(Nk):
        for j in range(Nk):
            _, V, R = c.solve_k([KX[i, j], KY[i, j]], band + 1, return_vecs=True)
            vecs[i, j] = V[:, band]

    def ov(u, v):
        return np.vdot(u, R @ v)
    F = np.zeros((Nk - 1, Nk - 1))
    for i in range(Nk - 1):
        for j in range(Nk - 1):
            U1 = ov(vecs[i, j], vecs[i + 1, j])
            U2 = ov(vecs[i + 1, j], vecs[i + 1, j + 1])
            U3 = ov(vecs[i + 1, j + 1], vecs[i, j + 1])
            U4 = ov(vecs[i, j + 1], vecs[i, j])
            prod = (U1/abs(U1))*(U2/abs(U2))*(U3/abs(U3))*(U4/abs(U4))
            F[i, j] = np.angle(prod)
    Cv = F.sum() / (2 * np.pi)
    return KX[:-1, :-1], KY[:-1, :-1], F, Cv


def valley_figure(rA=0.38, rB=0.16, save=True):
    c = ep.ElasticCrystal(rA, rB, nmax=5)
    b1, b2 = c.lat.b1, c.lat.b2
    K = (2 * b1 + b2) / 3.0
    Kp = (b1 + 2 * b2) / 3.0
    XK, YK, FK, CvK = berry_patch(c, K, LOWER_BAND)
    XKp, YKp, FKp, CvKp = berry_patch(c, Kp, LOWER_BAND)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.7))
    norm = (2 * np.pi / 22) ** 2  # patch cell area factor (approximate)
    vmax = max(np.abs(FK).max(), np.abs(FKp).max()) / norm
    for ax, (X, Y, F, Cv, name) in zip(axes, [
            (XK, YK, FK, CvK, "K valley"), (XKp, YKp, FKp, CvKp, "K' valley")]):
        sc = ax.pcolormesh(X, Y, F / norm, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                           shading="auto")
        ax.set_title(f"{name}:  valley Chern C_v = {Cv:+.2f}", fontsize=10.5)
        ax.set_xlabel("$k_x$"); ax.set_ylabel("$k_y$"); ax.set_aspect("equal")
        plt.colorbar(sc, ax=ax, fraction=.046)
    fig.suptitle("Elastic valley-Hall topology — Berry curvature concentrated with "
                 f"OPPOSITE sign at K / K' (C_v≈{CvK:+.2f} / {CvKp:+.2f}; total ≈ 0).\n"
                 "The correct invariant for the broken-inversion, T-symmetric class.",
                 fontsize=10.5)
    fig.tight_layout()
    if save:
        out = os.path.join(FIGDIR, "fig13_elastic_valley_chern.png")
        fig.savefig(out, dpi=130)
        print("saved", out)
    print(f"valley Chern: C_v(K)={CvK:+.3f}  C_v(K')={CvKp:+.3f}  sum={CvK+CvKp:+.3f}")
    return CvK, CvKp


def run():
    bands_figure()
    valley_figure()


if __name__ == "__main__":
    run()
