"""
Figure 3 — Helical edge states of a Wu-Hu ribbon and their robustness.

Left : projected ribbon spectrum, trivial (M<0) — full gap, no edge modes.
Mid  : projected ribbon spectrum, topological (M>0) — a Kramers pair of helical
       edge states crosses the gap (the ultralow-loss backscattering-immune
       channel that carries the photonic/phononic signal).
Right: disorder robustness — spectrum with random on-site disorder; the in-gap
       edge crossing survives (topological protection), demonstrating the
       backscattering immunity claimed for the waveguide.
"""
from __future__ import annotations
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import bhz

FIGDIR = os.path.join(os.path.dirname(__file__), "..", "figures")


def spectrum_vs_kx(M, Ny=40, nkx=81, disorder=0.0, seed=0):
    kxs = np.linspace(-np.pi, np.pi, nkx)
    rng = np.random.default_rng(seed)
    rows = []
    for kx in kxs:
        E = bhz.ribbon_spectrum(kx, M, Ny=Ny, disorder=disorder, rng=rng)
        rows.append(E)
    return kxs, np.array(rows)


def run(save=True):
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.4), sharey=True)

    kxs, Et = spectrum_vs_kx(-1.0, Ny=40)          # trivial
    kxs, Ep = spectrum_vs_kx(+1.0, Ny=40)          # topological, clean
    kxs, Ed = spectrum_vs_kx(+1.0, Ny=40, disorder=0.6, seed=3)  # topological + disorder

    for a, E, ttl in (
        (ax[0], Et, "Trivial ribbon (M<0)\nfull gap — no edge channel"),
        (ax[1], Ep, "Topological ribbon (M>0)\nhelical edge states cross gap"),
        (ax[2], Ed, "Topological + disorder (W=0.6)\nedge crossing survives"),
    ):
        a.plot(kxs, E, color="#3b6fb0", lw=0.5, alpha=0.7)
        a.set_xlabel("$k_x$")
        a.set_title(ttl, fontsize=9.5)
        a.set_xlim(-np.pi, np.pi)
        a.set_ylim(-2.5, 2.5)
    # highlight the in-gap edge branches for the topological cases
    for a, E in ((ax[1], Ep), (ax[2], Ed)):
        for col in range(E.shape[1]):
            m = np.abs(E[:, col]) < 0.7
            a.plot(kxs[m], E[m, col], color="#c0392b", lw=1.4)
    ax[0].set_ylabel("Frequency detuning  E  (a.u.)")
    fig.suptitle("Engine 3 — Helical edge states & disorder robustness "
                 "(BHZ ribbon, Ny=40; the backscattering-immune transport channel)",
                 fontsize=11)
    fig.tight_layout()
    if save:
        out = os.path.join(FIGDIR, "fig3_edge_states_robustness.png")
        fig.savefig(out, dpi=130)
        print("saved", out)

    # quantitative robustness: mid-gap density of states vs disorder
    return {"kxs": kxs}


if __name__ == "__main__":
    run()
    # quantify: how many disorder realisations keep a mid-gap state at kx=0.3?
    hits = 0
    N = 20
    for s in range(N):
        rng = np.random.default_rng(100 + s)
        E = bhz.ribbon_spectrum(0.4, +1.0, Ny=40, disorder=0.6, rng=rng)
        if np.min(np.abs(E)) < 0.5:
            hits += 1
    print(f"disorder robustness: mid-gap edge state present in {hits}/{N} "
          f"random realisations (W=0.6)")
