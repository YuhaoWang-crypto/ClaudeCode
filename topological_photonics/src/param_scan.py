"""
Figure 4 — Geometry parameter optimisation / phase diagram.

Sweeps the cluster radius R/a and reports, from the real PWE calculation:
  * the Gamma gap size (band gap of the photonic crystal), and
  * the topological phase (band-inversion diagnosis).
The gap closes at the critical R = a/3 (double Dirac cone) and reopens as a
topological gap for R > a/3. This is the "basic geometry parameter scan
optimisation" deliverable: it locates the phase boundary and the R that
maximises the protected gap (largest robust operating bandwidth).
"""
from __future__ import annotations
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import geometry as geo
import symmetry as sym
from bands import gamma_gap

FIGDIR = os.path.join(os.path.dirname(__file__), "..", "figures")


def run(Rvals=None, nmax=7, save=True):
    if Rvals is None:
        Rvals = np.linspace(0.28, 0.385, 22)
    gaps, phases = [], []
    for R in Rvals:
        pc = geo.wu_hu_crystal(R, nmax=nmax)
        g = gamma_gap(pc)
        gaps.append(g)
        if g > 0.008:
            phases.append(1 if sym.diagnose(pc)["inverted"] else 0)
        else:
            phases.append(-1)  # gapless (Dirac)
    gaps = np.array(gaps)
    phases = np.array(phases)

    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    triv = phases == 0
    topo = phases == 1
    ax.plot(Rvals[triv], gaps[triv], "o-", color="#2c7fb8", label="trivial gap", ms=5)
    ax.plot(Rvals[topo], gaps[topo], "s-", color="#c0392b", label="topological gap", ms=5)
    ax.axvline(1 / 3, color="k", ls="--", lw=1, label="critical R=a/3 (Dirac)")
    ax.fill_betweenx([0, gaps.max() * 1.1], 1 / 3, Rvals.max(),
                     color="#c0392b", alpha=.06)
    ax.fill_betweenx([0, gaps.max() * 1.1], Rvals.min(), 1 / 3,
                     color="#2c7fb8", alpha=.06)
    ax.set_xlabel("cluster radius  R / a")
    ax.set_ylabel(r"$\Gamma$ band gap  $\Delta(\omega a/2\pi c)$")
    ax.set_title("Geometry scan: photonic gap vs cluster radius\n"
                 "(gap closes at R=a/3, reopens topological for R>a/3)", fontsize=10.5)
    ax.set_ylim(0, gaps.max() * 1.1)
    ax.legend(fontsize=9)
    fig.tight_layout()
    if save:
        out = os.path.join(FIGDIR, "fig4_param_scan_phase.png")
        fig.savefig(out, dpi=130)
        print("saved", out)

    # optimum: largest topological gap
    if topo.any():
        i = np.argmax(np.where(topo, gaps, -1))
        print(f"optimal topological operating point: R/a={Rvals[i]:.3f}, "
              f"gap={gaps[i]:.4f} (omega a/2pi c)")
    return Rvals, gaps, phases


if __name__ == "__main__":
    run()
