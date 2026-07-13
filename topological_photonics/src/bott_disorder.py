"""
Figure 9 — Disorder topology via the real-space spin Bott index.

Two disorder sweeps (ensemble-averaged over random realisations), each a genuine
real-space invariant that needs no Brillouin zone:

  * M=+1 (clean topological): spin-Bott = 1 stays quantised, then collapses to 0
    at strong disorder — the disorder-driven destruction of topology.
  * M=-0.6 (clean trivial): spin-Bott jumps 0 -> 1 at intermediate disorder — a
    Topological Anderson Insulator (disorder-INDUCED topology), the mechanism of
    the cited photonic topological-Anderson insulator [3].

This replaces the naive "k-space Berry curvature survives disorder" reasoning
with the correct real-space diagnostic.
"""
from __future__ import annotations
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import bott

FIGDIR = os.path.join(os.path.dirname(__file__), "..", "figures")


def run(L=14, nreal=8, save=True):
    Ws = np.linspace(0, 9, 19)
    W1, m1, s1 = bott.sweep_disorder(L=L, M=1.0, Ws=Ws, nreal=nreal, seed=2)
    W2, m2, s2 = bott.sweep_disorder(L=L, M=-0.6, Ws=Ws, nreal=nreal, seed=3)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.errorbar(W1, m1, yerr=s1, fmt="o-", color="#1b3b6f", capsize=3,
                label="M=+1 (clean topological)")
    ax.errorbar(W2, m2, yerr=s2, fmt="s-", color="#c0392b", capsize=3,
                label="M=−0.6 (clean trivial → TAI)")
    ax.axhline(1, color="0.7", ls=":", lw=1)
    ax.axhline(0, color="0.7", ls=":", lw=1)
    ax.annotate("topology destroyed\nby strong disorder", (7.5, 0.7),
                fontsize=8.5, color="#1b3b6f", ha="center")
    ax.annotate("disorder-INDUCED\ntopology (TAI)", (6.4, 0.35),
                fontsize=8.5, color="#c0392b", ha="center")
    ax.set_xlabel("Anderson disorder strength  W")
    ax.set_ylabel("spin Bott index  $B_s$  (real space, no BZ)")
    ax.set_title("Disorder topology via the real-space spin Bott index\n"
                 f"(L={L} torus, {nreal} realisations, error bars = ensemble σ)",
                 fontsize=10.5)
    ax.set_ylim(-0.25, 1.25)
    ax.legend(fontsize=9, loc="center left")
    fig.tight_layout()
    if save:
        out = os.path.join(FIGDIR, "fig9_spin_bott_disorder.png")
        fig.savefig(out, dpi=130)
        print("saved", out)
    return (W1, m1), (W2, m2)


if __name__ == "__main__":
    run()
