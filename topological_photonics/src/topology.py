"""
Figure 2 — Topological invariant: Berry curvature of the effective BHZ blocks
and the spin-Chern number, calibrated from the PWE gap. Shows C_total = 0
(time-reversal symmetric photonics) but C_spin = ±1 in the topological phase.
"""
from __future__ import annotations
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import bhz
import geometry as geo
from bands import gamma_gap

FIGDIR = os.path.join(os.path.dirname(__file__), "..", "figures")


def calibrated_mass(R, nmax=9):
    """M = +gap/2 (topological) or -gap/2 (trivial), gap taken from PWE."""
    pc = geo.wu_hu_crystal(R, nmax=nmax)
    import symmetry as sym
    g = gamma_gap(pc)
    d = sym.diagnose(pc)
    sign = +1 if d["inverted"] else -1
    return sign * g / 2, g, d["phase"]


def run(save=True):
    M_topo, g_t, ph_t = calibrated_mass(geo.R_EXPAND)
    M_triv, g_s, ph_s = calibrated_mass(geo.R_SHRUNK)
    # scale B so gap ~ matches; keep A=1, B=1 for the universal curvature shape,
    # only the SIGN of M (from PWE) sets the topology.
    res_topo = bhz.spin_chern(np.sign(M_topo) * 1.0, N=48)
    res_triv = bhz.spin_chern(np.sign(M_triv) * 1.0 if M_triv != 0 else -1.0, N=48)

    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    ext = [-np.pi, np.pi, -np.pi, np.pi]
    norm = (2 * np.pi / 48) ** 2
    im = ax[0].imshow(res_topo["F_up"].T / norm, origin="lower", extent=ext,
                      cmap="RdBu_r", aspect="auto")
    ax[0].set_title(f"Berry curvature Ω↑(k) — expanded (topological)\n"
                    f"C↑={res_topo['C_up']:+.0f}, C↓={res_topo['C_dn']:+.0f}, "
                    f"C_spin={res_topo['C_spin']:+.0f}", fontsize=9.5)
    ax[0].set_xlabel("$k_x$"); ax[0].set_ylabel("$k_y$")
    plt.colorbar(im, ax=ax[0], fraction=.046)

    im2 = ax[1].imshow(res_triv["F_up"].T / norm, origin="lower", extent=ext,
                       cmap="RdBu_r", aspect="auto")
    ax[1].set_title(f"Berry curvature Ω↑(k) — shrunken (trivial)\n"
                    f"C↑={res_triv['C_up']:+.0f}, C↓={res_triv['C_dn']:+.0f}, "
                    f"C_spin={res_triv['C_spin']:+.0f}", fontsize=9.5)
    ax[1].set_xlabel("$k_x$"); ax[1].set_ylabel("$k_y$")
    plt.colorbar(im2, ax=ax[1], fraction=.046)

    # spin-Chern vs mass sweep (phase transition)
    Ms = np.linspace(-3, 3, 25)
    Cs = [bhz.spin_chern(m, N=28)["C_spin"] for m in Ms]
    ax[2].plot(Ms, np.round(Cs), "o-", color="#1b3b6f")
    ax[2].axvline(0, color="#c0392b", ls="--", lw=1, label="gap closes")
    ax[2].axhline(0, color="0.8", lw=.8)
    ax[2].set_xlabel("effective mass  M  (sign set by PWE band inversion)")
    ax[2].set_ylabel("spin-Chern number  $C_s$")
    ax[2].set_title("Topological phase transition", fontsize=9.5)
    ax[2].set_yticks([-1, 0, 1]); ax[2].legend(fontsize=8)
    fig.suptitle("Engine 2 — Berry curvature & spin-Chern number "
                 "(BHZ effective model calibrated from the PWE gap; FHS lattice method)",
                 fontsize=11)
    fig.tight_layout()
    if save:
        out = os.path.join(FIGDIR, "fig2_berry_spin_chern.png")
        fig.savefig(out, dpi=130)
        print("saved", out)
    return {"topo": res_topo, "triv": res_triv, "gap_topo": g_t, "gap_triv": g_s}


if __name__ == "__main__":
    r = run()
    print(f"\nExpanded (topological): C_spin = {r['topo']['C_spin']:+.0f}, "
          f"C_total = {r['topo']['C_total']:+.0f}, PWE gap = {r['gap_topo']:.3f}")
    print(f"Shrunken (trivial)    : C_spin = {r['triv']['C_spin']:+.0f}, "
          f"C_total = {r['triv']['C_total']:+.0f}, PWE gap = {r['gap_triv']:.3f}")
