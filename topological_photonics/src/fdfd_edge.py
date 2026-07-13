"""
Figure 5 — REAL-ROD edge states (FDFD supercell), replacing the BHZ ribbon.

Left  : projected band structure E(k_x) of a silicon-rod strip with a
        topological|trivial domain wall; states coloured by their localisation
        at the wall. Edge branches traverse the bulk gap (helical channel).
Right : |E_z|² of one in-gap edge mode — it hugs the domain wall, decaying into
        both bulk crystals. This is the ab-initio version of the transport
        channel (no effective model; ⚠️ removed).
"""
from __future__ import annotations
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import fdfd

FIGDIR = os.path.join(os.path.dirname(__file__), "..", "figures")


def run(Nx=48, Ncells_y=6, nkx=24, nstates=44, save=True):
    eps, h, Lx, Ly = fdfd.build_epsilon(Ncells_y, R_top=0.36, R_bot=0.30, Nx=Nx)
    Nxg, Nyg = eps.shape
    kxs = np.linspace(0.05, np.pi - 0.05, nkx)

    allf, allloc = [], []
    edge_vec = None
    for kx in kxs:
        f, V = fdfd.strip_spectrum(kx, eps, h, Lx, nstates=nstates)
        loc = np.array([fdfd.edge_localization(V[:, i], Nxg, Nyg) for i in range(len(f))])
        allf.append(f)
        allloc.append(loc)
        if edge_vec is None and abs(kx - kxs[nkx // 2]) < 1e-9:
            pass
    # grab a representative edge mode at mid-kx for the field plot
    kx_mid = kxs[nkx // 2]
    fmid, Vmid = fdfd.strip_spectrum(kx_mid, eps, h, Lx, nstates=nstates)
    locmid = np.array([fdfd.edge_localization(Vmid[:, i], Nxg, Nyg) for i in range(len(fmid))])
    cand = [i for i, fr in enumerate(fmid) if 0.44 < fr < 0.485]
    ie = max(cand, key=lambda i: locmid[i])
    edge_field = np.abs(Vmid[:, ie].reshape(Nxg, Nyg)) ** 2

    fig, ax = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={"width_ratios": [1.4, 1]})

    for kx, f, loc in zip(kxs, allf, allloc):
        sc = ax[0].scatter(np.full_like(f, kx), f, c=loc, cmap="turbo",
                           s=14, vmin=0, vmax=1, edgecolors="none")
    ax[0].axhspan(0.436, 0.487, color="0.85", alpha=.5, zorder=0, label="bulk gap")
    ax[0].set_xlabel("$k_x$  (units of 1/a)")
    ax[0].set_ylabel(r"Frequency  $\omega a / 2\pi c$")
    ax[0].set_ylim(0.40, 0.52)
    ax[0].set_title("Real-rod strip: projected bands\n(colour = domain-wall localisation)",
                    fontsize=10.5)
    cb = plt.colorbar(sc, ax=ax[0], fraction=.046)
    cb.set_label("edge localisation")
    ax[0].legend(loc="upper right", fontsize=8)

    im = ax[1].imshow(edge_field.T, origin="lower", cmap="magma",
                      extent=[0, Lx, 0, Ly], aspect="auto")
    # mark the domain wall (centre and boundary)
    ax[1].axhline(Ly / 2, color="cyan", ls="--", lw=1, alpha=.7)
    ax[1].set_title(f"Edge mode |E_z|² at $k_x$={kx_mid:.2f}, f={fmid[ie]:.3f}\n"
                    "localised at the topological|trivial wall", fontsize=10.5)
    ax[1].set_xlabel("x / a"); ax[1].set_ylabel("y / a")
    plt.colorbar(im, ax=ax[1], fraction=.046)

    fig.suptitle("Engine 3b (complete package) — ab-initio FDFD edge states of the "
                 "silicon-rod crystal (⚠️ effective model removed)", fontsize=11)
    fig.tight_layout()
    if save:
        out = os.path.join(FIGDIR, "fig5_fdfd_real_edge.png")
        fig.savefig(out, dpi=130)
        print("saved", out)
    return kxs, allf, allloc


if __name__ == "__main__":
    run()
