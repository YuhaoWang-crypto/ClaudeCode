"""
Figure 10 — Loss budget: "ultralow-loss" is NOT implied by topology.

Topological protection removes BACKSCATTERING loss (bends/defects — fig7), but
does nothing about material ABSORPTION, radiation, or interface/roughness loss,
which set the real floor. We quantify the absorption floor rigorously from the
ACTUAL FDFD edge-mode profile (fig5): its confinement factor Γ_rod (fraction of
modal energy inside the rods) feeds the standard confinement-factor modal-loss
formula

        α_power = (2π/λ0) · Γ_rod · κ_rod ,   κ = Im(n) = ε''/(2√ε') ,

for three rod materials at telecom (λ0 = a/f, a≈0.71 µm so λ0≈1.55 µm):

  * low-loss dielectric (Si/SiO2)  ε''≈1e-3   — the ultralow-loss baseline
  * SiO2-shelled Bi2Se3 (thin)     ε''≈0.5    — absorption reduced by the shell
  * bare Bi2Se3                    ε''≈8      — the strongly absorbing limit

Optical constants are REPRESENTATIVE (illustrative) — the message (absorption,
set by the Bi2Se3 fraction, decides ultralow-loss, not topology) is the
deliverable; dropping in measured Bi2Se3/SiO2 constants is a data step.
"""
from __future__ import annotations
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import fdfd

FIGDIR = os.path.join(os.path.dirname(__file__), "..", "figures")
FREQ = 0.46
A_UM = FREQ * 1.55            # a in microns (=> lambda0 = a/FREQ = 1.55 um)
EPS_ROD_RE = 11.7            # real part of the rod permittivity (Si-like host)


def edge_mode_confinement(Nx=48, Ncells_y=6):
    """Confinement factor Γ_rod = (∫_rod ε|E|²)/(∫ ε|E|²) for the FDFD edge mode."""
    eps, h, Lx, Ly = fdfd.build_epsilon(Ncells_y, R_top=0.36, R_bot=0.30, Nx=Nx)
    Nxg, Nyg = eps.shape
    f, V = fdfd.strip_spectrum(1.6, eps, h, Lx, nstates=44)
    cand = [i for i, fr in enumerate(f) if 0.44 < fr < 0.485]
    ie = max(cand, key=lambda i: fdfd.edge_localization(V[:, i], Nxg, Nyg))
    E2 = np.abs(V[:, ie].reshape(Nxg, Nyg)) ** 2
    energy = eps * E2
    rod = eps > 1.0
    return energy[rod].sum() / energy.sum()


def modal_loss(gamma, epp, eps_re=EPS_ROD_RE):
    """Power attenuation α (1/a) and dB/cm from confinement-factor formula."""
    kappa = epp / (2 * np.sqrt(eps_re))          # Im(n)
    alpha_per_a = (2 * np.pi * FREQ) * gamma * kappa   # 1/a (power)
    a_cm = A_UM * 1e-4
    db_cm = 8.686 * alpha_per_a / a_cm           # 8.686 = 10/ln10
    Lp_um = (1.0 / alpha_per_a) * A_UM if alpha_per_a > 0 else np.inf
    return alpha_per_a, db_cm, Lp_um


def run(save=True):
    gamma = edge_mode_confinement()
    cases = [
        ("low-loss\nSi/SiO2\nε''≈1e-3", 1e-3, "#1b3b6f"),
        ("SiO2-shelled\nBi2Se3\nε''≈0.5", 0.5, "#2c7fb8"),
        ("bare Bi2Se3\nε''≈8", 8.0, "#c0392b"),
    ]
    db, Lp = [], []
    for _, epp, _ in cases:
        _, d, lp = modal_loss(gamma, epp)
        db.append(d); Lp.append(lp)

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
    names = [c[0] for c in cases]; cols = [c[2] for c in cases]
    b1 = ax[0].bar(names, db, color=cols)
    ax[0].set_yscale("log")
    for bar, v in zip(b1, db):
        ax[0].text(bar.get_x() + bar.get_width() / 2, v, f"{v:.2g}",
                   ha="center", va="bottom", fontsize=8.5)
    ax[0].set_ylabel("absorption α  (dB/cm)")
    ax[0].set_title(f"Absorption floor from the real edge mode\n(Γ_rod={gamma:.2f} from FDFD, "
                    f"a≈{A_UM:.2f} µm)", fontsize=10)

    b2 = ax[1].bar(names, [min(l, 1e5) for l in Lp], color=cols)
    ax[1].set_yscale("log")
    for bar, v in zip(b2, Lp):
        txt = (">1 cm" if not np.isfinite(v) or v >= 1e5
               else f"{v:.1f} µm" if v < 10 else f"{v:.0f} µm")
        ax[1].text(bar.get_x() + bar.get_width() / 2, min(v, 1e5), txt,
                   ha="center", va="bottom", fontsize=8.5)
    ax[1].set_ylabel("propagation length  L_p  (µm)")
    ax[1].set_title("Propagation length (1/e amplitude)\n"
                    "topology cannot lengthen an absorbing mode", fontsize=10)

    fig.suptitle("Engine 6 — Loss budget: the Bi2Se3 fraction (absorption), not "
                 "topology, sets the ultralow-loss ceiling", fontsize=11)
    fig.tight_layout()
    if save:
        out = os.path.join(FIGDIR, "fig10_loss_budget.png")
        fig.savefig(out, dpi=130)
        print("saved", out)
    print(f"Γ_rod = {gamma:.3f}")
    for (n, epp, _), d, lp in zip(cases, db, Lp):
        lps = f"{lp:.1f} µm" if np.isfinite(lp) and lp < 1e5 else ">1 cm"
        print(f"  {n.replace(chr(10),' '):32s} ε''={epp:<6} α={d:8.1f} dB/cm  L_p={lps}")
    return {"gamma": gamma, "db": db, "Lp": Lp}


if __name__ == "__main__":
    run()
