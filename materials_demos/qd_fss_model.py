"""
Demo #1 — Exciton Fine-Structure Splitting (FSS) of self-assembled quantum
dots under uniaxial stress.

Analytical reproduction of the model in:
    M. Gong, W. Zhang, G.-C. Guo, L. He,
    "Exciton Polarization, Fine-Structure Splitting, and the Asymmetry of
     Quantum Dots under Uniaxial Stress",
    Phys. Rev. Lett. 106, 227401 (2011).

WHAT THIS IS
------------
A faithful, end-to-end reproduction of the *analytical* two-level bright-exciton
model (their Eq. (2)) and the resulting relations between:
    - applied uniaxial stress  p
    - exciton polarization angle  theta
    - fine-structure splitting  FSS = |E3 - E4|

This needs NO DFT engine — it is a 2x2 Hamiltonian diagonalised with numpy, so
it runs to completion in this lightweight environment and is a genuine
"can-we-do-theory" demo. It is NOT a supercell atomistic-pseudopotential
calculation (that is what the paper uses to *confirm* the analytical theory).

MODEL
-----
In the two bright-state basis, the excitonic Hamiltonian under stress p is

    H(p) = [[ Ebar + d + a3*p ,   k + b*p        ],
            [ k + b*p          ,  Ebar - d + a4*p ]]

Writing alpha = (a3 - a4)/2 (differential deformation potential) and dropping
the common diagonal shift (it does not affect the splitting), the splitting is

    FSS(p) = 2 * sqrt[ (d + alpha*p)^2 + (k + b*p)^2 ]

and the in-plane exciton polarization angle obeys

    tan(2*theta) = (k + b*p) / (d + alpha*p).

The stress that minimises the FSS (the "critical stress") and the FSS lower
bound follow in closed form:

    p_c      = -(alpha*d + b*k) / (alpha^2 + b^2)
    FSS_min  =  2 * |alpha*k - b*d| / sqrt(alpha^2 + b^2)

Two limiting cases reproduce the paper's Fig. 2:
    * b = 0      (stress // [100]):  p_c = -d/alpha,  FSS_min = 2|k|
    * alpha = 0  (stress // [110]):  p_c = -k/b,      FSS_min = 2|d|
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# --- units: energies in micro-eV (ueV), stress p in arbitrary "load" units ---

def fss(p, d, k, alpha, b):
    """Fine-structure splitting FSS(p) = 2 sqrt[(d+alpha p)^2 + (k+b p)^2]."""
    return 2.0 * np.sqrt((d + alpha * p) ** 2 + (k + b * p) ** 2)


def pol_angle_deg(p, d, k, alpha, b):
    """In-plane polarization angle theta(p) in degrees, tan(2 theta)=(k+bp)/(d+alpha p)."""
    return 0.5 * np.degrees(np.arctan2(k + b * p, d + alpha * p))


def critical_stress(d, k, alpha, b):
    return -(alpha * d + b * k) / (alpha ** 2 + b ** 2)


def fss_min(d, k, alpha, b):
    return 2.0 * abs(alpha * k - b * d) / np.sqrt(alpha ** 2 + b ** 2)


def make_figure(outfile):
    p = np.linspace(-6, 6, 1200)

    # A representative low-symmetry (C1) dot: both intrinsic terms non-zero,
    # and both stress couplings non-zero.
    d, k, alpha, b = 30.0, 18.0, 12.0, 7.0  # ueV, ueV, ueV/load, ueV/load

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))

    # (a) General C1 dot: FSS vs stress, showing the avoided-crossing minimum
    ax = axes[0, 0]
    y = fss(p, d, k, alpha, b)
    pc = critical_stress(d, k, alpha, b)
    ymin = fss_min(d, k, alpha, b)
    ax.plot(p, y, color="#2c6fbb", lw=2)
    ax.axvline(pc, ls="--", color="0.5", lw=1)
    ax.plot([pc], [ymin], "o", color="#d1495b", zorder=5)
    ax.annotate(f"critical stress\n$p_c$={pc:.2f}\nFSS$_{{min}}$={ymin:.1f} $\\mu$eV",
                xy=(pc, ymin), xytext=(pc + 1.2, ymin + 40),
                fontsize=9, color="#d1495b",
                arrowprops=dict(arrowstyle="->", color="#d1495b"))
    ax.set_title("(a) General $C_1$ dot — FSS is minimised, not zeroed", fontsize=10)
    ax.set_xlabel("uniaxial stress  $p$  (load units)")
    ax.set_ylabel("FSS  (μeV)")
    ax.grid(alpha=0.3)

    # (b) Two special cases from the paper's Fig. 2
    ax = axes[0, 1]
    y_b0 = fss(p, d, k, alpha, 0.0)      # stress // [100]
    y_a0 = fss(p, d, k, 0.0, b)          # stress // [110]
    ax.plot(p, y_b0, color="#2c6fbb", lw=2,
            label=r"$\beta=0$ (// [100]): min $=2|\kappa|$=%.0f" % (2 * abs(k)))
    ax.plot(p, y_a0, color="#e08e0b", lw=2,
            label=r"$\alpha=0$ (// [110]): min $=2|\delta|$=%.0f" % (2 * abs(d)))
    ax.axhline(2 * abs(k), ls=":", color="#2c6fbb", lw=1)
    ax.axhline(2 * abs(d), ls=":", color="#e08e0b", lw=1)
    ax.set_title("(b) Limiting cases reproduce PRL Fig. 2", fontsize=10)
    ax.set_xlabel("uniaxial stress  $p$")
    ax.set_ylabel("FSS  (μeV)")
    ax.legend(fontsize=8, loc="upper center")
    ax.grid(alpha=0.3)

    # (c) Polarization angle rotates through the critical stress
    ax = axes[1, 0]
    th = pol_angle_deg(p, d, k, alpha, b)
    ax.plot(p, th, color="#2a9d8f", lw=2)
    ax.axvline(pc, ls="--", color="0.5", lw=1)
    ax.set_title("(c) Exciton polarization angle $\\theta(p)$", fontsize=10)
    ax.set_xlabel("uniaxial stress  $p$")
    ax.set_ylabel("polarization angle  θ  (deg)")
    ax.grid(alpha=0.3)

    # (d) The paper's central message: which dots can be stress-tuned to ZERO FSS.
    #     FSS_min = 2|alpha*k - beta*d|/sqrt(alpha^2+b^2) vanishes exactly when
    #     kappa/delta = beta/alpha. Sweeping kappa shows a special kappa* at which
    #     the achievable lower bound touches zero; away from it a residual remains.
    ax = axes[1, 1]
    kk = np.linspace(0, 40, 300)
    lb = np.array([fss_min(d, kv, alpha, b) for kv in kk])
    k_star = b * d / alpha  # condition alpha*k = beta*d  ->  FSS_min = 0
    ax.plot(kk, lb, color="#8d5bbd", lw=2)
    ax.fill_between(kk, 0, lb, alpha=0.12, color="#8d5bbd")
    ax.axvline(k_star, ls="--", color="0.5", lw=1)
    ax.plot([k_star], [0], "o", color="#d1495b", zorder=5)
    ax.set_title("(d) Which dots can be tuned to zero FSS", fontsize=10)
    ax.set_xlabel(r"intrinsic coupling $\kappa$  (asymmetry proxy, μeV)")
    ax.set_ylabel("achievable FSS$_{min}$  (μeV)")
    ax.annotate(r"$\kappa^*=\beta\delta/\alpha$: FSS $\to 0$"
                "\n(dot tunable to a perfect\nentangled-photon emitter)",
                xy=(k_star, 0), xytext=(k_star + 2, 12),
                fontsize=8, color="#d1495b",
                arrowprops=dict(arrowstyle="->", color="#d1495b"))
    ax.grid(alpha=0.3)

    fig.suptitle("Demo #1 · Exciton FSS of quantum dots under uniaxial stress\n"
                 "analytical reproduction of Gong et al., PRL 106, 227401 (2011)",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(outfile, dpi=150)
    print(f"saved figure -> {outfile}")

    # console summary
    print("\n--- analytical results (general C1 dot) ---")
    print(f"  parameters: delta={d}, kappa={k}, alpha={alpha}, beta={b}  (ueV / ueV·load^-1)")
    print(f"  critical stress p_c   = {pc:.4f}")
    print(f"  FSS minimum           = {ymin:.4f} ueV")
    print(f"  FSS at zero stress    = {fss(0.0, d, k, alpha, b):.4f} ueV")
    print("  special cases:")
    print(f"    beta=0  (//[100]) -> p_c={critical_stress(d,k,alpha,0.0):.3f}, FSS_min={fss_min(d,k,alpha,1e-9):.3f} (=2|k|={2*abs(k):.1f})")
    print(f"    alpha=0 (//[110]) -> p_c={critical_stress(d,k,0.0,b):.3f}, FSS_min={fss_min(d,k,1e-9,b):.3f} (=2|d|={2*abs(d):.1f})")


if __name__ == "__main__":
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    make_figure(os.path.join(here, "qd_fss_demo.png"))
