"""
LF6 — Dielectric spectroscopy / membrane capacitance as a cell-state readout.

No open raw single-cell DEP dataset was available, so the demo uses the
standard single-shell (Maxwell–Wagner) model, parameterised ONLY with
published measurements (PMIDs below), and first checks that the model
reproduces the measured crossover frequencies.

Physics: a cell = conducting cytoplasm inside a thin capacitive membrane.
    Clausius–Mossotti factor  K(f) = (ε*_p − ε*_m)/(ε*_p + 2ε*_m)
    first DEP crossover        f_x0 ≈ σ_m / (√2 π r C_m)
so f_x0 reports radius r and specific membrane capacitance C_m (membrane
folding / microvilli / blebbing), NOT gene expression directly.

Published inputs
  Lannin 2016 Biomicrofluidics (PMID 27990211): PANC-1 ± TGF-β+EGF (EMT),
      BxPC3, white blood cells; σ_m = 70 mS/m; measured f_x0 and C_m.
  Pethig 2002 Electrophoresis (PMID 12210259): human T cells, PMA+ionomycin
      24 h: C_m 13.49 ± 4.72 -> 10.62 ± 5.13 mF/m².
  Pethig & Talary 2007 IET Nanobiotechnol (PMID 17500582): Jurkat, etoposide
      6 h (apoptosis): C_m 13.34 ± 2.88 -> 10.49 ± 4.00 mF/m².
"""
import numpy as np
from scipy.stats import norm
from . import _data

EPS0 = 8.854e-12

LANNIN = {  # name: (radius µm, C_m mF/m², σ_cyto S/m, measured f_x0 kHz)
    "PANC-1": (12.1, 20.2, 0.476, 64),
    "PANC-1 + EMT": (12.7, 20.1, 0.536, 54),
    "BxPC3": (8.1, None, None, 80),
    "white blood cells": (5.6, 9.8, 0.721, 271),
}
SIGMA_M = 0.070  # S/m, Lannin 2016


def fx0(r_um, Cm_mF, sigma_m=SIGMA_M):
    return sigma_m / (np.sqrt(2) * np.pi * r_um * 1e-6 * Cm_mF * 1e-3)


def cm_factor(f, r_um, Cm_mF, sigma_c=0.5, eps_c=60, sigma_m=SIGMA_M, eps_m=78, d=5e-9):
    """Single-shell model, exact (membrane of thickness d, ε_mem = C_m·d)."""
    w = 2 * np.pi * np.asarray(f)
    r = r_um * 1e-6
    e_mem = Cm_mF * 1e-3 * d / EPS0
    em = e_mem * EPS0 + 0j                         # membrane conductance neglected
    ec = eps_c * EPS0 + sigma_c / (1j * w)
    emed = eps_m * EPS0 + sigma_m / (1j * w)
    g = (r / (r - d)) ** 3
    a = (ec - em) / (ec + 2 * em)
    ep = em * (g + 2 * a) / (g - a)
    return (ep - emed) / (ep + 2 * emed)


def separability(m1, s1, m2, s2):
    d = abs(m1 - m2) / np.sqrt((s1 ** 2 + s2 ** 2) / 2)
    auc = norm.cdf(d / np.sqrt(2))
    n = int(np.ceil(2 * (norm.ppf(0.975) + norm.ppf(0.8)) ** 2 / d ** 2))
    return d, auc, n


def run(make_fig=True):
    print("=" * 72); print("LF6  Dielectric state readout (single-shell model + literature)"); print("=" * 72)
    print("[check] model f_x0 = σ_m/(√2·π·r·C_m) vs measured (Lannin 2016, σ_m = 70 mS/m):")
    chk = {}
    for k, (r, C, sc, meas) in LANNIN.items():
        if C is None:
            continue
        # crossover from the full single-shell model (root of Re K)
        f = np.logspace(3, 7, 4000)
        ReK = cm_factor(f, r, C, sigma_c=sc).real
        i = np.where(np.diff(np.sign(ReK)) > 0)[0][0]
        full = f[i] / 1e3
        chk[k] = (fx0(r, C) / 1e3, full, meas)
        print(f"    {k:18s} r={r:4.1f} µm C_m={C:4.1f} mF/m²: approx {fx0(r, C)/1e3:5.0f} kHz, "
              f"full model {full:5.0f} kHz, measured {meas} kHz")

    r0, C0, _, f0 = LANNIN["PANC-1"]; r1, C1, _, f1 = LANNIN["PANC-1 + EMT"]
    dln_r = np.log(r1 / r0); dln_C = np.log(C1 / C0)
    print(f"\n[EMT] PANC-1 -> EMT: Δln f_x0 = −(Δln r + Δln C_m) = −({dln_r:+.3f} {dln_C:+.3f})")
    print(f"      -> {abs(dln_r)/(abs(dln_r)+abs(dln_C)):.0%} of the predicted shift is cell SIZE, "
          f"{abs(dln_C)/(abs(dln_r)+abs(dln_C)):.0%} membrane capacitance (C_m 20.2 -> 20.1 mF/m²)")

    print("\n[state changes] single-cell separability from published mean ± SD of C_m:")
    states = {"T-cell activation (PMA+iono 24 h)": (13.49, 4.72, 10.62, 5.13),
              "Jurkat apoptosis (etoposide 6 h)": (13.34, 2.88, 10.49, 4.00)}
    sep = {}
    for k, (m1, s1, m2, s2) in states.items():
        d, auc, n = separability(m1, s1, m2, s2)
        sep[k] = (d, auc, n)
        print(f"    {k:34s} Cohen d = {d:.2f}; single-cell AUC ≈ {auc:.2f}; "
              f"cells/group for population test (80% power) ≈ {n}")
        r = 5.0  # µm, assumed lymphocyte-like radius for illustration
        print(f"      at r = {r} µm: f_x0 {fx0(r, m1)/1e3:.0f} -> {fx0(r, m2)/1e3:.0f} kHz")

    if make_fig:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        f = np.logspace(3, 8, 600)
        fig, ax = plt.subplots(1, 3, figsize=(16, 4.4))
        for k, (r, C, sc, meas) in LANNIN.items():
            if C is None:
                continue
            ax[0].semilogx(f, cm_factor(f, r, C, sigma_c=sc).real, label=f"{k} (meas. {meas} kHz)")
        ax[0].axhline(0, c="k", lw=0.5); ax[0].set_xlabel("frequency (Hz)")
        ax[0].set_ylabel("Re K(f)"); ax[0].legend(fontsize=8)
        ax[0].set_title("(a) EMT: spectra nearly overlap\n(PANC-1 data, Lannin 2016)")
        ks = list(chk)
        ax[1].scatter([chk[k][2] for k in ks], [chk[k][1] for k in ks])
        for k in ks:
            ax[1].annotate(k, (chk[k][2], chk[k][1]), fontsize=8)
        lim = [30, 400]; ax[1].plot(lim, lim, "k:"); ax[1].set_xscale("log"); ax[1].set_yscale("log")
        ax[1].set_xlabel("measured f_x0 (kHz)"); ax[1].set_ylabel("single-shell model (kHz)")
        ax[1].set_title("(b) model reproduces measured crossovers")
        x = np.linspace(0, 30, 400)
        for k, (m1, s1, m2, s2) in states.items():
            l, = ax[2].plot(x, norm.pdf(x, m1, s1), label=f"{k}: before")
            ax[2].plot(x, norm.pdf(x, m2, s2), "--", c=l.get_color(), label="after")
        ax[2].set_xlabel("specific membrane capacitance C_m (mF/m²)")
        ax[2].set_title("(c) single-cell distributions overlap\n(mean ± SD from Pethig 2002/2007)")
        ax[2].legend(fontsize=7)
        fig.tight_layout(); fig.savefig(_data.fig("lf6_dielectric.png"), dpi=130)
        plt.close(fig)
    return dict(check=chk, emt=(dln_r, dln_C), sep=sep)


if __name__ == "__main__":
    run()
