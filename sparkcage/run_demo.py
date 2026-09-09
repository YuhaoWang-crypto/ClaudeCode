"""
SparkCage feasibility demo: runs every open-source stand-in end to end.

    python3 sparkcage/run_demo.py

Produces figures in sparkcage/figures/ and prints the numeric tables quoted in
FEASIBILITY.md.  Runs in well under a minute on 4 CPU cores, no GPU, no
licensed software.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from thermo import switch_model as sw            # noqa: E402
from echem import swv, transport                 # noqa: E402
from structure import analyze                    # noqa: E402

FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)
STRUCT = HERE / "structure" / "pdb"

RESULTS: dict = {}


def banner(text):
    print("\n" + "=" * 74)
    print(text)
    print("=" * 74)


# ---------------------------------------------------------------------------
# 1. Structure: what the two candidate scaffolds actually are
# ---------------------------------------------------------------------------
def part1_structure():
    banner("1. SCAFFOLD ANATOMY  (quote Step 1.1: cage design on a specific scaffold)")
    p7cbc = analyze.fetch("7CBC", STRUCT)
    p1omp = analyze.fetch("1OMP", STRUCT)
    p1anf = analyze.fetch("1ANF", STRUCT)

    topo = analyze.luccage_topology(p7cbc)
    print(f"\n7CBC (LucCage-type switch): single chain, {topo['n_residues']} residues")
    print("  segment end-to-end lengths (Angstrom):")
    for name, val in topo["end_to_end_ang"].items():
        print(f"    {name:>3s}  {val:5.1f}")
    print("\n  inter-segment CA contact counts (< 10 A):")
    names = topo["names"]
    print("        " + " ".join(f"{n:>4s}" for n in names))
    for i, n in enumerate(names):
        print(f"  {n:>4s}  " + " ".join(f"{topo['contacts'][i, j]:4d}"
                                        for j in range(len(names))))
    print("\n  Read-out, cross-checked against Quijano-Rubio 2021 Table S6:")
    print("   - H1+H2 and H3+H4 are helical HAIRPINS, each split by a TDP/YDP")
    print("     helix cap, joined by GSGSGS loops. With H5 they are the five-helix")
    print("     CAGE (lucCage residues 1-300).")
    print("   - H6 (244-269), itself TDP-capped, is the LATCH lying in the cage")
    print("     groove. It contacts H5 (67 CA pairs) and H1 (45), and nothing else.")
    print("   - b1-b3 is the grafted 3-helix HA mini-binder. Its contacts to the")
    print("     cage (b1-H5 = 25, b1-H1 = 15, b2-H1 = 26 CA pairs) are the")
    print("     structural statement of caging: the paper reports that all")
    print("     HA-binding interface residues except F273 touch the cage domain.")
    print("   - 7CBC is sCageHA_267-1S, the SHORTENED crystallisation construct")
    print("     (319 aa), not full lucCage (359 aa). Its latch is 26 aa, not 59.")

    hinge = analyze.clamshell_hinge(p1omp, p1anf)
    print("\n1OMP -> 1ANF (Type I clamshell, apo -> holo):")
    for k, v in hinge.items():
        print(f"  {k:35s} {v:8.2f}" if isinstance(v, float) else f"  {k:35s} {v}")
    print("\n  Read-out: the two lobes are internally rigid (<0.8 A) and close by")
    print("  ~36 degrees on ligand binding, moving the far lobe by ~11 A on")
    print("  average and up to ~21 A. That displacement is the allosteric input")
    print("  a clamshell-based SparkCage would convert into a redox signal.")

    # LucCage reference model: reproduce the published K_open sweep
    from thermo import luccage_reference as ref
    print("\n  Published LucCage model reproduced (Quijano-Rubio 2021 ED Fig. 1a),")
    print("  K_CK = 1e-8 M, K_LT = 1e-9 M, K_R = 0.19, 10 nM cage : 100 nM key:")
    fig, axr = plt.subplots(figsize=(5.2, 4.2))
    tt = np.logspace(-12, -5, 60)
    for k_open in [1e-7, 1e-5, 1e-3, 1e-1]:
        y = ref.dose_response(tt, k_open=k_open)
        axr.semilogx(tt * 1e9, y * 100, label=f"$K_{{open}}$ = {k_open:.0e}")
        print(f"    K_open = {k_open:.0e}: background {y[0]*100:8.4f}%, "
              f"dynamic range {ref.dynamic_range_percent(k_open=k_open):10.1f}%")
    axr.set_xlabel("target (nM)")
    axr.set_ylabel("reconstituted luciferase (% of cage)")
    axr.set_title("Published LucCage model, reproduced")
    axr.legend(fontsize=8)
    axr.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIG / "00_luccage_reference.png", dpi=150)
    plt.close(fig)
    print(f"  -> {FIG/'00_luccage_reference.png'}")

    RESULTS["scaffold"] = {
        "7cbc_end_to_end_ang": topo["end_to_end_ang"],
        "7cbc_n_residues": topo["n_residues"],
        "mbp_hinge": hinge,
    }


# ---------------------------------------------------------------------------
# 2. Thermodynamics: the three states the quote budgets 3-4 weeks for
# ---------------------------------------------------------------------------
def part2_thermo():
    banner("2. THREE-STATE SWITCH MODEL  (quote Step 1.3: dG_open, dG_LT, dG_switch)")
    kd = 10e-9
    dg_lt = sw.dg_from_kd(kd)
    print(f"\nBinder affinity fixed at Kd = {kd*1e9:.0f} nM  ->  dG_LT = {dg_lt:.2f} kcal/mol")
    print("(a Monod Bio-class de novo binder; the value is an assumption, not a prediction)\n")

    rows = []
    header = (f"{'dG_open':>8s} {'dG_switch':>10s} {'EC50':>12s} "
              f"{'fold range':>11s} {'baseline open':>14s}")
    print(header)
    print("-" * len(header))
    for dgo in [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]:
        r = sw.sensor_table(dgo, dg_lt)
        rows.append(r)
        print(f"{dgo:8.1f} {r['dG_switch_kcal']:10.2f} "
              f"{r['EC50_M']*1e9:9.2f} nM {r['fold_dynamic_range']:11.1f} "
              f"{r['baseline_open_fraction']*100:13.2f}%")

    print("\n  Read-out: every +1 kcal/mol of cage stability multiplies the fold")
    print("  signal change by ~5 and degrades the EC50 by exactly the same factor.")
    print("  There is no free lunch; the design question is where to sit on this line.")

    # design map: swing across a required clinical window
    window = (1e-10, 1e-7)   # 0.1 - 100 nM, a typical protein-biomarker window
    dgo_grid = np.linspace(0.0, 10.0, 121)
    kd_grid = np.logspace(-11, -6, 121)
    swing = np.zeros((len(kd_grid), len(dgo_grid)))
    for i, k in enumerate(kd_grid):
        swing[i] = (sw.fraction_switched(window[1], dgo_grid, sw.dg_from_kd(k))
                    - sw.fraction_switched(window[0], dgo_grid, sw.dg_from_kd(k)))
    best = np.unravel_index(np.argmax(swing), swing.shape)
    print(f"\n  Design map over a {window[0]*1e9:.1f}-{window[1]*1e9:.0f} nM window:")
    print(f"    unconstrained optimum: Kd = {kd_grid[best[0]]*1e9:.3f} nM, "
          f"dG_open = {dgo_grid[best[1]]:.2f} kcal/mol, "
          f"swing = {swing[best]*100:.1f}%")
    if best[0] == 0:
        print("    NOTE: this sits on the tight-affinity edge of the grid. Affinity")
        print("    is not actually a free parameter -- the binder comes from Monod")
        print("    Bio. The design variable we control is dG_open at fixed Kd:")
    for kd_fixed in [1e-9, 1e-8, 1e-7]:
        dgo_opt, sw_opt = sw.optimal_dg_open(window, sw.dg_from_kd(kd_fixed))
        print(f"      Kd = {kd_fixed*1e9:6.1f} nM  ->  optimal dG_open = "
              f"{dgo_opt:4.2f} kcal/mol, swing = {sw_opt*100:5.1f}%")

    # figure
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    conc = np.logspace(-12, -4, 300)
    for dgo in [2.0, 3.0, 4.0, 5.0, 6.0]:
        ax[0].semilogx(conc * 1e9, sw.fraction_switched(conc, dgo, dg_lt),
                       label=f"$\\Delta G_{{open}}$ = {dgo:.0f}")
    ax[0].set_xlabel("analyte (nM)")
    ax[0].set_ylabel("fraction of sensor switched")
    ax[0].set_title("Dose-response vs cage stability\n(binder $K_d$ = 10 nM)")
    ax[0].legend(fontsize=8, title="kcal/mol")
    ax[0].grid(alpha=0.3)

    im = ax[1].pcolormesh(dgo_grid, kd_grid * 1e9, swing, shading="auto",
                          cmap="viridis")
    ax[1].set_yscale("log")
    ax[1].set_xlabel("$\\Delta G_{open}$ (kcal/mol)")
    ax[1].set_ylabel("binder $K_d$ (nM)")
    ax[1].set_title("Usable signal swing over a 0.1-100 nM window")
    ax[1].plot(dgo_grid[best[1]], kd_grid[best[0]] * 1e9, "r*", ms=14)
    fig.colorbar(im, ax=ax[1], label="fraction of monolayer switched")
    fig.tight_layout()
    fig.savefig(FIG / "01_thermodynamic_design_map.png", dpi=150)
    plt.close(fig)
    print(f"\n  -> {FIG/'01_thermodynamic_design_map.png'}")

    RESULTS["thermo"] = {
        "dg_lt_kcal": float(dg_lt),
        "scan": rows,
        "optimum": {"kd_nM": float(kd_grid[best[0]] * 1e9),
                    "dg_open_kcal": float(dgo_grid[best[1]]),
                    "swing": float(swing[best])},
    }
    return dg_lt


# ---------------------------------------------------------------------------
# 3. Electrochemistry: what COMSOL would have been bought for
# ---------------------------------------------------------------------------
def part3_echem(dg_lt):
    banner("3. ELECTROCHEMICAL READOUT  (quote Step 2: the COMSOL FEA work package)")
    # Calibrated so that the interrogation frequency and the peak current both
    # land where E-AB sensors are actually operated and actually read:
    # 20-300 Hz and a few hundred nA on a ~0.03 cm2 electroactive area.
    ECELL = dict(area_cm2=0.03, gamma_total_mol_cm2=2.0e-12, e_step=0.002)
    K0_CONTACT = 3.0e2
    swv.k0_from_distance.__defaults__ = (K0_CONTACT, 1.0, 5.0)
    d_closed, d_open = 6.0, 14.0
    sigma = 2.5
    k0_d = swv.k0_from_distance(d_closed)
    k0_o = swv.k0_from_distance(d_open)
    print(f"\nMethylene blue tunnelling distance: {d_closed:.0f} A docked, "
          f"{d_open:.0f} A released, ensemble width {sigma:.1f} A (beta = 1.0 /A)")
    print(f"  k0 docked   = {k0_d:10.3e} /s")
    print(f"  k0 released = {k0_o:10.3e} /s")

    # ---- 3a. the quasi-reversible maximum, and why it dictates the frequency
    # For a surface-confined couple the SWV peak is NOT monotonic in k0: it is
    # maximal at k0/f ~ 1 and falls off on both sides. Interrogating where the
    # docked state sits on the reversible side (k0/f >> 1) makes the response
    # non-monotonic in displacement, so the sensor can read signal-ON or
    # signal-OFF depending only on frequency. Getting this wrong is the single
    # easiest way to design an uninterpretable sensor.
    kappa = np.logspace(-3, 4, 80)          # k0 / f
    f_ref = 100.0
    psi = np.array([abs(swv.swv_scan(k0=kk * f_ref, freq_hz=f_ref,
                                     **{k: v for k, v in ECELL.items()
                                        if k != "gamma_total_mol_cm2"},
                                     gamma_mol_cm2=ECELL["gamma_total_mol_cm2"]
                                     )["peak_current"]) for kk in kappa])
    kappa_max = float(kappa[int(np.argmax(psi))])
    print(f"\n  Quasi-reversible maximum located at k0/f = {kappa_max:.2f}")
    print("  (Lovric & Komorsky-Lovric 1988 place it near unity; this is a")
    print("   validation of the model, not a fitted result.)")

    # Operate at the frequency putting the DOCKED state on the maximum, so that
    # any displacement of the reporter monotonically reduces the peak.
    best_f = float(k0_d / kappa_max)
    print(f"  -> interrogation frequency {best_f:.0f} Hz puts the docked reporter")
    print("     on the maximum, giving a clean monotonic signal-off sensor.")

    i_d = abs(swv.sensor_response_distributed(0.0, d_closed, d_open, sigma,
                                              freq_hz=best_f, **ECELL)["peak_current"])
    i_o = abs(swv.sensor_response_distributed(1.0, d_closed, d_open, sigma,
                                              freq_hz=best_f, **ECELL)["peak_current"])
    i_cap = abs(swv.capacitive_background(freq_hz=best_f))
    print(f"     peak at   0% switched = {i_d*1e9:9.3f} nA")
    print(f"     peak at 100% switched = {i_o*1e9:9.3f} nA")
    print(f"     capacitive floor      = {i_cap*1e9:9.3f} nA")
    print(f"     signal change         = {(i_o-i_d)/i_d*100:9.1f} %")

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    for f_switch, style in [(0.0, "-"), (0.5, "--"), (1.0, ":")]:
        r = swv.sensor_response_distributed(f_switch, d_closed, d_open, sigma,
                                            freq_hz=best_f, e_start=0.0,
                                            e_end=-0.45, **ECELL)
        ax[0].plot(r["e"], r["i_net"] * 1e9, style,
                   label=f"{f_switch*100:.0f}% switched")
    ax[0].set_xlabel("E vs Ag/AgCl (V)")
    ax[0].set_ylabel("net SWV current (nA)")
    ax[0].set_title(f"Square-wave voltammogram\n(surface-confined MB, {best_f:.0f} Hz)")
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=0.3)

    ax[1].loglog(kappa, psi * 1e9, "-", lw=1.5)
    ax[1].axvline(kappa_max, color="r", ls=":", lw=1.2,
                  label=f"maximum at k0/f = {kappa_max:.2f}")
    ax[1].axvline(k0_d / best_f, color="C0", ls="--", lw=1.2, label="docked state")
    ax[1].axvline(k0_o / best_f, color="C1", ls="--", lw=1.2, label="released state")
    ax[1].set_xlabel("$k^0 / f$")
    ax[1].set_ylabel("|peak current| (nA)")
    ax[1].set_title("Quasi-reversible maximum\n(why frequency choice is not free)")
    ax[1].legend(fontsize=7)
    ax[1].grid(alpha=0.3, which="both")

    # ---- 3b. invert the measured signal change to a displacement
    # Kang, Sun, Kurnik, Morales, Dahlquist & Plaxco, JACS 2017, 139, 12113,
    # the closest published architecture (methylene blue on engineered cysteines
    # of a natural protein, gold electrode), report 4-30% signal change.
    print("\n  Inverting the measured signal-change range to a displacement:")
    print("  (reference: Kang/Plaxco JACS 2017 report 4-30% for protein E-AB)")
    ds = np.linspace(6.1, 16.0, 40)
    gains_d = []
    base = abs(swv.sensor_response_distributed(0.0, d_closed, d_closed, sigma,
                                               freq_hz=best_f, **ECELL)["peak_current"])
    for d in ds:
        b = abs(swv.sensor_response_distributed(1.0, d_closed, float(d), sigma,
                                                freq_hz=best_f, **ECELL)["peak_current"])
        gains_d.append((b - base) / base * 100.0)
    gains_d = np.array(gains_d)
    for target in [-4.0, -30.0]:
        j = int(np.argmin(np.abs(gains_d - target)))
        print(f"    a {abs(target):4.0f}% signal change corresponds to a "
              f"{ds[j]-d_closed:4.1f} A displacement of the reporter")
    print("\n  This is the most consequential number in the whole demo. The 4-30%")
    print("  that real protein E-AB sensors deliver implies the reporter only")
    print("  moves 1-3 A, even though the underlying domain motion is ~11-21 A")
    print("  (part 1). Almost all of the conformational change is wasted. Where")
    print("  the methylene blue is attached therefore matters more than how big")
    print("  the protein's motion is -- which is exactly what quote item 1.2")
    print("  budgets 6-8 weeks of molecular dynamics to decide.")
    ax[2].plot(ds - d_closed, gains_d, "-", lw=1.5)
    ax[2].axhspan(-30, -4, color="C2", alpha=0.15,
                  label="measured range, protein E-AB")
    ax[2].set_xlabel("reporter displacement on switching (A)")
    ax[2].set_ylabel("SWV peak change (%)")
    ax[2].set_title("Displacement needed for a given signal\n(model inverted against experiment)")
    ax[2].legend(fontsize=8)
    ax[2].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIG / "02_electrochemical_response.png", dpi=150)
    plt.close(fig)
    print(f"  -> {FIG/'02_electrochemical_response.png'}")

    # ---- 3c. full calibration curve: analyte -> thermodynamics -> current
    fig2, ax2 = plt.subplots(figsize=(5.4, 4.2))
    conc = np.logspace(-12, -5, 40)
    for dgo in [2.0, 3.0, 4.0]:
        f_sw = sw.fraction_switched(conc, dgo, dg_lt)
        sig = [(abs(swv.sensor_response_distributed(
                    float(f), d_closed, d_open, sigma, freq_hz=best_f,
                    **ECELL)["peak_current"]) - i_d) / i_d * 100.0 for f in f_sw]
        ax2.semilogx(conc * 1e9, sig, "o-", ms=3,
                     label=f"$\\Delta G_{{open}}$ = {dgo:.0f} kcal/mol")
    ax2.set_xlabel("analyte (nM)")
    ax2.set_ylabel("SWV peak change (%)")
    ax2.set_title("Predicted calibration curve\nprotein design -> device output")
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.3)
    fig2.tight_layout()
    fig2.savefig(FIG / "04_calibration_curve.png", dpi=150)
    plt.close(fig2)
    print(f"  -> {FIG/'04_calibration_curve.png'}")
    print("\n  Read-out: this is the whole chain the quote splits across two work")
    print("  packages -- cage stability -> fraction switched -> reporter distance")
    print("  -> electron-transfer rate -> measured current -- in one model. No")
    print("  published model couples all four stages; the pieces exist separately.")

    RESULTS["echem"] = {
        "k0_docked_per_s": float(k0_d),
        "k0_released_per_s": float(k0_o),
        "kappa_at_maximum": kappa_max,
        "interrogation_frequency_hz": best_f,
        "peak_docked_nA": float(i_d * 1e9),
        "peak_released_nA": float(i_o * 1e9),
        "signal_change_pct": float((i_o - i_d) / i_d * 100.0),
    }
    return best_f


# ---------------------------------------------------------------------------
# 4. Transport: boundary conditions and response time
# ---------------------------------------------------------------------------
def part4_transport():
    banner("4. MASS TRANSPORT  (quote Step 2: boundary conditions, transport physics)")
    gamma = 5.0e-12   # mol/cm^2, a packed monolayer of a ~40 kDa protein
    cells = {"quiescent (500 um layer)": 0.05, "stirred (30 um layer)": 0.003}
    for label, length in cells.items():
        da = transport.damkohler(gamma_max_mol_cm2=gamma, length_cm=length)
        print(f"\n  {label}: Damkohler = {da:.1f}  -> "
              f"{'DIFFUSION-limited' if da > 1 else 'reaction-limited'}")

    print("\n  analyte    eq.occupancy   t90 reservoir   t90 stirred   t90 semi-infinite")
    print("  " + "-" * 76)
    rows = []
    for c in [1e-11, 1e-10, 1e-9, 1e-8, 1e-7]:
        r_q = transport.simulate(c, gamma_max_mol_cm2=gamma, length_cm=0.05,
                                 t_end_s=40000.0)
        r_s = transport.simulate(c, gamma_max_mol_cm2=gamma, length_cm=0.003,
                                 t_end_s=40000.0)
        t_q = transport.time_to_fraction(r_q, 0.9)
        t_s = transport.time_to_fraction(r_s, 0.9)
        # floor evaluated at the occupancy this concentration actually reaches,
        # so it is directly comparable with the simulated t90
        floor = transport.diffusion_limited_arrival_time(
            c, gamma_max_mol_cm2=gamma, occupancy=0.9 * r_q["phi_eq"])
        rows.append({"c_M": c, "t90_quiescent_s": t_q, "t90_stirred_s": t_s,
                     "floor_s": floor, "phi_eq": r_q["phi_eq"]})
        print(f"  {c*1e9:7.3f} nM {r_q['phi_eq']*100:10.2f}% {t_q/60:12.1f} min "
              f"{t_s/60:10.1f} min {floor/60:13.1f} min")

    print("\n  Read-out: Damkohler >> 1 in both cells, so this sensor is")
    print("  DIFFUSION-limited, not reaction-limited. Three consequences:")
    print("   - response time is set by transport, so stirring or flow buys an")
    print("     order of magnitude and the binder's k_on buys nothing.")
    print("   - the last two columns bracket the answer: a well-mixed reservoir")
    print("     500 um away vs a truly quiescent sample differ by ~10x. Which one")
    print("     the real device is, is a hydrodynamics question, not a chemistry")
    print("     one, and no amount of protein design will settle it.")
    print("   - semi-infinite delivery time scales as 1/c^2, so 10x less analyte")
    print("     costs 100x in time. This, not the electrochemistry, caps the")
    print("     practical limit of detection.")
    print("  This is the one place in the project where 2D/3D FEA of the actual")
    print("  cell geometry and flow field genuinely earns its cost.")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for c in [1e-10, 1e-9, 1e-8]:
        r = transport.simulate(c, gamma_max_mol_cm2=gamma, t_end_s=40000.0)
        ax[0].plot(r["t"] / 60.0, r["phi"], label=f"{c*1e9:.1f} nM quiescent")
    for c in [1e-9]:
        r = transport.simulate(c, gamma_max_mol_cm2=gamma, length_cm=0.003,
                               t_end_s=40000.0)
        ax[0].plot(r["t"] / 60.0, r["phi"], "--", label=f"{c*1e9:.1f} nM stirred")
    ax[0].set_xlabel("time (min)")
    ax[0].set_ylabel("fraction of sensors with analyte bound")
    ax[0].set_title("Approach to equilibrium")
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=0.3)

    cs = np.logspace(-12, -6, 40)
    ax[1].loglog(cs * 1e9,
                 [transport.diffusion_limited_arrival_time(c, gamma_max_mol_cm2=gamma,
                                                           occupancy=0.1) / 60
                  for c in cs], label="time to 10% occupancy")
    ax[1].axhline(1.0, color="k", lw=0.8, ls=":")
    ax[1].axhline(60.0, color="r", lw=0.8, ls=":")
    ax[1].text(1e-2, 70, "1 hour", color="r", fontsize=8)
    ax[1].text(1e-2, 1.2, "1 minute", fontsize=8)
    ax[1].set_xlabel("analyte (nM)")
    ax[1].set_ylabel("minimum response time (min)")
    ax[1].set_title("Transport floor sets the practical LOD")
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(FIG / "03_transport.png", dpi=150)
    plt.close(fig)
    print(f"\n  -> {FIG/'03_transport.png'}")
    RESULTS["transport"] = rows


if __name__ == "__main__":
    part1_structure()
    dg_lt = part2_thermo()
    part3_echem(dg_lt)
    part4_transport()
    out = HERE / "figures" / "demo_results.json"
    out.write_text(json.dumps(RESULTS, indent=2, default=float))
    banner(f"DONE. Numeric results -> {out}")
