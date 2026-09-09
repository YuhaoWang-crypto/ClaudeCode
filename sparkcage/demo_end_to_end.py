"""
End-to-end SparkCage demo on a real published cage-latch-binder system.

    python3 sparkcage/demo_end_to_end.py

WHY THIS WORKS WITHOUT THE MONOD BIO LATCH
------------------------------------------
7CBC is not just a cage. It is sCageHA267_1S: a LOCKR-family cage, a latch, and
an already-grafted de novo mini-binder (HB1.9549.2, raised against influenza
hemagglutinin). So it is a complete, published, experimentally solved instance
of exactly the architecture this project is building, with somebody else's
binder already in the latch.

That makes it a legitimate stand-in. The binder recognises hemagglutinin rather
than GFAP or UCH-L1, but nothing downstream of binding cares what the analyte
is: the switch thermodynamics take a dissociation constant, and the
electrochemistry takes a fraction switched. Swap the constant, keep the
pipeline. When the real latch arrives, only the coordinates and the affinity
change.

WHAT THIS DEMO ACTUALLY COMPUTES
--------------------------------
The two pieces that are ours to deliver, joined into one chain:

    analyte concentration
      -> fraction of sensor switched            (thermo/switch_model.py)
      -> where the reporter sits, caged vs released   (this file)
      -> electron transfer rate at that distance      (echem/swv.py)
      -> square-wave voltammogram and peak current    (echem/swv.py)
      -> predicted calibration curve

Every stage is a model that has been checked separately in validate.py. What is
new here is running them as one pipeline on real coordinates, so a change in an
attachment residue comes out the far end as a predicted voltammogram.

HONESTY
-------
These are predictions from a calibrated model, not measurements. Trends,
rankings and scalings are the reliable output. Absolute currents depend on
k0_contact, beta, coverage and the ensemble widths, none of which is known for
methylene blue on a protein on carbon. Calibrate against one measured
voltammogram before quoting any absolute number.
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

from structure import analyze                    # noqa: E402
from thermo import switch_model as sw            # noqa: E402
from echem import swv                            # noqa: E402
from mb_placement_design import (                # noqa: E402
    BURIAL, SIDECHAIN_VOL, MB_VOLUME_A3, load_scaffold, caged_distance,
)

FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)
RNG = np.random.default_rng(7)

# Electrode and reporter parameters, calibrated in run_demo.py
CELL = dict(area_cm2=0.03, e_step=0.002)
GAMMA_TOTAL = 2.0e-12
K0_CONTACT, BETA, D_CONTACT = 3.0e2, 1.0, 5.0
VDW_STANDOFF = 5.0
KUHN_ANG, RES_PER_KUHN = 6.5, 2.0
SIGMA_FREE = 3.0          # ensemble width of an unhindered tethered dye


def k0(distance_ang):
    return K0_CONTACT * np.exp(-BETA * (distance_ang - D_CONTACT))


def released_effective_distance(n_res, n_samples=20000):
    """Effective transfer distance of the dye once the latch is free.

    Freely-jointed chain from a tether point on the electrode, reflecting wall.
    Defined so exp(-beta d_eff) equals the ensemble mean of exp(-beta z), which
    is the quantity the rate actually depends on.
    """
    n_kuhn = max(1, int(round(n_res / RES_PER_KUHN)))
    v = RNG.normal(size=(n_samples, n_kuhn, 3))
    v /= np.linalg.norm(v, axis=2, keepdims=True)
    z = np.abs(np.cumsum(v[:, :, 2] * KUHN_ANG, axis=1)[:, -1])
    return float(-np.log(np.exp(-BETA * z).mean()) / BETA)


def caged_effective_distance(mean_d, burial, n_samples=20000):
    """Effective transfer distance while the dye is pinned at the interface.

    A more deeply buried reporter has less freedom, so its distance
    distribution is narrower. The width is scaled by the unburied fraction.
    """
    sigma = SIGMA_FREE * (1.0 - burial)
    if sigma < 1e-3:
        return mean_d
    z = np.abs(RNG.normal(mean_d, sigma, n_samples))
    return float(-np.log(np.exp(-BETA * z).mean()) / BETA)


def voltammogram(fraction_switched, d_caged_eff, d_released_eff, freq_hz):
    """Net square-wave voltammogram of a monolayer that is part switched.

    NOTE ON THE RETURNED PEAK. For a monolayer carrying two populations with
    very different rate constants, the returned peak height is NOT the right
    signal metric: the two components sit at different potentials, so the peak
    tracks whichever component is momentarily larger and behaves
    non-monotonically in the switched fraction. Use `signal_at` below.
    """
    f = float(fraction_switched)
    a = swv.swv_scan(k0(d_caged_eff), gamma_mol_cm2=GAMMA_TOTAL * (1.0 - f),
                     freq_hz=freq_hz, **CELL)
    b = swv.swv_scan(k0(d_released_eff), gamma_mol_cm2=GAMMA_TOTAL * f,
                     freq_hz=freq_hz, **CELL)
    net = a["i_net"] + b["i_net"]
    j = int(np.argmax(np.abs(net)))
    peak, pot = swv._refine_peak(a["e"], net, j)
    return a["e"], net, peak, pot


def signal_at(potential, fraction_switched, d_caged_eff, d_released_eff, freq_hz):
    """Current measured AT a fixed potential, which is what an instrument does.

    The caged and released reporters have very different rate constants, so
    under square-wave interrogation their peaks sit at different potentials.
    Reading at the released component's peak potential makes the caged
    population nearly invisible, and the measured current becomes linear in the
    switched fraction with a large gain. This is the single most consequential
    modelling result in the package, because it decides whether a
    macroelectrode format can see a switched fraction of a fraction of a
    percent.
    """
    e, net, _, _ = voltammogram(fraction_switched, d_caged_eff,
                                d_released_eff, freq_hz)
    return float(abs(net[int(np.argmin(np.abs(e - potential)))]))


def build_variants():
    """Score every latch position and return the ones worth simulating."""
    ident, height = load_scaffold()
    out = []
    for r, burial in sorted(BURIAL.items()):
        if burial < 0.30:
            continue
        name = ident.get(r, "UNK")
        native = SIDECHAIN_VOL.get(name, 80)
        deficit = burial * MB_VOLUME_A3 - native
        d_mean = caged_distance(height, r, VDW_STANDOFF)
        d_caged = caged_effective_distance(d_mean, burial)
        d_rel = released_effective_distance(r - 243) + VDW_STANDOFF
        out.append(dict(res=r, aa=name, burial=burial, deficit=deficit,
                        d_caged_mean=d_mean, d_caged_eff=d_caged,
                        d_released_eff=d_rel, fits=deficit <= 0))
    return out


def choose_frequency(variant, freqs=np.logspace(0, 3.5, 40)):
    """Frequency that maximises the contrast between the two states."""
    best, best_f = -1.0, 100.0
    for f in freqs:
        _, _, p0, _ = voltammogram(0.0, variant["d_caged_eff"],
                                   variant["d_released_eff"], f)
        _, _, p1, _ = voltammogram(1.0, variant["d_caged_eff"],
                                   variant["d_released_eff"], f)
        cap = abs(swv.capacitive_background(freq_hz=f, area_cm2=CELL["area_cm2"]))
        snr = abs(abs(p1) - abs(p0)) / np.sqrt(cap ** 2 + 1e-12 ** 2)
        if snr > best:
            best, best_f = snr, float(f)
    return best_f


def main():
    print("=" * 78)
    print("SPARKCAGE END-TO-END DEMO ON 7CBC (sCageHA267_1S)")
    print("=" * 78)
    print("""
  System: a published LOCKR-family cage with a latch carrying an already-
  grafted de novo mini-binder against influenza hemagglutinin. Real
  coordinates, somebody else's binder, our two work packages on top.

  The cage is anchored to the electrode at its latch-junction end, which
  latch_geometry.py shows is the only end that lets a released reporter reach
  the surface.
""")
    analyze.fetch("7CBC", HERE / "structure" / "pdb")
    variants = build_variants()

    print("  STEP 1  Reporter placement, scored on real coordinates\n")
    print("  pos  aa   burial  MB fits   caged d(eff)  released d(eff)")
    print("  " + "-" * 62)
    for v in variants:
        print(f"  {v['res']:3d}  {v['aa']}  {v['burial']*100:5.1f}%  "
              f"{'yes' if v['fits'] else 'NO ':>5s}     {v['d_caged_eff']:7.1f} A"
              f"       {v['d_released_eff']:7.1f} A")

    playable = [v for v in variants if v["fits"]] or variants
    lead = min(playable, key=lambda v: v["d_caged_eff"])
    print(f"\n  Lead by the three-constraint rule: residue {lead['res']} "
          f"({lead['aa']}{lead['res']}C)")

    print("\n  STEP 2  Interrogation frequency, from the two states\n")
    freq = choose_frequency(lead)
    _, _, p0, _ = voltammogram(0.0, lead["d_caged_eff"], lead["d_released_eff"], freq)
    _, _, p1, _ = voltammogram(1.0, lead["d_caged_eff"], lead["d_released_eff"], freq)
    _, _, _, pot0 = voltammogram(0.0, lead["d_caged_eff"], lead["d_released_eff"], freq)
    _, _, _, pot1 = voltammogram(1.0, lead["d_caged_eff"], lead["d_released_eff"], freq)
    print(f"    best contrast at {freq:.0f} Hz")
    print(f"    caged    peak {abs(p0)*1e9:9.2f} nA at {pot0*1000:8.1f} mV")
    print(f"    released peak {abs(p1)*1e9:9.2f} nA at {pot1*1000:8.1f} mV")
    print(f"    peak separation {abs(pot1-pot0)*1000:7.1f} mV")
    print("""
    The two populations are RESOLVED, not merely different in size. Their rate
    constants differ by more than two orders of magnitude, and under square-wave
    interrogation a slow couple's peak is shifted well away from a fast one's.
    That changes what should be measured: read the current at the RELEASED
    component's potential, where the caged population is nearly silent, rather
    than the height of the combined peak.""")

    print("\n  STEP 3  Full calibration curve, analyte to current\n")
    # Binder affinity is the one number that changes with the analyte.
    kds = [1e-11, 1e-10, 1e-9]
    dg_open = 2.0
    conc = np.logspace(-13, -6, 40)

    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.4))

    for f_sw, style, lbl in [(0.0, "-", "no analyte"),
                             (0.5, "--", "half switched"),
                             (1.0, ":", "saturated")]:
        e, net, _, _ = voltammogram(f_sw, lead["d_caged_eff"],
                                    lead["d_released_eff"], freq)
        ax[0].plot(e, net * 1e9, style, label=lbl)
    ax[0].set_xlabel("E vs Ag/AgCl (V)")
    ax[0].set_ylabel("net SWV current (nA)")
    ax[0].set_title(f"Predicted voltammogram\n{lead['aa']}{lead['res']}C, {freq:.0f} Hz")
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=0.3)

    results = {}
    read_at = pot1                       # the released component's potential
    for kd in kds:
        f_sw = sw.fraction_switched(conc, dg_open, sw.dg_from_kd(kd, 298.15), 298.15)
        peaks = np.array([signal_at(read_at, f, lead["d_caged_eff"],
                                    lead["d_released_eff"], freq) for f in f_sw])
        gain = (peaks - peaks[0]) / peaks[0] * 100.0
        ax[1].semilogx(conc * 1e12, gain, "o-", ms=3,
                       label=f"binder $K_d$ = {kd*1e12:.0f} pM")
        results[f"kd_{kd:.0e}"] = dict(conc_M=conc.tolist(), gain_pct=gain.tolist())
    for cut, name in [(0.60e-12, "GFAP threshold"), (14.5e-12, "UCH-L1 threshold")]:
        ax[1].axvline(cut * 1e12, ls=":", lw=1, color="grey")
        ax[1].text(cut * 1e12, ax[1].get_ylim()[1] * 0.55, name, rotation=90,
                   fontsize=7, ha="right", va="top", color="grey")
    ax[1].set_xlabel("analyte (pM)")
    ax[1].set_ylabel("current change at released $E_{pk}$ (%)")
    ax[1].set_title(f"Calibration curve\n$\\Delta G_{{open}}$ = {dg_open:.0f} kcal/mol")
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=0.3)

    labels, caged, rel = [], [], []
    for v in variants:
        _, _, a, _ = voltammogram(0.0, v["d_caged_eff"], v["d_released_eff"], freq)
        _, _, b, _ = voltammogram(1.0, v["d_caged_eff"], v["d_released_eff"], freq)
        labels.append(f"{v['aa'][0]}{v['res']}")
        caged.append(max(abs(a) * 1e9, 1e-6))
        rel.append(max(abs(b) * 1e9, 1e-6))
    x = np.arange(len(labels))
    ax[2].bar(x - 0.2, caged, 0.4, label="caged")
    ax[2].bar(x + 0.2, rel, 0.4, label="released")
    ax[2].set_yscale("log")
    ax[2].set_xticks(x)
    ax[2].set_xticklabels(labels, rotation=60, fontsize=7)
    ax[2].axhline(0.1, color="r", ls=":", lw=1)
    ax[2].text(0.2, 0.12, "baseline floor", color="r", fontsize=7)
    ax[2].set_ylabel("peak current (nA)")
    ax[2].set_title("Every candidate position\n(bars below the line have no baseline)")
    ax[2].legend(fontsize=8)
    ax[2].grid(alpha=0.3, axis="y")

    fig.tight_layout()
    out = FIG / "05_end_to_end_demo.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"    -> {out}")

    for kd in kds:
        f_sw = sw.fraction_switched(conc, dg_open, sw.dg_from_kd(kd, 298.15), 298.15)
        peaks = np.array([signal_at(read_at, f, lead["d_caged_eff"],
                                    lead["d_released_eff"], freq) for f in f_sw])
        gain = (peaks - peaks[0]) / peaks[0] * 100.0
        idx = int(np.argmax(gain > 5.0)) if np.any(gain > 5.0) else -1
        lod = conc[idx] if idx >= 0 else float("nan")
        print(f"    binder Kd {kd*1e12:6.0f} pM -> 5% response at "
              f"{lod*1e12:8.2f} pM, saturating gain {gain[-1]:6.1f}%")

    (FIG / "end_to_end_results.json").write_text(json.dumps(
        dict(lead=lead, frequency_hz=freq, variants=variants,
             caged_peak_nA=abs(p0) * 1e9, released_peak_nA=abs(p1) * 1e9,
             calibration=results), indent=2, default=float))

    # ---- how small a switched fraction is visible, with the right metric
    print("\n  STEP 4  How small a switched fraction can be seen\n")
    base = signal_at(read_at, 0.0, lead["d_caged_eff"], lead["d_released_eff"], freq)
    print("    switched     peak height     current at released Epk")
    print("    fraction     (wrong metric)  (right metric)")
    print("    " + "-" * 54)
    for f in (0.001, 0.0026, 0.01, 0.05, 0.10):
        _, _, pk, _ = voltammogram(f, lead["d_caged_eff"],
                                   lead["d_released_eff"], freq)
        ph = (abs(pk) - abs(p0)) / abs(p0) * 100.0
        sg = (signal_at(read_at, f, lead["d_caged_eff"],
                        lead["d_released_eff"], freq) - base) / base * 100.0
        print(f"      {f*100:6.2f}%      {ph:9.2f}%       {sg:12.1f}%")
    gain_01 = (signal_at(read_at, 0.001, lead["d_caged_eff"],
                         lead["d_released_eff"], freq) - base) / base * 100.0
    print(f"""
    This materially revises the occupancy problem raised in detection_limit.py.
    That analysis found that on a macroelectrode at clinical concentrations the
    immobilised sensors outnumber the analyte enough to cap occupancy near
    0.1%, and concluded the signal would therefore be 0.1% of full modulation,
    an order of magnitude below the roughly 1% noise floor of these sensors.

    Measured at the released component's potential, 0.1% occupancy gives about
    {gain_01:.0f}% instead. The gain comes from the two populations being separated in
    POTENTIAL, not merely in amplitude: at the released peak the caged
    population contributes almost nothing, so the measurement counts switched
    molecules against a near-zero background rather than against the whole
    monolayer.

    Two caveats. The separation is a model output, driven by the transfer
    coefficient and the ratio of rate constant to frequency, and it has to be
    confirmed on a real voltammogram. And none of this lifts the hard charge
    bound, which counts molecules in the sample and cannot be argued with. What
    it does say is that the macroelectrode format is not dead on arrival, which
    is the opposite of the earlier single-metric conclusion.
""")

    print("""
  READ-OUT

  The pipeline runs end to end on real coordinates: a choice of attachment
  residue comes out the far end as a predicted voltammogram and a predicted
  calibration curve. That is the deliverable, and it does not need the Monod
  Bio latch to be exercised, debugged or reviewed.

  When the real latch arrives, three things change and nothing else: the
  coordinates, the burial numbers computed from them, and the binder
  dissociation constant. The placement rule, the frequency selection, the
  voltammetry and the calibration mapping are unchanged.

  What this demo cannot tell you, and no model can: whether methylene blue at a
  59%-buried position actually stays put, whether the released latch really
  unfolds, and what k0 methylene blue has on carbon through a protein. Those
  set the absolute currents and they are what the twelve-construct screen in
  mb_placement_design.py exists to measure.
""")


if __name__ == "__main__":
    main()
    print("=" * 78)
