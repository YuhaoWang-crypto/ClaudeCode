"""Triage: can conformational dynamics show up in pulsed-voltammetry current at all?

Run this BEFORE building routes or fitting.  For each frequency and plausible effective rate it reports
how far the 1-D Smoluchowski route (R1) departs from the fast-equilibrium route (R0) at the given motion
scale, at 100x and 1000x slower motion, and with the chain frozen (static disorder) - all as a fraction of
the LARGEST net-current range over the tested frequencies (so a vanishing current never inflates it).

    # ferrocene (one electron, oxidation scan), heights from cg_md.py, D calibrated to simulated tau
    python timescale_check.py --z heights.npy --radius 0.30 --chem E1 --alpha 0.5 --E0 0.20 \
        --scan=-0.05,0.45 --freqs 50,500,5000 --k-eff 100,1000,10000 --tau-conf 3e-8
    # methylene blue (apparent concerted 2e-, reduction scan)
    python timescale_check.py --z heights.npy --radius 0.35 --chem C0 --E0=-0.26 --scan=0,-0.45

--z: 1-D .npy of reporter-centre heights above the monolayer (nm) from an equilibrium sampler. Never use
a mean height: the rate average needs the full distribution.
--tau-conf: correlation time of the reporter height measured in your trajectories (s). If given, D is
rescaled so R1's slowest relaxation equals it (Stokes-Einstein of the bare reporter ignores chain friction).
Rule of thumb seen in practice: |R1 - R0| / max signal ~ 5 * k_eff * tau_conf, saturating at the
frozen-chain value.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_backend as cb  # noqa: E402
import swv_backend as sb  # noqa: E402

KB = 1.380649e-23


def net_curves(kin, chem, p, wfs):
    return [o["net"] for o in cb.simulate(kin, chem, p, wfs)]


def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    ap.add_argument("--z", help=".npy of heights (nm)")
    ap.add_argument("--synthetic", action="store_true", help="gamma-distributed placeholder heights")
    ap.add_argument("--radius", type=float, default=0.35, help="reporter radius = contact height z_min (nm)")
    ap.add_argument("--eta", type=float, default=0.89e-3, help="viscosity (Pa s)")
    ap.add_argument("--T", type=float, default=298.15)
    ap.add_argument("--beta", type=float, default=5.0, help="distance decay (1/nm), an assumption")
    ap.add_argument("--chem", default="C0", choices=["E1", "C0"], help="E1 one-electron (ferrocene); C0 apparent 2e- (MB)")
    ap.add_argument("--alpha", type=float, default=0.5)
    ap.add_argument("--E0", type=float, default=-0.26, help="formal potential (V), same reference as --scan")
    ap.add_argument("--scan", default="0,-0.45", help="start,end potentials (V); end < start = reduction scan. Write --scan=-0.05,0.45 when it starts with a minus sign")
    ap.add_argument("--freqs", default="20,100,500")
    ap.add_argument("--amp", type=float, default=0.025)
    ap.add_argument("--k-eff", default="10,100,1000", help="effective rates k<w> to test (1/s)")
    ap.add_argument("--tau-conf", type=float, default=None, help="measured height correlation time (s)")
    ap.add_argument("--cells", type=int, default=40)
    a = ap.parse_args()
    if a.z:
        z = np.load(a.z).ravel()
    elif a.synthetic:
        z = a.radius + np.random.default_rng(0).gamma(3.0, 1.0, 200000)
        print("NOTE: synthetic heights (placeholder only)")
    else:
        ap.error("give --z heights.npy or --synthetic")
    e0, e1 = (float(x) for x in a.scan.split(","))
    D = KB * a.T / (6 * np.pi * a.eta * a.radius * 1e-9) * 1e18          # nm^2/s, bare-reporter proxy
    kin1 = sb.smoluchowski_kinetics(z, a.cells, D, a.beta, a.radius)
    lam = np.sort(np.linalg.eigvalsh(kin1.S))[::-1]
    t_relax = -1 / lam[1]
    if a.tau_conf:
        kin1 = kin1.scaled(t_relax / a.tau_conf)
        print(f"D rescaled by {t_relax / a.tau_conf:.3g} so the R1 relaxation equals the measured {a.tau_conf:.3g} s")
        t_relax = a.tau_conf
    wmean = float(np.mean(np.exp(-a.beta * (z - a.radius))))
    kin0 = sb.fast_eq_kinetics([1.0], [wmean])
    static = sb.Kinetics(kin1.pi, kin1.w, None)
    chem = cb.Chemistry(a.chem)
    freqs = [float(f) for f in a.freqs.split(",")]
    print(f"heights: n={z.size}, mean {z.mean():.2f} nm, P(z < z_min+0.2) = {np.mean(z < a.radius + 0.2):.4f}")
    print(f"<exp(-beta(z-z_min))> = {wmean:.4g} (beta {a.beta} /nm); slowest conformational relaxation {t_relax:.3g} s")
    for f in freqs:
        print(f"  {f:g} Hz: half-pulse {0.5/f:.3g} s -> separation {0.5/f/t_relax:.3g}x")
    wfs = [sb.Waveform(f, a.amp, e_start=e0, e_end=e1) for f in freqs]
    print(f"\n|R - R0| as % of the largest net-current range over {a.freqs} Hz ({a.chem}, alpha {a.alpha}, E0 {a.E0} V):")
    print("  k_eff(1/s)   s=1        s=1e-2     s=1e-3     frozen chain   curves with <0.1% of max signal")
    for ke in [float(k) for k in a.k_eff.split(",")]:
        p = dict(k=ke / wmean, alpha=a.alpha, E0=a.E0)
        r0 = net_curves(kin0, chem, p, wfs)
        big = max(np.ptp(x) for x in r0)
        weak = [f"{f:g}Hz" for f, x in zip(freqs, r0) if np.ptp(x) < 1e-3 * big]
        row = []
        for kin in (kin1, kin1.scaled(1e-2), kin1.scaled(1e-3), static):
            r = net_curves(kin, chem, p, wfs)
            row.append(100 * max(np.abs(x - y).max() for x, y in zip(r, r0)) / big)
        print(f"  {ke:<10g}  {row[0]:<9.3g}  {row[1]:<9.3g}  {row[2]:<9.3g}  {row[3]:<13.3g}  {', '.join(weak) or '-'}")
    print("\nReading: compare the s = 1 column with the data's noise floor (digitization or replicate scatter).\n"
          "If it is far below, R1/R2 cannot beat R0 at these frequencies: say so first and put effort into\n"
          "chemistry, rate heterogeneity and identifiability.  A large 'frozen chain' value with a tiny s = 1 value\n"
          "means heights are heterogeneous but exchange washes it out.  Effects that appear only at s <= 1e-2 need\n"
          "motion 100-1000x slower than modelled (e.g. sticky reporter-surface states the CG model cannot produce).")


if __name__ == "__main__":
    main()
