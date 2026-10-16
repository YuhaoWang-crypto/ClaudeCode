"""Triage: can conformational dynamics show up in pulsed-voltammetry current at all?

Run this BEFORE building routes or fitting.  It compares the slowest conformational relaxation of the
1-D Smoluchowski model (R1) with the half-pulse length, and directly computes how far R1 departs from the
fast-equilibrium route R0 for a plausible effective rate, at the Stokes-Einstein motion scale (s = 1) and
at 1000x slower motion (s = 1e-3).

    python timescale_check.py --z heights.npy --radius 0.35 --freqs 20,100,500 --k-eff 10,100,1000
    python timescale_check.py --synthetic --radius 0.35           # gamma-distributed placeholder heights

--z: 1-D .npy of reporter-centre heights above the monolayer plane (nm) from an equilibrium sampler
(CG Brownian dynamics, MD, or any ensemble).  Never use a mean height: the rate average needs the
full distribution.
"""
import argparse

import numpy as np

import swv_backend as sb

KB = 1.380649e-23


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--z", help=".npy of heights (nm)")
    ap.add_argument("--synthetic", action="store_true")
    ap.add_argument("--radius", type=float, default=0.35, help="reporter radius = contact height z_min (nm)")
    ap.add_argument("--eta", type=float, default=0.89e-3, help="viscosity (Pa s)")
    ap.add_argument("--T", type=float, default=298.15)
    ap.add_argument("--beta", type=float, default=5.0, help="distance decay (1/nm), an assumption")
    ap.add_argument("--freqs", default="20,100,500")
    ap.add_argument("--amp", type=float, default=0.025)
    ap.add_argument("--k-eff", default="10,100,1000", help="effective rates k_c<w> to test (1/s)")
    ap.add_argument("--cells", type=int, default=40)
    a = ap.parse_args()
    if a.z:
        z = np.load(a.z).ravel()
    else:
        z = a.radius + np.random.default_rng(0).gamma(3.0, 1.0, 200000)
        print("NOTE: synthetic heights (placeholder only)")
    D = KB * a.T / (6 * np.pi * a.eta * a.radius * 1e-9) * 1e18          # nm^2/s, free-draining proxy
    kin1 = sb.smoluchowski_kinetics(z, a.cells, D, a.beta, a.radius)
    wmean = float(np.mean(np.exp(-a.beta * (z - a.radius))))
    kin0 = sb.fast_eq_kinetics([1.0], [wmean])
    lam = np.sort(np.linalg.eigvalsh(kin1.S))[::-1]
    t_relax = -1 / lam[1]
    freqs = [float(f) for f in a.freqs.split(",")]
    print(f"heights: n={z.size}, mean {z.mean():.2f} nm, P(z < z_min+0.2) = {np.mean(z < a.radius + 0.2):.4f}")
    print(f"D (Stokes-Einstein, r={a.radius} nm) = {D:.3g} nm^2/s;  <exp(-beta(z-z_min))> = {wmean:.4g}")
    print(f"slowest conformational relaxation (R1, {a.cells} cells): {t_relax:.3g} s")
    for f in freqs:
        print(f"  {f:g} Hz: half-pulse {0.5/f:.3g} s  ->  separation {0.5/f/t_relax:.3g}x")
    wfs = [sb.Waveform(f, a.amp) for f in freqs]
    print("max |R1 - R0| / range of net current:")
    for ke in [float(k) for k in a.k_eff.split(",")]:
        kc = ke / wmean
        r0 = sb.simulate(kin0, wfs, kc, 0.4, -0.26)
        row = []
        for s in (1.0, 1e-3):
            r1 = sb.simulate(kin1.scaled(s), wfs, kc, 0.4, -0.26)
            d = max(np.abs(x["net"] - y["net"]).max() / max(np.ptp(y["net"]), 1e-300) for x, y in zip(r1, r0))
            row.append(d)
        print(f"  k_eff {ke:g}/s (k_c {kc:.3g}/s): s=1 -> {100*row[0]:.3g} %,  s=1e-3 -> {100*row[1]:.3g} %")
    print("\nReading: if the separation is >> 1e3 and |R1 - R0| at s = 1 is below the data's noise floor, R1/R2\n"
          "cannot outperform R0 on these frequencies; say so up front and spend effort on the chemistry instead.\n"
          "A large difference only at s = 1e-3 means dynamics matter only if motion is ~1000x slower than\n"
          "Stokes-Einstein; such fits will push s to its lower bound - check identifiability before believing them.")


if __name__ == "__main__":
    main()
