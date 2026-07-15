#!/usr/bin/env python3
"""
Toy NEMD desalination -- a dependency-light (numpy only) kinetic model that
stands in for the real LAMMPS/DeePMD run so the *analysis pipeline* and the
*physical picture* can be validated end-to-end in seconds on a laptop.

What it is
----------
A 1D-along-z overdamped Langevin ("Brownian dynamics") model of N water
"particles" and M ion pairs in a feed reservoir, pushed toward a membrane at
z = 0 by a constant pressure-like force. The membrane is a barrier at z = 0
with a pore modelled as a position-dependent transmission:

    * water sees a low, crossable free-energy barrier      -> high flux
    * ions see a high barrier (dehydration penalty)        -> mostly rejected

This reproduces the qualitative RO result -- water permeates, salt is
rejected -- WITHOUT any DFT or force field. It is deliberately NOT quantitative;
its job is to (a) demonstrate the full feed->membrane->permeate workflow and
(b) prove that analyze_flux_rejection.py recovers the correct numbers.

What it validates
-----------------
1. Geometry/counting: the analyzer's net-crossing count matches an independent
   ground-truth counter tracked inside the integrator.
2. Physics sign: higher pressure -> higher water flux; higher ion barrier ->
   higher rejection (monotonic trends asserted).
3. Unit conversion: molecular flux -> experimental permeance runs without error.

Run
---
    python toy_nemd.py --plot ../figures/toy_nemd.png
"""
from __future__ import annotations
import argparse
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "04_nemd"))
import analyze_flux_rejection as afr  # noqa: E402


def transmission_prob(pressure, barrier, k=1.6):
    """Logistic pore-transmission probability.

    Applied pressure lowers the effective free-energy barrier a particle must
    overcome to cross the pore; the crossing probability is a sigmoid of
    (k*pressure - barrier). Water has a small `barrier`, ions a large one
    (hydration-shell stripping penalty), so water passes and salt is rejected.
    This is the standard kinetic pore-gate picture used in reduced desalination
    models -- it makes the flux/rejection trends monotonic by construction,
    which is exactly what we want when validating the *analysis* code.
    """
    return 1.0 / (1.0 + np.exp(-(k * pressure - barrier)))


def run_bd(n_water=200, n_ion_pairs=20, pressure=1.5, water_barrier=0.5,
           ion_barrier=6.0, kT=1.0, gamma=1.0, dt=0.01, n_steps=4000,
           feed_z0=-14.0, z_absorb=8.0, record_every=5, seed=1):
    """Drift-diffusion of particles toward a gated membrane at z=0.

    Feed particles drift toward the membrane under a pressure-proportional
    force plus thermal noise. When a particle first reaches z=0 it attempts to
    cross: it transmits with probability `transmission_prob(pressure, barrier)`
    (then continues to the permeate side and is absorbed at z_absorb), or it is
    reflected back into the feed. Absorbing the permeate side means every
    transmitted particle contributes exactly one clean forward crossing, giving
    an unambiguous ground truth to test the analyzer against.

    Returns dict with z_water, z_ions trajectories, times, and the online
    ground-truth net water crossings.
    """
    rng = np.random.default_rng(seed)
    n_ion = 2 * n_ion_pairs
    p_water = transmission_prob(pressure, water_barrier)
    p_ion = transmission_prob(pressure, ion_barrier)

    def _init(n):
        return feed_z0 + rng.uniform(0.0, 8.0, size=n)

    zw, zi = _init(n_water), _init(n_ion)
    active_w = np.ones(n_water, dtype=bool)   # False once absorbed on permeate
    active_i = np.ones(n_ion, dtype=bool)
    gated_w = np.zeros(n_water, dtype=bool)   # True once the pore was attempted
    gated_i = np.zeros(n_ion, dtype=bool)
    trans_w = np.zeros(n_water, dtype=bool)   # True if the attempt succeeded
    trans_i = np.zeros(n_ion, dtype=bool)

    noise_amp = np.sqrt(2.0 * kT * gamma * dt)
    frames_w, frames_i, times = [], [], []
    gt_net = 0
    prev_side_w = (zw > 0.0).astype(int)

    def _step(z, active, gated, transmitted, p_transmit):
        # drift toward +z under pressure, plus diffusion (only active particles)
        dz = (pressure / gamma) * dt + noise_amp / gamma * rng.normal(size=z.size)
        z_new = np.where(active, z + dz, z)

        # feed particles reaching the membrane for the first time attempt the
        # pore exactly once (single attempt -> a rejected ion is not given
        # infinite retries, which would drive rejection to zero over long runs)
        reached = active & (~gated) & (z < 0.0) & (z_new >= 0.0)
        u = rng.random(size=z.size)
        gated |= reached
        transmitted |= reached & (u < p_transmit)

        # anything not (yet) transmitted is held on the feed side: rejected
        # ions and not-yet-arrived particles can never leak into the permeate.
        held = active & (~transmitted) & (z_new >= 0.0)
        z_new[held] = -0.1
        # far feed wall (the pushing piston)
        z_new = np.where(z_new < feed_z0 - 2.0, feed_z0 - 2.0, z_new)
        # absorb transmitted particles once they are well into the permeate
        absorbed = active & transmitted & (z_new > z_absorb)
        z_new[absorbed] = z_absorb + 2.0
        active[absorbed] = False
        return z_new

    for step in range(n_steps):
        zw = _step(zw, active_w, gated_w, trans_w, p_water)
        zi = _step(zi, active_i, gated_i, trans_i, p_ion)

        side_w = (zw > 0.0).astype(int)
        gt_net += int(((prev_side_w == 0) & (side_w == 1)).sum())
        gt_net -= int(((prev_side_w == 1) & (side_w == 0)).sum())
        prev_side_w = side_w

        if step % record_every == 0:
            frames_w.append(zw.copy())
            frames_i.append(zi.copy())
            times.append(step * dt)

    return {
        "z_water": np.array(frames_w),
        "z_ions": np.array(frames_i),
        "times_ns": np.array(times),
        "gt_net_water": gt_net,
        "n_ion_feed_initial": n_ion,
        "p_water": p_water,
        "p_ion": p_ion,
    }


def analyze(traj, area_A2=150.0, pressure_MPa=100.0):
    return afr.summarize(
        traj["z_water"], traj["z_ions"], z_membrane=0.0,
        times_ns=traj["times_ns"], area_A2=area_A2,
        n_ions_feed_initial=traj["n_ion_feed_initial"],
        pressure_MPa=pressure_MPa,
    )


def validate():
    """Run the assertions that make this a real test, not just a demo."""
    print("== Toy NEMD validation ==")
    ok = True

    # --- Test 1: analyzer counting matches the online ground truth ----------
    traj = run_bd(seed=1)
    summ = analyze(traj)
    net_analyzer = summ["net_water_permeated"]
    gt = traj["gt_net_water"]
    # the analyzer subsamples frames (record_every), so it can miss crossings
    # that happen and reverse between samples; require agreement within 10%.
    rel = abs(net_analyzer - gt) / max(1, gt)
    print(f"[1] crossings  analyzer={net_analyzer}  ground_truth={gt}  "
          f"rel_err={rel:.2%}  -> {'PASS' if rel < 0.10 else 'FAIL'}")
    ok &= rel < 0.10

    # --- Test 2: pressure monotonically increases water throughput ----------
    # Asserted on the count-based net_permeated (robust to the finite-feed
    # saturation of this toy). The analyzer's polyfit `water_flux_per_ns` is the
    # right estimator for a steady-state piston run where the curve is linear.
    perms = []
    for p in (0.5, 1.5, 3.0):
        s = analyze(run_bd(pressure=p, seed=2))
        perms.append(s["net_water_permeated"])
    mono = perms[0] < perms[1] < perms[2]
    print(f"[2] water permeated vs pressure {perms}  "
          f"-> {'PASS' if mono else 'FAIL'}")
    ok &= mono

    # --- Test 3: higher ion barrier increases rejection ---------------------
    r_lo = analyze(run_bd(ion_barrier=2.0, seed=3))["salt_rejection"]
    r_hi = analyze(run_bd(ion_barrier=8.0, seed=3))["salt_rejection"]
    trend = r_hi >= r_lo
    print(f"[3] rejection  low_barrier={r_lo:.2f}  high_barrier={r_hi:.2f}  "
          f"-> {'PASS' if trend else 'FAIL'}")
    ok &= trend

    # --- Test 4: rejection is selective (water passes, ions mostly don't) ---
    base = analyze(run_bd(seed=4))
    selective = base["net_water_permeated"] > 0 and base["salt_rejection"] > 0.5
    print(f"[4] selectivity  water_permeated={base['net_water_permeated']}  "
          f"rejection={base['salt_rejection']:.2f}  "
          f"-> {'PASS' if selective else 'FAIL'}")
    ok &= selective

    # --- Test 5: permeance unit conversion produces a finite number ---------
    finite = np.isfinite(base["water_permeance_L_cm2_day_MPa"])
    print(f"[5] permeance = {base['water_permeance_L_cm2_day_MPa']} "
          f"L/(cm^2.day.MPa)  -> {'PASS' if finite else 'FAIL'}")
    ok &= bool(finite)

    print(f"\nRESULT: {'ALL PASS' if ok else 'FAILURES PRESENT'}")
    return ok, base


def make_plot(path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pressures = [0.5, 1.5, 3.0]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))

    # panel 1: cumulative water permeated vs time for 3 pressures
    for p in pressures:
        tr = run_bd(pressure=p, seed=10)
        s = analyze(tr)
        axes[0].plot(tr["times_ns"], s["_cumulative_water"],
                     label=f"P={p}  (p_pass={tr['p_water']:.2f}, "
                           f"total={s['net_water_permeated']})")
    axes[0].set_xlabel("time (ns, toy units)")
    axes[0].set_ylabel("net water permeated")
    axes[0].set_title("Water throughput increases with pressure")
    axes[0].legend(fontsize=8)

    # panel 2: rejection vs ion barrier
    barriers = np.linspace(1.0, 9.0, 9)
    rej = [analyze(run_bd(ion_barrier=b, seed=11))["salt_rejection"]
           for b in barriers]
    axes[1].plot(barriers, rej, "o-", color="crimson")
    axes[1].set_xlabel("ion barrier (dehydration penalty, toy units)")
    axes[1].set_ylabel("salt rejection R")
    axes[1].set_ylim(0, 1.05)
    axes[1].set_title("Rejection increases with ion barrier")

    # panel 3: final z-distribution of water vs ions (selectivity picture)
    tr = run_bd(seed=12)
    axes[2].hist(tr["z_water"][-1], bins=30, alpha=0.6, label="water",
                 color="steelblue")
    axes[2].hist(tr["z_ions"][-1], bins=30, alpha=0.6, label="ions",
                 color="darkorange")
    axes[2].axvline(0.0, color="k", ls="--", lw=1, label="membrane")
    axes[2].set_xlabel("z (feed < 0 < permeate)")
    axes[2].set_ylabel("count")
    axes[2].set_title("Water crosses, ions stay in feed")
    axes[2].legend(fontsize=8)

    fig.suptitle("Toy NEMD desalination -- validates the analysis pipeline "
                 "(NOT a quantitative prediction)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fig.savefig(path, dpi=130)
    print(f"Wrote figure: {path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--plot", metavar="PATH", default=None,
                    help="Write a 3-panel validation figure to PATH")
    ap.add_argument("--no-validate", action="store_true")
    args = ap.parse_args()

    ok = True
    if not args.no_validate:
        ok, _ = validate()
    if args.plot:
        make_plot(args.plot)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
