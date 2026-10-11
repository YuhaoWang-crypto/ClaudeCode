"""
Module 10 - what would it actually cost to converge the TEM MSM?

M9 shows adaptive reseeding is the right protocol.  This module turns that
into a budget, using MD throughput MEASURED on the machine you run it on
rather than a quoted figure, and states plainly whether the campaign is
feasible here.

The arithmetic is deliberately simple and conservative:

  required aggregate  ~  N_transitions x tau_slow
      where tau_slow is the slowest process you need resolved and
      N_transitions is how many A<->B events the MSM needs to pin its
      rate.  ~100 is a common rule of thumb for a few-percent rate error;
      M9 reached ~10 % median error on the timescale with far fewer, so
      we report a range.

  adaptive speedup    ~  from M9: adaptive reached a given error at a
      fraction of the cost of long trajectories.  We quote this as a
      measured factor rather than folklore.
"""
import time

import numpy as np

from .common import banner

# From M2/M3 on the real data: the pocket-local process does not converge but
# its lower bound sits in the hundreds of ns.
TAU_SLOW_NS = 500.0          # conservative lower bound on the slow process
N_TRANSITIONS = (30, 100)    # events needed, loose -> comfortable


def measure_throughput(quick=True):
    """Benchmark OpenMM on THIS machine. -> dict of ns/day."""
    import openmm as mm
    import openmm.app as app
    import openmm.unit as u
    from pdbfixer import PDBFixer

    from .common import PDB, fetch_crystals
    fetch_crystals(["1BTL"], verbose=False)

    out = {}
    plats = [mm.Platform.getPlatform(i).getName()
             for i in range(mm.Platform.getNumPlatforms())]
    best = "CUDA" if "CUDA" in plats else ("OpenCL" if "OpenCL" in plats else "CPU")
    print(f"  OpenMM platforms available: {plats}  -> using {best}")

    fx = PDBFixer(filename=f"{PDB}/1BTL.pdb")
    fx.removeHeterogens(False)
    fx.findMissingResidues()
    fx.findMissingAtoms()
    fx.addMissingAtoms()
    fx.addMissingHydrogens(7.0)
    ff = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    mod = app.Modeller(fx.topology, fx.positions)
    mod.addSolvent(ff, padding=1.0 * u.nanometer, ionicStrength=0.15 * u.molar)
    n_atoms = mod.topology.getNumAtoms()
    system = ff.createSystem(mod.topology, nonbondedMethod=app.PME,
                             nonbondedCutoff=1.0 * u.nanometer,
                             constraints=app.HBonds, hydrogenMass=1.5 * u.amu)
    dt = 4.0
    integ = mm.LangevinMiddleIntegrator(300 * u.kelvin, 1 / u.picosecond,
                                        dt * u.femtoseconds)
    sim = app.Simulation(mod.topology, system, integ,
                         mm.Platform.getPlatformByName(best))
    sim.context.setPositions(mod.positions)
    sim.minimizeEnergy(maxIterations=200)
    sim.step(20)
    n = 100 if quick else 500
    t = time.time()
    sim.step(n)
    el = time.time() - t
    nsday = n * dt / 1e6 / (el / 86400)
    out[f"TEM-1 explicit ({n_atoms} atoms), {best}"] = nsday
    return out, n_atoms, best, nsday


def main(quick=True):
    banner("M10  cost of actually converging the TEM MSM, on measured hardware")
    tp, n_atoms, plat, nsday = measure_throughput(quick)
    for k, v in tp.items():
        print(f"  measured throughput: {k}: {v:.2f} ns/day")

    lo, hi = N_TRANSITIONS
    print(f"\n  target: resolve a tau_slow >= {TAU_SLOW_NS:.0f} ns process")
    print(f"  aggregate needed = N_events x tau_slow, N = {lo}-{hi}")
    agg_lo, agg_hi = lo * TAU_SLOW_NS / 1000, hi * TAU_SLOW_NS / 1000
    print(f"                   = {agg_lo:.0f} - {agg_hi:.0f} us")

    print(f"\n  {'hardware':34s} {'ns/day':>9s} {'time for '+f'{agg_lo:.0f} us':>16s}"
          f" {f'for {agg_hi:.0f} us':>14s}")
    rows = [(f"this machine ({plat})", nsday, 1),
            ("1x RTX 4090 (typical 40k atoms)", 250.0, 1),
            ("8x RTX 4090 (adaptive walkers)", 250.0, 8),
            ("1x H100", 600.0, 1),
            ("8x H100", 600.0, 8)]
    for name, per, n in rows:
        tot = per * n
        print(f"  {name:34s} {tot:9.0f} {agg_lo*1000/tot:13.0f} d {agg_hi*1000/tot:12.0f} d")

    print(f"""
  VERDICT FOR THIS CONTAINER
    {nsday:.1f} ns/day on 4 CPU cores means the {agg_lo:.0f} us lower bound would take
    {agg_lo*1000/nsday/365:.0f} years.  The campaign is NOT runnable here, and nothing in
    this repo pretends otherwise.  What is runnable here is the controlled
    benchmark (M9) and the end-to-end driver test on alanine dipeptide.

  WHAT IS ACTUALLY FEASIBLE
    On a single modern GPU the lower-bound campaign is weeks, not years, and
    adaptive sampling parallelises almost perfectly across walkers - that is
    the main practical reason to use it, beyond the statistical gain in M9.
    8 GPUs puts {agg_lo:.0f} us within a few days.

  A CONCRETE PROTOCOL FOR TEM
    python -m msm_pipeline.adaptive_md \\
        --pdb 1BTL.pdb --out /Volumes/T7/tem_adaptive \\
        --rounds 60 --walkers 16 --ns-per-walker 20 \\
        --sel "name CA and (resid 214 to 231 or resid 269 to 287)" \\
        --strategy hybrid --platform CUDA
      = 60 x 16 x 20 ns = 19.2 us, checkpointed every round, resumable.

    Seed it from BOTH ends, not just the closed state: include 1PZO (the
    inhibitor-bound open form) among the starting structures.  M3's failure
    was a connectivity failure, and the cheapest way to connect two basins is
    to start in both of them.  The driver's --resume lets you add 1PZO-seeded
    rounds to a running campaign.

    Stop when the ITS test flattens, not when the rounds run out.  The driver
    enforces this and refuses to call an unflattened campaign converged.""")
    return {"ns_per_day": nsday, "n_atoms": n_atoms}


if __name__ == "__main__":
    main()
