#!/usr/bin/env python3
"""
Production adaptive-sampling driver for real MD (OpenMM).

Round-based, checkpointed, resumable.  Point it at a PDB and it will:

    round r:  run N short walkers in parallel from stored configurations
              append frames + coordinates to the dataset on disk
              featurise -> TICA -> cluster -> MSM -> ITS convergence check
              pick N new seed configurations (see msm_pipeline/adaptive.py)
    ... until the ITS test passes or the round budget is spent.

Each round writes `round_XX/` with trajectories and a `state.npz`, so you can
kill it and resume with --resume.  That matters: these campaigns run for days.

  python -m msm_pipeline.adaptive_md \
      --pdb /Volumes/T7/tem/1BTL.pdb \
      --out /Volumes/T7/tem/adaptive \
      --rounds 40 --walkers 16 --ns-per-walker 20 \
      --sel "name CA and resid 214 to 231 290 to 300" \
      --strategy hybrid --platform CUDA

IMPORTANT - read before launching a long campaign:

  * Adaptive data is NOT Boltzmann-weighted.  You oversample the barrier on
    purpose.  The reversible MLE used here is still asymptotically correct,
    because transition counts are conditional on the starting state, but the
    raw frame histogram is NOT the equilibrium distribution.  Never read
    populations off the frame counts; read them off msm.stationary_distribution.

  * The ITS test is the stopping criterion, not the round budget.  A campaign
    that exhausts its rounds without the ITS flattening has not converged,
    and the driver says so rather than quietly reporting numbers.

  * Walkers are independent processes, so --walkers should be about the number
    of GPUs (CUDA) or physical cores (CPU) you can actually feed.
"""
import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np


# --------------------------------------------------------------- MD worker --
def _run_walker(args):
    """One short walker.  Runs in its own process; returns (coords, path)."""
    (sysxml, top_pdb, positions, seed, steps, report_every, out_npz,
     platform, temperature_K, timestep_fs) = args
    import openmm as mm
    import openmm.app as app
    import openmm.unit as u

    with open(sysxml) as fh:
        system = mm.XmlSerializer.deserialize(fh.read())
    pdb = app.PDBFile(top_pdb)
    integ = mm.LangevinMiddleIntegrator(temperature_K * u.kelvin,
                                        1.0 / u.picosecond,
                                        timestep_fs * u.femtoseconds)
    integ.setRandomNumberSeed(int(seed))
    plat = mm.Platform.getPlatformByName(platform)
    props = {"Threads": "1"} if platform == "CPU" else {}
    sim = app.Simulation(pdb.topology, system, integ, plat, props)
    sim.context.setPositions(positions * u.nanometer)
    sim.context.setVelocitiesToTemperature(temperature_K * u.kelvin, int(seed))

    frames = []
    n_rep = max(1, steps // report_every)
    for _ in range(n_rep):
        sim.step(report_every)
        st = sim.context.getState(getPositions=True)
        frames.append(st.getPositions(asNumpy=True).value_in_unit(u.nanometer))
    arr = np.asarray(frames, dtype=np.float32)
    np.savez_compressed(out_npz, xyz=arr)
    return out_npz, arr.shape


# ------------------------------------------------------------ system setup --
def build_system(pdb_path, work, implicit=False, padding=1.0, temperature=300.0,
                 keep_het=False):
    """Fix, solvate, minimise.  -> (system.xml, topology.pdb, positions)."""
    import openmm as mm
    import openmm.app as app
    import openmm.unit as u
    from pdbfixer import PDBFixer

    sysxml = os.path.join(work, "system.xml")
    toppdb = os.path.join(work, "topology.pdb")
    eqnpy = os.path.join(work, "equilibrated.npy")
    if all(os.path.exists(p) for p in (sysxml, toppdb, eqnpy)):
        print("  reusing prepared system")
        return sysxml, toppdb, np.load(eqnpy)

    print(f"  preparing {pdb_path}")
    fx = PDBFixer(filename=pdb_path)
    if not keep_het:
        # strips ligands and crystallographic waters - but ALSO strips ACE/NME
        # caps, which are HETATM records.  Use --keep-het for capped peptides.
        fx.removeHeterogens(False)
    fx.findMissingResidues()
    fx.findMissingAtoms()
    fx.addMissingAtoms()
    fx.addMissingHydrogens(7.0)

    if implicit:
        ff = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")
        top, pos = fx.topology, fx.positions
        system = ff.createSystem(top, nonbondedMethod=app.CutoffNonPeriodic,
                                 nonbondedCutoff=2.0 * u.nanometer,
                                 constraints=app.HBonds, hydrogenMass=1.5 * u.amu)
    else:
        ff = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
        mod = app.Modeller(fx.topology, fx.positions)
        mod.addSolvent(ff, padding=padding * u.nanometer, ionicStrength=0.15 * u.molar)
        top, pos = mod.topology, mod.positions
        system = ff.createSystem(top, nonbondedMethod=app.PME,
                                 nonbondedCutoff=1.0 * u.nanometer,
                                 constraints=app.HBonds, hydrogenMass=1.5 * u.amu)
        system.addForce(mm.MonteCarloBarostat(1 * u.bar, temperature * u.kelvin, 25))
    print(f"  system: {top.getNumAtoms()} atoms")

    integ = mm.LangevinMiddleIntegrator(temperature * u.kelvin, 1 / u.picosecond,
                                        4 * u.femtoseconds)
    sim = app.Simulation(top, system, integ)
    sim.context.setPositions(pos)
    print("  minimising + equilibrating 100 ps ...")
    sim.minimizeEnergy(maxIterations=2000)
    sim.context.setVelocitiesToTemperature(temperature * u.kelvin)
    sim.step(25000)
    st = sim.context.getState(getPositions=True)
    eq = st.getPositions(asNumpy=True).value_in_unit(u.nanometer)

    with open(sysxml, "w") as fh:
        fh.write(mm.XmlSerializer.serialize(system))
    with open(toppdb, "w") as fh:
        app.PDBFile.writeFile(top, st.getPositions(), fh)
    np.save(eqnpy, eq)
    return sysxml, toppdb, eq


# ------------------------------------------------------------- the MSM part --
def analyse(all_xyz, sel, tica_lag, dim, n_micro, lags, seed=0):
    """-> (msm, active, dtraj_pooled, its_series) on the selected atoms."""
    from deeptime.clustering import KMeans
    from deeptime.decomposition import TICA
    from deeptime.markov import TransitionCountEstimator
    from deeptime.markov.msm import MaximumLikelihoodMSM

    subs = [x[:, sel] * 10.0 for x in all_xyz]          # nm -> angstrom
    n = subs[0].shape[1]
    I, J = np.triu_indices(n, 3)
    feats = [np.linalg.norm(s[:, I] - s[:, J], axis=-1).astype(np.float32) for s in subs]

    usable_feats = [f for f in feats if len(f) > tica_lag]
    if not usable_feats:
        raise RuntimeError(
            f"no trajectory longer than the TICA lag ({tica_lag} frames); "
            f"longest is {max(len(f) for f in feats)}.  Increase "
            f"--ns-per-walker, decrease --frame-ps, or lower --tica-lag.")
    tica = TICA(lagtime=tica_lag, dim=dim, scaling="kinetic_map")
    for f in usable_feats:
        tica.partial_fit(f)
    model = tica.fetch_model()
    Y = [model.transform(f) for f in feats]

    km = (KMeans(n_clusters=min(n_micro, max(2, sum(len(y) for y in Y) // 20)),
                 max_iter=200, init_strategy="kmeans++", fixed_seed=seed,
                 n_jobs=max(1, (os.cpu_count() or 2) // 2))
          .fit(np.concatenate([y[::2] for y in Y])).fetch_model())
    dtrajs = [km.transform(y) for y in Y]

    its = []
    msm = active = None
    for lag in lags:
        usable = [d for d in dtrajs if len(d) > lag]
        if not usable:
            its.append(np.nan)
            continue
        try:
            cm = (TransitionCountEstimator(lag, "sliding").fit(usable)
                  .fetch_model().submodel_largest(directed=True))
            m = MaximumLikelihoodMSM(reversible=True).fit(cm).fetch_model()
            # deeptime returns timescales already in units of frames
            its.append(float(m.timescales(1)[0]))
            msm, active = m, cm.state_symbols
        except Exception:
            its.append(np.nan)
    return msm, active, np.concatenate(dtrajs), its, km, Y


def main(argv=None):
    from .adaptive import STRATEGIES, its_converged, pick_configurations, \
        select_seed_states

    ap = argparse.ArgumentParser(
        description="Adaptive-sampling MD driver (OpenMM + deeptime)",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--rounds", type=int, default=30)
    ap.add_argument("--walkers", type=int, default=8)
    ap.add_argument("--ns-per-walker", type=float, default=10.0)
    ap.add_argument("--frame-ps", type=float, default=10.0)
    ap.add_argument("--sel", default="name CA",
                    help="MDTraj selection defining the feature space")
    ap.add_argument("--strategy", default="hybrid", choices=STRATEGIES)
    ap.add_argument("--platform", default="CPU", choices=["CPU", "CUDA", "OpenCL", "HIP"])
    ap.add_argument("--implicit", action="store_true")
    ap.add_argument("--keep-het", action="store_true",
                    help="do not strip HETATM records; needed for ACE/NME-capped "
                         "peptides, whose caps are HETATMs")
    ap.add_argument("--temperature", type=float, default=300.0)
    ap.add_argument("--timestep-fs", type=float, default=4.0)
    ap.add_argument("--tica-lag", type=int, default=5, help="frames")
    ap.add_argument("--dim", type=int, default=3)
    ap.add_argument("--n-micro", type=int, default=100)
    ap.add_argument("--lags", type=int, nargs="+", default=[1, 2, 5, 10, 20])
    ap.add_argument("--resume", action="store_true")
    a = ap.parse_args(argv)

    os.makedirs(a.out, exist_ok=True)
    t_start = time.time()
    sysxml, toppdb, eq = build_system(a.pdb, a.out, a.implicit,
                                      temperature=a.temperature,
                                      keep_het=a.keep_het)

    import mdtraj as md
    sel = md.load(toppdb).topology.select(a.sel)
    if len(sel) < 4:
        sys.exit(f"selection '{a.sel}' matched {len(sel)} atoms; need >= 4")
    print(f"  feature selection: {len(sel)} atoms")

    steps = int(a.ns_per_walker * 1e6 / a.timestep_fs)   # ns -> fs -> steps
    report = max(1, int(a.frame_ps * 1000 / a.timestep_fs))
    print(f"  {a.rounds} rounds x {a.walkers} walkers x {a.ns_per_walker} ns"
          f"  =  {a.rounds*a.walkers*a.ns_per_walker/1000:.2f} us aggregate")

    all_xyz, start_round = [], 0
    if a.resume:
        for r in range(a.rounds):
            d = os.path.join(a.out, f"round_{r:03d}")
            fs = sorted(f for f in os.listdir(d) if f.endswith(".npz")) \
                if os.path.isdir(d) else []
            if not fs:
                break
            for f in fs:
                all_xyz.append(np.load(os.path.join(d, f))["xyz"])
            start_round = r + 1
        print(f"  resumed: {start_round} rounds, {len(all_xyz)} walkers on disk")

    seeds_xyz = [eq] * a.walkers
    rng = np.random.default_rng(0)
    history = []

    for r in range(start_round, a.rounds):
        rd = os.path.join(a.out, f"round_{r:03d}")
        os.makedirs(rd, exist_ok=True)
        jobs = [(sysxml, toppdb, seeds_xyz[k], 10_000 * (r + 1) + k, steps, report,
                 os.path.join(rd, f"w{k:02d}.npz"), a.platform, a.temperature,
                 a.timestep_fs) for k in range(a.walkers)]
        t0 = time.time()
        nproc = a.walkers if a.platform == "CPU" else min(a.walkers, 8)
        with ProcessPoolExecutor(max_workers=max(1, min(nproc, os.cpu_count() or 1))) as ex:
            for _, shape in ex.map(_run_walker, jobs):
                pass
        for k in range(a.walkers):
            all_xyz.append(np.load(os.path.join(rd, f"w{k:02d}.npz"))["xyz"])
        wall = time.time() - t0
        done_ns = (r + 1 - start_round) * a.walkers * a.ns_per_walker

        msm, active, dall, its, km, Y = analyse(all_xyz, sel, a.tica_lag, a.dim,
                                                a.n_micro, a.lags)
        ok, drift = its_converged(its)
        print(f"  round {r+1:3d}/{a.rounds}  {wall/60:6.1f} min  "
              f"{done_ns:7.1f} ns total  states {len(np.unique(dall)):4d}  "
              f"ITS1 {its[-1] if np.isfinite(its[-1]) else float('nan'):9.2f} frames  "
              f"drift {drift:5.2f}  {'CONVERGED' if ok else ''}")
        history.append({"round": r, "ns": done_ns, "its": its, "drift": drift,
                        "converged": bool(ok), "n_states": int(len(np.unique(dall)))})
        json.dump(history, open(os.path.join(a.out, "history.json"), "w"), indent=1)

        if ok and r > 3:
            print("  ITS test passed - stopping early.")
            break

        sel_states = select_seed_states(a.walkers, dall, rng, strategy=a.strategy,
                                        msm=msm, active=active,
                                        round_index=r, n_rounds=a.rounds)
        cfg = pick_configurations(sel_states, dall, rng)
        pooled = np.concatenate(all_xyz, 0)
        seeds_xyz = [pooled[i] for i in cfg]

    print(f"\n  total wall time {(time.time()-t_start)/3600:.2f} h")
    print(f"  history -> {os.path.join(a.out, 'history.json')}")
    if history and not history[-1]["converged"]:
        print("  NOT CONVERGED: the ITS test never flattened.  Do not quote rates;")
        print("  extend the campaign (--resume) or widen --lags.")


if __name__ == "__main__":
    main()
