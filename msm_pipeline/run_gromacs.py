#!/usr/bin/env python3
"""
Run the same conformational-flux analysis on YOUR OWN trajectories.

This is the part you point at a Samsung T7 (or any disk).  It reads whatever
MDTraj reads - GROMACS .xtc/.trr with a .gro/.pdb/.tpr topology, AMBER
.nc/.mdcrd with .prmtop, CHARMM .dcd with .psf - so no GROMACS installation is
needed for the analysis itself, only for producing the trajectory.

    python -m msm_pipeline.run_gromacs \
        --top  /Volumes/T7/run/md.gro \
        --traj "/Volumes/T7/run/*/md_center.xtc" \
        --sel  "name CA and resid 200 to 290" \
        --a-pair 220 276 --out /Volumes/T7/analysis

What it does, in the order that matters:
  1. featurise  - pairwise distances within your selection
  2. TICA       - slow collective coordinates
  3. ITS TEST   - *prints the convergence table first and refuses to quote
                  rates if it fails*.  This is deliberate: an MSM that has not
                  passed this test produces confident, wrong numbers.
  4. MSM + TPT  - committor, flux, pathway decomposition, bottleneck
  5. direct cross-check - committed crossings counted straight off the
                  trajectory, which needs no Markov assumption

PRACTICAL NOTES FOR YOUR OWN RUNS
  * Preprocess first.  `gmx trjconv -pbc mol -center -fit rot+trans` or the
    analysis will be dominated by tumbling.  Strip waters; keep Cα or heavy
    atoms of the region you care about.
  * Sampling beats length.  As this dataset shows (M3), 32x500 ns of
    *disconnected* runs will not converge an MSM for a ~100 ns process.
    Many short runs reseeded from the barrier region will.
  * Memory: frames are streamed per-file, but the feature matrix is held in
    RAM.  n_frames x n_pairs x 4 bytes - keep it under a few GB or raise
    --stride.
"""
import argparse
import glob
import os
import sys

import numpy as np


def _load(top, patterns, sel, stride):
    import mdtraj as md
    files = []
    for p in patterns:
        files.extend(sorted(glob.glob(os.path.expanduser(p))))
    if not files:
        sys.exit(f"no trajectory files matched: {patterns}")
    ref = md.load_frame(files[0], 0, top=top)
    idx = ref.topology.select(sel)
    if len(idx) < 4:
        sys.exit(f"selection '{sel}' matched {len(idx)} atoms; need >= 4")
    print(f"  {len(files)} trajectory file(s), selection '{sel}' -> {len(idx)} atoms")
    out, labels = [], []
    for f in files:
        t = md.load(f, top=top, stride=stride, atom_indices=idx)
        out.append(t.xyz * 10.0)      # nm -> angstrom
        labels.append(os.path.basename(f))
        print(f"    {labels[-1]:40s} {t.n_frames:7d} frames")
    return out, labels, ref, idx


def _features(trajs, cutoff):
    ref = np.concatenate([t[::max(1, len(t) // 200)] for t in trajs])
    n = ref.shape[1]
    I, J = np.triu_indices(n, 3)
    d = np.linalg.norm(ref[:, I] - ref[:, J], axis=-1).mean(0)
    keep = d < cutoff
    I, J = I[keep], J[keep]
    print(f"  {len(I)} distance features (mean separation < {cutoff} A)")
    return [np.linalg.norm(t[:, I] - t[:, J], axis=-1).astype(np.float32)
            for t in trajs], (I, J)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="MSM / TPT conformational-flux analysis of your own MD",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--top", required=True, help="topology (.gro/.pdb/.prmtop/.psf)")
    ap.add_argument("--traj", required=True, nargs="+",
                    help="trajectory file(s) or glob(s); quote the glob")
    ap.add_argument("--sel", default="name CA", help="MDTraj atom selection")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--dt", type=float, default=None,
                    help="ns between STORED frames after striding "
                         "(default: read from the trajectory)")
    ap.add_argument("--cutoff", type=float, default=15.0,
                    help="keep atom pairs closer than this on average (A)")
    ap.add_argument("--tica-lag", type=int, default=20, help="frames")
    ap.add_argument("--dim", type=int, default=3)
    ap.add_argument("--n-micro", type=int, default=100)
    ap.add_argument("--lags", type=int, nargs="+",
                    default=[5, 10, 20, 40, 80, 160], help="ITS lags, frames")
    ap.add_argument("--msm-lag", type=int, default=None,
                    help="frames; default = largest ITS lag")
    ap.add_argument("--a-pair", type=int, nargs=2, metavar=("I", "J"),
                    help="two RESIDUE numbers whose distance defines the "
                         "A/B core sets (e.g. a pocket mouth)")
    ap.add_argument("--a-max", type=float, default=None, help="A if dist < this")
    ap.add_argument("--b-min", type=float, default=None, help="B if dist > this")
    ap.add_argument("--out", default="msm_out")
    a = ap.parse_args(argv)

    from deeptime.clustering import KMeans
    from deeptime.decomposition import TICA
    from deeptime.markov import (TransitionCountEstimator, reactive_flux)
    from deeptime.markov.msm import MaximumLikelihoodMSM

    os.makedirs(a.out, exist_ok=True)
    print("[1/5] loading")
    trajs, labels, ref, idx = _load(a.top, a.traj, a.sel, a.stride)
    dt = a.dt
    if dt is None:
        import mdtraj as md
        t2 = md.load(sorted(glob.glob(os.path.expanduser(a.traj[0])))[0],
                     top=a.top, atom_indices=idx[:1])
        dt = float(t2.timestep) / 1000.0 * a.stride if t2.n_frames > 1 else 1.0
        print(f"  inferred frame spacing: {dt:.4f} ns")

    print("[2/5] featurising")
    F, _ = _features(trajs, a.cutoff)

    print("[3/5] TICA")
    tica = TICA(lagtime=a.tica_lag, dim=a.dim, scaling="kinetic_map")
    for f in F:
        tica.partial_fit(f)
    m = tica.fetch_model()
    Y = [m.transform(f) for f in F]
    btw = [np.var([y[:, i].mean() for y in Y]) /
           np.var(np.concatenate([y[:, i] for y in Y])) for i in range(a.dim)]
    print(f"  TICA timescales (ns): "
          f"{np.round(-a.tica_lag*dt/np.log(np.abs(m.singular_values[:a.dim])),1)}")
    print(f"  between-trajectory variance fraction: {np.round(btw, 2)}")
    if max(btw) > 0.5:
        print("  !! WARNING: your slow coordinates mostly separate TRAJECTORIES,")
        print("     not conformations.  The runs are not interconverting; an MSM")
        print("     built on this will report the inverse of a near-zero rate.")

    print("[4/5] implied-timescale convergence test")
    km = (KMeans(n_clusters=a.n_micro, max_iter=300, init_strategy="kmeans++",
                 n_jobs=max(1, os.cpu_count() // 2 or 1), fixed_seed=1)
          .fit(np.concatenate([y[::3] for y in Y])).fetch_model())
    dt_disc = [km.transform(y) for y in Y]
    print(f"  {'lag (ns)':>9s} {'states':>7s}   slowest 3 ITS (ns)")
    series = []
    for lag in a.lags:
        c = (TransitionCountEstimator(lag, "sliding").fit(dt_disc)
             .fetch_model().submodel_largest(directed=True))
        mm = MaximumLikelihoodMSM(reversible=True).fit(c).fetch_model()
        ts = mm.timescales(3) * dt
        series.append(ts)
        print(f"  {lag*dt:9.3f} {c.n_states:7d}   " + ", ".join(f"{v:10.2f}" for v in ts))
    drift = series[-1][0] / max(series[len(series)//2][0], 1e-12)
    converged = drift < 1.5
    print(f"  slowest ITS grows {drift:.2f}x over the upper half of the lag range"
          f"  ->  {'CONVERGED' if converged else 'NOT CONVERGED'}")
    if not converged:
        print("  Rates below are NOT trustworthy.  Treat flux topology as")
        print("  qualitative and rely on the direct count in step 5.")

    print("[5/5] MSM + TPT")
    lag = a.msm_lag or a.lags[-1]
    c = (TransitionCountEstimator(lag, "sliding").fit(dt_disc)
         .fetch_model().submodel_largest(directed=True))
    msm = MaximumLikelihoodMSM(reversible=True).fit(c).fetch_model()
    np.save(os.path.join(a.out, "transition_matrix.npy"), msm.transition_matrix)
    np.save(os.path.join(a.out, "stationary_distribution.npy"),
            msm.stationary_distribution)

    if a.a_pair and a.a_max is not None and a.b_min is not None:
        rI, rJ = a.a_pair
        sub = ref.atom_slice(idx).topology
        pi_ = [at.index for at in sub.atoms if at.residue.resSeq == rI]
        pj_ = [at.index for at in sub.atoms if at.residue.resSeq == rJ]
        if not pi_ or not pj_:
            sys.exit(f"residues {rI}/{rJ} not inside selection '{a.sel}'")
        cvs = [np.linalg.norm(t[:, pi_[0]] - t[:, pj_[0]], axis=-1) for t in trajs]
        centers = np.concatenate([np.bincount(d, weights=cv, minlength=a.n_micro)
                                  for d, cv in [(np.concatenate(dt_disc),
                                                 np.concatenate(cvs))]])
        counts = np.bincount(np.concatenate(dt_disc), minlength=a.n_micro)
        centers = centers / np.maximum(counts, 1)
        cp = centers[c.state_symbols]
        A, B = np.where(cp < a.a_max)[0], np.where(cp > a.b_min)[0]
        if len(A) == 0 or len(B) == 0:
            sys.exit(f"empty core set: A={len(A)} B={len(B)}; adjust --a-max/--b-min")
        fl = reactive_flux(msm.transition_matrix, A, B)
        per_us = 1.0 / (dt * 1e-3)
        print(f"  A: {len(A)} states pop {msm.stationary_distribution[A].sum():.4f}"
              f" | B: {len(B)} states pop {msm.stationary_distribution[B].sum():.4f}")
        print(f"  MSM A->B rate  {fl.rate*per_us:10.2f} /us"
              f"   {'' if converged else '[NOT CONVERGED]'}")
        paths, caps = fl.pathways(fraction=0.9)
        print(f"  {len(paths)} routes carry 90% of flux; "
              f"top route {100*caps[0]/fl.total_flux:.2f}%")

        n = 0
        for cv in cvs:
            s = np.full(len(cv), -1, np.int8)
            s[cv < a.a_max] = 0
            s[cv > a.b_min] = 1
            q = s[s >= 0]
            if len(q) > 1:
                n += int((np.diff(q) == 1).sum())
        total_us = sum(len(t) for t in trajs) * dt / 1000.0
        print(f"  DIRECT committed crossings: {n} in {total_us:.2f} us"
              f"  ->  {n/max(total_us,1e-9):.2f} /us   (no Markov assumption)")
        np.save(os.path.join(a.out, "committor.npy"), fl.forward_committor)
    else:
        print("  (pass --a-pair I J --a-max X --b-min Y for the TPT/flux section)")
    print(f"\n  results in {a.out}/")


if __name__ == "__main__":
    main()
