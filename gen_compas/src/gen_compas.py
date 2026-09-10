"""Gen-COMPAS: generative committor-guided path sampling.

Reimplementation of the framework of Tang, Pandey et al., arXiv:2510.24979v1,
applied to the paper's own proof-of-concept system (NANMA in vacuum).

The loop, following Fig. 1 of the paper:

  (0) very short unbiased MD in each metastable state -> seed dataset
  repeat:
  (1) train a denoising diffusion model on the accumulated conformations
  (2) generate intermediates connecting the two states (latent interpolation)
  (3) learn the committor q(x) in conformational space from the accumulated
      commitment outcomes; keep generated targets with q ~ 1/2
  (4) targeted MD from A and from B onto each target -> physical structures
  (5) unbiased shooting from those structures -> new data + outcomes
  (6) feed everything back into (1)

No collective variable is used anywhere in the loop.  (phi, psi) are computed
only to define the two metastable states and, afterwards, to plot results.
"""

import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402
import mdops  # noqa: E402
import models  # noqa: E402
from openmm import unit as openmm_unit  # noqa: E402
from build_system import build_coordinates, PDB_ORDER  # noqa: E402


# ---------------------------------------------------------------------------
# Data container
# ---------------------------------------------------------------------------


class Store:
    """Accumulated unbiased sampling."""

    def __init__(self):
        self.coords = []      # list of (n_frames, 22, 3)
        self.labels = []      # list of (n_frames,) per-frame committor labels
        self.outcome = []     # per trajectory: 0, 1 or -1
        self.source = []      # per trajectory: str
        self.iteration = []   # per trajectory: int

    def add(self, coords, outcome, source, iteration, labels=None):
        coords = np.asarray(coords, dtype=np.float32)
        if labels is None:
            labels = np.full(len(coords), int(outcome), dtype=np.int8)
        self.coords.append(coords)
        self.labels.append(np.asarray(labels, dtype=np.int8))
        self.outcome.append(int(outcome))
        self.source.append(source)
        self.iteration.append(int(iteration))

    def n_frames(self):
        return sum(len(c) for c in self.coords)

    def labelled_frames(self, stride=3, basin_cap=400, rng=None):
        """Frames whose own future reaches a core, for committor fitting.

        Each such frame contributes one Bernoulli sample of q(x); this is the
        definition of the committor, so the fit is a Monte-Carlo estimator of
        it.  Deep-basin frames are numerous and redundant, so they are capped.
        """
        rng = rng or np.random.default_rng(0)
        xs, ys, groups = [], [], []
        for j, (c, lab, s) in enumerate(zip(self.coords, self.labels,
                                            self.source)):
            sel = np.arange(0, len(c), stride)
            sel = sel[lab[sel] >= 0]
            if len(sel) == 0:
                continue
            if s.startswith("seed") and len(sel) > basin_cap:
                sel = np.sort(rng.choice(sel, basin_cap, replace=False))
            xs.append(c[sel])
            ys.append(lab[sel].astype(np.float32))
            groups.append(np.full(len(sel), j))
        if not xs:
            return None, None, None
        return (np.concatenate(xs), np.concatenate(ys), np.concatenate(groups))

    def core_frames(self, label):
        """Frames actually inside a macrostate, used as interpolation ends."""
        out = []
        for c in self.coords:
            phi, _ = common.phi_psi(c)
            m = common.which_basin(phi) == label
            if m.any():
                out.append(c[m])
        return np.concatenate(out) if out else np.zeros((0, 22, 3), np.float32)

    def all_frames(self, stride=1):
        return np.concatenate([c[::stride] for c in self.coords])

    def save(self, path):
        np.savez_compressed(
            path,
            coords=np.concatenate(self.coords),
            labels=np.concatenate(self.labels),
            lengths=np.array([len(c) for c in self.coords]),
            outcome=np.array(self.outcome),
            source=np.array(self.source),
            iteration=np.array(self.iteration),
        )


# ---------------------------------------------------------------------------
# Seeding
# ---------------------------------------------------------------------------


def find_basin_minimum(engine, label, verbose=True):
    """Locate the deepest potential-energy minimum inside a macrostate.

    A single idealised Z-matrix geometry is not good enough to start from: for
    NANMA it relaxes into a side minimum ~11 kcal/mol above the bottom of the
    basin, and releasing that much energy into 22 atoms is itself enough to
    push the molecule over the 7 kcal/mol barrier within a few hundred
    picoseconds.  Instead the whole basin is scanned on a coarse (phi, psi)
    grid, every point is minimised, and the lowest minimum that is still
    inside the macrostate is returned.  Only the state definition is used.
    """
    phis = np.arange(-165.0, 180.0, 30.0)
    psis = np.arange(-165.0, 180.0, 30.0)
    best = None
    for phi in phis:
        for psi in psis:
            coords = build_coordinates(float(phi), float(psi))
            xyz = np.array([coords[f"{n}_{r}"] for n, r in PDB_ORDER]) / 10.0
            engine.set_positions(xyz)
            try:
                engine.sim.minimizeEnergy(maxIterations=4000)
            except Exception:
                engine.recover()
                continue
            st = engine.sim.context.getState(getEnergy=True, getPositions=True)
            e = st.getPotentialEnergy().value_in_unit(
                openmm_unit.kilocalorie_per_mole)
            x = st.getPositions(asNumpy=True).value_in_unit(
                openmm_unit.nanometer)
            ph, ps = common.phi_psi(x)
            if int(common.which_basin(np.array(ph))) != label:
                continue
            if best is None or e < best[0]:
                best = (e, np.array(x), float(ph), float(ps))
    if best is None:
        raise RuntimeError(f"no minimum found inside state {label}")
    if verbose:
        print(f"    basin minimum: E = {best[0]:.2f} kcal/mol at "
              f"phi={best[2]:.0f}, psi={best[3]:.0f}", flush=True)
    return best[1]


def make_state_structure(engine, label, warmup_ps=200.0, tries=8,
                         verbose=True):
    """Basin minimum, then equilibrate, then check it is still in the basin."""
    x_min = find_basin_minimum(engine, label, verbose=verbose)
    for t in range(tries):
        engine.set_positions(x_min)
        engine.randomize_velocities()
        engine.step(int(warmup_ps / 0.002))
        x = engine.positions()
        ph, ps = common.phi_psi(x)
        if int(common.which_basin(np.array(ph))) == label:
            if verbose:
                print(f"    equilibrated {warmup_ps:.0f} ps, "
                      f"phi={ph:.0f} psi={ps:.0f}"
                      + (f" (attempt {t + 1})" if t else ""), flush=True)
            return x
        if verbose:
            print(f"    equilibration left the basin (phi={ph:.0f}); "
                  f"redrawing velocities", flush=True)
    raise RuntimeError(f"could not equilibrate inside state {label}")


def seed(engine, store, ns_per_state=1.0, save_ps=0.2, verbose=True):
    """Very short unbiased MD in each metastable state (paper: 1-2 ns)."""
    starts = {}
    for label, name in enumerate([common.CORE_A["name"], common.CORE_B["name"]]):
        if verbose:
            print(f"  [seed {name}]", flush=True)
        x0 = make_state_structure(engine, label, verbose=verbose)
        engine.set_positions(x0)
        save_every = int(save_ps / 0.002)
        n_saves = int(ns_per_state * 1000 / save_ps)
        frames = []
        for _ in range(n_saves):
            engine.step(save_every)
            frames.append(engine.positions())
        frames = np.asarray(frames, dtype=np.float32)
        ph, ps = common.phi_psi(frames)
        core = common.which_basin(ph)
        frac = float(np.mean(core == label))
        if verbose:
            print(f"    {ns_per_state} ns unbiased: "
                  f"{frac * 100:.0f}% of frames in the target core, "
                  f"phi in [{ph.min():.0f}, {ph.max():.0f}]", flush=True)
        # Keep only in-core frames, split into time-contiguous runs: a masked
        # trajectory is not a trajectory, and gluing the pieces together would
        # invent transitions that never happened.
        inside = core == label
        kept = 0
        for piece in _contiguous(inside):
            if len(piece) < 5:
                continue
            store.add(frames[piece], label, f"seed_{label}", 0)
            kept += len(piece)
        starts[label] = frames[inside]
        if verbose and kept < inside.sum():
            print(f"    (kept {kept} of {int(inside.sum())} in-core frames "
                  f"after splitting into contiguous runs)", flush=True)
    return starts


def _contiguous(mask):
    """Index arrays for each maximal run of True in a boolean mask."""
    idx = np.where(mask)[0]
    if len(idx) == 0:
        return []
    breaks = np.where(np.diff(idx) > 1)[0] + 1
    return np.split(idx, breaks)


# ---------------------------------------------------------------------------
# One Gen-COMPAS iteration
# ---------------------------------------------------------------------------


def iteration(engine, store, it, cfg, rng, log, sep_store):
    t0 = time.time()
    ns_before = engine.ns_used()

    # ---- (1) train the generative model on everything sampled so far ------
    diff = None
    if cfg["generator"] == "diffusion":
        all_x = store.all_frames(stride=cfg["diff_stride"])
        feats = common.featurize(all_x)
        diff = models.DiffusionModel(feats.shape[1], n_steps=cfg["ddpm_steps"],
                                     hidden=cfg["ddpm_hidden"], seed=it)
        diff.fit(feats, steps=cfg["ddpm_steps_train"], batch=256,
                 lr=2e-4, verbose=cfg["verbose"], seed=it)

    # ---- (3a) learn the committor from the accumulated outcomes -----------
    lx, ly, lg = store.labelled_frames(stride=cfg["committor_stride"], rng=rng)
    lf = common.featurize(lx)
    qnet = models.train_committor(lf, ly, dim=lf.shape[1],
                                  steps=cfg["committor_steps_train"], seed=it,
                                  verbose=cfg["verbose"])

    # ---- (2) generate intermediates ---------------------------------------
    a_frames = store.core_frames(0)
    b_frames = store.core_frames(1)
    m = cfg["n_pairs"]
    ia = rng.choice(len(a_frames), m)
    ib = rng.choice(len(b_frames), m)
    fa = common.featurize(a_frames[ia])
    fb = common.featurize(b_frames[ib])
    lambdas = np.linspace(0.15, 0.85, cfg["n_lambda"])
    if diff is not None:
        gen = diff.interpolate(fa, fb, lambdas, n_infer=cfg["ddim_steps"])
    else:
        # Ablation: straight-line interpolation in the same feature space,
        # i.e. the generative model removed and nothing else changed.
        gen = np.concatenate([(1.0 - lam) * fa + lam * fb for lam in lambdas])
    n_physical = int(sum(mdops.target_is_physical(g.reshape(-1, 3))
                         for g in gen))

    # ---- (3b) committor filter: keep targets on the separatrix ------------
    qgen = qnet.predict(gen)
    band = cfg["q_band"]
    mask = np.abs(qgen - 0.5) < band
    if mask.sum() < cfg["n_targets"]:
        order = np.argsort(np.abs(qgen - 0.5))
        idx = order[: cfg["n_targets"]]
    else:
        cand = np.where(mask)[0]
        idx = rng.choice(cand, cfg["n_targets"], replace=False)
    targets = gen[idx]
    q_targets = qgen[idx]

    # ---- (4) targeted MD from A and from B onto each target ---------------
    n_heavy = len(common.heavy_atoms())
    heavy_idx = common.heavy_atoms()
    tmd_products = []
    tmd_rmsd = []
    n_refined = 0
    for k, tgt in enumerate(targets):
        tgt_xyz = tgt.reshape(n_heavy, 3)
        for side, pool in ((0, a_frames), (1, b_frames)):
            start = pool[rng.integers(len(pool))]
            res = mdops.targeted_md(engine, start, tgt_xyz,
                                    n_steps=cfg["tmd_steps"],
                                    k=cfg["tmd_k"],
                                    relax_steps=cfg["tmd_relax"])
            if res is None:
                continue
            tmd_rmsd.append(res["rmsd"])
            # The steering carries the molecule across the barrier, so the
            # separatrix lies between two consecutive path points.  Pick the
            # path structures whose predicted committor is closest to 1/2 --
            # a one-dimensional bracketing problem, which the committor can
            # solve reliably even while it is still poorly calibrated in the
            # full conformational space.
            path = res["path"]
            qp = qnet.predict(common.featurize(path))

            # The crossing occupies about one steering window out of sixty, so
            # the two path points nearest q = 1/2 straddle it rather than sit
            # on it.  Where consecutive points bracket 1/2, steer slowly from
            # one to the other: the same machinery, applied over the width of
            # a single window, resolves the top of the barrier.
            cross = np.where((qp[:-1] - 0.5) * (qp[1:] - 0.5) < 0)[0]
            used_refined = False
            if len(cross):
                j = int(cross[len(cross) // 2])
                ref = mdops.targeted_md(
                    engine, path[j], path[j + 1][heavy_idx],
                    n_steps=cfg["refine_steps"], k=cfg["tmd_k"],
                    relax_steps=1, n_windows=cfg["refine_windows"],
                    n_path=cfg["refine_windows"])
                if ref is not None:
                    rpath = ref["path"]
                    rq = qnet.predict(common.featurize(rpath))
                    for jj in np.argsort(np.abs(rq - 0.5))[
                            : cfg["n_path_points"]]:
                        tmd_products.append((rpath[jj], side, k, float(rq[jj])))
                    n_refined += 1
                    used_refined = True
            if not used_refined:
                for j in np.argsort(np.abs(qp - 0.5))[: cfg["n_path_points"]]:
                    tmd_products.append((path[j], side, k, float(qp[j])))

    # ---- (5) unbiased shooting from the separatrix ------------------------
    n_sep = 0
    shot_q = []
    shot_qpred = []
    commit_times = []
    shot_pe = []
    sep_points = []
    for p, (x, side, k, q_pred) in enumerate(tmd_products):
        shot_pe.append(engine.potential_energy(x))
        outs = []
        for s in range(cfg["n_shots"]):
            r = mdops.shoot(engine, x, total_ps=cfg["shoot_total_ps"],
                            save_ps=cfg["save_ps"],
                            seed=int(rng.integers(1, 2 ** 30)))
            if r is None:
                continue
            store.add(r["coords"], r["outcome"], "shoot", it,
                      labels=r["labels"])
            if r["outcome"] >= 0:
                outs.append(r["outcome"])
                if np.isfinite(r["commit_ps"]):
                    commit_times.append(r["commit_ps"])
        if outs:
            qe = float(np.mean(outs))
            shot_q.append(qe)
            shot_qpred.append(q_pred)
            if 0.2 <= qe <= 0.8:
                n_sep += 1
                sep_points.append((np.asarray(x, dtype=np.float32), qe,
                                   len(outs)))

    ns_used = engine.ns_used() - ns_before
    rec = dict(
        iteration=it,
        n_targets=int(len(targets)),
        n_generated=int(len(gen)),
        frac_physical=float(n_physical / max(1, len(gen))),
        q_targets_mean=float(np.mean(q_targets)),
        n_tmd=len(tmd_products),
        n_refined=int(n_refined),
        tmd_rmsd_nm=float(np.mean(tmd_rmsd)) if tmd_rmsd else float("nan"),
        n_shot_points=len(shot_q),
        n_separatrix_points=n_sep,
        frac_separatrix=float(n_sep / max(1, len(shot_q))),
        empirical_q_mean=float(np.mean(shot_q)) if shot_q else float("nan"),
        committor_mae=float(np.mean(np.abs(np.array(shot_qpred)
                                           - np.array(shot_q))))
        if shot_q else float("nan"),
        commit_ps_median=float(np.median(commit_times))
        if commit_times else float("nan"),
        shot_pe_median=float(np.median(shot_pe)) if shot_pe else float("nan"),
        frames=store.n_frames(),
        ns_this_iteration=ns_used,
        ns_cumulative=engine.ns_used(),
        wall_s=time.time() - t0,
    )
    log.append(rec)
    print(f"  [iter {it}] gen={rec['n_generated']} "
          f"physical={rec['frac_physical'] * 100:.0f}% "
          f"targets={rec['n_targets']} tmd={rec['n_tmd']} "
          f"refined={rec['n_refined']} "
          f"shot_pts={rec['n_shot_points']} "
          f"on-separatrix={rec['frac_separatrix'] * 100:.0f}% "
          f"<q_emp>={rec['empirical_q_mean']:.2f} "
          f"|q_pred-q_emp|={rec['committor_mae']:.2f} "
          f"t_commit={rec['commit_ps_median']:.1f}ps "
          f"E={rec['shot_pe_median']:.0f} "
          f"frames={rec['frames']} "
          f"ns={rec['ns_this_iteration']:.2f} (cum {rec['ns_cumulative']:.2f}) "
          f"[{rec['wall_s']:.0f}s]", flush=True)
    sep_store.extend(sep_points)
    return diff, qnet


# ---------------------------------------------------------------------------


DEFAULT_CFG = dict(
    seed_ns=1.0,
    n_iterations=5,
    n_pairs=24,
    n_lambda=7,
    n_targets=10,
    q_band=0.15,
    n_shots=3,
    # The phi crossing in NANMA is quasi-ballistic: commitment happens
    # within a few hundred femtoseconds, so a 20 ps shot would spend 95%
    # of its cost watching the molecule sit in a basin.  Shots are short
    # and finely saved; the seed runs, which sample the basins, are not.
    shoot_total_ps=6.0,
    save_ps=0.1,
    seed_save_ps=0.2,
    tmd_steps=6000,
    tmd_relax=100,
    n_path_points=2,
    refine_steps=4000,
    refine_windows=40,
    tmd_k=20000.0,
    generator="diffusion",
    ddpm_steps=400,
    ddpm_hidden=384,
    ddpm_steps_train=8000,
    ddim_steps=40,
    diff_stride=2,
    committor_stride=3,
    committor_steps_train=4000,
    verbose=False,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iterations", type=int, default=DEFAULT_CFG["n_iterations"])
    ap.add_argument("--seed-ns", type=float, default=DEFAULT_CFG["seed_ns"])
    ap.add_argument("--targets", type=int, default=DEFAULT_CFG["n_targets"])
    ap.add_argument("--shots", type=int, default=DEFAULT_CFG["n_shots"])
    ap.add_argument("--rng", type=int, default=0)
    ap.add_argument("--tag", type=str, default="run")
    ap.add_argument("--generator", choices=["diffusion", "linear"],
                    default="diffusion",
                    help="linear = ablation with the generative model removed")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    cfg = dict(DEFAULT_CFG)
    cfg.update(n_iterations=args.iterations, seed_ns=args.seed_ns,
               n_targets=args.targets, n_shots=args.shots,
               generator=args.generator, verbose=args.verbose)

    rng = np.random.default_rng(args.rng)
    engine = mdops.Engine(threads=1, seed=args.rng + 1)
    store = Store()
    log = []
    sep_store = []

    print(f"[gen-compas] tag={args.tag} seed={args.rng} "
          f"generator={args.generator}", flush=True)
    t_all = time.time()
    seed(engine, store, ns_per_state=cfg["seed_ns"],
         save_ps=cfg["seed_save_ps"])
    print(f"  [seed] cumulative MD = {engine.ns_used():.2f} ns", flush=True)

    for it in range(cfg["n_iterations"]):
        diff, qnet = iteration(engine, store, it, cfg, rng, log, sep_store)

    out = os.path.join(common.RESULTS, f"gencompas_{args.tag}")
    store.save(out + "_store.npz")
    if sep_store:
        np.savez_compressed(out + "_separatrix.npz",
                            coords=np.stack([s0 for s0, _, _ in sep_store]),
                            q_emp=np.array([q for _, q, _ in sep_store]),
                            n_shots=np.array([n for _, _, n in sep_store]))
    with open(out + "_log.json", "w") as fh:
        json.dump(dict(cfg={k: v for k, v in cfg.items()},
                       log=log,
                       total_ns=engine.ns_used(),
                       wall_s=time.time() - t_all), fh, indent=2)

    import torch
    torch.save(qnet.state_dict(), out + "_committor.pt")
    print(f"[gen-compas] done: {engine.ns_used():.2f} ns of MD, "
          f"{store.n_frames()} frames, {time.time() - t_all:.0f} s wall",
          flush=True)


if __name__ == "__main__":
    main()
