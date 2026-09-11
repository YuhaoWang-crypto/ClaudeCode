"""Downstream analysis of a Gen-COMPAS run.

Produces, from the accumulated unbiased sampling alone:
  (1) the transition-state ensemble,
  (2) the committor projected onto interpretable CVs defined after the fact,
  (3) committor-consistent transition pathways,
  (4) the free-energy landscape,
and compares (4) with the metadynamics reference.
"""

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402
import models  # noqa: E402
import msm as msmlib  # noqa: E402

KT = common.KT_KCAL


def load_store(tag):
    d = np.load(os.path.join(common.RESULTS, f"gencompas_{tag}_store.npz"),
                allow_pickle=True)
    lengths = d["lengths"]
    coords = d["coords"]
    segs = np.split(coords, np.cumsum(lengths)[:-1])
    labels = d["labels"] if "labels" in d.files else None
    lsegs = (np.split(labels, np.cumsum(lengths)[:-1])
             if labels is not None else
             [np.full(L, int(o), np.int8)
              for L, o in zip(lengths, d["outcome"])])
    return dict(segs=segs, lengths=lengths, coords=coords, lsegs=lsegs,
                outcome=d["outcome"], source=d["source"],
                iteration=d["iteration"])


def _jsonable(o):
    if hasattr(o, "tolist"):
        return o.tolist()
    if isinstance(o, (np.integer, np.floating, np.bool_)):
        return o.item()
    return str(o)


# Save intervals actually used by the sampler, in picoseconds.
SAVE_PS = {"seed": 0.2, "shoot": 0.1}
DT_PS = 0.2  # common time step the Markov model is built on


def uniform_dt_frames(store):
    """Frame indices, and per-trajectory lengths, at one common time step.

    Trajectories recorded more finely than DT_PS are strided down.  Returns
    (indices into store["coords"], lengths of each resulting trajectory).
    """
    idx, lengths = [], []
    start = 0
    for L, src in zip(store["lengths"], store["source"]):
        name = str(src).split("_")[0]
        dt = SAVE_PS.get(name, DT_PS)
        stride = max(1, int(round(DT_PS / dt)))
        sel = np.arange(start, start + L, stride)
        if len(sel) >= 2:
            idx.append(sel)
            lengths.append(len(sel))
        start += L
    return np.concatenate(idx), np.array(lengths)


def committor_training_set(store, rng_seed=0, core_ratio=1.0):
    """Frames to fit the committor on.

    Frames already inside a core carry no information about where the
    separatrix is -- their committor is 0 or 1 by definition -- and there are
    an order of magnitude more of them than transition-region frames, because
    a shot spends most of its length in a basin after committing.  Training on
    all of them makes the fit a hard classifier that saturates at 0 and 1 and
    never places a separatrix at all.  So every labelled transition-region
    frame is kept, and core frames are subsampled to a comparable number,
    enough to impose the boundary conditions.
    """
    rng = np.random.default_rng(rng_seed)
    trans_x, trans_y, core_x, core_y = [], [], [], []
    trans_g, core_g = [], []
    for gi, (c, lab) in enumerate(zip(store["segs"], store["lsegs"])):
        phi, _ = common.phi_psi(c)
        in_core = common.which_core(phi) >= 0
        ok = lab >= 0
        t = ok & ~in_core
        k = ok & in_core
        if t.any():
            trans_x.append(c[t])
            trans_y.append(lab[t].astype(float))
            trans_g.append(np.full(int(t.sum()), gi))
        if k.any():
            core_x.append(c[k])
            core_y.append(lab[k].astype(float))
            core_g.append(np.full(int(k.sum()), gi))
    tx = np.concatenate(trans_x) if trans_x else np.zeros((0, 22, 3), np.float32)
    ty = np.concatenate(trans_y) if trans_y else np.zeros(0)
    tg = np.concatenate(trans_g) if trans_g else np.zeros(0, int)
    cx = np.concatenate(core_x) if core_x else np.zeros((0, 22, 3), np.float32)
    cy = np.concatenate(core_y) if core_y else np.zeros(0)
    cg = np.concatenate(core_g) if core_g else np.zeros(0, int)
    n_core = int(min(len(cx), max(200, core_ratio * len(tx))))
    if len(cx) > n_core:
        pick = rng.choice(len(cx), n_core, replace=False)
        cx, cy, cg = cx[pick], cy[pick], cg[pick]
    return (np.concatenate([tx, cx]), np.concatenate([ty, cy]),
            len(tx), len(cx), np.concatenate([tg, cg]))


GRID = [(128, 1e-3, 4000), (64, 1e-2, 3000), (32, 3e-2, 2000),
        (16, 1e-1, 1500), (8, 3e-1, 800)]


def fit_committor(store, rng_seed=0, verbose=True):
    """Fit the committor, choosing capacity by held-out validation.

    Each frame contributes one Bernoulli sample, and frames within a
    trajectory are strongly correlated, so the effective sample size is far
    smaller than the frame count.  An unregularised network simply memorises
    the binary labels: its bin-averaged predictions match the observed
    frequencies perfectly while every individual prediction is 0 or 1, so it
    places no separatrix at all.  Capacity and weight decay are therefore
    chosen on trajectories held out of the fit, and the winner is refitted on
    everything.
    """
    xr, y, n_trans, n_core, groups = committor_training_set(store, rng_seed)
    x = common.featurize(xr)
    if verbose:
        print(f"[analysis] committor training set: {n_trans} transition-region "
              f"frames + {n_core} core frames")

    # Holding out whole trajectories is not enough: the three shots fired from
    # one shooting point are different trajectories but almost the same
    # configuration, so a memorising model transfers straight across the
    # split and wins on held-out loss.  Whole regions of configuration space
    # are held out instead, by clustering the frames and reserving clusters.
    rng = np.random.default_rng(rng_seed)
    _, cl = msmlib.cluster(x, n_clusters=min(40, max(4, len(x) // 40)),
                           seed=rng_seed)
    uniq = np.unique(cl)
    val_g = set(rng.choice(uniq, max(1, len(uniq) // 5),
                           replace=False).tolist())
    val = np.array([c in val_g for c in cl])
    if val.sum() < 20 or (~val).sum() < 20 or len(np.unique(y[val])) < 2:
        val = np.zeros(len(y), bool)

    best = None
    for hidden, wd, steps in GRID:
        if val.any():
            m = models.train_committor(x[~val], y[~val], dim=x.shape[1],
                                       steps=steps, hidden=hidden, n_blocks=1,
                                       weight_decay=wd, seed=rng_seed)
            p = np.clip(m.predict(x[val]), 1e-6, 1 - 1e-6)
            loss = float(-(y[val] * np.log(p)
                           + (1 - y[val]) * np.log(1 - p)).mean())
        else:
            loss = float("inf")
        if verbose:
            print(f"    hidden={hidden:4d} weight_decay={wd:<6} "
                  f"held-out log-loss = {loss:.4f}")
        if best is None or loss < best[0]:
            best = (loss, hidden, wd, steps)

    _, hidden, wd, steps = best
    if verbose:
        print(f"[analysis] committor: hidden={hidden}, weight_decay={wd}")
    model = models.train_committor(x, y, dim=x.shape[1], steps=steps,
                                   hidden=hidden, n_blocks=1,
                                   weight_decay=wd, seed=rng_seed)
    if verbose:
        calibration_table(model, xr, y, n_trans)
    return model, len(y)


def calibration_table(model, xr, y, n_trans):
    """Predicted committor against the frequency actually observed."""
    phi, _ = common.phi_psi(xr)
    q = model.predict(common.featurize(xr[:n_trans]))
    print("    committor calibration in the transition region:")
    print("      phi window      n   observed   predicted")
    for lo, hi in [(-45, -20), (-20, -10), (-10, 0), (0, 10), (10, 20),
                   (20, 55)]:
        m = (phi[:n_trans] > lo) & (phi[:n_trans] < hi)
        if m.sum():
            print(f"      [{lo:4d},{hi:4d})  {int(m.sum()):5d}   "
                  f"{y[:n_trans][m].mean():8.2f}   {q[m].mean():9.2f}")


def build_msm(store, n_clusters=150, lag=4, seed=0):
    # A Markov model needs one time step.  Seed runs are saved every 0.2 ps
    # and shots every 0.1 ps, so the shots are strided down to the common
    # 0.2 ps before any transition is counted; mixing them would count a lag
    # of different physical length in different trajectories.
    idx, lengths = uniform_dt_frames(store)
    feats = common.featurize(store["coords"][idx])
    centers, labels = msmlib.cluster(feats, n_clusters=n_clusters, seed=seed)
    lags = [1, 2, 3, 4, 6, 8, 12, 16]
    its = msmlib.implied_timescales(None, lags, labels, lengths)
    cm = msmlib.count_matrix(labels, lengths, lag)
    keep = msmlib.largest_connected_set(cm)
    t, pi_sub = msmlib.reversible_mle(cm[np.ix_(keep, keep)])
    pi = np.zeros(labels.max() + 1)
    pi[keep] = pi_sub
    mask = np.isin(labels, keep)
    w_sub = np.zeros(len(labels))
    w_sub[mask] = msmlib.frame_weights(labels[mask], pi)
    # map the weights back onto the full frame list
    w = np.zeros(len(store["coords"]))
    w[idx] = w_sub
    # macrostate membership of the (connected) microstates, from the frames
    phi_all, _ = common.phi_psi(store["coords"][idx])
    micro_state = np.full(labels.max() + 1, -1)
    for k in range(labels.max() + 1):
        sel = labels == k
        if not sel.any():
            continue
        c = common.which_basin(phi_all[sel])
        vals, cnt = np.unique(c, return_counts=True)
        micro_state[k] = vals[cnt.argmax()]
    sub_a = np.where(micro_state[keep] == 0)[0]
    sub_b = np.where(micro_state[keep] == 1)[0]
    print(f"[analysis] Markov model: {labels.max() + 1} microstates, "
          f"{len(keep)} in the largest connected set, "
          f"{100 * mask.mean():.0f}% of frames retained")
    print(f"           microstates assigned to A / B / neither: "
          f"{int((micro_state[keep] == 0).sum())} / "
          f"{int((micro_state[keep] == 1).sum())} / "
          f"{int((micro_state[keep] == -1).sum())}")
    if len(sub_a) == 0 or len(sub_b) == 0:
        print("           WARNING: a macrostate is missing from the connected "
              "set; free energies and passage times are not defined")
    lag_ps = lag * 0.2  # uniform frame spacing, see uniform_dt_frames
    mfpt_ab = msmlib.mfpt(t, pi_sub, sub_a, sub_b, lag_ps)
    mfpt_ba = msmlib.mfpt(t, pi_sub, sub_b, sub_a, lag_ps)
    return dict(feats=feats, labels=labels, centers=centers, pi=pi,
                weights=w, mask=mask, its=its, lags=lags, lag=lag, T=t,
                keep=keep, micro_state=micro_state, idx=idx,
                dt_ps=0.2, mfpt_ab=mfpt_ab, mfpt_ba=mfpt_ba)


def fel_from_weights(phi, psi, w, nbins=36):
    edges = np.linspace(-180.0, 180.0, nbins + 1)
    h, _, _ = np.histogram2d(phi, psi, bins=[edges, edges], weights=w)
    p = h / h.sum()
    with np.errstate(divide="ignore"):
        f = -KT * np.log(p)
    f -= np.nanmin(f[np.isfinite(f)])
    return f, edges


def metad_fes(nbins=36):
    from read_metad import load_fes
    fes, centers, _ = load_fes()
    if fes is None:
        return None, None
    # fes is indexed [psi, phi]; transpose to [phi, psi] and coarse-grain
    f = fes.T
    k = f.shape[0] // nbins
    if k > 1:
        p = np.exp(-f / KT)
        p = p.reshape(nbins, k, nbins, k).sum(axis=(1, 3))
        f = -KT * np.log(p)
    f -= f.min()
    edges = np.linspace(-180.0, 180.0, nbins + 1)
    return f, edges


def state_free_energies(f, edges):
    """Delta G between the two macrostates from a (phi, psi) landscape.

    Returns nan (not +-inf) when one of the states carries no weight, which
    happens if the Markov model's connected set does not reach it.
    """
    centers = 0.5 * (edges[1:] + edges[:-1])
    p = np.exp(-f / KT)
    p[~np.isfinite(p)] = 0.0
    core = common.which_basin(centers)
    pa = p[core == 0].sum()
    pb = p[core == 1].sum()
    if pa <= 0 or pb <= 0:
        return float("nan")
    return -KT * np.log(pb / pa)


def free_energy_along_phi(f, edges):
    """Marginalise a (phi, psi) landscape onto phi."""
    p = np.exp(-f / KT)
    p[~np.isfinite(p)] = 0.0
    pp = p.sum(1)
    with np.errstate(divide="ignore"):
        fphi = -KT * np.log(pp)
    return fphi - np.nanmin(fphi[np.isfinite(fphi)])


def barrier_height(f, edges, side=0):
    """Height of the phi ~ 0 saddle above each basin minimum, on the free
    energy marginalised onto phi."""
    centers = 0.5 * (edges[1:] + edges[:-1])
    fphi = free_energy_along_phi(f, edges)
    core = common.which_basin(centers)
    sel = (centers > -25) & (centers < 25)
    if not np.any(np.isfinite(fphi[sel])):
        return dict(saddle=float("nan"), from_A=float("nan"),
                    from_B=float("nan"))
    saddle = float(np.nanmin(fphi[sel]))
    fa = float(np.nanmin(fphi[core == 0])) if np.any(core == 0) else np.nan
    fb = float(np.nanmin(fphi[core == 1])) if np.any(core == 1) else np.nan
    return dict(saddle=saddle, from_A=saddle - fa, from_B=saddle - fb)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="main")
    ap.add_argument("--clusters", type=int, default=150)
    ap.add_argument("--lag", type=int, default=10)
    args = ap.parse_args()

    store = load_store(args.tag)
    phi, psi = common.phi_psi(store["coords"])
    print(f"[analysis] {len(store['segs'])} trajectories, "
          f"{len(store['coords'])} frames")

    # ---- committor -------------------------------------------------------
    qnet, n_lab = fit_committor(store)
    q = qnet.predict(common.featurize(store["coords"]))
    print(f"[analysis] committor trained on {n_lab} labelled frames")

    # ---- MSM / free energy -----------------------------------------------
    m = build_msm(store, n_clusters=args.clusters, lag=args.lag)
    print("[analysis] implied timescales (ps) vs lag (ps):")
    for lag, ts in zip(m["lags"], m["its"]):
        print(f"    lag={lag * 0.2:5.2f}   " +
              "  ".join(f"{t * 0.2:8.2f}" for t in ts))

    print(f"[analysis] MSM mean first passage time  A->B = "
          f"{m['mfpt_ab'] / 1000:.1f} ns,  B->A = {m['mfpt_ba'] / 1000:.1f} ns")

    w = m["weights"]
    fel, edges = fel_from_weights(phi, psi, w)
    ref, _ = metad_fes()

    dg = state_free_energies(fel, edges)
    print(f"\n[analysis] Gen-COMPAS  dG(B-A) = {dg:+.2f} kcal/mol")
    bar = barrier_height(fel, edges)
    print(f"[analysis] Gen-COMPAS  phi~0 barrier: "
          f"from A = {bar['from_A']:.2f}, from B = {bar['from_B']:.2f} kcal/mol")
    if ref is not None:
        dgr = state_free_energies(ref, edges)
        barr = barrier_height(ref, edges)
        print(f"[analysis] metadynamics dG(B-A) = {dgr:+.2f} kcal/mol")
        print(f"[analysis] metadynamics phi~0 barrier: "
              f"from A = {barr['from_A']:.2f}, from B = {barr['from_B']:.2f} kcal/mol")
        both = np.isfinite(fel) & np.isfinite(ref) & (ref < 8.0)
        rmse = float(np.sqrt(np.mean((fel[both] - ref[both]) ** 2)))
        print(f"[analysis] RMSE over the region below 8 kcal/mol "
              f"({both.sum()} bins) = {rmse:.2f} kcal/mol")

    # ---- transition-state ensemble ---------------------------------------
    tse = np.abs(q - 0.5) < 0.05
    print(f"\n[analysis] TSE: {tse.sum()} frames with |q-0.5|<0.05")
    if tse.sum():
        print(f"    phi  {np.percentile(phi[tse], [10, 50, 90]).round(0)}")
        print(f"    psi  {np.percentile(psi[tse], [10, 50, 90]).round(0)}")

    # ---- is sampling actually concentrating on the barrier? ---------------
    frame_iter = np.concatenate([np.full(L, it) for L, it in
                                 zip(store["lengths"], store["iteration"])])
    print("\n[analysis] sampling concentrated on the barrier region")
    print("    iteration   frames   in 0.1<q<0.9   fraction")
    per_iter = []
    for it in sorted(set(frame_iter.tolist())):
        sel = frame_iter == it
        in_barrier = sel & (q > 0.1) & (q < 0.9)
        per_iter.append(dict(iteration=int(it), frames=int(sel.sum()),
                             barrier=int(in_barrier.sum()),
                             fraction=float(in_barrier.sum()
                                            / max(1, sel.sum()))))
        print(f"    {it:9d}   {sel.sum():6d}   {in_barrier.sum():12d}   "
              f"{100 * in_barrier.sum() / max(1, sel.sum()):6.1f}%")

    # write representative structures for inspection
    write_representatives(store, q, phi, args.tag)

    # ---- committor-consistent pathways -----------------------------------
    paths = committor_paths(store, q, w, phi, psi)

    out = dict(
        tag=args.tag,
        n_traj=int(len(store["segs"])),
        n_frames=int(len(store["coords"])),
        dG_gencompas=float(dg),
        mfpt_ab_ns=float(m["mfpt_ab"]) / 1000.0,
        mfpt_ba_ns=float(m["mfpt_ba"]) / 1000.0,
        barrier_gencompas=bar,
        tse_n=int(tse.sum()),
        tse_phi=[float(v) for v in np.percentile(phi[tse], [10, 50, 90])]
        if tse.sum() else [],
        tse_psi=[float(v) for v in np.percentile(psi[tse], [10, 50, 90])]
        if tse.sum() else [],
        paths=paths,
        per_iteration=per_iter,
    )
    if ref is not None:
        out.update(dG_metad=float(dgr), barrier_metad=barr, fel_rmse=rmse)
    with open(os.path.join(common.RESULTS, f"analysis_{args.tag}.json"), "w") as fh:
        json.dump(out, fh, indent=2, default=_jsonable)
    np.savez_compressed(
        os.path.join(common.RESULTS, f"analysis_{args.tag}.npz"),
        fel=fel, edges=edges, ref=ref if ref is not None else np.zeros(0),
        phi=phi, psi=psi, q=q, weights=w, tse=tse,
        its=m["its"], lags=np.array(m["lags"]),
    )
    print(f"[analysis] wrote results/analysis_{args.tag}.*")


def write_representatives(store, q, phi, tag, n_each=10):
    """Multi-model PDB of reactant, transition-state and product structures."""
    from openmm.app import PDBFile
    from openmm import unit as u

    top, _, _, _ = common.get_system()
    groups = {
        "reactant": np.abs(q - 0.0) < 0.02,
        "transition": np.abs(q - 0.5) < 0.05,
        "product": np.abs(q - 1.0) < 0.02,
    }
    rng = np.random.default_rng(0)
    for name, sel in groups.items():
        idx = np.where(sel)[0]
        if len(idx) == 0:
            continue
        idx = rng.choice(idx, min(n_each, len(idx)), replace=False)
        path = os.path.join(common.RESULTS, f"{tag}_{name}.pdb")
        with open(path, "w") as fh:
            for k, i in enumerate(idx):
                PDBFile.writeModel(top, store["coords"][i] * u.nanometer,
                                   fh, modelIndex=k + 1)
        print(f"    wrote {os.path.basename(path)} "
              f"({len(idx)} structures, phi = "
              f"{np.percentile(phi[idx], [10, 90]).round(0)})")


def committor_paths(store, q, w, phi, psi, n_bins=16):
    """Committor-consistent strings.

    Transition-region frames are split into channels by k-means *in feature
    space* (no CV), then within each channel the equilibrium-weighted mean
    conformation is computed in each committor slab.  Projecting those means
    onto (phi, psi) gives the pathway.
    """
    sel = (q > 0.03) & (q < 0.97)
    if sel.sum() < 50:
        return []
    f = common.featurize(store["coords"][sel])
    centers, labels = msmlib.cluster(f, n_clusters=2, seed=1)
    qs, ws = q[sel], w[sel]
    ph, ps = phi[sel], psi[sel]
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    paths = []
    for ch in range(labels.max() + 1):
        m = labels == ch
        pts = []
        for i in range(n_bins):
            b = m & (qs >= edges[i]) & (qs < edges[i + 1])
            if b.sum() < 3:
                continue
            ww = ws[b]
            ww = ww / ww.sum() if ww.sum() > 0 else np.ones(b.sum()) / b.sum()
            # circular means
            def cmean(a):
                r = np.deg2rad(a)
                return float(np.degrees(np.arctan2((ww * np.sin(r)).sum(),
                                                   (ww * np.cos(r)).sum())))
            pts.append(dict(q=float(0.5 * (edges[i] + edges[i + 1])),
                            phi=cmean(ph[b]), psi=cmean(ps[b]),
                            n=int(b.sum())))
        if len(pts) >= 4:
            paths.append(dict(channel=int(ch), n_frames=int(m.sum()),
                              points=pts))
    print(f"\n[analysis] {len(paths)} committor-consistent pathway(s)")
    for p in paths:
        mid = [x for x in p["points"] if abs(x["q"] - 0.5) < 0.1]
        if mid:
            print(f"    channel {p['channel']}: {p['n_frames']} frames, "
                  f"crossing at phi={mid[0]['phi']:.0f} psi={mid[0]['psi']:.0f}")
    return paths


if __name__ == "__main__":
    main()
