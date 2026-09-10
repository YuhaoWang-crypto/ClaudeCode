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


def fit_committor(store, rng_seed=0, steps=8000):
    """Refit the committor on the complete accumulated dataset."""
    xs, ys = [], []
    rng = np.random.default_rng(rng_seed)
    for c, lab, s in zip(store["segs"], store["lsegs"], store["source"]):
        sel = np.arange(0, len(c), 3)
        sel = sel[lab[sel] >= 0]
        if len(sel) == 0:
            continue
        if str(s).startswith("seed") and len(sel) > 400:
            sel = np.sort(rng.choice(sel, 400, replace=False))
        xs.append(c[sel])
        ys.append(lab[sel].astype(float))
    x = common.featurize(np.concatenate(xs))
    y = np.concatenate(ys)
    return models.train_committor(x, y, dim=x.shape[1], steps=steps,
                                  seed=rng_seed), len(y)


def build_msm(store, n_clusters=150, lag=10, seed=0):
    feats = common.featurize(store["coords"])
    centers, labels = msmlib.cluster(feats, n_clusters=n_clusters, seed=seed)
    lags = [1, 2, 4, 6, 8, 10, 14, 20]
    its = msmlib.implied_timescales(None, lags, labels, store["lengths"])
    cm = msmlib.count_matrix(labels, store["lengths"], lag)
    keep = msmlib.largest_connected_set(cm)
    t, pi_sub = msmlib.reversible_mle(cm[np.ix_(keep, keep)])
    pi = np.zeros(labels.max() + 1)
    pi[keep] = pi_sub
    mask = np.isin(labels, keep)
    w = np.zeros(len(labels))
    w[mask] = msmlib.frame_weights(labels[mask], pi)
    # macrostate membership of the (connected) microstates, from the frames
    phi_all, _ = common.phi_psi(store["coords"])
    micro_state = np.full(labels.max() + 1, -1)
    for k in range(labels.max() + 1):
        sel = labels == k
        if not sel.any():
            continue
        c = common.which_core(phi_all[sel])
        vals, cnt = np.unique(c, return_counts=True)
        micro_state[k] = vals[cnt.argmax()]
    sub_a = np.where(micro_state[keep] == 0)[0]
    sub_b = np.where(micro_state[keep] == 1)[0]
    lag_ps = lag * 0.2
    mfpt_ab = msmlib.mfpt(t, pi_sub, sub_a, sub_b, lag_ps)
    mfpt_ba = msmlib.mfpt(t, pi_sub, sub_b, sub_a, lag_ps)
    return dict(feats=feats, labels=labels, centers=centers, pi=pi,
                weights=w, mask=mask, its=its, lags=lags, lag=lag, T=t,
                keep=keep, micro_state=micro_state,
                mfpt_ab=mfpt_ab, mfpt_ba=mfpt_ba)


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
    """Delta G between the two macrostates from a (phi, psi) landscape."""
    centers = 0.5 * (edges[1:] + edges[:-1])
    p = np.exp(-f / KT)
    p[~np.isfinite(p)] = 0.0
    core = common.which_core(centers)
    pa = p[core == 0].sum()
    pb = p[core == 1].sum()
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
    core = common.which_core(centers)
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
        print(f"    lag={lag * 0.2:5.1f}   " +
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
    )
    if ref is not None:
        out.update(dG_metad=float(dgr), barrier_metad=barr, fel_rmse=rmse)
    with open(os.path.join(common.RESULTS, f"analysis_{args.tag}.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    np.savez_compressed(
        os.path.join(common.RESULTS, f"analysis_{args.tag}.npz"),
        fel=fel, edges=edges, ref=ref if ref is not None else np.zeros(0),
        phi=phi, psi=psi, q=q, weights=w, tse=tse,
        its=m["its"], lags=np.array(m["lags"]),
    )
    print(f"[analysis] wrote results/analysis_{args.tag}.*")


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
