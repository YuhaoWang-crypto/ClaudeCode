"""
Two-stage / hybrid retrieval on a REAL 1,000,000-cell Census library.

Library  1,000,000 primary blood cells from every Census blood dataset except
         the 3 held-out labs, mapped onto the identity tree; at most 150k
         cells per dataset; sampled as runs of 5,000 consecutive cells
         (random single-cell access from S3 is ~50x slower than runs).
Queries  300 random cells from each held-out lab: COMBAT (10x 5' v1),
         SLE_1.2M (10x 3' v2), Wilk (Seq-Well).
Genes    2,000 highest-variance genes chosen on the NON-held-out part of
         the earlier 97k subset; log1p(count / raw_sum * 1e4), raw_sum =
         Census per-cell total over all genes.

Steps (each resumable, everything cached under data_cache/c1m/):
    python3 -m dcsc.real.census_1m select
    python3 -m dcsc.real.census_1m fetch      # ~1M cells x 2000 genes from S3
    python3 -m dcsc.real.census_1m prepare    # z-score, 32-program basis, dense float16 memmap
    python3 -m dcsc.real.census_1m evaluate   # -> dcsc/real/results_census_1m.json
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .census_fetch import CENSUS_VERSION, DATASETS, LEAF2LIN, map_cell_type

D = "data_cache/c1m"
HELD_OUT = {"ebc2e1ff-c8f9-466a-acf4-9d291afaf8b3": "COMBAT (10x 5' v1)",
            "218acb0f-9f2f-4f76-b90b-15a4b7c7f629": "SLE_1.2M (10x 3' v2)",
            "456e8b9b-f872-488b-871d-94534090a865": "Wilk (Seq-Well)"}
N_LIB, CAP, RUN, N_Q, N_HVG, K = 1_000_000, 150_000, 5_000, 300, 2000, 32
OBS_COLS = ["soma_joinid", "raw_sum", "dataset_id", "assay", "cell_type", "disease", "donor_id"]


# ---------------------------------------------------------------------------
def select(seed=0):
    rng = np.random.default_rng(seed)
    obs = pd.read_parquet("data_cache/census_blood_obs.parquet")
    obs["leaf"] = obs.cell_type.astype(str).map(map_cell_type)
    obs = obs[obs.leaf.notna()]
    lib = obs[~obs.dataset_id.isin(HELD_OUT)]
    runs = []
    for ds, g in lib.groupby("dataset_id", observed=True):
        ids = np.sort(g.soma_joinid.values)
        blocks = [ids[i:i + RUN] for i in range(0, len(ids), RUN)]
        order = rng.permutation(len(blocks))
        take, n = [], 0
        for b in order:
            if n >= CAP:
                break
            take.append(blocks[b])
            n += len(blocks[b])
        runs += take
    rng.shuffle(runs)
    ids = np.concatenate(runs)
    if len(ids) > N_LIB:      # trim whole runs, then the last partial run
        cum = np.cumsum([len(r) for r in runs])
        k = np.searchsorted(cum, N_LIB)
        ids = np.concatenate(runs[:k] + [runs[k][:N_LIB - (cum[k - 1] if k else 0)]])
    ids = np.sort(ids)
    q = []
    for ds in HELD_OUT:
        q.append(rng.choice(obs.soma_joinid.values[obs.dataset_id == ds], N_Q, replace=False))
    os.makedirs(D, exist_ok=True)
    np.save(f"{D}/lib_ids.npy", ids)
    np.save(f"{D}/query_ids.npy", np.sort(np.concatenate(q)))
    # HVGs from the non-held-out part of the 97k subset
    import anndata as ad
    a = ad.read_h5ad("data_cache/census_pbmc.h5ad")
    a = a[~a.obs.dataset_id.isin(HELD_OUT).values]
    X = a.X.tocsr().astype(np.float32)
    lib_sz = np.asarray(X.sum(1)).ravel()
    X = sp.diags(1e4 / np.maximum(lib_sz, 1)) @ X
    X.data = np.log1p(X.data)
    mu = np.asarray(X.mean(0)).ravel()
    var = np.asarray(X.multiply(X).mean(0)).ravel() - mu ** 2
    hv = np.argsort(-var)[:N_HVG]
    pd.DataFrame({"feature_id": a.var.feature_id.values[hv]}).to_parquet(f"{D}/hvg.parquet")
    n_ds = lib[lib.soma_joinid.isin(ids)].dataset_id.nunique()
    print(f"library {len(ids):,} cells from {n_ds} datasets; queries {N_Q * len(HELD_OUT)}")


# ---------------------------------------------------------------------------
def _fetch(census, ids, var_ids, order):
    import cellxgene_census as cc
    a = cc.get_anndata(census, "Homo sapiens", obs_coords=ids, var_coords=var_ids,
                       obs_column_names=OBS_COLS)
    # align columns to the HVG order by feature_id
    pos = pd.Series(np.arange(a.n_vars), index=a.var.feature_id.values)
    X = a.X.tocsr()[:, pos.loc[order].values].astype(np.float32)
    X = sp.diags(1e4 / a.obs.raw_sum.values.astype(np.float64)) @ X
    X.data = np.log1p(X.data)
    return X.astype(np.float32).tocsr(), a.obs.reset_index(drop=True)


def fetch(chunk=50_000):
    import cellxgene_census as cc
    ids = np.load(f"{D}/lib_ids.npy")
    hv = pd.read_parquet(f"{D}/hvg.parquet").feature_id.values
    with cc.open_soma(census_version=CENSUS_VERSION) as census:
        var = census["census_data"]["homo_sapiens"].ms["RNA"].var.read(
            column_names=["soma_joinid", "feature_id"]).concat().to_pandas()
        var_ids = np.sort(var.soma_joinid.values[var.feature_id.isin(hv)])
        if not os.path.exists(f"{D}/query.npz"):
            Xq, oq = _fetch(census, np.load(f"{D}/query_ids.npy"), var_ids, hv)
            sp.save_npz(f"{D}/query.npz", Xq)
            oq.to_parquet(f"{D}/query_obs.parquet")
            print("queries fetched", Xq.shape, flush=True)
        for i, s in enumerate(range(0, len(ids), chunk)):
            f = f"{D}/lib_{i:03d}.npz"
            if os.path.exists(f):
                continue
            t = time.time()
            for attempt in range(4):
                try:
                    X, o = _fetch(census, ids[s:s + chunk], var_ids, hv)
                    break
                except Exception as e:          # transient S3 errors
                    print("retry", attempt, e, flush=True)
                    time.sleep(2 ** attempt * 5)
            else:
                raise RuntimeError(f"chunk {i} failed")
            sp.save_npz(f, X)
            o.to_parquet(f"{D}/lib_{i:03d}_obs.parquet")
            print(f"chunk {i}: {X.shape[0]:,} cells, nnz {X.nnz:,}, {time.time() - t:.0f}s", flush=True)


# ---------------------------------------------------------------------------
def _lib_chunks():
    i = 0
    while os.path.exists(f"{D}/lib_{i:03d}.npz"):
        yield sp.load_npz(f"{D}/lib_{i:03d}.npz"), pd.read_parquet(f"{D}/lib_{i:03d}_obs.parquet")
        i += 1


def prepare(seed=0):
    """z-score + 32-program basis from a 50k-cell library sample; write the
    dense z-scored library as a float16 memmap (for raw-space tests only)."""
    rng = np.random.default_rng(seed)
    obs, n = [], 0
    s1 = np.zeros(N_HVG)
    s2 = np.zeros(N_HVG)
    samp = []
    for X, o in _lib_chunks():
        obs.append(o)
        s1 += np.asarray(X.sum(0)).ravel()
        s2 += np.asarray(X.multiply(X).sum(0)).ravel()
        n += X.shape[0]
        samp.append(X[rng.choice(X.shape[0], min(2500, X.shape[0]), replace=False)])
    mu = s1 / n
    sd = np.sqrt(s2 / n - mu ** 2) + 1e-3
    S = np.clip((sp.vstack(samp).toarray() - mu) / sd, -10, 10).astype(np.float32)
    _, _, vt = np.linalg.svd(S, full_matrices=False)
    np.savez(f"{D}/model.npz", mu=mu, sd=sd, W=vt[:K].T.astype(np.float32), sample=S)
    obs = pd.concat(obs, ignore_index=True)
    obs["leaf"] = obs.cell_type.astype(str).map(map_cell_type)
    obs.to_parquet(f"{D}/lib_obs_all.parquet")
    F = np.lib.format.open_memmap(f"{D}/lib_dense_f16.npy", mode="w+", dtype=np.float16, shape=(n, N_HVG))
    r = 0
    for X, _ in _lib_chunks():
        F[r:r + X.shape[0]] = np.clip((X.toarray() - mu) / sd, -10, 10).astype(np.float16)
        r += X.shape[0]
    F.flush()
    print(f"prepared {n:,} cells")


# ---------------------------------------------------------------------------
def _queries(mu, sd):
    Xq = sp.load_npz(f"{D}/query.npz")
    oq = pd.read_parquet(f"{D}/query_obs.parquet")
    oq["leaf"] = oq.cell_type.astype(str).map(map_cell_type)
    return np.clip((Xq.toarray() - mu) / sd, -10, 10).astype(np.float32), oq


def _chunks_dense(F, chunk=50_000):
    for s in range(0, len(F), chunk):
        yield s, np.asarray(F[s:s + chunk], dtype=np.float32)


def _time(fn, reps=10):
    fn()
    t = time.perf_counter()
    for _ in range(reps):
        fn()
    return (time.perf_counter() - t) / reps


def evaluate(out="dcsc/real/results_census_1m.json"):
    import faiss
    from ..twostage import HybridIndex
    m = np.load(f"{D}/model.npz")
    mu, sd, W, S = m["mu"], m["sd"], m["W"], m["sample"]
    F = np.load(f"{D}/lib_dense_f16.npy", mmap_mode="r")
    obs = pd.read_parquet(f"{D}/lib_obs_all.parquet")
    Fq, oq = _queries(mu, sd)
    N = len(F)
    leaf = obs.leaf.values
    sick = (obs.disease.values != "normal")
    q_leaf, q_sick = oq.leaf.values, (oq.disease.values != "normal")
    q_ds = oq.dataset_id.map(HELD_OUT).values
    has_leaf = pd.notna(q_leaf)
    print(f"library {N:,} x {F.shape[1]}, queries {len(Fq)} ({has_leaf.sum()} with a mapped cell type)")
    res = dict(N=N, n_datasets=int(obs.dataset_id.nunique()), n_donors=int(obs.donor_id.nunique()),
               assays=obs.assay.value_counts().to_dict(),
               leaf_mix=obs.leaf.value_counts(normalize=True).round(4).to_dict(), rows=[], timing=[])

    def bio(nbs, name, extra=None):
        nbs = np.asarray(nbs)
        r = dict(method=name,
                 same_leaf=float(np.mean([(leaf[n] == ql).mean() for n, ql, h in zip(nbs, q_leaf, has_leaf) if h])),
                 same_lineage=float(np.mean([np.mean([LEAF2LIN.get(x) == LEAF2LIN[ql] for x in leaf[n]])
                                             for n, ql, h in zip(nbs, q_leaf, has_leaf) if h])),
                 same_status=float(np.mean([(sick[n] == qs).mean() for n, qs in zip(nbs, q_sick)])),
                 n_source_datasets=float(np.mean([len(set(obs.dataset_id.values[n])) for n in nbs])))
        for lab in HELD_OUT.values():
            sel = (q_ds == lab) & has_leaf
            r[f"same_leaf | {lab}"] = float(np.mean([(leaf[n] == ql).mean() for n, ql in zip(nbs[sel], q_leaf[sel])]))
        r.update(extra or {})
        res["rows"].append(r)
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
        return nbs

    # ---- exact references (chunked brute force, multi-threaded BLAS) -------
    t0 = time.time()
    Zq = Fq @ W
    Zqn = Zq / np.linalg.norm(Zq, axis=1, keepdims=True)
    Fqn = Fq / np.linalg.norm(Fq, axis=1, keepdims=True)
    best_raw = [np.full((len(Fq), 10), -np.inf), np.zeros((len(Fq), 10), np.int64)]
    Z = np.zeros((N, K), np.float32)
    for s, X in _chunks_dense(F):
        Z[s:s + len(X)] = X @ W
        cos = (Fqn @ X.T) / (np.linalg.norm(X, axis=1) + 1e-12)
        _merge_topk(best_raw, cos, s)
    raw_nb = best_raw[1][np.arange(len(Fq))[:, None], np.argsort(-best_raw[0], 1)]
    print(f"exact raw scan over {N:,}: {time.time() - t0:.0f}s total", flush=True)
    bio(raw_nb, "exact raw-gene scan (2000 genes, float)", dict(bytes=N_HVG * 4))
    Zn = Z / (np.linalg.norm(Z, axis=1, keepdims=True) + 1e-12)
    prog_nb = np.array([np.argsort(-(Zn @ z))[:10] for z in Zqn])
    bio(prog_nb, "exact program scan (32 PCs, float)", dict(bytes=K * 4))

    # ---- indices: built with all cores, queried and timed on ONE core -------
    from threadpoolctl import threadpool_limits
    builds = {"IVF-PQ m=8 program": dict(pq_m=8, coarse="ivfpq", coarse_space="program", nlist=4096),
              "IVF-PQ m=32 raw": dict(pq_m=32, coarse="ivfpq", coarse_space="raw", nlist=4096)}
    ref = None
    for name, kw in builds.items():
        t0 = time.time()
        faiss.omp_set_num_threads(os.cpu_count())
        ix = HybridIndex(S, K=K, basis=W, add_init=False, **kw)
        for s, X in _chunks_dense(F):
            ix.add(X)
        len(ix)
        build = time.time() - t0
        faiss.omp_set_num_threads(1)
        lim = threadpool_limits(1)
        if ref is None:   # reference = exhaustive int8 program-code scan
            ref = np.array([ix.program_scan(q, 10) for q in Fq])
            t = _time(lambda: ix.program_scan(Fq[0], 10))
            bio(ref, "int8 program-code scan (exhaustive, reference)",
                dict(bytes=K, recall=1.0, ms=t * 1e3))
            res["timing"].append(dict(method="int8 program-code scan", ms=t * 1e3, N=N))
            # PQ alone on raw genes, no rescoring: needs the raw index; done below
        for nprobe in (8, 32, 128):
            for M in (200, 1000):
                got = np.array([_pad10(ix.search(q, 10, M=M, nprobe=nprobe)) for q in Fq])
                t = _time(lambda: ix.search(Fq[0], 10, M=M, nprobe=nprobe))
                rec = float(np.mean([len(set(a) & set(b)) / 10 for a, b in zip(got, ref)]))
                bio(got, f"{name} -> int8 program code, nprobe={nprobe}, M={M}",
                    dict(bytes=ix.bytes_per_item, recall=rec, ms=t * 1e3, build_s=round(build)))
        if "raw" in name:
            got = np.array([_pad10(ix.search(q, 10, M=10, nprobe=32, rescore="none")) for q in Fq])
            t = _time(lambda: ix.search(Fq[0], 10, M=10, nprobe=32, rescore="none"))
            rec = float(np.mean([len(set(a) & set(b)) / 10 for a, b in zip(got, ref)]))
            bio(got, "IVF-PQ m=32 raw only (no rescoring), nprobe=32",
                dict(bytes=ix.pq_m, recall=rec, ms=t * 1e3, build_s=round(build)))
        lim.unregister()
        del ix
    rnd = np.random.default_rng(0).integers(0, N, (len(Fq), 10))
    bio(rnd, "random library cells")
    with open(out, "w") as f:
        json.dump(res, f, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))


def _pad10(a):
    a = np.asarray(a)
    return np.concatenate([a, np.full(10 - len(a), a[0] if len(a) else 0)]) if len(a) < 10 else a[:10]


def _merge_topk(best, cos, offset, k=10):
    idx = np.argpartition(-cos, k - 1, axis=1)[:, :k]
    val = np.take_along_axis(cos, idx, 1)
    allv = np.hstack([best[0], val])
    alli = np.hstack([best[1], idx + offset])
    o = np.argsort(-allv, 1)[:, :k]
    best[0][:] = np.take_along_axis(allv, o, 1)
    best[1][:] = np.take_along_axis(alli, o, 1)


if __name__ == "__main__":
    {"select": select, "fetch": fetch, "prepare": prepare, "evaluate": evaluate}[sys.argv[1]]()
