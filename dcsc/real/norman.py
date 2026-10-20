"""
Real test of brain-hole 2 (Δ-state + residual pursuit) and brain-hole 3
(partial-view query) on Norman et al., Science 2019 (CRISPRa in K562,
scPerturb release). The ground truth is real: each double perturbation's two
genes are known, and both single perturbations were measured separately.

  library  = Δ (perturbed − control pseudobulk) of every single perturbation
  query    = Δ of a double perturbation A+B whose singles are in the library
  question = can residual pursuit name {A, B} from the library alone?
  bonus    = does the unexplained-residual flag track real non-additivity
             (genetic interaction), measured as |Δ_AB − Δ_A − Δ_B| / |Δ_AB| ?

Program space: the first K right-singular vectors of the library Δ matrix
(data-driven perturbation-response "programs", the PROGENy idea), so the
code is 2 bits × K and needs no curated gene sets.

Search designs compared:
  full      residual pursuit on the 2000-gene Δ (no compression)
  code      residual pursuit on the K=32 ternary code only
  two-stage code Hamming shortlist of M entries -> full-gene pursuit inside it
and recovery is scored exactly, up to near-duplicates (cos > 0.9), and up to
functional clusters of the library (average linkage, cos > 0.5).

    python3 -m dcsc.real.norman
"""
import json
import os

import anndata as ad
import numpy as np
import scipy.sparse as sp
from scipy.stats import spearmanr

from .. import codec as C
from .. import delta as DL

PATH = "data_cache/NormanWeissman2019_filtered.h5ad"
K = 32
N_HVG = 2000


def load():
    a = ad.read_h5ad(PATH)
    X = a.X.tocsc() if not sp.isspmatrix_csc(a.X) else a.X
    lib = np.asarray(X.sum(1)).ravel()
    mean = np.asarray(X.sum(0)).ravel() / X.shape[0]
    keep = np.where(mean > 0.05)[0]
    X = X[:, keep].tocsr()
    X = sp.diags(1e4 / lib) @ X
    X.data = np.log1p(X.data)
    genes = a.var_names[keep]
    pert = a.obs.perturbation.astype(str).values
    gem = a.obs.gemgroup.values
    return X.astype(np.float32), genes, pert, gem


def pseudobulk(X, idx):
    return np.asarray(X[idx].mean(0)).ravel()


def compute_deltas(cache="data_cache/norman_deltas.npz"):
    """All Δ vectors (z-units, HVG) the analyses need; cached after first run."""
    if os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        return {k: z[k] for k in z.files}
    rng = np.random.default_rng(0)
    X, genes, pert, gem = load()
    ctrl = np.where(pert == "control")[0]
    groups = {p: np.where(pert == p)[0] for p in np.unique(pert) if p != "control"}
    singles = sorted(p for p in groups if "_" not in p)
    doubles = sorted(p for p in groups if "_" in p and all(g in groups for g in p.split("_")))
    mu_c = pseudobulk(X, ctrl)
    Xc = X[ctrl]
    var_c = np.asarray(Xc.multiply(Xc).mean(0)).ravel() - mu_c ** 2
    D_single = np.array([pseudobulk(X, groups[p]) - mu_c for p in singles])
    score = (D_single ** 2).mean(0) / (var_c + 1e-3)
    hv = np.argsort(-score)[:N_HVG]
    sd = np.sqrt(var_c[hv] + 1e-3)
    D_double, H1, H2, S1, S2 = [], [], [], [], []
    for p in doubles:
        idx = groups[p]
        D_double.append((pseudobulk(X, idx) - mu_c)[hv] / sd)
        h = rng.permutation(idx)
        H1.append((pseudobulk(X, h[: len(h) // 2]) - mu_c)[hv] / sd)
        H2.append((pseudobulk(X, h[len(h) // 2:]) - mu_c)[hv] / sd)
    for p in singles:   # split halves of singles: library half vs held-out half
        h = rng.permutation(groups[p])
        S1.append((pseudobulk(X, h[: len(h) // 2]) - mu_c)[hv] / sd)
        S2.append((pseudobulk(X, h[len(h) // 2:]) - mu_c)[hv] / sd)
    out = dict(singles=np.array(singles), doubles=np.array(doubles), genes=np.array(genes[hv]),
               Ls=D_single[:, hv] / sd, Dd=np.array(D_double), H1=np.array(H1), H2=np.array(H2),
               S1=np.array(S1), S2=np.array(S2),
               n_single=np.array([len(groups[p]) for p in singles]),
               n_double=np.array([len(groups[p]) for p in doubles]))
    np.savez_compressed(cache, **out)
    return out


def _pad(Z):
    out = np.zeros((len(Z), 32))
    out[:, :Z.shape[1]] = Z
    return out


def _clusters(Ln, thr=0.5):
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform
    Dm = np.clip(1 - Ln @ Ln.T, 0, 2)
    np.fill_diagonal(Dm, 0)
    return fcluster(linkage(squareform(Dm, checks=False), "average"), 1 - thr, "distance")


def _match(got, truth, rel):
    """got/truth: sets of library indices; rel: boolean relation matrix."""
    return len(got) == len(truth) and all(any(rel[g, t] for g in got) for t in truth)


def _two_stage(d, Ls, code_L, code_q, M):
    from ..index import mod_distance
    ham = mod_distance(code_L, code_q)
    short = np.argsort(ham, kind="stable")[:M]
    r = DL.residual_pursuit(d, Ls[short], None)
    return {int(short[i]) for i in r["stimuli"]}, r, set(short.tolist())


def _design(W, prior_sd, n_cells, k):
    """Greedy A-optimal gene panel for the K program coordinates."""
    cov = np.diag(prior_sd ** 2)
    chosen = []
    for _ in range(k):
        best, bv, bc = None, np.inf, None
        for g in range(W.shape[0]):
            if g in chosen:
                continue
            w = W[g]
            cw = cov @ w
            new = cov - np.outer(cw, cw) / (1.0 / n_cells + w @ cw)
            v = np.trace(new)
            if v < bv:
                best, bv, bc = g, v, new
        chosen.append(best)
        cov = bc
    return np.array(chosen)


def main(out_json="dcsc/real/results_norman.json"):
    d = compute_deltas()
    singles, doubles = list(d["singles"]), list(d["doubles"])
    Ls, Dd = d["Ls"], d["Dd"]
    Ln = Ls / np.linalg.norm(Ls, axis=1, keepdims=True)
    dup = (Ln @ Ln.T) > 0.9
    cl = _clusters(Ln)
    same_cl = cl[:, None] == cl[None, :]
    noise_s = np.linalg.norm(d["S1"] - d["S2"], axis=1) / 2 * np.sqrt(2)  # noise of a full-n Δ
    strength = np.linalg.norm(Ls, axis=1)

    _, s, vt = np.linalg.svd(Ls, full_matrices=False)
    W = vt[:K].T
    Lz = Ls @ W
    DL.TAU_DELTA = float(np.quantile(np.abs(Lz), 0.8))
    code_L = DL.signature(_pad(Lz))

    designs = ["full 2000 genes", "code only (K=32, 8 B)", "two-stage M=8", "two-stage M=16",
               "two-stage M=32", "top-2 cosine (baseline)"]
    score = {k: dict(exact=0, dup=0, cluster=0, one=0) for k in designs}
    shortlist_recall = {8: 0, 16: 0, 32: 0}
    rows = []
    for i, p in enumerate(doubles):
        a, b = p.split("_")
        t = {singles.index(a), singles.index(b)}
        q = Dd[i]
        zq = q @ W
        code_q = DL.signature(_pad(zq[None]))[0]
        res = {}
        r_full = DL.residual_pursuit(q, Ls, None)
        res["full 2000 genes"] = set(r_full["stimuli"])
        res["code only (K=32, 8 B)"] = set(DL.residual_pursuit(_pad(zq[None])[0], _pad(Lz), code_L,
                                                               prefilter=16)["stimuli"])
        for M in (8, 16, 32):
            got, _, short = _two_stage(q, Ls, code_L, code_q, M)
            res[f"two-stage M={M}"] = got
            shortlist_recall[M] += t <= short
        res["top-2 cosine (baseline)"] = set(np.argsort(-(Ln @ q))[:2].tolist())
        for k, got in res.items():
            score[k]["exact"] += got == t
            score[k]["dup"] += _match(got, t, dup | np.eye(len(Ls), dtype=bool))
            score[k]["cluster"] += _match(got, t, same_cl)
            score[k]["one"] += len(got & t) > 0
        ia, ib = sorted(t)
        nonadd = np.linalg.norm(q - Ls[ia] - Ls[ib]) / np.linalg.norm(q)
        noise_q = np.linalg.norm(d["H1"][i] - d["H2"][i]) / 2 / np.linalg.norm(q)
        weak = min(t, key=lambda j: strength[j] / noise_s[j])
        rows.append(dict(double=p, nonadd=float(nonadd), noise=float(noise_q),
                         detectable=bool(strength[weak] > 3 * noise_s[weak]),
                         unexplained=float(r_full["unexplained"]),
                         full_ok=bool(_match(res["full 2000 genes"], t, dup | np.eye(len(Ls), dtype=bool))),
                         got=[singles[j] for j in r_full["stimuli"]],
                         dose=[round(float(x), 2) for x in r_full["dose"]]))
    n = len(doubles)
    table = {k: {m: v / n for m, v in sc.items()} for k, sc in score.items()}
    det = np.array([r["detectable"] for r in rows])
    ok = np.array([r["full_ok"] for r in rows])
    na = np.array([r["nonadd"] for r in rows])
    un = np.array([r["unexplained"] for r in rows])
    rho, pv = spearmanr(na, un)
    noise = np.array([r["noise"] for r in rows])
    # excess residual beyond the query's own sampling noise (split-half estimate)
    excess = np.sqrt(np.maximum(un ** 2 - noise ** 2, 0))
    flagged = excess > 0.35
    hi = na > np.quantile(na, 2 / 3)
    from sklearn.metrics import roc_auc_score
    auc_raw, auc_excess = roc_auc_score(hi, un), roc_auc_score(hi, excess)

    # partial view: query only k genes; library restricted to the same genes
    n_med = int(np.median(d["n_double"]))
    panels = {"12 genes (A-optimal)": _design(W, Lz.std(0), n_med, 12),
              "12 genes (top-variance)": np.argsort(-(Ls ** 2).mean(0))[:12],
              "50 genes (A-optimal)": _design(W, Lz.std(0), n_med, 50),
              "200 genes (A-optimal)": _design(W, Lz.std(0), n_med, 200),
              "all 2000 genes": np.arange(Ls.shape[1])}
    partial = []
    for name, g in panels.items():
        sc = dict(exact=0, dup=0, cluster=0, one=0)
        for i, p in enumerate(doubles):
            a, b = p.split("_")
            t = {singles.index(a), singles.index(b)}
            got = set(DL.residual_pursuit(Dd[i][g], Ls[:, g], None)["stimuli"])
            sc["exact"] += got == t
            sc["dup"] += _match(got, t, dup | np.eye(len(Ls), dtype=bool))
            sc["cluster"] += _match(got, t, same_cl)
            sc["one"] += len(got & t) > 0
        partial.append(dict(panel=name, **{k: v / n for k, v in sc.items()}))

    summary = dict(
        n_singles=len(singles), n_doubles=n, n_clusters=int(cl.max()),
        near_duplicate_singles=int((dup.sum(1) > 1).sum()),
        var_explained_K32=float((s[:K] ** 2).sum() / (s ** 2).sum()),
        single_effect_over_noise_quartiles=np.quantile(strength / noise_s, [.25, .5, .75]).round(2).tolist(),
        detectable_doubles=int(det.sum()),
        search=table,
        shortlist_recall={M: v / n for M, v in shortlist_recall.items()},
        full_on_detectable=float(ok[det].mean()), full_on_undetectable=float(ok[~det].mean()),
        nonadditivity=dict(spearman_rho=float(rho), p=float(pv),
                           auroc_top_third_raw=float(auc_raw), auroc_top_third_excess=float(auc_excess),
                           share_flagged=float(flagged.mean()),
                           flagged_among_top_third_nonadditive=float(flagged[hi].mean()),
                           flagged_among_rest=float(flagged[~hi].mean()),
                           median_split_half_noise=float(np.median([r["noise"] for r in rows]))),
        partial_view=partial,
    )
    with open(out_json, "w") as f:
        json.dump(dict(summary=summary, per_double=rows), f, indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
