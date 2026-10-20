"""
Real test of Δ-state transfer across identities and dose read-out on
sci-Plex3 (Srivatsan et al., Science 2020; scPerturb release): 188 drugs ×
4 doses (10 nM–10 µM) × 3 cancer cell lines (A549, MCF7, K562), 24 h.

  library  Δ of each drug at 10 µM in line A, replicate 1
  queries  (a) same line, replicate 2           — independent cells, baseline
           (b) another line, replicate 2         — cross-identity transfer
  scored   exact drug, same annotated target, same pathway (top-1 / top-5)
  dose     NNLS coefficient of the drug's own library entry fitted to the
           query at 10 nM … 10 µM: does it rise with the true dose?

Only "active" drugs are scored: replicate-1 and replicate-2 Δ at 10 µM agree
(cosine > 0.3 on the reproducibly-moved genes) in every line involved; an
inactive drug has no state change to retrieve.

    python3 -m dcsc.real.sciplex
"""
import json
import os

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.optimize import nnls
from scipy.stats import spearmanr

from .. import delta as DL
from ..index import mod_distance

PATH = "data_cache/SrivatsanTrapnell2020_sciplex3.h5ad"
CACHE = "data_cache/sciplex_pseudobulk.npz"
N_GENES = 6000
N_HVG = 2000
K = 32
ACTIVE_COS = 0.3


def pseudobulks():
    if os.path.exists(CACHE):
        z = np.load(CACHE, allow_pickle=True)
        return {k: z[k] for k in z.files}
    a = ad.read_h5ad(PATH, backed="r")
    o = a.obs
    o = o[o.time.astype(float) == 24.0]
    rows = np.where(a.obs_names.isin(o.index))[0]
    key = (o.cell_line.astype(str) + "|" + o.perturbation.astype(str) + "|" +
           o.dose_value.astype(float).astype(int).astype(str) + "|" + o.replicate.astype(str)).values
    groups, gid = np.unique(key, return_inverse=True)
    # gene filter from a sample of cells
    samp = np.sort(np.random.default_rng(0).choice(rows, 50_000, replace=False))
    Xs = a.X[samp]
    det = np.asarray((Xs > 0).mean(0)).ravel()
    genes = np.argsort(-det)[:N_GENES]
    S = np.zeros((len(groups), N_GENES))
    Q = np.zeros((len(groups), N_GENES))
    n = np.bincount(gid, minlength=len(groups)).astype(float)
    for s in range(0, len(rows), 40_000):
        r = rows[s:s + 40_000]
        X = a.X[r]
        lib = np.asarray(X.sum(1)).ravel()
        X = sp.diags(1e4 / np.maximum(lib, 1)) @ X[:, genes]
        X.data = np.log1p(X.data)
        G = sp.csr_matrix((np.ones(len(r)), (gid[s:s + 40_000], np.arange(len(r)))),
                          shape=(len(groups), len(r)))
        S += (G @ X).toarray()
        Q += (G @ X.multiply(X)).toarray()
        print(f"  {s + len(r):,}/{len(rows):,}")
    meta = o.drop_duplicates("perturbation").set_index("perturbation")[["target", "pathway", "pathway_level_1"]]
    out = dict(groups=groups, n=n, mean=S / n[:, None], sq=Q / n[:, None],
               genes=np.array(a.var_names[genes]),
               meta_drug=np.array(meta.index), meta_target=meta.target.astype(str).values,
               meta_pathway=meta.pathway.astype(str).values)
    np.savez_compressed(CACHE, **out)
    return out


def main(out_json="dcsc/real/results_sciplex.json"):
    pb = pseudobulks()
    g = pd.DataFrame([x.split("|") for x in pb["groups"]], columns=["line", "drug", "dose", "rep"])
    g["dose"] = g.dose.astype(int)
    g["n"] = pb["n"]
    target = dict(zip(pb["meta_drug"], pb["meta_target"]))
    pathway = dict(zip(pb["meta_drug"], pb["meta_pathway"]))
    lines = sorted(g.line.unique())
    mean, sq = pb["mean"], pb["sq"]

    def ctrl(line, rep=None):
        m = (g.line == line) & (g.drug == "control") & ((g.rep == rep) if rep else True)
        w = g.n[m].values
        mu = (mean[m.values] * w[:, None]).sum(0) / w.sum()
        var = (sq[m.values] * w[:, None]).sum(0) / w.sum() - mu ** 2
        return mu, np.sqrt(var + 1e-3)

    def delta(line, drug, dose, rep):
        m = ((g.line == line) & (g.drug == drug) & (g.dose == dose) & (g.rep == rep)).values
        if m.sum() != 1 or g.n[m].values[0] < 15:
            return None
        mu, sd = ctrl(line)
        return (mean[m][0] - mu) / sd

    drugs = sorted(d for d in g.drug.unique() if d != "control")
    D1 = {l: {} for l in lines}
    D2 = {l: {} for l in lines}
    for l in lines:
        for d in drugs:
            a, b = delta(l, d, 10000, "rep1"), delta(l, d, 10000, "rep2")
            if a is not None and b is not None:
                D1[l][d], D2[l][d] = a, b
    # HVG = genes with the most REPRODUCIBLE drug-induced change (rep1·rep2 is
    # unbiased for signal energy, noise cancels), then active = replicate
    # cosine of the 10 µM Δ on those genes > 0.3
    rep_energy = np.mean([D1[l][d] * D2[l][d] for l in lines for d in D1[l]], 0)
    hv = np.argsort(-rep_energy)[:N_HVG]
    def rcos(l, d):
        a, b = D1[l][d][hv], D2[l][d][hv]
        return a @ b / (np.linalg.norm(a) * np.linalg.norm(b))
    active = {l: {d for d in D1[l] if rcos(l, d) > ACTIVE_COS} for l in lines}
    print({l: len(v) for l, v in active.items()})

    results = []
    for la in lines:
        lib_drugs = sorted(active[la])
        L = np.array([D1[la][d][hv] for d in lib_drugs])
        Ln = L / np.linalg.norm(L, axis=1, keepdims=True)
        _, s, vt = np.linalg.svd(L, full_matrices=False)
        k = min(K, len(L))
        W = vt[:k].T
        Lz = L @ W
        DL.TAU_DELTA = float(np.quantile(np.abs(Lz), 0.8))
        pad = lambda Z: np.hstack([Z, np.zeros((len(Z), 32 - Z.shape[1]))])
        sig = DL.signature(pad(Lz))
        for lb in lines:
            qs = [d for d in lib_drugs if d in active[lb] and d in D2[lb]]
            sc = {m: dict(drug1=0, drug5=0, target1=0, pathway1=0) for m in ["full cosine", "code (8 B)"]}
            for d in qs:
                q = D2[lb][d][hv]
                rank_full = np.argsort(-(Ln @ q))
                rank_code = np.argsort(mod_distance(sig, DL.signature(pad((q @ W)[None]))[0]) -
                                       1e-3 * (Ln @ q), kind="stable")   # Hamming, cosine tie-break
                for m, rk in [("full cosine", rank_full), ("code (8 B)", rank_code)]:
                    top = [lib_drugs[i] for i in rk[:5]]
                    sc[m]["drug1"] += top[0] == d
                    sc[m]["drug5"] += d in top
                    sc[m]["target1"] += target[top[0]] == target[d]
                    sc[m]["pathway1"] += pathway[top[0]] == pathway[d]
            nq = max(len(qs), 1)
            # chance level of the target / pathway match among library drugs
            ch_t = np.mean([np.mean([target[x] == target[d] for x in lib_drugs]) for d in qs]) if qs else 0
            for m, v in sc.items():
                results.append(dict(library=la, query=lb, method=m, n_queries=len(qs), n_library=len(lib_drugs),
                                    **{k: x / nq for k, x in v.items()}, chance_drug1=1 / len(lib_drugs),
                                    chance_target1=float(ch_t)))

    # dose: fit the drug's own 10 µM library entry (rep1) to rep2 queries at every dose
    dose_rows = []
    for la in lines:
        for lb in lines:
            rhos, coefs = [], {10: [], 100: [], 1000: [], 10000: []}
            for d in active[la] & active[lb]:
                lib = D1[la][d][hv]
                cs = []
                for dose in (10, 100, 1000, 10000):
                    q = delta(lb, d, dose, "rep2")
                    if q is None:
                        break
                    c = nnls(lib[:, None], q[hv])[0][0]
                    cs.append(c)
                    coefs[dose].append(c)
                if len(cs) == 4:
                    if np.ptp(cs) > 0:
                        rhos.append(spearmanr([10, 100, 1000, 10000], cs)[0])
                    else:
                        rhos.append(0.0)   # flat response counts as not monotone
            dose_rows.append(dict(library=la, query=lb, n_drugs=len(rhos),
                                  median_spearman=float(np.nanmedian(rhos)) if rhos else None,
                                  share_monotone_up=float(np.mean(np.array(rhos) > 0.7)) if rhos else None,
                                  median_coef={k: float(np.median(v)) for k, v in coefs.items() if v}))
    out = dict(active=dict((l, len(v)) for l, v in active.items()), retrieval=results, dose=dose_rows)
    with open(out_json, "w") as f:
        json.dump(out, f, indent=1)
    df = pd.DataFrame(results)
    df["cross"] = df["library"] != df["query"]
    print(df.round(2).to_string())
    print(df.groupby(["cross", "method"])[["drug1", "drug5", "target1", "pathway1", "chance_drug1",
                                             "chance_target1"]].mean().round(3))
    print(pd.DataFrame(dose_rows).to_string())


if __name__ == "__main__":
    main()
