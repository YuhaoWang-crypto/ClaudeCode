"""
Donor-level disease-state retrieval on the real 1M-cell Census library.

Single-cell neighbours are dominated by cell type (REAL_RESULTS §5.5: same
disease status only 0.54 vs 0.47 random). Here every DONOR becomes one
compact profile and we ask: do a donor's nearest donors — from OTHER labs —
share its disease status?

Donor profile variants (donors with >= 100 mapped cells):
  comp       cell-type composition (9 fractions, sqrt-transformed)
  program    mean int8 program code (32 dims) inside each of 6 major cell
             types -> 192 dims  (what the DCSC code gives you for free)
  comp+prog  both, each block scaled to unit total variance
  raw        mean z-scored expression of 2000 genes in each of 6 cell types
             (12,000 dims, the uncompressed baseline)
  ISG score  1-D: Kang IFN-β loadings x CD14-monocyte raw profile (marker
             baseline, not a retrieval method)

Lab (batch) handling:
  none          profiles as they are
  lab-centred   subtract the mean of the SAME lab's healthy donors (needs
                >= 3 healthy controls per lab; leave-one-out for healthy
                queries). Labs without controls are dropped from the library.
  label-free    subtract the mean of all the lab's donors (no labels needed)

Evaluation: leave-one-lab-out. Every library lab with >= 3 sick and >= 3
healthy donors is held out in turn; its donors query the donors of all other
labs (k = 10, cosine). Plus 3 external labs never in the 1M library
(COMBAT, SLE_1.2M, Wilk; from the 97k subset).
  AUROC     fraction of sick neighbours as a score for "query is sick"
  COVID@10  for COVID queries: share of COVID-19 neighbours (vs base rate)
  ceiling   same, but searching only the query's OWN lab (other donors)

    python3 -m dcsc.real.donor_search
"""
import json
import os

import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.metrics import roc_auc_score

from .census_1m import D, HELD_OUT
from .census_fetch import map_cell_type

LEAVES = ["CD14_Mono", "CD16_Mono", "cDC", "pDC", "CD4_T", "CD8_T", "NK", "B", "Plasma"]
MAJOR = ["CD14_Mono", "CD16_Mono", "CD4_T", "CD8_T", "NK", "B"]
MIN_CELLS, MIN_LEAF, K_NN = 100, 10, 10


def _int8(Z, scale):
    return np.clip(np.round(Z / scale), -127, 127) * scale


def _aggregate(Z, Fchunks, leaf, donor_idx, n_donor):
    """Per donor x leaf: cell counts, mean int8 program code, mean raw profile."""
    li = pd.Series(leaf).map({l: i for i, l in enumerate(LEAVES)}).fillna(-1).astype(int).values
    ok = (li >= 0) & (donor_idx >= 0)
    g = np.where(ok, donor_idx * len(LEAVES) + li, -1)
    G = n_donor * len(LEAVES)
    cnt = np.bincount(g[ok], minlength=G).astype(float)
    zs = np.zeros((G, Z.shape[1]))
    for j in range(Z.shape[1]):
        zs[:, j] = np.bincount(g[ok], weights=Z[ok, j], minlength=G)
    major = np.isin(li, [LEAVES.index(m) for m in MAJOR])
    gm = np.where(ok & major, donor_idx * len(MAJOR) + pd.Series(leaf).map(
        {l: i for i, l in enumerate(MAJOR)}).fillna(0).astype(int).values, -1)
    raw = None
    for s, X in Fchunks:
        gg = gm[s:s + len(X)]
        sel = gg >= 0
        A = sp.csr_matrix((np.ones(sel.sum()), (gg[sel], np.where(sel)[0])), shape=(n_donor * len(MAJOR), len(X)))
        part = A @ X
        raw = part if raw is None else raw + part
    return cnt.reshape(n_donor, len(LEAVES)), zs.reshape(n_donor, len(LEAVES), -1), raw


def build_profiles():
    cache = f"{D}/donor_profiles.npz"
    if os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        return {k: z[k] for k in z.files}
    m = np.load(f"{D}/model.npz")
    mu, sd, W, S = m["mu"], m["sd"], m["W"], m["sample"]
    scale = np.abs(S @ W).max(0) / 127
    F = np.load(f"{D}/lib_dense_f16.npy", mmap_mode="r")
    obs = pd.read_parquet(f"{D}/lib_obs_all.parquet")
    obs["donor"] = obs.dataset_id.astype(str) + "|" + obs.donor_id.astype(str)
    Z = np.zeros((len(F), W.shape[1]), np.float32)
    for s in range(0, len(F), 50_000):
        Z[s:s + 50_000] = np.asarray(F[s:s + 50_000], np.float32) @ W
    Zq = _int8(Z, scale)

    def chunks(Fm, n):
        for s in range(0, n, 50_000):
            yield s, np.asarray(Fm[s:s + 50_000], np.float32)

    lib = _profiles_from(obs, Zq, lambda: chunks(F, len(F)))
    # external labs from the 97k subset (never in the 1M library)
    import anndata as ad
    a = ad.read_h5ad("data_cache/census_pbmc.h5ad")
    a = a[a.obs.dataset_id.isin(HELD_OUT).values]
    hv = pd.read_parquet(f"{D}/hvg.parquet").feature_id.values
    X = a.X.tocsr().astype(np.float32)
    libsz = np.asarray(X.sum(1)).ravel()          # ~ raw_sum (genes in < 50 cells were dropped)
    pos = pd.Series(np.arange(a.n_vars), index=a.var.feature_id.values)
    X = sp.diags(1e4 / np.maximum(libsz, 1)) @ X[:, pos.reindex(hv).fillna(-1).astype(int).clip(0).values]
    X = X.toarray()
    X[:, pos.reindex(hv).isna().values] = 0
    Fe = np.clip((np.log1p(X) - mu) / sd, -10, 10).astype(np.float32)
    oe = a.obs.reset_index(drop=True)[["dataset_id", "donor_id", "disease", "cell_type"]].copy()
    oe["leaf"] = oe.cell_type.astype(str).map(map_cell_type)
    oe["donor"] = oe.dataset_id.astype(str) + "|" + oe.donor_id.astype(str)
    ext = _profiles_from(oe, _int8(Fe @ W, scale), lambda: iter([(0, Fe)]))
    ext["external"] = np.ones(len(ext["donor"]), bool)
    lib["external"] = np.zeros(len(lib["donor"]), bool)
    out = {k: np.concatenate([lib[k], ext[k]]) for k in lib}
    # ISG baseline: Kang IFN loadings on the HVGs, applied to the CD14-monocyte raw profile
    from .kang import ifn_loadings
    w = ifn_loadings()
    wi = np.where(np.abs(w.pooled) > 3 * w.pooled_se, w.pooled, 0)
    wv = pd.Series(wi, index=w.index).reindex(hv).fillna(0).values * sd   # back to z-units
    out["isg"] = out["raw"][:, MAJOR.index("CD14_Mono"), :] @ wv
    np.savez_compressed(cache, **out)
    return out


def _profiles_from(obs, Zq, chunk_fn):
    donors = obs.donor.value_counts()
    keep = donors.index[donors >= MIN_CELLS]
    didx = pd.Series(np.arange(len(keep)), index=keep)
    donor_idx = obs.donor.map(didx).fillna(-1).astype(int).values
    cnt, zs, raw = _aggregate(Zq, chunk_fn(), obs.leaf.values, donor_idx, len(keep))
    raw = np.asarray(raw.todense() if sp.issparse(raw) else raw).reshape(len(keep), len(MAJOR), -1)
    mcnt = cnt[:, [LEAVES.index(m) for m in MAJOR]]
    raw = raw / np.maximum(mcnt, 1)[:, :, None]
    prog = zs[:, [LEAVES.index(m) for m in MAJOR], :] / np.maximum(mcnt, 1)[:, :, None]
    present = mcnt >= MIN_LEAF
    meta = obs.groupby("donor", observed=True).agg(dataset=("dataset_id", "first"),
                                                     disease=("disease", lambda s: s.value_counts().index[0]))
    meta = meta.loc[keep]
    return dict(donor=np.array(keep), dataset=meta.dataset.astype(str).values,
                disease=meta.disease.astype(str).values, comp=cnt / cnt.sum(1, keepdims=True),
                prog=prog.astype(np.float32), raw=raw.astype(np.float32), present=present)


# ---------------------------------------------------------------------------
def features(P, kind):
    comp = np.sqrt(P["comp"])
    prog = np.where(P["present"][:, :, None], P["prog"], np.nan).reshape(len(comp), -1)
    raw = np.where(P["present"][:, :, None], P["raw"], np.nan).reshape(len(comp), -1)
    blocks = {"comp": [comp], "program": [prog], "comp+prog": [comp, prog], "raw": [raw]}[kind]
    return blocks


def center(blocks, P, mode):
    """Per-lab centring; NaN (missing cell type) -> 0 after centring."""
    ds, sick = P["dataset"], P["disease"] != "normal"
    out, keep = [], np.ones(len(ds), bool)
    for B in blocks:
        B = B.copy()
        gmean = np.nanmean(B, 0)
        if mode != "none":
            for d in np.unique(ds):
                m = ds == d
                ref = m & ~sick if mode == "lab-centred" else m
                if mode == "lab-centred" and ref.sum() < 3:
                    keep[m] = False
                    continue
                R = B[ref]
                mu = np.nanmean(R, 0)
                if mode == "lab-centred":
                    # leave-one-out for the controls themselves
                    n = np.sum(~np.isnan(R), 0)
                    tot = np.nansum(R, 0)
                    for i in np.where(ref)[0]:
                        xi = np.nan_to_num(B[i])
                        has = ~np.isnan(B[i])
                        B[i] = B[i] - np.where(n - has > 0, (tot - xi * has) / np.maximum(n - has, 1), mu)
                    B[m & ~ref] = B[m & ~ref] - mu
                else:
                    B[m] = B[m] - mu
        else:
            B = B - gmean
        B = np.nan_to_num(B)
        B /= np.sqrt((B ** 2).sum(1).mean()) + 1e-9          # each block: unit mean energy
        out.append(B)
    return np.hstack(out), keep


def knn_eval(X, P, keep, within=False, metric="cosine"):
    ds, dis = P["dataset"], P["disease"]
    sick = dis != "normal"
    ext = P["external"]
    Xn = X / (np.linalg.norm(X, axis=1, keepdims=True) + 1e-9)
    sq = (X ** 2).sum(1)
    rows = []
    labs = [d for d in np.unique(ds) if keep[ds == d].all()
            and (sick[ds == d]).sum() >= 3 and (~sick[ds == d]).sum() >= 3]
    for d in labs:
        q = np.where(ds == d)[0]
        lib = np.where((ds == d) if within else ((ds != d) & keep & ~ext))[0]
        score, covid_share, cscore = [], [], []
        for i in q:
            cand = lib[lib != i]
            if metric == "cosine":
                s = Xn[cand] @ Xn[i]
            else:                                   # euclidean: rank by -distance²
                s = -(sq[cand] - 2 * X[cand] @ X[i])
            nb = cand[np.argsort(-s)[:K_NN]]
            score.append(sick[nb].mean())
            cscore.append((dis[nb] == "COVID-19").mean())
            if dis[i] == "COVID-19":
                covid_share.append((dis[nb] == "COVID-19").mean())
        cv = np.isin(dis[q], ["COVID-19", "normal"])
        y = dis[q][cv] == "COVID-19"
        cauc = roc_auc_score(y, np.array(cscore)[cv]) if y.sum() >= 3 and (~y).sum() >= 3 else np.nan
        rows.append(dict(lab=d, external=bool(ext[q].all()), n=len(q), n_sick=int(sick[q].sum()),
                         covid_vs_healthy_auroc=cauc, n_covid=int(y.sum()),
                         auroc=roc_auc_score(sick[q], score), covid_at10=np.mean(covid_share) if covid_share else np.nan,
                         covid_base=float((dis[lib] == "COVID-19").mean())))
    return pd.DataFrame(rows)


def main(out="dcsc/real/results_donor_search.json"):
    P = build_profiles()
    sick = P["disease"] != "normal"
    print(f"{len(P['donor'])} donors ({P['external'].sum()} external), {sick.sum()} sick, "
          f"{(P['disease'] == 'COVID-19').sum()} COVID-19, {len(np.unique(P['dataset']))} labs")
    res = dict(n_donors=int(len(P["donor"])), n_external=int(P["external"].sum()), rows=[], per_lab=[])
    # 1-D marker baseline, scored directly (no retrieval)
    for mode in ["none", "label-free", "lab-centred"]:
        Xi, keep = center([P["isg"][:, None]], P, mode)
        x = Xi[:, 0]
        within = []
        for d in np.unique(P["dataset"]):
            m = (P["dataset"] == d) & keep
            if sick[m].sum() >= 3 and (~sick[m]).sum() >= 3:
                within.append(roc_auc_score(sick[m], x[m]))
        lib = keep & ~P["external"]
        ext = keep & P["external"]
        res["rows"].append(dict(profile="ISG score (1-D marker, no retrieval)", batch=mode,
                                within_lab_ceiling=float(np.median(within)),
                                pooled_auroc_library=float(roc_auc_score(sick[lib], x[lib])),
                                pooled_auroc_external=float(roc_auc_score(sick[ext], x[ext]))))
        print(res["rows"][-1], flush=True)
    for kind in ["comp", "program", "comp+prog", "raw", "lab-only control"]:
        for mode in ["none", "label-free", "lab-centred"]:
            if kind == "lab-only control":
                if mode != "none":
                    continue
                # every donor = its lab's mean program profile (+ tiny jitter):
                # retrieval can only match LABS, so this measures lab confounding
                X0, keep = center(features(P, "comp+prog"), P, "none")
                X = np.zeros_like(X0)
                for d in np.unique(P["dataset"]):
                    m = P["dataset"] == d
                    X[m] = X0[m].mean(0)
                X += np.random.default_rng(0).normal(0, 1e-6, X.shape)
            else:
                X, keep = center(features(P, kind), P, mode)
            for metric in ["cosine", "euclidean"]:
                t = knn_eval(X, P, keep, metric=metric)
                ceil = knn_eval(X, P, keep, within=True, metric=metric)
                r = dict(profile=kind, batch=mode, metric=metric, dims=X.shape[1],
                         auroc_library_labs=float(t[~t.external].auroc.median()),
                         auroc_external_labs=float(t[t.external].auroc.median()) if t.external.any() else None,
                         covid_at10=float(t.covid_at10.mean()), covid_base=float(t.covid_base.mean()),
                         within_lab_ceiling=float(ceil.auroc.median()), n_query_labs=int(len(t)),
                         covid_auroc_library=float(t[~t.external].covid_vs_healthy_auroc.median()),
                         covid_auroc_external=float(t[t.external].covid_vs_healthy_auroc.median()),
                         n_covid_labs=int(t.covid_vs_healthy_auroc.notna().sum()),
                         covid_within_lab_ceiling=float(ceil.covid_vs_healthy_auroc.median()))
                for lab, nm in HELD_OUT.items():
                    x = t[t.lab == lab]
                    r[f"auroc | {nm}"] = float(x.auroc.iloc[0]) if len(x) else None
                res["rows"].append(r)
                for _, x in t.iterrows():
                    res["per_lab"].append(dict(profile=kind, batch=mode, metric=metric, **x.to_dict()))
                print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
    with open(out, "w") as f:
        json.dump(res, f, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))


if __name__ == "__main__":
    main()
