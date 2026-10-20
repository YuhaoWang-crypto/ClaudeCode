"""
Real test of brain-holes 1 and 3 on a CELLxGENE Census subset
(dcsc/real/census_fetch.py: 8 blood datasets, 5 assays, healthy vs COVID-19 /
SLE donors). Every number is leave-one-dataset-out (LODO): the held-out lab /
assay contributes nothing to centroids, technical axes or panel design.

A. Identity tree   hierarchical nearest-centroid (lineage -> leaf) on genes
                   chosen from the training datasets; truth = author labels.
B. IFN program     ONE real calibrated program: Kang 2018 in-vitro IFN-β
                   (another lab, another chemistry). Estimated three ways:
                     naive  univariate projection on the top-50 ISGs
                     solo   weighted regression on all genes, IFN column only
                     joint  IFN column + 6 technical "Env" axes learned from
                            how HEALTHY donors differ between training datasets
                   Truth proxy: COVID-19 and SLE carry a type-I IFN signature
                   (well established); healthy donors do not.
                     within-dataset AUROC  disease vs healthy (ceiling)
                     pooled AUROC          all datasets mixed: batch hurts here
                     healthy eta²          share of healthy-donor variance
                                           explained by dataset (leakage)
C. Partial views   clinical 6-gene IFN score genes, A-optimal 12 genes, and an
                   8-marker "flow panel" for identity with abstention.

Not testable here: the Memory block (needs chromatin / re-challenge data).

    python3 -m dcsc.real.census
"""
import json

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.metrics import roc_auc_score

from .census_fetch import LEAF2LIN, TREE
from .kang import ifn_loadings

PATH = "data_cache/census_pbmc.h5ad"
N_ENV = 6
TAU = 0.25          # IFN+ bit: a quarter of a full in-vitro IFN-β stimulation
IFN6 = ["IFI27", "IFI44L", "IFIT1", "ISG15", "RSAD2", "SIGLEC1"]   # Rice et al. 2013
FLOW8 = ["CD3E", "CD4", "CD8A", "CD14", "FCGR3A", "CD19", "NCAM1", "HLA-DRA"]
LEAVES = [l for ls in TREE.values() for l in ls]


def load():
    a = ad.read_h5ad(PATH)
    W = ifn_loadings()
    keep = a.var.feature_id.isin(W.index).values
    a = a[:, keep].copy()
    X = a.X.tocsr().astype(np.float32)
    lib = np.asarray(X.sum(1)).ravel()
    X = sp.diags(1e4 / np.maximum(lib, 1)) @ X
    X.data = np.log1p(X.data)
    w = W.loc[a.var.feature_id.values]
    obs = a.obs.reset_index(drop=True)
    obs["sick"] = (obs.disease != "normal").astype(int)
    return X.tocsr(), obs, np.array(a.var_names), w


def group_mean(X, mask):
    return np.asarray(X[np.where(mask)[0]].mean(0)).ravel()


def corr_argmax(Xd, C):
    Xc = Xd - Xd.mean(1, keepdims=True)
    Cc = C - C.mean(1, keepdims=True)
    Xc /= np.linalg.norm(Xc, axis=1, keepdims=True) + 1e-9
    Cc /= np.linalg.norm(Cc, axis=1, keepdims=True) + 1e-9
    return np.argmax(Xc @ Cc.T, 1)


def posterior_matrix(W, D, prior_sd):
    F = W.T @ (W / D[:, None])
    return np.linalg.solve(F + np.diag(1 / prior_sd ** 2), (W / D[:, None]).T)


def a_optimal(W, D, prior_sd, k, n_cells, target=0, cand=None):
    cov = np.diag(prior_sd ** 2)
    cand = np.arange(W.shape[0]) if cand is None else cand
    chosen = []
    for _ in range(k):
        best, bv, bc = None, np.inf, None
        for g in cand:
            if g in chosen:
                continue
            w = W[g]
            cw = cov @ w
            new = cov - np.outer(cw, cw) / (D[g] / n_cells + w @ cw)
            if new[target, target] < bv:
                best, bv, bc = g, new[target, target], new
        chosen.append(best)
        cov = bc
    return np.array(chosen)


def lodo(X, obs, genes, wdf):
    datasets = obs.dataset.unique()
    gidx = {g: i for i, g in enumerate(genes)}
    # keep only loadings that clear 3 s.e. across Kang donors (as in codec.fit):
    # thousands of near-zero noisy loadings otherwise dilute the estimate
    w_ifn = np.where(np.abs(wdf.pooled.values) > 3 * wdf.pooled_se.values, wdf.pooled.values, 0.0)
    print(f"IFN loadings kept: {(w_ifn != 0).sum()} of {len(w_ifn)} genes")
    cell_rows, id_rows, flow_rows = [], [], []
    for test in datasets:
        tr, te = (obs.dataset != test).values, (obs.dataset == test).values
        healthy = tr & (obs.sick == 0).values
        # ---- reference centroids: mean of per-dataset healthy means --------
        cent, dsoff = {}, []
        for leaf in LEAVES:
            per = [group_mean(X, healthy & (obs.dataset == d).values & (obs.leaf == leaf).values)
                   for d in obs.dataset[tr].unique()
                   if (healthy & (obs.dataset == d).values & (obs.leaf == leaf).values).sum() >= 20]
            if per:
                cent[leaf] = np.mean(per, 0)
                dsoff += [p - cent[leaf] for p in per]
        leaves = list(cent)
        Cm = np.array([cent[l] for l in leaves])
        # noise: pooled within-(dataset, leaf) variance of healthy training cells
        sq = np.zeros(X.shape[1])
        mu2 = np.zeros(X.shape[1])
        n = 0
        for d in obs.dataset[tr].unique():
            for leaf in leaves:
                m = healthy & (obs.dataset == d).values & (obs.leaf == leaf).values
                if m.sum() >= 20:
                    Xi = X[np.where(m)[0]]
                    mu = np.asarray(Xi.mean(0)).ravel()
                    sq += np.asarray(Xi.multiply(Xi).sum(0)).ravel() - m.sum() * mu ** 2
                    n += m.sum()
        D = sq / n + 0.01
        # ---- A. identity ------------------------------------------------------
        mk = np.argsort(-(Cm.var(0) / D))[:800]
        lin_names = [l for l in TREE if any(x in leaves for x in TREE[l])]
        Lm = np.array([np.mean([cent[x] for x in TREE[l] if x in cent], 0) for l in lin_names])
        Xt = X[np.where(te)[0]]
        Xd = Xt[:, mk].toarray()
        lin_hat = np.array(lin_names)[corr_argmax(Xd, Lm[:, mk])]
        leaf_hat = np.empty(len(Xd), object)
        for l in lin_names:
            sel = lin_hat == l
            cand = [x for x in TREE[l] if x in cent]
            if sel.any():
                leaf_hat[sel] = np.array(cand)[corr_argmax(Xd[sel], np.array([cent[x] for x in cand])[:, mk])]
        ot = obs[te].reset_index(drop=True)
        id_rows.append(dict(dataset=test, n=int(te.sum()),
                            lineage_acc=float(np.mean(lin_hat == ot.lineage.values)),
                            leaf_acc=float(np.mean(leaf_hat == ot.leaf.values))))
        # ---- B. IFN program --------------------------------------------------
        U, s, Vt = np.linalg.svd(np.array(dsoff), full_matrices=False)
        W_env = (Vt[:N_ENV] * (s[:N_ENV, None] / np.sqrt(len(dsoff)))).T     # prior sd 1
        W_joint = np.column_stack([w_ifn, W_env])
        # noise for a ~50-cell donor pseudobulk, plus the between-dataset floor
        # (same rule as codec.LinearGaussianView: D/n + D_sys)
        D_sys = np.var(np.array(dsoff), 0)
        D_eff = D / 50 + D_sys
        M_joint = posterior_matrix(W_joint, D, np.ones(1 + N_ENV))[0]
        M_solo = posterior_matrix(w_ifn[:, None], D, np.ones(1))[0]
        M_joint_f = posterior_matrix(W_joint, D_eff, np.ones(1 + N_ENV))[0]
        M_solo_f = posterior_matrix(w_ifn[:, None], D_eff, np.ones(1))[0]
        top = np.argsort(-w_ifn)[:50]
        M_naive = np.zeros(len(w_ifn))
        M_naive[top] = w_ifn[top] / (w_ifn[top] @ w_ifn[top])
        # panels (designed on training data only)
        ifn6 = np.array([gidx[g] for g in IFN6 if g in gidx])
        a12 = a_optimal(W_joint, D + 50 * D_sys, np.ones(1 + N_ENV), 12, n_cells=50,
                        cand=np.argsort(-np.abs(W_joint).max(1))[:600])
        M_ifn6 = posterior_matrix(W_joint[ifn6], D_eff[ifn6], np.ones(1 + N_ENV))[0]
        M_a12 = posterior_matrix(W_joint[a12], D_eff[a12], np.ones(1 + N_ENV))[0]
        Mfull = {"naive": M_naive, "solo": M_solo, "joint": M_joint,
                 "solo + batch floor": M_solo_f, "joint + batch floor": M_joint_f}
        Mpart = {"6-gene IFN score genes": (ifn6, M_ifn6), "12 genes A-optimal": (a12, M_a12)}
        leaf_idx = np.array([leaves.index(l) for l in leaf_hat])
        res = {}
        for k, M in Mfull.items():
            res[k] = np.asarray(Xt @ M).ravel() - (Cm @ M)[leaf_idx]
        for k, (g, M) in Mpart.items():
            res[k] = np.asarray(Xt[:, g] @ M).ravel() - (Cm[:, g] @ M)[leaf_idx]
        df = ot[["dataset", "donor_id", "disease", "sick", "leaf"]].copy()
        for k, v in res.items():
            df[k] = v
        cell_rows.append(df)
        # ---- C2. flow-like 8-marker identity panel, Gaussian with abstention --
        f = np.array([gidx[g] for g in FLOW8 if g in gidx])
        mu = Cm[:, f]
        var = np.array([[np.asarray(X[np.where(healthy & (obs.leaf == l).values)[0]][:, [j]]
                                    .toarray()).var() + 0.05 for j in f] for l in leaves])
        Xf = Xt[:, f].toarray()
        ll = -0.5 * (((Xf[:, None, :] - mu[None]) ** 2) / var[None] + np.log(var[None])).sum(2)
        post = np.exp(ll - ll.max(1, keepdims=True))
        post /= post.sum(1, keepdims=True)
        pred = np.array(leaves)[post.argmax(1)]
        lin_post = np.array([post[:, [i for i, x in enumerate(leaves) if LEAF2LIN[x] == l]].sum(1)
                             for l in lin_names]).T
        for thr in (0.0, 0.9):
            acc_leaf = post.max(1) >= thr
            acc_lin = lin_post.max(1) >= thr
            flow_rows.append(dict(dataset=test, thr=thr,
                                  leaf_cov=float(acc_leaf.mean()),
                                  leaf_acc=float(np.mean(pred[acc_leaf] == ot.leaf.values[acc_leaf])),
                                  lin_cov=float(acc_lin.mean()),
                                  lin_acc=float(np.mean(np.array(lin_names)[lin_post.argmax(1)][acc_lin]
                                                        == ot.lineage.values[acc_lin]))))
        print(id_rows[-1])
    return pd.concat(cell_rows), pd.DataFrame(id_rows), pd.DataFrame(flow_rows)


def evaluate(cells, methods):
    don = cells.groupby(["dataset", "donor_id", "leaf", "disease", "sick"], observed=True)
    cnt = don.size()
    d = don[methods].mean()[cnt >= 10].reset_index()
    out = {}
    for m in methods:
        within = []
        for (_, leaf), g in d.groupby(["dataset", "leaf"], observed=True):
            if g.sick.sum() >= 3 and (g.sick == 0).sum() >= 3:
                within.append(roc_auc_score(g.sick, g[m]))
        pooled = [roc_auc_score(g.sick, g[m]) for _, g in d.groupby("leaf", observed=True) if g.sick.nunique() == 2]
        h = d[d.sick == 0]
        eta = []
        for _, g in h.groupby("leaf", observed=True):
            if g.dataset.nunique() >= 3:
                tot = ((g[m] - g[m].mean()) ** 2).sum()
                bet = sum(len(x) * (x[m].mean() - g[m].mean()) ** 2 for _, x in g.groupby("dataset", observed=True))
                eta.append(bet / tot)
        pos = d[m] > TAU
        out[m] = dict(within_auroc=float(np.median(within)), n_strata=len(within),
                      pooled_auroc=float(np.median(pooled)), healthy_eta2=float(np.mean(eta)),
                      ifn_pos_sick=float(pos[d.sick == 1].mean()), ifn_pos_healthy=float(pos[d.sick == 0].mean()),
                      corr_with_joint=float(np.corrcoef(d[m], d["joint"])[0, 1]))
    by_dis = {}
    for dis in ["COVID-19", "systemic lupus erythematosus"]:
        sub = d[(d.disease == dis) | (d.sick == 0)]
        sub = sub[sub.dataset.isin(d.dataset[d.disease == dis].unique())]
        aucs = [roc_auc_score(g.sick, g["joint"]) for _, g in sub.groupby(["dataset", "leaf"], observed=True)
                if g.sick.sum() >= 3 and (g.sick == 0).sum() >= 3]
        by_dis[dis] = dict(within_auroc=float(np.median(aucs)), n_strata=len(aucs),
                           median_ifn_sick=float(sub[sub.sick == 1]["joint"].median()),
                           median_ifn_healthy=float(sub[sub.sick == 0]["joint"].median()))
    return out, by_dis, d


def main(out_json="dcsc/real/results_census.json"):
    X, obs, genes, wdf = load()
    print(f"{X.shape[0]:,} cells x {X.shape[1]:,} genes; {obs.donor_id.nunique()} donors; "
          f"{obs.dataset.nunique()} datasets")
    cells, ident, flow = lodo(X, obs, genes, wdf)
    methods = ["naive", "solo", "joint", "solo + batch floor", "joint + batch floor",
               "6-gene IFN score genes", "12 genes A-optimal"]
    ifn, by_dis, donors = evaluate(cells, methods)
    summary = dict(
        n_cells=int(X.shape[0]), n_genes=int(X.shape[1]), n_donors=int(obs.donor_id.nunique()),
        datasets=obs.groupby("dataset", observed=True).agg(cells=("leaf", "size"), assay=("assay", "first"),
                                            donors=("donor_id", "nunique")).reset_index().to_dict("records"),
        identity=ident.to_dict("records"),
        identity_mean=dict(lineage=float(ident.lineage_acc.mean()), leaf=float(ident.leaf_acc.mean())),
        flow_panel=flow.groupby("thr")[["leaf_cov", "leaf_acc", "lin_cov", "lin_acc"]].mean()
                        .reset_index().to_dict("records"),
        ifn=ifn, ifn_by_disease=by_dis, n_donor_leaf=int(len(donors)),
    )
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=1, default=str)
    print(json.dumps({k: summary[k] for k in ["identity_mean", "flow_panel", "ifn", "ifn_by_disease"]},
                     indent=1))


if __name__ == "__main__":
    main()
