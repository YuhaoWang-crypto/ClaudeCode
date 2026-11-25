"""CytoTRACE 2 (published pretrained weights, no splicing information) on the unseen queries.
Stratified subsample (<=300 cells per annotated type) to keep the CPU run short.
"""
import os, numpy as np, pandas as pd, anndata as ad
from cytotrace2_py.cytotrace2_py import cytotrace2
from common import *

rng = np.random.default_rng(0)
for name in ["q_triana_wta", "q_ts_bm"]:
    q = ad.read_h5ad(f"{WORK}/{name}_allgenes.h5ad")
    keep = np.concatenate([rng.choice(ix, min(300, len(ix)), replace=False)
                           for ix in q.obs.groupby("fine", observed=True).indices.values()])
    q = q[np.sort(keep)].copy()
    q.var_names = q.var["feature_name"].astype(str).str.split("_ENSG").str[0].values
    q = q[:, ~q.var_names.duplicated() & ~q.var_names.isin(["NA", "NaN", "nan", "NULL", "null", "N/A"])].copy()
    q = q[:, np.asarray((q.X > 0).sum(0)).ravel() >= 5].copy()
    d = f"{WORK}/ct2_{name}"; os.makedirs(d, exist_ok=True)
    df = pd.DataFrame(q.X.T.toarray().astype(int), index=q.var_names, columns=q.obs_names)
    df.to_csv(f"{d}/expr.txt", sep="\t")
    res = cytotrace2(f"{d}/expr.txt", species="human", output_dir=d, disable_plotting=True, max_cores=4)
    res = res.reindex(q.obs_names)
    out = q.obs[["fine", "coarse"]].join(res)
    out["n_genes_detected"] = np.asarray((q.X > 0).sum(1)).ravel()   # CytoTRACE-1-style baseline
    out["potency_truth"] = out.fine.map(FINE2POTENCY).fillna(0).astype(int)
    out.to_csv(f"{OUT}/cytotrace2_{name}.csv")
    print(name, out.shape, list(res.columns))
