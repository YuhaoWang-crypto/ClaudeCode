"""Build the leak-free reference and the two unseen queries on a common gene space.

Reference  = HCA bone-marrow atlas minus the Tabula Sapiens study minus the held-out lineages.
Query A    = Triana 2021 (different lab, different platform: BD Rhapsody WTA).
Query B    = Tabula Sapiens bone marrow, 10x cells (study removed from the reference).
HVGs are chosen on the reference only; genes a query lacks are zero-filled.
"""
import os, json
import numpy as np, scipy.sparse as sp, anndata as ad, scanpy as sc
from common import *

os.makedirs(WORK, exist_ok=True); os.makedirs(OUT, exist_ok=True)

def load(name):
    a = ad.read_h5ad(f"{DATA}/{name}.h5ad")
    X = a.raw.X if a.raw is not None else a.X
    var = a.raw.var if a.raw is not None else a.var
    a = ad.AnnData(sp.csr_matrix(X, dtype=np.float32), obs=a.obs.copy(), var=var.copy())
    assert (a.X.data[:100000] % 1 == 0).all(), "expected raw counts"
    a.obs["fine"] = a.obs["cell_type"].astype(str)
    a.obs["coarse"] = a.obs["fine"].map(FINE2COARSE)
    assert a.obs["coarse"].notna().all(), a.obs.loc[a.obs.coarse.isna(), "fine"].unique()
    return a

ref = load("ref_hca_bm")
n0 = ref.n_obs
ref = ref[ref.obs["Dataset"] != HELD_OUT_STUDY].copy()
n1 = ref.n_obs
ref_full = ref
ref = ref[~ref.obs["coarse"].isin(HELD_OUT_LINEAGES)].copy()
ref.obs["batch"] = ref.obs["donor_id"].astype(str)
sc.pp.filter_genes(ref, min_cells=20)
sc.pp.highly_variable_genes(ref, n_top_genes=3000, flavor="seurat_v3", batch_key="batch", span=0.5)
hvg = ref.var_names[ref.var.highly_variable]
symbols = ref.var.loc[hvg, "feature_name"].astype(str)
ref[:, hvg].write(f"{WORK}/ref.h5ad")
# full-gene copy of the reference for CellTypist custom training (log1p-normalised later)
ref.write(f"{WORK}/ref_allgenes.h5ad")

info = {"ref_cells_original": n0, "ref_cells_after_study_removal": n1, "ref_cells_training": ref.n_obs,
        "ref_donors": ref.obs.batch.nunique(), "ref_fine_labels": ref.obs.fine.nunique(), "queries": {}}

def prep_query(name, keep=None):
    q = load(name)
    if keep is not None: q = q[keep(q.obs)].copy()
    q.obs["batch"] = name + "_" + q.obs["donor_id"].astype(str)
    q.write(f"{WORK}/{name}_allgenes.h5ad")
    present = np.isin(hvg, q.var_names)
    idx = q.var_names.get_indexer(hvg[present])
    X = sp.lil_matrix((q.n_obs, len(hvg)), dtype=np.float32).tocsc()
    X = sp.csc_matrix(q.X)[:, idx]
    full = sp.csr_matrix((q.n_obs, len(hvg)), dtype=np.float32).tolil()
    full = sp.hstack([X, sp.csr_matrix((q.n_obs, int((~present).sum())), dtype=np.float32)]).tocsc()
    order = np.argsort(np.concatenate([np.where(present)[0], np.where(~present)[0]]))
    qh = ad.AnnData(full[:, order].tocsr(), obs=q.obs.copy())
    qh.var_names = hvg
    qh.write(f"{WORK}/{name}.h5ad")
    info["queries"][name] = {"cells": q.n_obs, "hvg_missing_zero_filled": int((~present).sum()),
                             "donors": sorted(q.obs.donor_id.astype(str).unique()),
                             "coarse_counts": q.obs.coarse.value_counts().to_dict()}

prep_query("q_triana_wta")
prep_query("q_ts_bm", keep=lambda o: o["method"] == "10X")
overlap = set(ref.obs.donor_id.astype(str)) & {d for q in info["queries"].values() for d in q["donors"]}
assert not overlap, overlap
info["donor_overlap_ref_vs_query"] = sorted(overlap)
json.dump(info, open(f"{OUT}/design.json", "w"), indent=1)
print(json.dumps(info, indent=1))
