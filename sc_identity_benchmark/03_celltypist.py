"""CellTypist, two ways:
 (a) custom model trained only on our leak-free reference (same held-out design as scANVI);
 (b) the published pretrained Immune_All_Low weights (zero training here; training set may overlap queries).
"""
import os, re, numpy as np, pandas as pd, anndata as ad, scanpy as sc, celltypist
from celltypist import models
from common import *

def norm(a):
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a); return a

mpath = f"{WORK}/celltypist_custom.pkl"
if not os.path.exists(mpath):
    ref = norm(ad.read_h5ad(f"{WORK}/ref_allgenes.h5ad"))
    m = celltypist.train(ref, labels="fine", n_jobs=4, feature_selection=True, use_SGD=True,
                         mini_batch=True, check_expression=False)
    m.write(mpath)
models.download_models(model="Immune_All_Low.pkl")

RULES = [  # pretrained label -> coarse lineage
    (r"plasma cell|plasmablast", "Plasma"), (r"pDC", "pDC"),
    (r"HSC|MPP|CMP|GMP|ELP|MNP|Myeloid progenitor|Early lymphoid", "HSPC"),
    (r"erythro|MEMP|megakaryocyte|platelet|Early MK", "Ery/Mk"),
    (r"neutrophil|promyelocyte|myelocyte|granulocyte", "Neutrophil lineage"),
    (r"monocyte|macrophage|mono-mac|kupffer|microglia|Hofbauer", "Mono/Mac"), (r"DC", "cDC"),
    (r"mast|basophil", "Mast/Baso"),
    (r"pro-B|pre-B|pre-pro-B|B-cell lineage|Transitional B|Immature B", "B progenitor"),
    (r"B cell", "B mature"), (r"NK cells|ILC|Transitional NK", "NK"),
    (r"T cell|Treg|MAIT|thymocyte|Tcm|Tem|Trm|helper|cytotoxic|NKT|T\(|ETP|CD8a", "T"),
    (r"endothelial|fibroblast|epithelial", "Stroma"),
]
def to_coarse(lab):
    for pat, c in RULES:
        if re.search(pat, lab, flags=re.I): return c
    return "other"

for name in ["q_triana_wta", "q_ts_bm"]:
    q = norm(ad.read_h5ad(f"{WORK}/{name}_allgenes.h5ad"))
    out = q.obs[["fine", "coarse"]].copy()
    r = celltypist.annotate(q, model=mpath, majority_voting=False)
    out["ct_custom_fine"] = r.predicted_labels.predicted_labels.values
    out["ct_custom_coarse"] = out.ct_custom_fine.map(FINE2COARSE)
    out["ct_custom_maxp"] = r.probability_matrix.max(1).values
    qs = q.copy(); qs.var_names = q.var["feature_name"].astype(str).str.split("_ENSG").str[0].values
    qs.var_names_make_unique()
    r = celltypist.annotate(qs, model="Immune_All_Low.pkl", majority_voting=False)
    out["ct_pre_label"] = r.predicted_labels.predicted_labels.values
    out["ct_pre_coarse"] = out.ct_pre_label.map(to_coarse)
    out["ct_pre_maxp"] = r.probability_matrix.max(1).values
    print(name, "pretrained labels left unmapped:", out.loc[out.ct_pre_coarse == "other", "ct_pre_label"].value_counts().to_dict())
    out.to_csv(f"{OUT}/pred_celltypist_{name}.csv")
