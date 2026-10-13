"""
Replogle et al. 2022 (Cell 185:2559) genome-scale CRISPRi Perturb-seq,
figshare article 20029387 (CC BY 4.0), downloaded with curl into
figures/_data_cache/perturbseq/ (git-ignored):

  K562_essential_normalized_bulk_01.h5ad        80 MB  2,285 perturbations
  K562_gwps_normalized_bulk_01.h5ad            375 MB  ~11k perturbations
  rpe1_normalized_bulk_01.h5ad                  95 MB  RPE1 essential screen
  K562_essential_normalized_singlecell_01.h5ad  10.7 GB 310k cells x 8,563 genes

"normalized" = per-gem-group z-score against that gem group's non-targeting
cells, so a pseudo-bulk row is the mean shift of each gene (in control-cell SD
units) when the perturbed gene is knocked down.

`reduce_singlecell` streams the 10.7 GB dense matrix once and keeps:
  * per-perturbation means in two random halves of the cells (split-half
    reproducibility for every inferred edge), and
  * every non-targeting control cell (the unperturbed "equilibrium ensemble"
    whose fluctuations the fluctuation-response test needs).
"""
import os

import numpy as np

from netjac import DATA_CACHE

DIR = os.path.join(DATA_CACHE, "perturbseq")
URLS = {
    "K562_essential_bulk": "https://ndownloader.figshare.com/files/35780870",
    "K562_gwps_bulk": "https://ndownloader.figshare.com/files/35773217",
    "rpe1_bulk": "https://ndownloader.figshare.com/files/35775512",
    "K562_essential_sc": "https://ndownloader.figshare.com/files/35773075",
}
FILES = {
    "K562_essential_bulk": "K562_essential_normalized_bulk_01.h5ad",
    "K562_gwps_bulk": "K562_gwps_normalized_bulk_01.h5ad",
    "rpe1_bulk": "rpe1_normalized_bulk_01.h5ad",
    "K562_essential_sc": "K562_essential_normalized_singlecell_01.h5ad",
}


def path(key):
    fn = os.path.join(DIR, FILES[key])
    if not os.path.exists(fn):
        import subprocess
        os.makedirs(DIR, exist_ok=True)
        subprocess.run(["curl", "-sSL", "-C", "-", "-o", fn, URLS[key]], check=True)
    return fn


def load_bulk(key="K562_essential_bulk"):
    """Pseudo-bulk response matrix. Returns dict(X: perturbations x genes,
    pert_gene, pert_id (Ensembl), gene, gene_id, obs (DataFrame), control
    (bool, non-targeting rows))."""
    import anndata as ad
    a = ad.read_h5ad(path(key))
    idx = list(a.obs.index)
    parts = [s.split("_") for s in idx]
    pert_gene = np.array([p[1] for p in parts])
    pert_id = np.array([p[-1] for p in parts])
    return {"X": np.asarray(a.X, dtype=np.float32), "pert_gene": pert_gene,
            "pert_id": pert_id, "gene": np.array(a.var["gene_name"].astype(str)),
            "gene_id": np.array(a.var.index.astype(str)), "obs": a.obs,
            "control": pert_gene == "non-targeting", "key": key}


def reduce_singlecell(key="K562_essential_sc", chunk=4000, seed=0):
    """One pass over the single-cell matrix -> cached npz with split-half
    pseudo-bulk means and all control cells."""
    out = os.path.join(DIR, f"{key}_reduced.npz")
    if os.path.exists(out):
        return dict(np.load(out, allow_pickle=True))
    import h5py
    f = h5py.File(path(key), "r")
    X = f["X"]
    n, g = X.shape
    code = f["obs/gene_transcript"][:]
    cats = np.array([c.decode() if isinstance(c, bytes) else c
                     for c in f["obs/__categories/gene_transcript"][:]])
    gene_id = np.array([c.decode() if isinstance(c, bytes) else c
                        for c in f["var/gene_id"][:]])
    gn_cats = np.array([c.decode() if isinstance(c, bytes) else c
                        for c in f["var/__categories/gene_name"][:]])
    gene = gn_cats[f["var/gene_name"][:]]
    gem = f["obs/gem_group"][:]
    P = len(cats)
    rng = np.random.default_rng(seed)
    half = rng.integers(0, 2, n)
    ctrl_codes = [i for i, c in enumerate(cats) if "non-targeting" in c]
    is_ctrl = np.isin(code, ctrl_codes)
    S = np.zeros((P, 2, g))
    N = np.zeros((P, 2))
    ctrl = np.zeros((is_ctrl.sum(), g), dtype=np.float32)
    ctrl_gem = gem[is_ctrl]
    k = 0
    for s in range(0, n, chunk):
        e = min(n, s + chunk)
        B = X[s:e]
        c, h = code[s:e], half[s:e]
        for hh in (0, 1):
            m = h == hh
            np.add.at(S[:, hh, :], c[m], B[m])
            np.add.at(N[:, hh], c[m], 1)
        mc = is_ctrl[s:e]
        ctrl[k:k + mc.sum()] = B[mc]
        k += mc.sum()
    f.close()
    means = S / np.maximum(N, 1)[:, :, None]
    res = {"half_means": means.astype(np.float32), "half_n": N,
           "pert": cats, "gene": gene, "gene_id": gene_id,
           "ctrl_cells": ctrl, "ctrl_gem": ctrl_gem}
    np.savez(out, **res)
    return res
