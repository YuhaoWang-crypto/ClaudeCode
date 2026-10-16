"""
LINCS L1000 Phase 2 (GSE70138) level-5 signatures as a DRUG lever library.

The keratinocyte library in `netjac.data.keratinocyte` has nine measured
cytokine levers. This one has ~2,100 compounds, but in cancer / immortalised
cell lines, and only 978 of its 12,328 genes are MEASURED (the "landmark"
genes); the rest are inferred by the L1000 model. Both limits matter for any
conclusion drawn from it, so `signature_matrix` labels the gene space.

Files (NCBI FTP, open access), cached under figures/_data_cache/lincs_gse70138:
  GSE70138_Broad_LINCS_sig_info      one row per signature (compound, cell, dose, time)
  GSE70138_Broad_LINCS_pert_info     compound annotation
  GSE70138_Broad_LINCS_gene_info     gene id -> symbol, landmark flag
  GSE70138_Broad_LINCS_Level5_COMPZ  118,050 x 12,328 z-score matrix (5.9 GB)

The matrix is downloaded straight through gunzip (never stored compressed)
because the container's disk allowance does not fit both copies.
"""
import gzip
import os
import subprocess

import numpy as np

from netjac import DATA_CACHE

DIR = os.path.join(DATA_CACHE, "lincs_gse70138")
FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl"
_META = {"sig": "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz",
         "pert": "GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz",
         "gene": "GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz",
         "metrics": "GSE70138_Broad_LINCS_sig_metrics_2017-03-06.txt.gz"}
_MATRIX = "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz"
GCTX = os.path.join(DIR, "level5.gctx")


def _meta(key):
    import pandas as pd
    fn = os.path.join(DIR, _META[key])
    if not os.path.exists(fn):
        os.makedirs(DIR, exist_ok=True)
        subprocess.run(["curl", "-sS", "--max-time", "900", "-o", fn,
                        f"{FTP}/{_META[key]}"], check=True)
    return pd.read_csv(fn, sep="\t", dtype=str)


def sig_info():
    return _meta("sig")


def sig_metrics():
    return _meta("metrics")


def gene_info():
    """gene_info with pr_gene_symbol and the landmark flag pr_is_lm."""
    return _meta("gene")


def fetch_matrix():
    """Download and decompress the level-5 matrix if it is not there yet."""
    if os.path.exists(GCTX):
        return GCTX
    os.makedirs(DIR, exist_ok=True)
    part = GCTX + ".part"
    cmd = f'curl -sS --max-time 7000 "{FTP}/{_MATRIX}" | gunzip > "{part}"'
    subprocess.run(["bash", "-c", cmd], check=True)
    os.replace(part, GCTX)
    return GCTX


def signature_matrix(symbols=None, space="landmark", block=4000):
    """Read the level-5 z-scores for a set of genes.

    space   "landmark" = the 978 genes L1000 actually measures
            "all"      = all 12,328, i.e. measured plus model-inferred
    symbols restrict further to these gene symbols (after the space filter).

    Returns dict with
      Z       signatures x genes, float32
      genes   gene symbols of the columns
      sig_id  signature ids of the rows
      space   which space was used
    """
    import h5py
    fetch_matrix()
    G = gene_info()
    if space == "landmark":
        G = G[G.pr_is_lm == "1"]
    elif space != "all":
        raise ValueError("space must be 'landmark' or 'all'")
    if symbols is not None:
        G = G[G.pr_gene_symbol.isin(set(symbols))]
    with h5py.File(GCTX, "r") as f:
        rid = [x.decode() if isinstance(x, bytes) else str(x)
               for x in f["/0/META/ROW/id"][:]]
        cid = [x.decode() if isinstance(x, bytes) else str(x)
               for x in f["/0/META/COL/id"][:]]
        dset = f["/0/DATA/0/matrix"]
        # GCTX is written column-major, so h5py sees (n_col, n_row) =
        # (signatures, genes) here. Check instead of assuming.
        if dset.shape != (len(cid), len(rid)):
            raise ValueError(f"unexpected gctx layout {dset.shape} "
                             f"vs ({len(cid)}, {len(rid)})")
        pos = {g: i for i, g in enumerate(rid)}
        keep = [(s, pos[gid]) for s, gid in zip(G.pr_gene_symbol, G.pr_gene_id)
                if gid in pos]
        idx = np.array([k[1] for k in keep])
        order = np.argsort(idx)
        idx_sorted = idx[order]
        out = np.empty((dset.shape[0], len(idx_sorted)), dtype=np.float32)
        for a in range(0, dset.shape[0], block):
            b = min(a + block, dset.shape[0])
            out[a:b] = dset[a:b, :][:, idx_sorted]
    genes = np.array([keep[i][0] for i in order])
    return {"Z": out, "genes": genes, "sig_id": np.array(cid), "space": space}
