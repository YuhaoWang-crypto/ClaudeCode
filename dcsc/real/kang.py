"""
Real calibration contrast for one program: type-I interferon.

Kang et al., Nat Biotechnol 2018 (GSE96583, batch 2): PBMCs from 8 lupus
donors, control vs IFN-β (6 h), 10x 3' v2. Loading of the IFN program on a
gene = mean over donors of [stim − ctrl] pseudobulk (log1p CP10k), per Kang
cell type and pooled. Unit 1 = one full in-vitro IFN-β stimulation.

Kang cell types are mapped onto the Census identity leaves.
"""
import gzip
import os

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

DIR = "data_cache"
KANG2LEAF = {"CD14+ Monocytes": "CD14_Mono", "FCGR3A+ Monocytes": "CD16_Mono",
             "Dendritic cells": "cDC", "CD4 T cells": "CD4_T", "CD8 T cells": "CD8_T",
             "NK cells": "NK", "B cells": "B"}


def _read(mtx, bc):
    X = scipy.io.mmread(gzip.open(os.path.join(DIR, "kang", mtx))).T.tocsr()
    b = pd.read_csv(os.path.join(DIR, "kang", bc), header=None)[0].values
    return X, b


def ifn_loadings(cache=os.path.join(DIR, "kang_ifn.parquet")):
    """DataFrame indexed by Ensembl ID: columns = leaves + 'pooled' + 'pooled_se'."""
    if os.path.exists(cache):
        return pd.read_parquet(cache)
    genes = pd.read_csv(os.path.join(DIR, "GSE96583_batch2.genes.tsv.gz"), sep="\t", header=None)
    meta = pd.read_csv(os.path.join(DIR, "GSE96583_batch2.total.tsne.df.tsv.gz"), sep="\t", index_col=0)
    Xc, bc = _read("GSM2560248_2.1.mtx.gz", "GSM2560248_barcodes.tsv.gz")
    Xs, bs = _read("GSM2560249_2.2.mtx.gz", "GSM2560249_barcodes.tsv.gz")
    X = sp.vstack([Xc, Xs]).tocsr().astype(np.float32)
    m = meta.loc[np.concatenate([bc, bs])]
    lib = np.asarray(X.sum(1)).ravel()
    X = sp.diags(1e4 / lib) @ X
    X.data = np.log1p(X.data)
    single = (m.multiplets == "singlet").values
    cols, ses = {}, {}
    for kct, leaf in KANG2LEAF.items():
        diffs = []
        for ind in m.ind.unique():
            base = single & (m.cell == kct).values & (m.ind == ind).values
            c, s = base & (m.stim == "ctrl").values, base & (m.stim == "stim").values
            if c.sum() >= 10 and s.sum() >= 10:
                diffs.append(np.asarray(X[s].mean(0)).ravel() - np.asarray(X[c].mean(0)).ravel())
        cols[leaf] = np.mean(diffs, 0)
        ses[leaf] = np.std(diffs, 0, ddof=1) / np.sqrt(len(diffs))
    df = pd.DataFrame(cols, index=genes[0].values)
    se = pd.DataFrame(ses, index=genes[0].values)
    df["pooled"] = df[list(KANG2LEAF.values())].mean(1)
    # s.e. of the pooled loading across 8 donors x 7 cell types (independent-ish)
    df["pooled_se"] = np.sqrt((se ** 2).sum(1)) / len(KANG2LEAF)
    df = df[~df.index.duplicated()]
    df.to_parquet(cache)
    return df


if __name__ == "__main__":
    df = ifn_loadings()
    sym = pd.read_csv(os.path.join(DIR, "GSE96583_batch2.genes.tsv.gz"), sep="\t", header=None,
                      index_col=0)[1]
    top = df.pooled.sort_values(ascending=False).head(15)
    print(pd.DataFrame({"gene": sym.reindex(top.index).values, "lfc": top.values.round(2)}))
