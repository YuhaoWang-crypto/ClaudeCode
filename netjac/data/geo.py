"""
GEO series-matrix loader (NCBI FTP, open access): expression table plus a
tidy sample sheet built from the `!Sample_characteristics_ch1` lines, and
GPL probe -> gene symbol maps from the platform `.annot.gz` files.
Files are cached under figures/_data_cache/geo_<gse>/ and geo_<gpl>/.
"""
import gzip
import io
import os
import subprocess

import numpy as np

from netjac import DATA_CACHE

FTP = "https://ftp.ncbi.nlm.nih.gov/geo"


def _stub(acc):
    return acc[:-3] + "nnn" if len(acc) > 6 else acc[:3] + "nnn"


def _fetch(url, fn):
    if not os.path.exists(fn):
        os.makedirs(os.path.dirname(fn), exist_ok=True)
        subprocess.run(["curl", "-sS", "--max-time", "1800", "-o", fn + ".part", url], check=True)
        os.replace(fn + ".part", fn)
    return fn


def series_matrix(gse):
    """Returns (X genes-or-probes x samples DataFrame, meta DataFrame).
    meta has one column per characteristics key ('key: value' lines) plus
    'title'."""
    import pandas as pd
    fn = _fetch(f"{FTP}/series/{_stub(gse)}/{gse}/matrix/{gse}_series_matrix.txt.gz",
                os.path.join(DATA_CACHE, f"geo_{gse.lower()}", "matrix.txt.gz"))
    lines = gzip.open(fn, "rt", errors="replace").read().split("\n")
    meta, titles, ids, desc = {}, None, None, None
    for l in lines:
        if l.startswith("!Sample_title"):
            titles = [x.strip('"') for x in l.split("\t")[1:]]
        elif l.startswith("!Sample_description") and desc is None:
            desc = [x.strip('"') for x in l.split("\t")[1:]]
        elif l.startswith("!Sample_geo_accession"):
            ids = [x.strip('"') for x in l.split("\t")[1:]]
        elif l.startswith("!Sample_characteristics_ch1"):
            vals = [x.strip('"') for x in l.split("\t")[1:]]
            key = vals[0].split(":")[0].strip()
            meta[key] = [v.split(":", 1)[1].strip() if ":" in v else v for v in vals]
    i = lines.index("!series_matrix_table_begin")
    j = lines.index("!series_matrix_table_end")
    X = pd.read_csv(io.StringIO("\n".join(lines[i + 1:j])), sep="\t", index_col=0)
    M = pd.DataFrame(meta, index=ids)
    M["title"] = titles
    if desc is not None and len(desc) == len(M):
        M["description"] = desc
    return X[ids], M


def gpl_symbols(gpl):
    """probe ID -> gene symbol (first symbol of multi-mapping probes;
    probes without a symbol are dropped)."""
    import pandas as pd
    fn = _fetch(f"{FTP}/platforms/{_stub(gpl)}/{gpl}/annot/{gpl}.annot.gz",
                os.path.join(DATA_CACHE, f"geo_{gpl.lower()}", f"{gpl}.annot.gz"))
    lines = gzip.open(fn, "rt", errors="replace").read().split("\n")
    i = lines.index("!platform_table_begin")
    j = lines.index("!platform_table_end")
    A = pd.read_csv(io.StringIO("\n".join(lines[i + 1:j])), sep="\t", index_col=0,
                    usecols=[0, 2], dtype=str)
    s = A.iloc[:, 0].dropna()
    return s.str.split("///").str[0].str.strip()


def collapse_to_genes(X, sym):
    """Keep, for each gene, the probe with the highest mean expression."""
    X = X.loc[X.index.intersection(sym.index)]
    g = sym.loc[X.index]
    best = X.mean(axis=1).groupby(g).idxmax()
    out = X.loc[best.values]
    out.index = best.index
    return out


def quantile_normalise(X):
    """Column-wise quantile normalisation of a genes x samples DataFrame."""
    A = X.values
    order = np.argsort(A, axis=0)
    ref = np.sort(A, axis=0).mean(axis=1)
    out = np.empty_like(A, dtype=float)
    for k in range(A.shape[1]):
        out[order[:, k], k] = ref
    return X.__class__(out, index=X.index, columns=X.columns)
