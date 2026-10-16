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
    probes without a symbol are dropped).

    Uses the curated GEO `.annot.gz` when GEO has one. Affymetrix Gene ST /
    HTA platforms (e.g. GPL17586) have none, so fall back to the platform
    table itself and parse its `gene_assignment` column."""
    import pandas as pd
    cache = os.path.join(DATA_CACHE, f"geo_{gpl.lower()}", f"{gpl}.annot.gz")
    try:
        fn = _fetch(f"{FTP}/platforms/{_stub(gpl)}/{gpl}/annot/{gpl}.annot.gz", cache)
        with open(fn, "rb") as f:
            if f.read(2) != b"\x1f\x8b":
                raise OSError("not gzip: GEO served an error page")
    except Exception:
        if os.path.exists(cache):
            os.remove(cache)
        return _platform_table_symbols(gpl)
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


def _platform_table_symbols(gpl):
    """probe -> symbol from the GEO platform table's `gene_assignment` column.

    The full table is ~190 MB, so it is parsed once into a small cached map
    and then deleted; only the map is kept.
    """
    import pandas as pd
    small = os.path.join(DATA_CACHE, f"geo_{gpl.lower()}", f"{gpl}_sym.tsv.gz")
    if not os.path.exists(small):
        big = os.path.join(DATA_CACHE, f"geo_{gpl.lower()}", f"{gpl}_table.txt")
        _fetch("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
               f"?acc={gpl}&targ=self&form=text&view=full", big)
        rows, hdr, on = [], None, False
        with open(big, "rt", errors="replace") as f:
            for line in f:
                if line.startswith("!platform_table_begin"):
                    on = True
                    continue
                if line.startswith("!platform_table_end"):
                    break
                if not on:
                    continue
                p = line.rstrip("\n").split("\t")
                if hdr is None:
                    hdr = p
                    try:
                        gi = hdr.index("gene_assignment")
                    except ValueError:
                        raise ValueError(f"{gpl} platform table has no gene_assignment column")
                    continue
                if len(p) <= gi or not p[gi] or p[gi] == "---":
                    continue
                # "ENST... // SYMBOL // description // cytoband // entrez /// ..."
                first = p[gi].split("///")[0].split("//")
                if len(first) < 2:
                    continue
                sym = first[1].strip()
                if sym and sym != "---":
                    rows.append((p[0], sym))
        os.makedirs(os.path.dirname(small), exist_ok=True)
        pd.DataFrame(rows, columns=["ID", "symbol"]).to_csv(small, sep="\t", index=False)
        os.remove(big)
    return pd.read_csv(small, sep="\t", index_col=0)["symbol"]
