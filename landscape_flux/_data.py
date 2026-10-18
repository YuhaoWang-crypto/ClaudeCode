"""Download-and-cache helpers for the public datasets used by the demos.

Everything is fetched at run time into ./data_cache (git-ignored).
"""
import gzip
import os
import subprocess
import tarfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "data_cache")
FIG = os.path.join(ROOT, "figures")

URLS = {
    # scVelo endocrine pancreas (Bastidas-Ponce 2019), spliced/unspliced layers
    "pancreas.h5ad": "https://raw.githubusercontent.com/theislab/scvelo_notebooks/"
                     "master/data/Pancreas/endocrinogenesis_day15.h5ad",
    # Yeast metabolic cycle (Tu et al. 2005 Science), GEO GSE3431
    "GSE3431_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE3nnn/"
                                    "GSE3431/matrix/GSE3431_series_matrix.txt.gz",
    # 10x Genomics PBMC 3k
    "pbmc3k.tar.gz": "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/"
                     "pbmc3k_filtered_gene_bc_matrices.tar.gz",
    # Paul et al. 2015 myeloid progenitors (scanpy mirror)
    "paul15.h5": "https://falexwolf.de/data/paul15.h5",
}


def fetch(name):
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, name)
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        subprocess.run(["curl", "-sSL", "--max-time", "600", "-o", path, URLS[name]],
                       check=True)
    return path


def fig(name):
    os.makedirs(FIG, exist_ok=True)
    return os.path.join(FIG, name)


def load_pancreas():
    import anndata as ad
    return ad.read_h5ad(fetch("pancreas.h5ad"))


def load_paul15():
    import scanpy as sc
    sc.settings.datasetdir = CACHE
    fetch("paul15.h5")
    return sc.datasets.paul15()


def load_pbmc3k():
    """Return (counts cells×genes as scipy.sparse CSR, gene symbols list)."""
    import scipy.io
    path = fetch("pbmc3k.tar.gz")
    out = os.path.join(CACHE, "pbmc3k")
    if not os.path.isdir(out):
        with tarfile.open(path) as t:
            t.extractall(out, filter="data")
    d = os.path.join(out, "filtered_gene_bc_matrices", "hg19")
    X = scipy.io.mmread(os.path.join(d, "matrix.mtx")).T.tocsr()
    genes = [l.split("\t")[1].strip() for l in open(os.path.join(d, "genes.tsv"))]
    return X, genes


def load_ymc():
    """Yeast metabolic cycle: returns (probe ids, gene-by-time matrix, times in min)."""
    import numpy as np
    path = fetch("GSE3431_series_matrix.txt.gz")
    ids, rows = [], []
    with gzip.open(path, "rt") as f:
        inside = False
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                inside = True; next(f); continue
            if line.startswith("!series_matrix_table_end"):
                break
            if inside:
                p = line.rstrip("\n").split("\t")
                try:
                    vals = [float(x) for x in p[1:]]
                except ValueError:
                    continue
                ids.append(p[0].strip('"')); rows.append(vals)
    X = np.array(rows)
    t = np.arange(X.shape[1]) * 25.0  # Tu et al.: one sample every ~25 min
    return ids, X, t


# ------------------------------------------------------------------ BBBC021
BBBC021 = "https://data.broadinstitute.org/bbbc/BBBC021/"


class _HTTPRangeFile:
    """Read-only, seekable file object over HTTP range requests, so that
    zipfile can pull single members out of an ~850 MB plate archive."""

    def __init__(self, url):
        import urllib.request
        self.url, self.pos = url, 0
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=60) as r:
            self.size = int(r.headers["Content-Length"])

    def seekable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else (self.pos + off if whence == 1 else self.size + off)
        return self.pos

    def read(self, n=-1):
        import urllib.request
        if n is None or n < 0:
            n = self.size - self.pos
        if n == 0 or self.pos >= self.size:
            return b""
        end = min(self.pos + n, self.size) - 1
        req = urllib.request.Request(self.url, headers={"Range": f"bytes={self.pos}-{end}"})
        import time
        for attempt in range(6):
            try:
                with urllib.request.urlopen(req, timeout=120) as r:
                    data = r.read()
                break
            except OSError:
                if attempt == 5:
                    raise
                time.sleep(2 ** attempt)
        self.pos += len(data)
        return data


def bbbc021_metadata():
    import pandas as pd
    d = os.path.join(CACHE, "bbbc021"); os.makedirs(d, exist_ok=True)
    for f in ("BBBC021_v1_image.csv", "BBBC021_v1_moa.csv"):
        p = os.path.join(d, f)
        if not os.path.exists(p):
            subprocess.run(["curl", "-sSL", "-o", p, BBBC021 + f], check=True)
    im = pd.read_csv(os.path.join(d, "BBBC021_v1_image.csv"))
    moa = pd.read_csv(os.path.join(d, "BBBC021_v1_moa.csv"))
    return im.merge(moa, left_on=["Image_Metadata_Compound", "Image_Metadata_Concentration"],
                    right_on=["compound", "concentration"])


def bbbc021_dapi(plate_path, filenames):
    """Fetch DAPI TIFFs of one plate (e.g. 'Week1/Week1_22123') into the cache,
    downloading only the needed zip members.  Returns list of local paths."""
    import zipfile
    d = os.path.join(CACHE, "bbbc021", plate_path.split("/")[-1])
    os.makedirs(d, exist_ok=True)
    need = [f for f in filenames if not os.path.exists(os.path.join(d, f))]
    if need:
        zname = "BBBC021_v1_images_" + plate_path.split("/")[-1] + ".zip"
        with zipfile.ZipFile(_HTTPRangeFile(BBBC021 + zname)) as z:
            members = {os.path.basename(m): m for m in z.namelist()}
            for f in need:
                with z.open(members[f]) as src, open(os.path.join(d, f), "wb") as dst:
                    dst.write(src.read())
    return [os.path.join(d, f) for f in filenames]
