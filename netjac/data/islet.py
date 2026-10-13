"""
Human islet single-cell RNA-seq, Segerstolpe et al. 2016 (Cell Metab 24:593),
ArrayExpress E-MTAB-5061: 3,514 Smart-seq2 cells from 6 non-diabetic and 4
type-2-diabetic donors, with donor HbA1c and BMI.
Files are fetched from the EBI BioStudies FTP mirror into
figures/_data_cache/islet/ (git-ignored); the parsed beta-cell matrix is
cached as an npz.
"""
import os
import subprocess

import numpy as np

from netjac import DATA_CACHE

DIR = os.path.join(DATA_CACHE, "islet")
BASE = "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/061/E-MTAB-5061/Files/"
FILES = ["E-MTAB-5061.sdrf.txt", "pancreas_refseq_rpkms_counts_3514sc.txt"]


def _fetch():
    os.makedirs(DIR, exist_ok=True)
    for f in FILES:
        fn = os.path.join(DIR, f)
        if not os.path.exists(fn):
            subprocess.run(["curl", "-sS", "-o", fn, BASE + f], check=True)


def load(cell_type="beta cell"):
    """log2(RPKM+1) for one inferred cell type (quality 'OK' only).
    Returns dict(X cells x genes, genes, donor, t2d (bool), hba1c, bmi)."""
    out = os.path.join(DIR, f"{cell_type.replace(' ', '_')}.npz")
    if not os.path.exists(out):
        import pandas as pd
        _fetch()
        s = pd.read_csv(os.path.join(DIR, FILES[0]), sep="\t")
        s = s.drop_duplicates("Source Name").set_index("Source Name")
        with open(os.path.join(DIR, FILES[1])) as fh:
            head = fh.readline().rstrip("\n").split("\t")[1:]
            nc = len(head)
            keep = [i for i, c in enumerate(head)
                    if c in s.index
                    and s.loc[c, "Characteristics [inferred cell type]"] == cell_type
                    and s.loc[c, "Characteristics [single cell quality]"] == "OK"]
            genes, rows = [], []
            for ln in fh:
                f = ln.rstrip("\n").split("\t")
                genes.append(f[0])
                v = np.array(f[2:2 + nc], dtype=np.float32)
                rows.append(v[keep])
        X = np.log2(np.array(rows).T + 1)
        cells = [head[i] for i in keep]
        meta = s.loc[cells]
        hb = meta["Characteristics [clinical information]"].str.extract(
            r"HbA1c ([\d.]+)%")[0].astype(float).values
        np.savez_compressed(
            out, X=X.astype(np.float32), genes=np.array(genes, dtype=str),
            donor=np.array([str(x) for x in meta["Characteristics [individual]"]]),
            t2d=np.array([str(x) == "type II diabetes mellitus"
                          for x in meta["Characteristics [disease]"]]),
            hba1c=np.asarray(hb, dtype=float),
            bmi=np.asarray(meta["Characteristics [body mass index]"], dtype=float))
    z = np.load(out, allow_pickle=False)
    X, genes = z["X"], z["genes"]
    # collapse duplicate gene symbols (RefSeq isoforms): keep the max-mean row
    order = np.argsort(-X.mean(0))
    _, first = np.unique(genes[order], return_index=True)
    idx = np.sort(order[first])
    return {"X": X[:, idx], "genes": genes[idx], "donor": z["donor"],
            "t2d": z["t2d"], "hba1c": z["hba1c"], "bmi": z["bmi"]}


# ---------------------------------------------------------------- Xin 2016
XIN_DIR = os.path.join(DATA_CACHE, "islet_xin")
XIN = {"rpkm.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81608/"
                      "suppl/GSE81608_human_islets_rpkm.txt.gz",
       "series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/"
                               "GSE81608/matrix/GSE81608_series_matrix.txt.gz",
       "Homo_sapiens.gene_info.gz": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/"
                                    "GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"}


def load_xin(cell_type="beta"):
    """Xin et al. 2016 (Cell Metab 24:608), GEO GSE81608: 1,600 SMARTer
    islet cells from 12 non-diabetic and 6 T2D donors, with the authors'
    cell-subtype calls. log2(RPKM+1); Entrez ids mapped to symbols."""
    import gzip
    out = os.path.join(XIN_DIR, f"xin_{cell_type}.npz")
    if not os.path.exists(out):
        os.makedirs(XIN_DIR, exist_ok=True)
        for f, u in XIN.items():
            fn = os.path.join(XIN_DIR, f)
            if not os.path.exists(fn):
                subprocess.run(["curl", "-sS", "-o", fn, u], check=True)
        sym = {}
        with gzip.open(os.path.join(XIN_DIR, "Homo_sapiens.gene_info.gz"), "rt") as fh:
            next(fh)
            for ln in fh:
                f = ln.split("\t")
                sym[f[1]] = f[2]
        meta = {}
        with gzip.open(os.path.join(XIN_DIR, "series_matrix.txt.gz"), "rt") as fh:
            for ln in fh:
                if ln.startswith("!Sample_characteristics_ch1"):
                    v = [x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
                    k = v[0].split(":")[0]
                    meta[k] = [x.split(": ", 1)[-1] for x in v]
        keep = [i for i, c in enumerate(meta["cell subtype"]) if c == cell_type]
        genes, rows = [], []
        with gzip.open(os.path.join(XIN_DIR, "rpkm.txt.gz"), "rt") as fh:
            next(fh)
            for ln in fh:
                f = ln.rstrip("\n").split("\t")
                genes.append(sym.get(f[0], "ENTREZ" + f[0]))
                rows.append(np.array(f[1:], dtype=np.float32)[keep])
        X = np.log2(np.array(rows).T + 1)
        np.savez_compressed(
            out, X=X.astype(np.float32), genes=np.array(genes, dtype=str),
            donor=np.array([meta["donor id"][i] for i in keep], dtype=str),
            t2d=np.array([meta["condition"][i] == "T2D" for i in keep]),
            age=np.array([float(meta["age"][i]) for i in keep]))
    z = np.load(out)
    X, genes = z["X"], z["genes"]
    order = np.argsort(-X.mean(0))
    _, first = np.unique(genes[order], return_index=True)
    idx = np.sort(order[first])
    return {"X": X[:, idx], "genes": genes[idx], "donor": z["donor"],
            "t2d": z["t2d"], "age": z["age"]}
