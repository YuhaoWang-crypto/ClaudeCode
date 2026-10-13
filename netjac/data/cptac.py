"""
CPTAC global proteomics via the NCI Proteomic Data Commons GraphQL API
(https://pdc.cancer.gov/graphql, open access, `acceptDUA: true`).

`quantDataMatrix(data_type: "log2_ratio")` returns the harmonised CDAP
gene-level protein abundance (log2 ratio to the TMT pooled reference) for a
whole study: ~10-12k proteins x all aliquots. `biospecimenPerStudy` maps each
aliquot to Primary Tumor / Solid Tissue Normal. Tumour and normal tissue are
run in the same TMT plexes against the same reference, so the units match.

Cache: figures/_data_cache/cptac/<study>.npz (git-ignored).
"""
import json
import os
import subprocess
import time

import numpy as np

from netjac import DATA_CACHE

API = "https://pdc.cancer.gov/graphql"
DIR = os.path.join(DATA_CACHE, "cptac")

# CPTAC discovery cohorts that contain paired solid-tissue normals
STUDIES = {
    "ccRCC": "PDC000127",
    "LUAD": "PDC000153",
    "LSCC": "PDC000234",
    "HNSCC": "PDC000221",
    "PDAC": "PDC000270",
    "UCEC": "PDC000125",
    "HCC": "PDC000198",
    "COAD": "PDC000116",
}


def _gql(query, tries=4):
    body = json.dumps({"query": query})
    for k in range(tries):
        r = subprocess.run(["curl", "-sf", "--max-time", "600", "-X", "POST",
                            API, "-H", "Content-Type: application/json",
                            "-d", body], capture_output=True)
        if r.returncode == 0 and r.stdout:
            d = json.loads(r.stdout)
            if d.get("data"):
                return d["data"]
        time.sleep(2 ** (k + 1))
    raise RuntimeError("PDC query failed: " + query[:80])


# matching phosphoproteome studies (same cases, same TMT plexes)
PHOSPHO = {
    "ccRCC": "PDC000128",
    "LUAD": "PDC000149",
    "LSCC": "PDC000232",
    "HNSCC": "PDC000222",
    "PDAC": "PDC000271",
    "UCEC": "PDC000126",
    "HCC": "PDC000199",
    "COAD": "PDC000117",
}


def load(name, min_present=0.8, kind="proteome"):
    """Return dict(X: samples x proteins log2 ratio (NaN allowed), genes,
    aliquots, cases, sample_type) for one indication. Proteins missing in more
    than (1-min_present) of samples are dropped.
    kind="phospho": rows are phosphosites 'NP_xxx.v:s123' (CDAP site-level
    log2 ratios, not normalised to protein); see `load_phospho`."""
    sid = (PHOSPHO if kind == "phospho" else STUDIES).get(name, name)
    os.makedirs(DIR, exist_ok=True)
    fn = os.path.join(DIR, f"{sid}.npz")
    if not os.path.exists(fn):
        m = _gql(f'{{ quantDataMatrix(pdc_study_id: "{sid}", '
                 f'data_type: "log2_ratio" acceptDUA: true) }}')["quantDataMatrix"]
        bio = _gql(f'{{ biospecimenPerStudy(pdc_study_id: "{sid}" acceptDUA: true)'
                   ' { aliquot_id aliquot_submitter_id sample_type '
                   'case_submitter_id } }')["biospecimenPerStudy"]
        st = {b["aliquot_id"]: b for b in bio}
        st.update({b["aliquot_submitter_id"]: b for b in bio})
        head = m[0][1:]
        genes = np.array([r[0] for r in m[1:]])
        V = np.array([[np.nan if v in ("", None, "NA", "null") else float(v)
                       for v in r[1:]] for r in m[1:]], dtype=np.float32)
        al, typ, case = [], [], []
        for h in head:
            a, _, sub = h.partition(":")
            b = st.get(a) or st.get(sub) or {}
            al.append(sub or a)
            typ.append(b.get("sample_type", "unknown"))
            case.append(b.get("case_submitter_id", ""))
        np.savez_compressed(fn, V=V, genes=genes, aliquots=np.array(al),
                            sample_type=np.array(typ), cases=np.array(case))
    z = np.load(fn, allow_pickle=False)
    V, genes = z["V"].T, z["genes"]
    typ = z["sample_type"].astype(object)
    case = z["cases"].astype(object)
    # HBV-HCC (PDC000198) labels aliquots T<n> (tumour) / P<n> (paired
    # para-tumour liver) without a biospecimen sample type
    for i, a in enumerate(z["aliquots"]):
        if typ[i] == "unknown" and len(a) > 1 and a[0] in "TP" and a[1:].isdigit():
            typ[i] = "Primary Tumor" if a[0] == "T" else "Solid Tissue Normal"
            case[i] = a[1:]
    keep_s = np.isin(typ, ["Primary Tumor", "Solid Tissue Normal"])
    V, typ, al, case = V[keep_s], typ[keep_s], z["aliquots"][keep_s], case[keep_s]
    # duplicate aliquots of the same case/type: keep the first
    seen, idx = set(), []
    for i, (c, t) in enumerate(zip(case, typ)):
        if (c, t) not in seen:
            seen.add((c, t))
            idx.append(i)
    V, typ, al, case = V[idx], typ[idx], al[idx], case[idx]
    pres = np.isfinite(V).mean(0)
    keep_g = pres >= min_present
    return {"X": V[:, keep_g], "genes": genes[keep_g], "aliquots": al,
            "cases": case, "tumour": typ == "Primary Tumor", "name": name,
            "study": sid}


_NP2SYM = None


def _np2sym():
    """RefSeq protein accession (no version) -> gene symbol, from NCBI
    gene2refseq (human rows, streamed once) + Homo_sapiens.gene_info."""
    global _NP2SYM
    if _NP2SYM is None:
        import gzip
        fn = os.path.join(DIR, "np2geneid.tsv")
        if not os.path.exists(fn):
            os.makedirs(DIR, exist_ok=True)
            subprocess.run(
                "curl -s https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz | zcat | "
                "awk -F'\\t' '$1==9606 && $6!=\"-\" {split($6,a,\".\"); print a[1]\"\\t\"$2}' "
                f"| sort -u > {fn}", shell=True, check=True)
        info = os.path.join(DATA_CACHE, "islet_xin", "Homo_sapiens.gene_info.gz")
        if not os.path.exists(info):
            os.makedirs(os.path.dirname(info), exist_ok=True)
            subprocess.run(["curl", "-sS", "-o", info,
                            "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/"
                            "Mammalia/Homo_sapiens.gene_info.gz"], check=True)
        sym = {}
        with gzip.open(info, "rt") as fh:
            next(fh)
            for ln in fh:
                f = ln.split("\t")
                sym[f[1]] = f[2]
        _NP2SYM = {}
        for ln in open(fn):
            a, g = ln.rstrip("\n").split("\t")
            if g in sym:
                _NP2SYM[a] = sym[g]
    return _NP2SYM


def load_phospho(name, min_present=0.7):
    """Phosphosite matrix for one indication with sites renamed
    GENE_S123 (single sites) or GENE_S123_S125 (multiply phosphorylated
    peptides). Returns the `load` dict plus `site_gene`, `single` (bool)."""
    import re
    d = load(name, min_present=min_present, kind="phospho")
    m = _np2sym()
    names, gene, single = [], [], []
    for r in d["genes"]:
        acc, _, sites = r.partition(":")
        g = m.get(acc.split(".")[0], acc)
        ss = re.findall(r"([sty])(\d+)", sites)
        names.append(g + "_" + "_".join(f"{a.upper()}{b}" for a, b in ss))
        gene.append(g)
        single.append(len(ss) == 1)
    d["site"] = np.array(names)
    d["site_gene"] = np.array(gene)
    d["single"] = np.array(single)
    return d


def impute(X):
    """Per-protein median imputation (after the presence filter this touches
    < 20 % of each column); returns a copy."""
    X = np.array(X, dtype=float)
    med = np.nanmedian(X, axis=0)
    i, j = np.where(~np.isfinite(X))
    X[i, j] = med[j]
    return X
