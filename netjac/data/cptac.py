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


def load(name, min_present=0.8):
    """Return dict(X: samples x proteins log2 ratio (NaN allowed), genes,
    aliquots, cases, sample_type) for one indication. Proteins missing in more
    than (1-min_present) of samples are dropped."""
    sid = STUDIES.get(name, name)
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


def impute(X):
    """Per-protein median imputation (after the presence filter this touches
    < 20 % of each column); returns a copy."""
    X = np.array(X, dtype=float)
    med = np.nanmedian(X, axis=0)
    i, j = np.where(~np.isfinite(X))
    X[i, j] = med[j]
    return X
