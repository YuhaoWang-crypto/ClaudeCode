"""
TCGA via the GDC open-access API (https://api.gdc.cancer.gov).

Gene expression: STAR-Counts quantification files (GENCODE v36), one per
aliquot, ~4 MB each. We stream every file of a project and keep only the
requested genes (TPM), so a whole project costs ~2 min and ~1 MB on disk.
Tumour AND solid-tissue-normal samples come through the same pipeline, so they
are on identical units. (The GDC `gene_expression/values` endpoint is faster
but returns one primary-tumour column per case and no normals.)

Clinical: vital status, days to death / last follow-up, AJCC stage, grade,
age at diagnosis, from the GDC `cases` endpoint.

Cache: figures/_data_cache/tcga/ (git-ignored).
"""
import json
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from netjac import DATA_CACHE

API = "https://api.gdc.cancer.gov"
DIR = os.path.join(DATA_CACHE, "tcga")


def _get_json(endpoint, params, tries=4):
    import time
    import urllib.parse
    url = f"{API}/{endpoint}?" + urllib.parse.urlencode(params)
    for k in range(tries):
        r = subprocess.run(["curl", "-sf", "--max-time", "120", url],
                           capture_output=True)
        if r.returncode == 0:
            return json.loads(r.stdout)
        time.sleep(2 ** (k + 1))
    raise RuntimeError(f"GDC {endpoint} failed")


def list_expression_files(project):
    filt = {"op": "and", "content": [
        {"op": "in", "content": {"field": "cases.project.project_id",
                                 "value": [project]}},
        {"op": "in", "content": {"field": "data_type",
                                 "value": ["Gene Expression Quantification"]}},
        {"op": "in", "content": {"field": "analysis.workflow_type",
                                 "value": ["STAR - Counts"]}}]}
    d = _get_json("files", {
        "filters": json.dumps(filt), "size": 5000,
        "fields": "file_id,cases.case_id,cases.submitter_id,"
                  "cases.samples.sample_type,cases.samples.submitter_id"})
    rows = []
    for h in d["data"]["hits"]:
        c = h["cases"][0]
        s = c["samples"][0]
        rows.append({"file_id": h["file_id"], "case_id": c["case_id"],
                     "patient": c["submitter_id"],
                     "sample": s["submitter_id"],
                     "sample_type": s["sample_type"]})
    return rows


def _fetch_one(file_id, genes, tries=4):
    import time
    for k in range(tries):
        r = subprocess.run(["curl", "-sf", "--max-time", "180",
                            f"{API}/data/{file_id}"], capture_output=True)
        if r.returncode == 0 and r.stdout:
            out = {}
            for line in r.stdout.decode().splitlines():
                f = line.split("\t")
                if len(f) > 6 and f[1] in genes:
                    out[f[1]] = out.get(f[1], 0.0) + float(f[6])   # TPM
            return out
        time.sleep(2 ** (k + 1))
    raise RuntimeError(f"download failed: {file_id}")


def expression(project, genes, workers=8, force=False):
    """DataFrame: rows = aliquots (index = sample barcode), columns = genes
    (TPM) + metadata columns patient, case_id, sample_type."""
    import pandas as pd
    os.makedirs(DIR, exist_ok=True)
    path = os.path.join(DIR, f"{project}_expr.csv")
    genes = sorted(set(genes))
    if os.path.exists(path) and not force:
        df = pd.read_csv(path, index_col=0)
        if set(genes) <= set(df.columns):
            return df
    files = list_expression_files(project)
    gset = set(genes)
    with ThreadPoolExecutor(workers) as ex:
        vals = list(ex.map(lambda f: _fetch_one(f["file_id"], gset), files))
    rows = []
    for f, v in zip(files, vals):
        rows.append({**{g: v.get(g, np.nan) for g in genes},
                     "patient": f["patient"], "case_id": f["case_id"],
                     "sample_type": f["sample_type"], "sample": f["sample"]})
    df = pd.DataFrame(rows).set_index("sample")
    df.to_csv(path)
    return df


def clinical(project, force=False):
    import pandas as pd
    os.makedirs(DIR, exist_ok=True)
    path = os.path.join(DIR, f"{project}_clinical.csv")
    if os.path.exists(path) and not force:
        return pd.read_csv(path)
    filt = {"op": "in", "content": {"field": "project.project_id",
                                    "value": [project]}}
    d = _get_json("cases", {
        "filters": json.dumps(filt), "size": 5000,
        "fields": "submitter_id,demographic.vital_status,"
                  "demographic.days_to_death,diagnoses.days_to_last_follow_up,"
                  "diagnoses.ajcc_pathologic_stage,diagnoses.tumor_grade,"
                  "diagnoses.age_at_diagnosis,diagnoses.diagnosis_is_primary_disease"})
    rows = []
    for h in d["data"]["hits"]:
        dem = h.get("demographic", {}) or {}
        dx = h.get("diagnoses", []) or [{}]
        prim = [x for x in dx if x.get("diagnosis_is_primary_disease")] or dx
        x = prim[0]
        fu = [y.get("days_to_last_follow_up") for y in dx
              if y.get("days_to_last_follow_up") is not None]
        rows.append({"patient": h["submitter_id"],
                     "vital_status": dem.get("vital_status"),
                     "days_to_death": dem.get("days_to_death"),
                     "days_to_last_follow_up": max(fu) if fu else None,
                     "stage": x.get("ajcc_pathologic_stage"),
                     "grade": x.get("tumor_grade"),
                     "age_days": x.get("age_at_diagnosis")})
    df = pd.DataFrame(rows)
    dead = df["vital_status"] == "Dead"
    df["event"] = dead.astype(int)
    df["time"] = np.where(dead, df["days_to_death"],
                          df["days_to_last_follow_up"]).astype(float)
    df.to_csv(path, index=False)
    return df
