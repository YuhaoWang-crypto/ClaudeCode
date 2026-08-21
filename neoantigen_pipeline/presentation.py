"""Step 5a: HLA presentation prediction.

Primary backend is the IEDB cloud REST service (NetMHCpan-4.1 EL for class I,
NetMHCIIpan for class II) -- open, no licence, no local install. An optional
MHCflurry backend is used automatically when the package is importable, which
removes the network round-trip for large runs.

Everything is cached per (method, allele, length, peptide-batch) on disk, so a
re-run of the pipeline costs nothing.
"""

from __future__ import annotations

import hashlib
import io
import json
import os
import time
from typing import Dict, Iterable, List, Optional, Sequence

import pandas as pd
import requests

from .config import CACHE_DIR

IEDB_MHCI = "https://tools-cluster-interface.iedb.org/tools_api/mhci/"
IEDB_MHCII = "https://tools-cluster-interface.iedb.org/tools_api/mhcii/"

RANK_STRONG_I, RANK_WEAK_I = 0.5, 2.0
RANK_STRONG_II, RANK_WEAK_II = 2.0, 10.0


def classify(rank: Optional[float], mhc_class: str = "I") -> str:
    if rank is None or pd.isna(rank):
        return "none"
    s, w = (RANK_STRONG_I, RANK_WEAK_I) if mhc_class == "I" else (RANK_STRONG_II, RANK_WEAK_II)
    if rank <= s:
        return "strong"
    if rank <= w:
        return "weak"
    return "none"


def _batch_key(method: str, alleles: Sequence[str], length: int,
               peptides: Sequence[str]) -> str:
    h = hashlib.sha1()
    h.update(method.encode())
    h.update("|".join(sorted(alleles)).encode())
    h.update(str(length).encode())
    for p in peptides:
        h.update(p.encode())
    return f"pred_{method}_{length}_{h.hexdigest()[:16]}.json"


def _iedb_call(url: str, data: dict, tries: int = 4, timeout: int = 1200) -> str:
    last = None
    for i in range(tries):
        try:
            r = requests.post(url, data=data, timeout=timeout)
            r.raise_for_status()
            if r.text.strip().lower().startswith(("error", "<!doctype", "<html")):
                raise RuntimeError(r.text[:300])
            return r.text
        except Exception as exc:            # noqa: BLE001 - network retry
            last = exc
            time.sleep(3 * (i + 1))
    raise RuntimeError(f"IEDB call failed after {tries} tries: {last}")


def predict_iedb(peptides: Sequence[str], alleles: Sequence[str], length: int,
                 mhc_class: str = "I", method: Optional[str] = None,
                 batch_size: int = 500, cache: bool = True,
                 progress: bool = True) -> pd.DataFrame:
    """-> DataFrame[peptide, allele, score, percentile_rank, binder]."""
    method = method or ("netmhcpan_el" if mhc_class == "I" else "netmhciipan_el")
    url = IEDB_MHCI if mhc_class == "I" else IEDB_MHCII
    peptides = [p for p in dict.fromkeys(peptides) if p and len(p) == length
                and set(p) <= set("ACDEFGHIKLMNPQRSTVWY")]
    if not peptides or not alleles:
        return pd.DataFrame(columns=["peptide", "allele", "score", "percentile_rank", "binder"])

    os.makedirs(CACHE_DIR, exist_ok=True)
    frames = []
    n_batches = (len(peptides) + batch_size - 1) // batch_size
    for bi in range(n_batches):
        chunk = peptides[bi * batch_size:(bi + 1) * batch_size]
        key = os.path.join(CACHE_DIR, _batch_key(method, alleles, length, chunk))
        if cache and os.path.exists(key):
            frames.append(pd.DataFrame(json.load(open(key))))
            continue
        fasta = "\n".join(f">p{i}\n{p}" for i, p in enumerate(chunk))
        data = {"method": method, "sequence_text": fasta,
                "allele": ",".join(alleles),
                "length": ",".join([str(length)] * len(alleles))}
        text = _iedb_call(url, data)
        raw = pd.read_csv(io.StringIO(text), sep="\t")
        raw.columns = [c.strip().lower().replace(" ", "_") for c in raw.columns]
        pep_col = "peptide" if "peptide" in raw.columns else "sequence"
        rank_col = next((c for c in ("percentile_rank", "rank", "adjusted_rank")
                         if c in raw.columns), None)
        score_col = next((c for c in ("score", "ic50", "affinity")
                          if c in raw.columns), None)
        out = pd.DataFrame({
            "peptide": raw[pep_col].astype(str),
            "allele": raw["allele"].astype(str),
            "score": pd.to_numeric(raw[score_col], errors="coerce") if score_col else None,
            "percentile_rank": pd.to_numeric(raw[rank_col], errors="coerce") if rank_col else None,
        })
        out = out[out["peptide"].isin(set(chunk))]
        if cache:
            json.dump(out.to_dict("list"), open(key, "w"))
        frames.append(out)
        if progress:
            print(f"    IEDB {method} L={length}: batch {bi+1}/{n_batches} "
                  f"({len(chunk)} peptides x {len(alleles)} alleles)", flush=True)

    res = pd.concat(frames, ignore_index=True)
    res["binder"] = res["percentile_rank"].map(lambda r: classify(r, mhc_class))
    return res


def predict_mhcflurry(peptides: Sequence[str], alleles: Sequence[str]) -> Optional[pd.DataFrame]:
    """Optional local backend. Returns None when mhcflurry is not installed."""
    try:
        from mhcflurry import Class1PresentationPredictor
    except Exception:                       # noqa: BLE001 - optional dependency
        return None
    pred = Class1PresentationPredictor.load()
    df = pred.predict(peptides=list(dict.fromkeys(peptides)),
                      alleles={a: [a] for a in alleles}, verbose=0)
    return pd.DataFrame({
        "peptide": df["peptide"],
        "allele": df["best_allele"],
        "score": df["presentation_score"],
        "percentile_rank": df["presentation_percentile"],
        "binder": df["presentation_percentile"].map(lambda r: classify(r, "I")),
    })


def predict_all(pep_table: pd.DataFrame, alleles_i: Sequence[str],
                alleles_ii: Sequence[str] = (), backend: str = "iedb",
                batch_size: int = 500, include_wt: bool = True) -> pd.DataFrame:
    """Run every (length, class) group needed by a peptide table.

    Returns a long table [peptide, allele, mhc_class, percentile_rank, score, binder]
    covering both mutant and wild-type peptides (the WT calls are what
    agretopicity is computed from).
    """
    frames = []
    for mhc_class, alleles in (("I", list(alleles_i)), ("II", list(alleles_ii))):
        if not alleles:
            continue
        sub = pep_table[pep_table["mhc_class"] == mhc_class]
        for L, grp in sub.groupby("length"):
            peps = list(grp["mut_peptide"].dropna().unique())
            if include_wt:
                peps += list(grp["wt_peptide"].dropna().unique())
            peps = [p for p in dict.fromkeys(peps) if len(p) == L]
            if not peps:
                continue
            got = None
            if backend == "mhcflurry" and mhc_class == "I":
                got = predict_mhcflurry(peps, alleles)
            if got is None:
                got = predict_iedb(peps, alleles, int(L), mhc_class,
                                   batch_size=batch_size)
            got["mhc_class"] = mhc_class
            frames.append(got)
    if not frames:
        return pd.DataFrame(columns=["peptide", "allele", "mhc_class",
                                     "score", "percentile_rank", "binder"])
    return pd.concat(frames, ignore_index=True)


def join_predictions(pep_table: pd.DataFrame, preds: pd.DataFrame) -> pd.DataFrame:
    """Attach mutant and wild-type ranks to each (peptide, allele) row."""
    p = preds.rename(columns={"peptide": "mut_peptide",
                              "percentile_rank": "mut_rank",
                              "score": "mut_score"})
    merged = pep_table.merge(p[["mut_peptide", "allele", "mhc_class",
                                "mut_rank", "mut_score", "binder"]],
                             on=["mut_peptide", "mhc_class"], how="inner")
    w = preds.rename(columns={"peptide": "wt_peptide",
                              "percentile_rank": "wt_rank",
                              "score": "wt_score"})
    merged = merged.merge(w[["wt_peptide", "allele", "mhc_class", "wt_rank", "wt_score"]],
                          on=["wt_peptide", "allele", "mhc_class"], how="left")
    return merged
