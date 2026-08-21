"""Benchmark on the public TESLA mirror -- real labels, real negatives.

`benchmark.py` mines IEDB for positives and uses *unlabelled* TCGA peptides as
decoys, which makes every AUC a lower bound and leaves the true negative rate
unknown. This module removes that weakness where it can be removed: the TESLA
consortium's blinded neoantigen study tested peptide-HLA pairs against patient
T cells and published both the hits and the misses, so a peptide labelled 0
here was *assayed and found negative*, not merely never tested.

Data
----
`data/tesla_deepimmuno_public.csv` -- 522 peptide-HLA pairs, 35 experimentally
immunogenic, across 6 patients, with the per-model scores published alongside
the DeepImmuno evaluation (rf/ada/cnn classify+regress, the DeepImmuno
"immunogenic score", and an IEDB immunogenicity score). It is a processed
public mirror, NOT the full 608-pMHC supplemental table, and every number
computed from it inherits that scope.

Why this is the right benchmark to argue over
---------------------------------------------
* 35 positives in 522 rows is a 6.7% base rate. Accuracy is meaningless here;
  average precision (PR-AUC) and top-N recovery are the metrics that move.
* 6 patients. Metrics are reported per patient as well as pooled, because a
  single patient with many positives can carry a pooled number on its own.
* A vaccine picks a fixed budget per patient, so `top-N recovery per patient`
  at N = 20 and N = 34 is the operationally meaningful readout.
"""

from __future__ import annotations

import os
import re
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from . import features as F
from . import fetch, presentation, selfindex
from .benchmark import auc

DATA = os.path.join(os.path.dirname(__file__), "data", "tesla_deepimmuno_public.csv")

# Score columns published with the mirror (not produced by this package).
PUBLISHED = ["rf_classify", "ada_classify", "rf_regress", "ada_regress",
             "cnn_classify", "cnn_regress", "IEDB", "immunogenic score"]


def normalize_hla(value: str) -> str:
    """'HLA-A*0201' -> 'HLA-A*02:01' (IEDB / NetMHCpan four-digit form)."""
    v = str(value).strip().upper().replace(" ", "").replace("_", "")
    if re.match(r"^HLA-[A-C]\*\d{2}:\d{2}$", v):
        return v
    m = re.match(r"^HLA-([A-C])\*?(\d{2})(\d{2})$", v)
    if m:
        return f"HLA-{m.group(1)}*{m.group(2)}:{m.group(3)}"
    return v


def load(path: str = DATA, proteome: Optional[Dict[str, str]] = None,
         indexes: Optional[Dict[int, object]] = None,
         verbose: bool = True) -> pd.DataFrame:
    """Read the mirror, normalize alleles, and recover each peptide's wild-type
    counterpart from the reference proteome so agretopicity is computable.

    The mirror ships peptide + HLA + label only. The wild-type peptide is
    recovered exactly the way `benchmark.py` mines its positives: the self
    k-mer that is one substitution away. Peptides with no such counterpart keep
    `wt_peptide = NaN` and get the neutral agretopicity value -- they are not
    dropped, and the count is reported.
    """
    d = pd.read_csv(path)
    d = d.rename(columns={"peptide": "peptide", "HLA": "allele_raw",
                          "immunogenicity": "label"})
    d["allele"] = d["allele_raw"].map(normalize_hla)
    d["length"] = d["peptide"].str.len()
    d = d[d["peptide"].str.match(r"^[ACDEFGHIKLMNPQRSTVWY]+$", na=False)].copy()

    proteome = proteome if proteome is not None else fetch.human_proteome()
    indexes = indexes or selfindex.build_indexes(proteome, sorted(set(d["length"])))
    wt = []
    for pep, L in zip(d["peptide"], d["length"]):
        hits = indexes[int(L)].one_mismatch(pep)
        wt.append(hits[0] if hits else None)
    d["wt_peptide"] = wt
    if verbose:
        n = int(d["wt_peptide"].notna().sum())
        print(f"  TESLA mirror: {len(d)} pMHC, {int(d['label'].sum())} positive, "
              f"{d['patient'].nunique()} patients; wild-type counterpart "
              f"recovered for {n}/{len(d)}")
    return d


def score(d: pd.DataFrame, weights: Dict[str, float],
          iedb_reference: Optional[Sequence[str]] = None,
          verbose: bool = True) -> pd.DataFrame:
    """Run this package's peptide-intrinsic features over the mirror.

    Expression and clonality are unavailable for a TESLA pMHC pair and are left
    at zero rather than imputed -- so what is measured here is the peptide-
    intrinsic half of the score, the same half `benchmark.py` evaluates.
    """
    d = d.copy()
    frames = []
    for (allele, L), grp in d.groupby(["allele", "length"]):
        peps = list(dict.fromkeys(list(grp["peptide"]) + list(grp["wt_peptide"].dropna())))
        pr = presentation.predict_iedb(peps, [allele], int(L), "I", progress=False)
        frames.append(pr)
        if verbose:
            print(f"  predicted {allele} L={L}: {len(grp)} pMHC", flush=True)
    P = pd.concat(frames, ignore_index=True)
    rank = dict(zip(zip(P["peptide"], P["allele"]), P["percentile_rank"]))
    d["mut_rank"] = [rank.get((p, a), np.nan) for p, a in zip(d["peptide"], d["allele"])]
    d["wt_rank"] = [rank.get((w, a), np.nan) if isinstance(w, str) else np.nan
                    for w, a in zip(d["wt_peptide"], d["allele"])]

    if iedb_reference is None:
        rows = fetch.iedb_positive_epitopes(lengths=tuple(sorted(set(d["length"]))))
        iedb_reference = [r["linear_sequence"] for r in rows if r.get("linear_sequence")]

    d["feat_presentation"] = d["mut_rank"].map(F.f_presentation)
    d["feat_agretopicity"] = [F.f_agretopicity(m, w) for m, w in zip(d["mut_rank"], d["wt_rank"])]
    d["feat_dissimilarity"] = [F.f_dissimilarity(p, w) for p, w in zip(d["peptide"], d["wt_peptide"])]
    d["feat_hydrophobicity"] = d["peptide"].map(F.f_hydrophobicity)
    d["feat_tcr_prior"] = 0.0
    for L, grp in d.groupby("length"):
        q = set(grp["peptide"])
        ref = [r for r in iedb_reference if len(r) == L and r not in q]
        d.loc[grp.index, "feat_tcr_prior"] = F.tcr_prior_scores(list(grp["peptide"]), ref)

    intrinsic = {k: v for k, v in weights.items()
                 if k in ("presentation", "agretopicity", "dissimilarity",
                          "tcr_prior", "hydrophobicity")}
    tot = sum(intrinsic.values()) or 1.0
    intrinsic = {k: v / tot for k, v in intrinsic.items()}
    d["score_composite"] = sum(w * d["feat_" + k].fillna(0) for k, w in intrinsic.items())
    no_tcr = {k: v for k, v in intrinsic.items() if k != "tcr_prior"}
    t2 = sum(no_tcr.values()) or 1.0
    d["score_composite_no_tcr"] = sum((w / t2) * d["feat_" + k].fillna(0)
                                      for k, w in no_tcr.items())
    d["score_netmhcpan_only"] = d["feat_presentation"]
    return d


def average_precision(labels: Sequence[int], scores: Sequence[float]) -> float:
    """PR-AUC by the standard interpolation-free definition (same convention the
    uploaded package uses, so the two are directly comparable)."""
    pairs = sorted(zip(np.asarray(scores, dtype=float), np.asarray(labels)),
                   key=lambda t: -t[0])
    pos = int(sum(l for _, l in pairs))
    if not pos:
        return float("nan")
    tp = 0
    acc = 0.0
    for i, (_, label) in enumerate(pairs, start=1):
        if label:
            tp += 1
            acc += tp / i
    return acc / pos


def topn_recovery(d: pd.DataFrame, score_col: str, n: int = 20) -> Dict[str, float]:
    """Per-patient: how many true positives land in a budget of N picks."""
    hits = total = 0
    per = {}
    for patient, grp in d.groupby("patient"):
        top = grp.sort_values(score_col, ascending=False).head(n)
        h, t = int(top["label"].sum()), int(grp["label"].sum())
        per[str(patient)] = h
        hits += h
        total += t
    return {"hits": hits, "positives": total,
            "recall": hits / total if total else float("nan"), "per_patient": per}


def evaluate(d: pd.DataFrame, score_cols: Optional[Sequence[str]] = None,
             ns=(20, 34)) -> pd.DataFrame:
    """One row per score: pooled AP and AUC, plus per-patient top-N recovery."""
    cols = list(score_cols or (["score_netmhcpan_only", "score_composite_no_tcr",
                                "score_composite"]
                               + [c for c in PUBLISHED if c in d.columns]
                               + [c for c in d.columns if c.startswith("feat_")]))
    base = float(d["label"].mean())
    rows = []
    for c in cols:
        if c not in d.columns or d[c].nunique() <= 1:
            continue
        rec = {"score": c,
               "AP": round(average_precision(d["label"], d[c].fillna(0)), 4),
               "AP_vs_baseline": round(average_precision(d["label"], d[c].fillna(0)) / base, 2),
               "AUC": round(auc(d["label"], d[c].fillna(0)), 3)}
        for n in ns:
            r = topn_recovery(d, c, n)
            rec[f"top{n}_hits"] = r["hits"]
            rec[f"top{n}_recall"] = round(r["recall"], 3)
        rows.append(rec)
    out = pd.DataFrame(rows).sort_values("AP", ascending=False).reset_index(drop=True)
    out.attrs["baseline_AP"] = round(base, 4)
    return out


def per_patient_table(d: pd.DataFrame, score_col: str, n: int = 20) -> pd.DataFrame:
    rows = []
    for patient, grp in d.groupby("patient"):
        top = grp.sort_values(score_col, ascending=False).head(n)
        rows.append({"patient": patient, "pMHC": len(grp),
                     "positives": int(grp["label"].sum()),
                     f"hits_in_top{n}": int(top["label"].sum())})
    return pd.DataFrame(rows)
