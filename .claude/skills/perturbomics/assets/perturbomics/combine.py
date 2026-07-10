"""Combination analysis across modalities — the payoff of one representation.

Given a *disease* (or unwanted-state) signature and a library of perturbation
signatures drawn from DIFFERENT modalities (drug profiles from the Repurposing
Hub / CMap, CRISPRko/i/a screen hits, Geneformer in-silico knockouts), find:

  1. single perturbations that REVERSE the disease signature (rank_reversers)
  2. drug + CRISPR PAIRS whose reversals are COMPLEMENTARY, i.e. together they
     cover more of the disease signature than either alone (best_combinations)

The complementarity model here is deliberately simple and transparent: a pair
is scored on how much of the disease signature the two jointly reverse, with a
penalty when they merely reverse the SAME genes (redundant) rather than
different ones (complementary). This is a screening heuristic to prioritise
wet-lab combinations, NOT a quantitative synergy prediction.

RIGOR: the reversal/coverage arithmetic is deterministic. "This pair is
synergistic in cells" is a HYPOTHESIS — validate with an actual combination
assay (e.g. Bliss/Loewe on a dose matrix).
"""

from __future__ import annotations

from itertools import combinations
from typing import Sequence

import numpy as np
import pandas as pd

from .signature import Signature
from .connectivity import weighted_connectivity_score


def rank_reversers(
    disease: Signature, library: Sequence[Signature], k: int = 50,
    weight: float = 1.0
) -> pd.DataFrame:
    """Rank library perturbations by how strongly they reverse ``disease``.

    Returns a DataFrame sorted by connectivity ascending (most negative =
    strongest reverser first), with modality and any target/MOA meta carried
    through so drug and CRISPR hits sit in one ranked table. RIGOROUS score;
    the therapeutic call is a hypothesis.
    """
    rows = []
    for sig in library:
        wtcs = weighted_connectivity_score(disease, sig, k, weight)
        rows.append({
            "name": sig.name,
            "modality": sig.modality,
            "connectivity": wtcs,
            "target": sig.meta.get("target", ""),
            "moa": sig.meta.get("moa", ""),
        })
    df = pd.DataFrame(rows).sort_values("connectivity").reset_index(drop=True)
    return df


def _reversal_vector(disease: Signature, pert: Signature,
                     genes: list[str]) -> np.ndarray:
    """Per-gene reversal: positive where ``pert`` pushes a disease gene back
    toward baseline (opposite sign, on the disease's own gene set)."""
    d = disease.align(genes).values
    p = pert.align(genes).values
    # reversal is the amount of -sign(d) motion the perturbation supplies,
    # capped at the disease magnitude (you cannot "over-correct" for credit).
    desired = -np.sign(d) * np.abs(d)          # target displacement
    supplied = p
    # component of supplied in the desired direction, clipped to [0, |d|]
    aligned = np.sign(desired) * supplied
    return np.clip(aligned, 0.0, np.abs(d))


def combination_score(
    disease: Signature, pert_a: Signature, pert_b: Signature,
    k_genes: int = 100
) -> dict:
    """Score a two-perturbation combination against a disease signature.

    Model
    -----
    * Work on the top ``k_genes`` |disease| genes (the state to fix).
    * Each perturbation supplies a per-gene reversal vector (capped at the
      disease magnitude — no credit for overshoot).
    * The COMBINATION reversal at each gene is the elementwise max of the two
      (whichever partner fixes that gene better), so covering DIFFERENT genes
      adds up while covering the SAME gene does not double-count.
    * coverage      = combined reversal / total disease magnitude in [0,1]
    * complementarity = combined coverage - max(single coverages) >= 0
      (how much the partner adds beyond the better single agent)

    Returns a dict of the coverages + a ``combo_score`` = coverage weighted by
    complementarity, which is what ``best_combinations`` ranks on.

    RIGOROUS arithmetic; SYNERGY remains a hypothesis (see module docstring).
    """
    genes = disease.top(k_genes, "abs")
    total = np.abs(disease.align(genes).values).sum()
    if total == 0:
        return {"combo_score": 0.0, "coverage": 0.0, "complementarity": 0.0,
                "cov_a": 0.0, "cov_b": 0.0}

    ra = _reversal_vector(disease, pert_a, genes)
    rb = _reversal_vector(disease, pert_b, genes)
    combo = np.maximum(ra, rb)

    cov_a = ra.sum() / total
    cov_b = rb.sum() / total
    coverage = combo.sum() / total
    complementarity = max(0.0, coverage - max(cov_a, cov_b))
    combo_score = coverage * (1.0 + complementarity)
    return {
        "combo_score": float(combo_score),
        "coverage": float(coverage),
        "complementarity": float(complementarity),
        "cov_a": float(cov_a),
        "cov_b": float(cov_b),
    }


def best_combinations(
    disease: Signature, library: Sequence[Signature], k_genes: int = 100,
    require_cross_modality: bool = False, top_n: int = 20
) -> pd.DataFrame:
    """Rank all pairs in ``library`` by combination score against ``disease``.

    Parameters
    ----------
    require_cross_modality : bool
        If True, only score pairs whose ``modality`` differs (e.g. a drug +
        a CRISPR perturbation) — the classic "small molecule plus genetic
        target" combination question.
    top_n : int
        How many top pairs to return.

    Returns a ranked DataFrame. RIGOROUS ranking; hypotheses to validate.
    """
    rows = []
    for a, b in combinations(library, 2):
        if require_cross_modality and a.modality == b.modality:
            continue
        sc = combination_score(disease, a, b, k_genes)
        rows.append({
            "pert_a": a.name, "modality_a": a.modality,
            "pert_b": b.name, "modality_b": b.modality,
            **sc,
        })
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows).sort_values(
        "combo_score", ascending=False).reset_index(drop=True)
    return df.head(top_n)
