"""Connectivity scoring — how similar are two perturbation signatures?

This is the mathematical bridge that lets a drug signature, a CRISPR signature
and a disease signature all be compared on the same axis. The method is the
weighted connectivity score (WTCS) of Subramanian et al., Cell 2017 (the
"Next-Generation Connectivity Map"), built on a weighted running-sum
enrichment score (the GSEA statistic).

Interpretation convention used throughout:
    +1  the two perturbations do the SAME thing (mimic)
    -1  one REVERSES the other (a candidate corrector / therapy)
     0  unrelated

RIGOR: ``enrichment_score`` and ``weighted_connectivity_score`` are exact
deterministic statistics. Turning a score into "these share a mechanism" or
"this drug will treat this disease" is a HYPOTHESIS that needs orthogonal /
experimental validation.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd

from .signature import Signature


def enrichment_score(
    ranked: pd.Series, tag_set: Sequence[str], weight: float = 1.0
) -> float:
    """Weighted running-sum enrichment score (the GSEA ES). RIGOROUS.

    Parameters
    ----------
    ranked : Series
        The reference signature's scores, sorted high -> low (index = genes).
    tag_set : sequence of gene ids
        The query gene set (e.g. the up-tags of a drug signature).
    weight : float
        Exponent on |score| when walking the hit. ``weight=0`` recovers the
        classic Kolmogorov-Smirnov statistic; ``weight=1`` (default) is the
        weighted ES used by CMap.

    Returns
    -------
    float in [-1, 1] — the max-deviation running-sum enrichment score.
    """
    genes = list(ranked.index)
    hits = set(g for g in tag_set if g in ranked.index)
    n = len(genes)
    n_hit = len(hits)
    if n_hit == 0 or n_hit == n:
        return 0.0

    is_hit = np.array([g in hits for g in genes])
    corr = np.abs(ranked.values) ** weight
    hit_norm = corr[is_hit].sum()
    if hit_norm == 0:
        # degenerate: all hit weights zero -> unweighted
        corr = np.ones(n)
        hit_norm = float(n_hit)

    p_hit = np.where(is_hit, corr / hit_norm, 0.0)
    p_miss = np.where(is_hit, 0.0, 1.0 / (n - n_hit))
    running = np.cumsum(p_hit - p_miss)
    # signed max-deviation (ES): the extreme the walk reaches
    hi, lo = running.max(), running.min()
    return float(hi if hi >= -lo else lo)


def weighted_connectivity_score(
    query: Signature, reference: Signature, k: int = 50, weight: float = 1.0
) -> float:
    """WTCS between a query and a reference signature. RIGOROUS.

    Follows Subramanian 2017: derive up/down tag sets from the query, score
    each against the reference's ranked list, and combine

        WTCS = (ES_up - ES_down) / 2      if sign(ES_up) != sign(ES_down)
             = 0                          otherwise

    The zero-out rule keeps the score meaningful only when the query's up and
    down tags land at opposite ends of the reference — a coherent connection.

    Returns a value in [-1, 1]; positive = mimic, negative = reversal.
    """
    up, down = query.query_sets(k)
    ranked = reference.ranked_index()
    es_up = enrichment_score(ranked, up, weight)
    es_down = enrichment_score(ranked, down, weight)
    if np.sign(es_up) == np.sign(es_down):
        return 0.0
    return float((es_up - es_down) / 2.0)


def connectivity_matrix(
    signatures: Sequence[Signature], k: int = 50, weight: float = 1.0,
    symmetric: bool = True
) -> pd.DataFrame:
    """All-pairs WTCS. RIGOROUS.

    WTCS is not symmetric in general (query tags vs reference ranks). With
    ``symmetric=True`` (default) we report the mean of both directions, which
    is the standard way to get a usable distance for clustering.
    """
    names = [s.name for s in signatures]
    n = len(signatures)
    m = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                m[i, j] = 1.0
                continue
            m[i, j] = weighted_connectivity_score(
                signatures[i], signatures[j], k, weight)
    if symmetric:
        m = (m + m.T) / 2.0
    return pd.DataFrame(m, index=names, columns=names)


def normalized_connectivity(
    raw: pd.DataFrame, groups: dict[str, str] | None = None
) -> pd.DataFrame:
    """Normalize raw WTCS within perturbation-type/cell-line groups -> NCS.

    CMap normalizes each connection by the mean of its |scores| within the
    signature's own class (cell line x perturbation type), separately for the
    positive and negative sides, so that scores are comparable across a
    heterogeneous reference collection.

    ``groups`` maps a column (reference) name -> its group key. If omitted, one
    global group is used. RIGOROUS (deterministic rescaling); the resulting NCS
    is what you threshold, but the biological call remains a hypothesis.
    """
    if groups is None:
        groups = {c: "all" for c in raw.columns}
    ncs = raw.copy().astype(float)
    for grp in set(groups.values()):
        cols = [c for c in raw.columns if groups.get(c) == grp]
        block = raw[cols].values
        pos = block[block > 0]
        neg = block[block < 0]
        mu_pos = np.abs(pos).mean() if pos.size else 1.0
        mu_neg = np.abs(neg).mean() if neg.size else 1.0
        b = block.copy()
        b[block > 0] = block[block > 0] / (mu_pos or 1.0)
        b[block < 0] = block[block < 0] / (mu_neg or 1.0)
        ncs[cols] = b
    return ncs
