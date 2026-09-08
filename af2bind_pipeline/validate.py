"""Score AF2BIND predictions against ligand-contact ground truth.

Given a holo structure, the residues contacting the crystallographic ligand are
the labels. This is the self-check that tells you whether a run is sane before
you trust it on an apo or predicted structure. Metrics are computed here rather
than imported so the package stays numpy-only.
"""

from __future__ import annotations

import numpy as np


def roc_auc(labels: np.ndarray, scores: np.ndarray) -> float:
    """AUC via the rank-sum identity; ties get averaged ranks."""
    labels = np.asarray(labels).astype(bool)
    scores = np.asarray(scores, dtype=float)
    n_pos, n_neg = int(labels.sum()), int((~labels).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    order = np.argsort(scores, kind="mergesort")
    ranks = np.empty(len(scores), dtype=float)
    sorted_scores = scores[order]
    i = 0
    while i < len(scores):
        j = i
        while j + 1 < len(scores) and sorted_scores[j + 1] == sorted_scores[i]:
            j += 1
        ranks[order[i : j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    return float((ranks[labels].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def average_precision(labels: np.ndarray, scores: np.ndarray) -> float:
    labels = np.asarray(labels).astype(bool)
    scores = np.asarray(scores, dtype=float)
    n_pos = int(labels.sum())
    if n_pos == 0:
        return float("nan")
    order = np.argsort(-scores, kind="mergesort")
    hits = labels[order]
    cum = np.cumsum(hits)
    precision_at_hit = cum[hits] / (np.flatnonzero(hits) + 1)
    return float(precision_at_hit.sum() / n_pos)


def top_n_metrics(labels: np.ndarray, scores: np.ndarray, n: int) -> dict:
    labels = np.asarray(labels).astype(bool)
    scores = np.asarray(scores, dtype=float)
    n = min(n, len(scores))
    order = np.argsort(-scores, kind="mergesort")[:n]
    tp = int(labels[order].sum())
    n_pos = int(labels.sum())
    base = n_pos / len(labels) if len(labels) else float("nan")
    precision = tp / n if n else float("nan")
    return {
        "n": n,
        "true_positives": tp,
        f"precision_at_{n}": round(precision, 4),
        f"recall_at_{n}": round(tp / n_pos, 4) if n_pos else float("nan"),
        "enrichment_over_random": round(precision / base, 3) if base else float("nan"),
    }


def evaluate(
    keys: list[tuple[str, int]],
    scores: np.ndarray,
    positives: set[tuple[str, int]],
    top_ns: tuple[int, ...] = (10, 15, 25),
) -> dict:
    """Full report for one target."""
    labels = np.array([k in positives for k in keys])
    n_pos = int(labels.sum())
    out = {
        "n_residues": len(keys),
        "n_contact_residues": n_pos,
        "positive_rate": round(n_pos / len(keys), 4) if keys else float("nan"),
        "roc_auc": round(roc_auc(labels, scores), 4),
        "average_precision": round(average_precision(labels, scores), 4),
        "top_n": {str(n): top_n_metrics(labels, scores, n) for n in top_ns},
    }
    if n_pos:
        ranks = np.argsort(np.argsort(-scores, kind="mergesort"), kind="mergesort")
        out["best_rank_of_a_true_site"] = int(ranks[labels].min()) + 1
        out["median_rank_of_true_sites"] = float(np.median(ranks[labels]) + 1)
    return out
