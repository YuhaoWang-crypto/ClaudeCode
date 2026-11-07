"""Small statistics helpers (kept scipy-free so the package needs only
numpy + biopython + matplotlib)."""

from __future__ import annotations

import numpy as np

__all__ = ["rankdata", "spearman", "pearson"]


def rankdata(x) -> np.ndarray:
    """Average ranks, ties shared (matches scipy.stats.rankdata default)."""
    a = np.asarray(x, dtype=np.float64)
    order = np.argsort(a, kind="mergesort")
    ranks = np.empty(len(a), dtype=np.float64)
    ranks[order] = np.arange(1, len(a) + 1, dtype=np.float64)
    # average tied ranks
    sorted_a = a[order]
    i = 0
    while i < len(a):
        j = i
        while j + 1 < len(a) and sorted_a[j + 1] == sorted_a[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = ranks[order[i:j + 1]].mean()
        i = j + 1
    return ranks


def pearson(x, y) -> float:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if len(x) < 2 or x.std() == 0 or y.std() == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def spearman(x, y) -> float:
    """Spearman rank correlation."""
    return pearson(rankdata(x), rankdata(y))
