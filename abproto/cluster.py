"""Backbone RMSD clustering to maximise pose diversity (plan step 5).

The paper does average-linkage agglomerative clustering at a 3 A backbone-RMSD
cutoff, then picks one representative per cluster. scipy gives average linkage;
without it we fall back to single-linkage via union-find (connected components
of the <cutoff graph), which is dependency-free and adequate for the scaffold.
"""

from __future__ import annotations

from typing import Callable

import numpy as np

from .geometry import kabsch_rmsd
from .types import Design


def pairwise_rmsd(designs: list[Design]) -> np.ndarray:
    """Symmetric (n, n) antibody-CA RMSD matrix (Kabsch-aligned)."""
    n = len(designs)
    cas = [d.antibody_ca() for d in designs]
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            if cas[i].shape == cas[j].shape:
                r = kabsch_rmsd(cas[i], cas[j])
            else:
                r = np.inf   # different length => different topology
            D[i, j] = D[j, i] = r
    return D


def cluster(designs: list[Design], cutoff: float = 3.0) -> np.ndarray:
    """Return a cluster-id array (len n). Uses scipy average linkage if present,
    else single-linkage union-find."""
    n = len(designs)
    if n == 0:
        return np.array([], dtype=int)
    if n == 1:
        return np.array([0])

    D = pairwise_rmsd(designs)
    try:
        return _agglomerative_average(D, cutoff)
    except ImportError:
        return _single_linkage(D, cutoff)


def _agglomerative_average(D: np.ndarray, cutoff: float) -> np.ndarray:
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform

    finite = np.where(np.isinf(D), D[np.isfinite(D)].max() * 10 if np.isfinite(D).any() else 1e6, D)
    condensed = squareform(finite, checks=False)
    Z = linkage(condensed, method="average")
    labels = fcluster(Z, t=cutoff, criterion="distance")
    return labels.astype(int) - 1


def _single_linkage(D: np.ndarray, cutoff: float) -> np.ndarray:
    n = D.shape[0]
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for i in range(n):
        for j in range(i + 1, n):
            if D[i, j] <= cutoff:
                union(i, j)

    roots = {}
    labels = np.empty(n, dtype=int)
    for i in range(n):
        r = find(i)
        if r not in roots:
            roots[r] = len(roots)
        labels[i] = roots[r]
    return labels


def pick_representatives(designs: list[Design], labels: np.ndarray,
                         score: Callable[[int], float]) -> list[int]:
    """One representative index per cluster: the highest-scoring member."""
    reps: dict[int, int] = {}
    for idx, lab in enumerate(labels):
        if lab not in reps or score(idx) > score(reps[lab]):
            reps[lab] = idx
    # Return representatives ordered by descending score.
    return sorted(reps.values(), key=lambda i: score(i), reverse=True)
