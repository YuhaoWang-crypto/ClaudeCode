"""Sequence-identity clustering for redundancy removal and homology splitting.

A dependency-free k-mer Jaccard proxy for sequence identity, clustered by
single-linkage union-find. In production replace with mmseqs2 easy-cluster
(--min-seq-id 0.9 for dedup, 0.4 for homology split). The interface (list of
sequences -> cluster labels) is identical so the swap is local.
"""

from __future__ import annotations

import numpy as np


def _kmers(seq: str, k: int = 5) -> set[str]:
    seq = seq.upper()
    if len(seq) < k:
        return {seq} if seq else set()
    return {seq[i : i + k] for i in range(len(seq) - k + 1)}


def jaccard(a: set[str], b: set[str]) -> float:
    if not a or not b:
        return 0.0
    inter = len(a & b)
    return inter / (len(a) + len(b) - inter)


def cluster_sequences(seqs: list[str], threshold: float = 0.5, k: int = 5) -> np.ndarray:
    """Single-linkage cluster labels; two seqs linked if Jaccard >= threshold.

    threshold ~0.9 behaves like strict dedup; ~0.3-0.5 like a homology cutoff.
    """
    n = len(seqs)
    if n == 0:
        return np.array([], dtype=int)
    kms = [_kmers(s, k) for s in seqs]
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(n):
        for j in range(i + 1, n):
            if jaccard(kms[i], kms[j]) >= threshold:
                ra, rb = find(i), find(j)
                if ra != rb:
                    parent[rb] = ra

    roots: dict[int, int] = {}
    labels = np.empty(n, dtype=int)
    for i in range(n):
        r = find(i)
        labels[i] = roots.setdefault(r, len(roots))
    return labels


def dedup_representatives(seqs: list[str], threshold: float = 0.9) -> list[int]:
    """Indices of one representative per near-identical cluster (first seen)."""
    labels = cluster_sequences(seqs, threshold=threshold)
    seen: dict[int, int] = {}
    for idx, lab in enumerate(labels):
        seen.setdefault(int(lab), idx)
    return sorted(seen.values())
