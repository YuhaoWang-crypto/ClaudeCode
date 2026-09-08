"""Group high-scoring residues into spatially distinct candidate pockets.

AF2BIND scores residues independently, so a raw ranked list mixes together
residues from several different sites. Single-linkage clustering on side-chain
centroids turns that list into "here are N pockets, ranked", which is what a
downstream docking or design decision actually needs.

This clustering step is an addition to the published method, not part of it.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .structure import Residue

BACKBONE = {"N", "CA", "C", "O", "OXT"}


@dataclass
class Pocket:
    rank: int
    residues: list[tuple[str, int, str]]  # (chain, resi, resn)
    scores: list[float]
    centroid: tuple[float, float, float]

    @property
    def size(self) -> int:
        return len(self.residues)

    @property
    def score_sum(self) -> float:
        return float(np.sum(self.scores))

    @property
    def score_mean(self) -> float:
        return float(np.mean(self.scores))

    @property
    def score_max(self) -> float:
        return float(np.max(self.scores))

    def pymol_selection(self, name: str | None = None) -> str:
        name = name or f"pocket_{self.rank}"
        by_chain: dict[str, list[int]] = {}
        for c, r, _ in self.residues:
            by_chain.setdefault(c, []).append(r)
        parts = [
            f"(chain {c} and resi {'+'.join(str(r) for r in sorted(v))})"
            for c, v in sorted(by_chain.items())
        ]
        return f"select {name}, " + " or ".join(parts)

    def as_dict(self) -> dict:
        return {
            "rank": self.rank,
            "size": self.size,
            "score_sum": round(self.score_sum, 4),
            "score_mean": round(self.score_mean, 4),
            "score_max": round(self.score_max, 4),
            "centroid": [round(v, 2) for v in self.centroid],
            "residues": [
                {"chain": c, "resi": r, "resn": n, "p_bind": round(s, 4)}
                for (c, r, n), s in zip(self.residues, self.scores)
            ],
            "pymol": self.pymol_selection(),
        }


def sidechain_centroid(res: Residue) -> np.ndarray:
    """Centroid of side-chain heavy atoms, falling back to CA for GLY/missing."""
    names = [
        n for n in res.atoms
        if n not in BACKBONE and not n.startswith("H")
    ]
    if not names:
        return res.atoms.get("CA", next(iter(res.atoms.values())))
    return np.mean([res.atoms[n] for n in names], axis=0)


def cluster_pockets(
    residues: list[Residue],
    scores: dict[tuple[str, int], float],
    threshold: float = 0.5,
    top_n: int | None = None,
    link_cutoff: float = 10.0,
    min_size: int = 3,
) -> list[Pocket]:
    """Single-linkage cluster the selected residues into pockets.

    Selection is `p_bind >= threshold`, or the `top_n` highest-scoring residues
    if `top_n` is given (which is the more robust choice, because absolute
    p(bind) is not calibrated across targets).
    """
    by_key = {r.key: r for r in residues}
    ranked = sorted(
        (k for k in scores if k in by_key), key=lambda k: -scores[k]
    )
    if top_n is not None:
        selected = ranked[:top_n]
    else:
        selected = [k for k in ranked if scores[k] >= threshold]
    if not selected:
        return []

    pts = np.stack([sidechain_centroid(by_key[k]) for k in selected])
    n = len(selected)

    # union-find over pairs closer than link_cutoff
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[rj] = ri

    d2 = ((pts[:, None, :] - pts[None, :, :]) ** 2).sum(-1)
    ii, jj = np.where(d2 <= link_cutoff * link_cutoff)
    for i, j in zip(ii, jj):
        if i < j:
            union(int(i), int(j))

    groups: dict[int, list[int]] = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)

    pockets = []
    for members in groups.values():
        if len(members) < min_size:
            continue
        keys = [selected[i] for i in members]
        pockets.append(
            {
                "keys": keys,
                "scores": [scores[k] for k in keys],
                "centroid": pts[members].mean(0),
            }
        )
    pockets.sort(key=lambda p: -float(np.sum(p["scores"])))

    out = []
    for rank, p in enumerate(pockets, start=1):
        order = np.argsort([-s for s in p["scores"]])
        keys = [p["keys"][i] for i in order]
        out.append(
            Pocket(
                rank=rank,
                residues=[(k[0], k[1], by_key[k].one_letter) for k in keys],
                scores=[p["scores"][i] for i in order],
                centroid=tuple(float(v) for v in p["centroid"]),
            )
        )
    return out
