"""
design_his_sites.py
-------------------
Turn ranked candidate regions into a small His-pair / His-triplet library for a
Zn(II)/Ni(II)/Co(II) switch.

We propose grafting His at positions whose backbone geometry can plausibly
present two (bidentate) or three (facial) imidazole nitrogens to a single metal:

  - Cα-Cα distance in a window (defaults: pair 5.5-9.0 A, triplet 5.0-9.5 A).
  - Side chains point the SAME way (toward a shared metal): the Cα->Cβ vectors
    should converge, i.e. each Cβ is closer to the partner than its Cα is.
  - Solvent exposure proxy: few CA neighbours within 10 A (buried positions can't
    host a surface metal site and shouldn't be mutated to His).

These are COARSE backbone heuristics to prioritise a wet-lab library, NOT a
substitute for rotamer-level metal-site modelling (Rosetta, MIB, HADDOCK metal,
or manual Zn-placement). Every suggestion should be checked with a rotamer build
before synthesis.

Direction of regulation (from the source paper) is geometry-dependent: the same
metal can ACTIVATE or INHIBIT depending on whether the clamp rigidifies the
catalytically-productive or the inactive conformer. So the library is meant to be
screened for kcat / Km / (kcat/Km) shifts under Zn2+, Ni2+, Co2+ -- not assumed.
"""
from __future__ import annotations
import itertools
import numpy as np


def _neighbor_counts(ca, cutoff=10.0):
    d = np.linalg.norm(ca[:, None, :] - ca[None, :, :], axis=2)
    return (d < cutoff).sum(1) - 1


def _convergent(ca_i, cb_i, ca_j, cb_j):
    """True if both side chains point toward each other (Cβ closer to partner)."""
    if cb_i is None or cb_j is None:
        return True  # GLY etc.: can't judge, don't exclude
    return (np.linalg.norm(cb_i - ca_j) < np.linalg.norm(ca_i - ca_j) and
            np.linalg.norm(cb_j - ca_i) < np.linalg.norm(ca_j - ca_i))


def design_library(result, regions, cfg, focus_regions=None, max_sites=40):
    """Enumerate His-pair and His-triplet candidates within/near top regions."""
    resnums = result["resnums"]
    ca = result["ref_ca"]
    index = {int(r): k for k, r in enumerate(resnums)}
    cb_lookup = result.get("cb_lookup", {})  # optional; filled by run_pipeline

    # residue pool: union of the focused candidate regions (default: top regions)
    regs = focus_regions if focus_regions is not None else regions
    pool = []
    for reg in regs:
        pool.extend(range(reg["start"], reg["end"] + 1))
    pool = [p for p in sorted(set(pool)) if p in index]

    nbr = _neighbor_counts(ca)
    exposed = {int(r): nbr[index[int(r)]] for r in pool}
    # prefer the more solvent-exposed half of the pool
    if exposed:
        med = np.median(list(exposed.values()))
        pool = [p for p in pool if exposed[p] <= med + 2]

    pmin, pmax = cfg["his_pair_dist_A"]
    tmin, tmax = cfg["his_triplet_dist_A"]

    def coords(p):
        k = index[p]
        return ca[k], cb_lookup.get(p)

    pairs = []
    for a, b in itertools.combinations(pool, 2):
        if abs(a - b) < 2:          # avoid i, i+1 (can't both chelate cleanly)
            continue
        ca_a, cb_a = coords(a)
        ca_b, cb_b = coords(b)
        dist = float(np.linalg.norm(ca_a - ca_b))
        if pmin <= dist <= pmax and _convergent(ca_a, cb_a, ca_b, cb_b):
            pairs.append({
                "type": "His-pair",
                "positions": [a, b],
                "ca_dist": round(dist, 2),
                "score": round(_pair_score(dist, exposed.get(a, 0), exposed.get(b, 0)), 3),
            })

    triplets = []
    for a, b, c in itertools.combinations(pool, 3):
        if min(abs(a - b), abs(b - c), abs(a - c)) < 2:
            continue
        cas = [coords(a)[0], coords(b)[0], coords(c)[0]]
        ds = [np.linalg.norm(cas[i] - cas[j]) for i, j in [(0, 1), (1, 2), (0, 2)]]
        if all(tmin <= d <= tmax for d in ds):
            triplets.append({
                "type": "His-triplet",
                "positions": [a, b, c],
                "ca_dists": [round(float(x), 2) for x in ds],
                "score": round(_triplet_score(ds), 3),
            })

    pairs.sort(key=lambda d: d["score"], reverse=True)
    triplets.sort(key=lambda d: d["score"], reverse=True)
    return {"pairs": pairs[:max_sites], "triplets": triplets[:max_sites]}


def _pair_score(dist, exp_a, exp_b):
    # best around ~7 A Cα-Cα; reward exposure (low neighbour count)
    geom = np.exp(-((dist - 7.0) ** 2) / (2 * 1.2 ** 2))
    expo = np.exp(-(exp_a + exp_b) / 30.0)
    return geom * expo


def _triplet_score(ds):
    ds = np.array(ds)
    geom = np.exp(-((ds - 7.0) ** 2 / (2 * 1.5 ** 2)).sum())
    balance = np.exp(-ds.std() / 1.5)   # reward equilateral-ish facial arrangement
    return float(geom * balance)
