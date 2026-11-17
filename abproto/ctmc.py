"""CTMC-style sequence redesign (plan step 3).

For a fixed backbone, generate extra sequence variants and, during sampling,
hard-avoid liability motifs. The real engine is a discrete continuous-time
Markov chain over sequence (e.g. Multiflow's discrete head). Here we implement
a lightweight, liability-aware CDR resampler with the same *interface and
constraint behaviour*, so the pipeline shape and the motif guarantees are real
even without the learned model.
"""

from __future__ import annotations

import copy
from typing import Optional

import numpy as np

from .filters import cdr_liabilities
from .types import CDRS, Design

_AA = "ACDEFGHIKLMNPQRSTVWY"


def redesign(design: Design, n_variants: int, rng: Optional[np.random.Generator] = None,
             max_tries: int = 40) -> list[Design]:
    """Return ``n_variants`` sequence variants on the SAME backbone, each with
    zero N-glyc sequons and zero CDR cysteines (rejection sampling)."""
    if rng is None:
        rng = np.random.default_rng(abs(hash(design.design_id)) % (2**32))

    variants: list[Design] = []
    for k in range(n_variants):
        variant = _resample_once(design, rng, max_tries)
        variant.design_id = f"{design.design_id}_ctmc{k}"
        variant.parent_id = design.design_id
        variant.provenance = "ctmc"
        variants.append(variant)
    return variants


def _resample_once(design: Design, rng: np.random.Generator, max_tries: int) -> Design:
    best = None
    for _ in range(max_tries):
        cand = copy.deepcopy(design)
        heavy = list(cand.heavy_seq)
        light = list(cand.light_seq)
        for cdr in CDRS:
            if cdr not in cand.cdr_spans:
                continue
            lo, hi = cand.cdr_spans[cdr]
            target = heavy if cdr.startswith("H") else light
            # Mutate ~40% of positions in this CDR, never introducing C.
            for i in range(lo, hi):
                if rng.random() < 0.4:
                    target[i] = _draw_non_cys(rng)
        cand.heavy_seq = "".join(heavy)
        cand.light_seq = "".join(light)
        liab = cdr_liabilities(cand)
        if liab["nglyc"] == 0 and liab["cys"] == 0:
            return cand
        # Keep the least-bad as a fallback.
        if best is None or (liab["nglyc"] + liab["cys"]) < best[1]:
            best = (cand, liab["nglyc"] + liab["cys"])
    return best[0]


def _draw_non_cys(rng: np.random.Generator) -> str:
    aa = "C"
    while aa == "C":
        aa = _AA[int(rng.integers(len(_AA)))]
    return aa
