"""CTMC-style sequence redesign for the pipeline (plan step 3).

This now delegates to the EXACT motif-constrained CTMC in ``abproto.model.ctmc``:
liability motifs are removed by zeroing the corresponding transitions in the
rate matrix, so every variant is guaranteed free of N-glyc sequons, CDR
cysteines, and (by default) deamidation/isomerization hotspots — not merely
unlikely as with rejection sampling.

The learned denoiser p_theta(x_1 | x_t, c) is swapped in for ``_parent_biased``
once the model is trained; here a parent-biased posterior keeps variants close
to the generated sequence while the constraint guarantees developability.
"""

from __future__ import annotations

import copy
from typing import Optional

import numpy as np

from .model.ctmc import AA, ConstrainedCTMC, seq_to_tokens, tokens_to_seq
from .types import CDRS, Design


def _region_mask(length: int, spans: list[tuple[int, int]]) -> np.ndarray:
    m = np.zeros(length, dtype=bool)
    for lo, hi in spans:
        m[lo:hi] = True
    return m


def _parent_biased(parent_tokens: np.ndarray, noise: float = 0.5):
    """Denoiser posterior concentrated on the parent residue (temperature via
    ``noise``), so CTMC variants are local edits of the parent sequence."""
    def denoiser(tokens, t, cond=None):
        L = len(parent_tokens)
        p = np.full((L, len(AA)), noise / len(AA))
        for i, tok in enumerate(parent_tokens):
            if 0 <= int(tok) < len(AA):
                p[i, int(tok)] += (1.0 - noise)
        p /= p.sum(axis=1, keepdims=True)
        return p
    return denoiser


def redesign(design: Design, n_variants: int, rng: Optional[np.random.Generator] = None,
             noise: float = 0.5, steps: int = 60) -> list[Design]:
    """Return ``n_variants`` sequence variants on the SAME backbone, each with
    guaranteed-zero liability motifs in the CDRs (exact rate-matrix constraint)."""
    if rng is None:
        rng = np.random.default_rng(abs(hash(design.design_id)) % (2**32))

    h_spans = [design.cdr_spans[c] for c in ("H1", "H2", "H3") if c in design.cdr_spans]
    l_spans = [design.cdr_spans[c] for c in ("L1", "L2", "L3") if c in design.cdr_spans]
    h_tokens = seq_to_tokens(design.heavy_seq)
    l_tokens = seq_to_tokens(design.light_seq)
    h_mask = _region_mask(len(h_tokens), h_spans)
    l_mask = _region_mask(len(l_tokens), l_spans)

    h_ctmc = ConstrainedCTMC(denoiser=_parent_biased(h_tokens, noise), steps=steps)
    l_ctmc = ConstrainedCTMC(denoiser=_parent_biased(l_tokens, noise), steps=steps)

    variants: list[Design] = []
    for k in range(n_variants):
        seed = int(rng.integers(2**31))
        new_h = h_ctmc.sample(len(h_tokens), h_mask, seed=seed, init_tokens=h_tokens.copy())
        new_l = l_ctmc.sample(len(l_tokens), l_mask, seed=seed + 1, init_tokens=l_tokens.copy())
        v = copy.deepcopy(design)
        v.heavy_seq = tokens_to_seq(new_h)
        v.light_seq = tokens_to_seq(new_l)
        v.design_id = f"{design.design_id}_ctmc{k}"
        v.parent_id = design.design_id
        v.provenance = "ctmc"
        variants.append(v)
    return variants
