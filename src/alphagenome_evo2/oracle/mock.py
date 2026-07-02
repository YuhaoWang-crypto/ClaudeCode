"""A deterministic, dependency-free stand-in for AlphaGenome.

The mock is *not* a biological model. It is a smooth, reproducible fitness
surface that behaves qualitatively like a cell-type-specific regulatory oracle:

* each :class:`CellContext` is assigned a set of "activating" transcription-factor
  motifs; more matches → higher predicted expression/accessibility/TF-binding;
* motifs are context-specific, so a sequence rich in context-A motifs and poor in
  context-B motifs reads as *A-active, B-silent* — exactly the specificity signal
  the loop is asked to optimise;
* dense canonical splice-site patterns raise the ``splice_anomaly`` penalty.

This lets the whole closed loop run and be tested with no GPUs or API keys, and
lets the demo show fitness climbing over rounds. Swap in
:class:`~alphagenome_evo2.oracle.alphagenome.AlphaGenomeOracle` for real biology.
"""

from __future__ import annotations

import hashlib
import math
from typing import Dict, List, Sequence as Seq

from ..types import CellContext, OraclePrediction, Sequence
from .base import Oracle

# Canonical-ish TF motifs used to define per-context "activation grammars".
_MOTIF_LIBRARY = [
    "TGACTCA",   # AP-1
    "GATAAG",    # GATA
    "CACGTG",    # E-box / MYC
    "TTGCGCAA",  # CTCF-like
    "AACCGGTT",  # synthetic palindrome
    "TGASTCA".replace("S", "C"),  # AP-1 variant
    "CCAAT",     # NF-Y / CAAT box
    "TATAAA",    # TATA box
    "GGGCGG",    # SP1 / GC box
    "TGANTCA".replace("N", "G"),
]

# Canonical splice donor/acceptor dinucleotides; clusters imply cryptic splicing.
_SPLICE_DONOR = "GT"
_SPLICE_ACCEPTOR = "AG"


def _stable_hash(text: str) -> int:
    return int(hashlib.sha256(text.encode()).hexdigest(), 16)


def _motifs_for_context(context: CellContext) -> List[str]:
    """Deterministically assign 3 activating motifs to a context by its name.

    Distinct contexts get largely disjoint motif sets, which is what creates the
    tension the specificity objective must resolve.
    """
    h = _stable_hash(context.ontology_term or context.name)
    picks: List[str] = []
    i = 0
    while len(picks) < 3:
        m = _MOTIF_LIBRARY[(h >> (i * 5)) % len(_MOTIF_LIBRARY)]
        if m not in picks:
            picks.append(m)
        i += 1
        if i > 32:  # pragma: no cover - defensive
            break
    return picks


def _count_overlapping(hay: str, needle: str) -> int:
    n = 0
    start = 0
    while True:
        idx = hay.find(needle, start)
        if idx == -1:
            return n
        n += 1
        start = idx + 1


def _saturating(x: float, k: float = 3.0) -> float:
    """Map a non-negative count-like value into (0, 1) with diminishing returns."""
    return 1.0 - math.exp(-x / k)


def _gc_content(seq: str) -> float:
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


class MockOracle(Oracle):
    """Deterministic motif-grammar oracle. Same input → same output, always."""

    def __init__(self, motif_weight: float = 1.0, gc_optimum: float = 0.55):
        self.motif_weight = motif_weight
        self.gc_optimum = gc_optimum

    def predict(
        self,
        sequence: Sequence,
        contexts: Seq[CellContext],
    ) -> Dict[str, OraclePrediction]:
        seq = sequence.seq
        length = max(1, len(seq))
        gc = _gc_content(seq)
        # A mild GC-content preference shared across contexts.
        gc_term = math.exp(-((gc - self.gc_optimum) ** 2) / (2 * 0.12 ** 2))

        splice = self._splice_anomaly(seq)

        out: Dict[str, OraclePrediction] = {}
        for ctx in contexts:
            motifs = _motifs_for_context(ctx)
            hits = sum(_count_overlapping(seq, m) for m in motifs)
            # Density-normalise so long sequences do not win by size alone.
            density = hits / (length / 1000.0)
            activation = _saturating(self.motif_weight * density)

            expression = max(0.0, min(1.0, activation * (0.7 + 0.3 * gc_term)))
            accessibility = max(0.0, min(1.0, 0.6 * activation + 0.4 * gc_term))
            tf_binding = _saturating(density, k=2.0)

            out[ctx.name] = OraclePrediction(
                context=ctx,
                expression=expression,
                accessibility=accessibility,
                tf_binding=tf_binding,
                splice_anomaly=splice,
                raw={
                    "motif_hits": float(hits),
                    "motif_density_per_kb": float(density),
                    "gc_content": float(gc),
                },
            )
        return out

    @staticmethod
    def _splice_anomaly(seq: str) -> float:
        """Penalise dense donor/acceptor dinucleotide pairs (cryptic splicing)."""
        donors = _count_overlapping(seq, _SPLICE_DONOR)
        acceptors = _count_overlapping(seq, _SPLICE_ACCEPTOR)
        pairs = min(donors, acceptors)
        density = pairs / (max(1, len(seq)) / 100.0)
        # Only unusually dense clustering should register as a strong penalty.
        return _saturating(max(0.0, density - 2.0), k=4.0)
