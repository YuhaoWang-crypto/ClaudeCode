"""Genome-wide mining of candidate regulatory elements (pipeline step 1).

Slides a window across supplied genomic sequence(s), scores every window with the
oracle, and returns the windows that are **active in the positive context** and
**silent in the negative contexts** — the natural enhancers/promoters that seed
the Evo 2 generation stage.

Input genomes are provided by the caller (e.g. read from FASTA); this module does
not fetch reference genomes itself, keeping it dependency-free and testable.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Sequence as Seq, Tuple

from .config import PipelineConfig
from .oracle.base import Oracle
from .scoring import score_candidate
from .types import Candidate, Sequence


@dataclass
class MiningConfig:
    window: int = 200
    stride: int = 100
    top_n: int = 32
    min_specificity_gap: float = 0.0  # require expr(pos) - max expr(neg) >= this


def _windows(seq: Sequence, window: int, stride: int) -> Iterable[Tuple[int, str]]:
    s = seq.seq
    if len(s) <= window:
        yield 0, s
        return
    for start in range(0, len(s) - window + 1, stride):
        yield start, s[start:start + window]


def mine_elements(
    genomes: Seq[Sequence],
    oracle: Oracle,
    config: PipelineConfig,
    mining: MiningConfig | None = None,
) -> List[Candidate]:
    """Scan ``genomes`` and return the top specificity-scoring windows.

    The returned candidates carry oracle predictions and fitness, ready to seed
    :class:`~alphagenome_evo2.loop.DesignLoop`.
    """
    mining = mining or MiningConfig()
    contexts = config.all_contexts()

    tiles: List[Sequence] = []
    for genome in genomes:
        for start, sub in _windows(genome, mining.window, mining.stride):
            tiles.append(
                Sequence(
                    id=f"{genome.id}:{start}-{start + len(sub)}",
                    seq=sub,
                    metadata={"origin": "mined", "source": genome.id, "start": start},
                )
            )

    evaluated = oracle.predict_batch(tiles, contexts)

    scored: List[Candidate] = []
    for cand in evaluated:
        fitness, terms = score_candidate(cand, config)
        if terms.get("specificity_gap", 0.0) < mining.min_specificity_gap:
            continue
        scored.append(cand.scored(fitness, terms))

    scored.sort(key=lambda c: c.fitness or float("-inf"), reverse=True)
    return scored[: mining.top_n]
