"""Core data structures shared across the AlphaGenome ↔ Evo 2 pipeline.

These types are deliberately backend-agnostic: an ``Oracle`` (AlphaGenome) and a
``Generator`` (Evo 2) exchange plain :class:`Sequence` / :class:`Candidate`
objects, so the mock and real adapters are fully interchangeable.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Mapping, Optional, Sequence as Seq

VALID_BASES = frozenset("ACGT")


def validate_dna(seq: str) -> str:
    """Return an upper-cased DNA string, raising on non-ACGT characters."""
    s = seq.strip().upper()
    if not s:
        raise ValueError("DNA sequence is empty")
    bad = set(s) - VALID_BASES
    if bad:
        raise ValueError(f"sequence contains non-ACGT characters: {sorted(bad)}")
    return s


@dataclass(frozen=True)
class CellContext:
    """A cell line / tissue context that AlphaGenome can predict for.

    ``ontology_term`` mirrors AlphaGenome's ontology inputs (e.g. an EFO/CL/UBERON
    identifier such as ``"EFO:0002067"`` for K562). ``name`` is a human label.
    """

    name: str
    ontology_term: Optional[str] = None

    def __str__(self) -> str:  # pragma: no cover - trivial
        return self.name


@dataclass(frozen=True)
class Sequence:
    """An immutable DNA sequence with a stable identifier and free-form metadata."""

    id: str
    seq: str
    metadata: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "seq", validate_dna(self.seq))

    def __len__(self) -> int:
        return len(self.seq)

    def with_metadata(self, **updates: object) -> "Sequence":
        merged = {**self.metadata, **updates}
        return replace(self, metadata=merged)


@dataclass(frozen=True)
class OraclePrediction:
    """AlphaGenome-style read-outs for a single sequence in a single context.

    All scores are normalised to roughly ``[0, 1]`` by the oracle so that
    downstream fitness functions can combine them without per-track calibration.
    """

    context: CellContext
    expression: float          # aggregate transcriptional output (RNA-seq/CAGE)
    accessibility: float       # chromatin accessibility (DNASE/ATAC)
    tf_binding: float          # aggregate transcription-factor occupancy (ChIP)
    splice_anomaly: float      # penalty in [0, 1]; higher = more cryptic splicing
    raw: Mapping[str, float] = field(default_factory=dict)


@dataclass(frozen=True)
class Candidate:
    """A design candidate: a sequence plus the predictions gathered for it."""

    sequence: Sequence
    predictions: Mapping[str, OraclePrediction] = field(default_factory=dict)
    fitness: Optional[float] = None
    fitness_terms: Mapping[str, float] = field(default_factory=dict)

    @property
    def id(self) -> str:
        return self.sequence.id

    def prediction_for(self, context: CellContext) -> Optional[OraclePrediction]:
        return self.predictions.get(context.name)

    def scored(
        self,
        fitness: float,
        terms: Optional[Mapping[str, float]] = None,
    ) -> "Candidate":
        return replace(self, fitness=fitness, fitness_terms=dict(terms or {}))


@dataclass
class RoundResult:
    """Snapshot of one optimisation round for reporting / convergence checks."""

    round_index: int
    n_generated: int
    n_evaluated: int
    best_fitness: float
    mean_elite_fitness: float
    best_candidate_id: str


def elite_stats(candidates: Seq[Candidate]) -> tuple[float, float]:
    """Return (best, mean) fitness over ``candidates`` (0.0 for an empty set)."""
    scored = [c.fitness for c in candidates if c.fitness is not None]
    if not scored:
        return 0.0, 0.0
    return max(scored), sum(scored) / len(scored)
