"""Configuration for the design loop.

A single :class:`PipelineConfig` drives the closed loop. Weights control the
cell-type-specificity fitness (reward on-target activity in the "positive"
context, penalise activity in "negative" contexts, off-target and splicing
anomalies).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List

from .types import CellContext


@dataclass
class FitnessWeights:
    """Coefficients for the specificity fitness function (see ``scoring.py``)."""

    on_target_expression: float = 1.0
    on_target_accessibility: float = 0.3
    off_target_expression: float = 1.0   # subtracted: activity in negative contexts
    splice_penalty: float = 0.75          # subtracted: cryptic-splicing risk
    tf_binding_bonus: float = 0.2


@dataclass
class PipelineConfig:
    """Everything the :class:`~alphagenome_evo2.loop.DesignLoop` needs.

    ``positive_context`` is the cell line where the element should be ACTIVE;
    ``negative_contexts`` are lines where it must stay SILENT.
    """

    positive_context: CellContext
    negative_contexts: List[CellContext] = field(default_factory=list)

    weights: FitnessWeights = field(default_factory=FitnessWeights)

    # Loop control.
    rounds: int = 8
    population_size: int = 32       # candidates kept as the breeding pool each round
    children_per_round: int = 64    # new sequences Evo 2 proposes each round
    elite_fraction: float = 0.25    # fraction of population carried over unchanged

    # Generation control (passed through to the generator).
    mutation_rate: float = 0.03
    temperature: float = 0.8

    # Convergence: stop early if best fitness improves by < ``tol`` for
    # ``patience`` consecutive rounds.
    tol: float = 1e-3
    patience: int = 3

    # Reproducibility.
    seed: int = 0

    def all_contexts(self) -> List[CellContext]:
        return [self.positive_context, *self.negative_contexts]

    def n_elites(self) -> int:
        return max(1, int(round(self.population_size * self.elite_fraction)))
