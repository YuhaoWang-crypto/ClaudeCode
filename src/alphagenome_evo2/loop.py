"""The AlphaGenome ↔ Evo 2 closed-loop optimiser (pipeline steps 2–4).

Implements the active-learning cycle from the project spec:

    seed pool ──▶ (2) Evo 2 generates variants ──▶ (3) AlphaGenome scores them
        ▲                                                     │
        └──────────── (4) select elites, iterate ◀────────────┘

The generator (Evo 2) and oracle (AlphaGenome) are injected, so the exact same
loop runs against the mock backends (offline, in tests) or the real models.
"""

from __future__ import annotations

import random
from dataclasses import dataclass, field
from typing import Callable, List, Optional, Sequence as Seq

from .config import PipelineConfig
from .generator.base import Generator
from .oracle.base import Oracle
from .scoring import score_candidate
from .selection import next_population
from .types import Candidate, RoundResult, Sequence, elite_stats


@dataclass
class LoopResult:
    """Everything a completed run produces."""

    best: Candidate
    population: List[Candidate]           # final breeding pool, fittest first
    history: List[RoundResult] = field(default_factory=list)
    converged: bool = False

    @property
    def best_fitness(self) -> float:
        return self.best.fitness if self.best.fitness is not None else float("-inf")


class DesignLoop:
    """Runs the generate → score → select cycle to convergence."""

    def __init__(
        self,
        oracle: Oracle,
        generator: Generator,
        config: PipelineConfig,
        on_round: Optional[Callable[[RoundResult], None]] = None,
    ):
        self.oracle = oracle
        self.generator = generator
        self.config = config
        self.on_round = on_round
        self._rng = random.Random(config.seed)

    # -- public API ---------------------------------------------------------
    def run(self, seeds: Seq[Sequence | Candidate]) -> LoopResult:
        """Optimise starting from ``seeds`` (raw sequences or scored candidates)."""
        population = self._score(self._as_sequences(seeds))
        population = self._trim(population, self.config.population_size)

        history: List[RoundResult] = []
        best = self._best(population)
        stale_rounds = 0
        converged = False

        for r in range(self.config.rounds):
            parents = [c.sequence for c in population]
            children_seqs = self.generator.generate(
                parents,
                self.config.children_per_round,
                temperature=self.config.temperature,
                mutation_rate=self.config.mutation_rate,
            )
            children = self._score(children_seqs)

            combined = population + children
            new_best = self._best(combined)

            population = next_population(
                combined,
                self.config.population_size,
                self.config.n_elites(),
                self._rng,
            )

            best_fit, mean_fit = elite_stats(population)
            record = RoundResult(
                round_index=r,
                n_generated=len(children_seqs),
                n_evaluated=len(children),
                best_fitness=best_fit,
                mean_elite_fitness=mean_fit,
                best_candidate_id=new_best.id if new_best else "",
            )
            history.append(record)
            if self.on_round:
                self.on_round(record)

            improvement = new_best.fitness - best.fitness if (new_best and best) else 0.0
            if new_best and (best is None or new_best.fitness > best.fitness):
                best = new_best
            if improvement < self.config.tol:
                stale_rounds += 1
            else:
                stale_rounds = 0
            if stale_rounds >= self.config.patience:
                converged = True
                break

        population = sorted(
            population, key=lambda c: c.fitness or float("-inf"), reverse=True
        )
        return LoopResult(
            best=best, population=population, history=history, converged=converged
        )

    # -- internals ----------------------------------------------------------
    def _as_sequences(self, seeds: Seq[Sequence | Candidate]) -> List[Sequence]:
        out: List[Sequence] = []
        for s in seeds:
            out.append(s.sequence if isinstance(s, Candidate) else s)
        if not out:
            raise ValueError("DesignLoop.run requires at least one seed sequence")
        return out

    def _score(self, sequences: Seq[Sequence]) -> List[Candidate]:
        contexts = self.config.all_contexts()
        evaluated = self.oracle.predict_batch(sequences, contexts)
        scored: List[Candidate] = []
        for cand in evaluated:
            fitness, terms = score_candidate(cand, self.config)
            scored.append(cand.scored(fitness, terms))
        return scored

    @staticmethod
    def _trim(candidates: List[Candidate], size: int) -> List[Candidate]:
        return sorted(
            candidates, key=lambda c: c.fitness or float("-inf"), reverse=True
        )[:size]

    @staticmethod
    def _best(candidates: Seq[Candidate]) -> Optional[Candidate]:
        scored = [c for c in candidates if c.fitness is not None]
        return max(scored, key=lambda c: c.fitness) if scored else None
