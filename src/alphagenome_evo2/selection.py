"""Selection strategies for the evolutionary loop."""

from __future__ import annotations

import random
from typing import List, Sequence as Seq

from .types import Candidate


def _sorted_by_fitness(candidates: Seq[Candidate]) -> List[Candidate]:
    return sorted(
        candidates,
        key=lambda c: (c.fitness if c.fitness is not None else float("-inf")),
        reverse=True,
    )


def top_k(candidates: Seq[Candidate], k: int) -> List[Candidate]:
    """Deterministic truncation selection: keep the ``k`` fittest candidates."""
    return _sorted_by_fitness(candidates)[: max(0, k)]


def tournament(
    candidates: Seq[Candidate],
    k: int,
    tournament_size: int = 3,
    rng: random.Random | None = None,
) -> List[Candidate]:
    """Stochastic tournament selection — preserves diversity better than top-k.

    Repeatedly samples ``tournament_size`` contenders and keeps the fittest until
    ``k`` winners are chosen (without replacement).
    """
    rng = rng or random.Random()
    pool = list(candidates)
    winners: List[Candidate] = []
    while pool and len(winners) < k:
        contenders = rng.sample(pool, min(tournament_size, len(pool)))
        best = max(contenders, key=lambda c: c.fitness if c.fitness is not None else float("-inf"))
        winners.append(best)
        pool.remove(best)
    return winners


def next_population(
    candidates: Seq[Candidate],
    population_size: int,
    n_elites: int,
    rng: random.Random,
    tournament_size: int = 3,
) -> List[Candidate]:
    """Form the next breeding pool: elites (kept verbatim) + tournament picks."""
    ranked = _sorted_by_fitness(candidates)
    elites = ranked[:n_elites]
    remaining_slots = max(0, population_size - len(elites))
    rest_pool = ranked[n_elites:]
    picks = tournament(rest_pool, remaining_slots, tournament_size, rng)
    # Backfill deterministically if the tournament pool was too small.
    if len(picks) < remaining_slots:
        already = {id(c) for c in elites} | {id(c) for c in picks}
        for c in ranked:
            if len(picks) >= remaining_slots:
                break
            if id(c) not in already:
                picks.append(c)
                already.add(id(c))
    return elites + picks
