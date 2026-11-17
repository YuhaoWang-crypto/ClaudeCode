"""Confidence-guided Beam Search = test-time compute scaling.

Instead of one-shot sampling, we keep a beam of the best candidates and spend
inference compute iteratively refining them. Each step:
  1. expand : every survivor spawns ``branching`` refined children
  2. score  : rank all candidates by the confidence objective
  3. prune  : keep the top ``beam_width``

This is the single highest-ROI mechanism in the plan; the generator's ``refine``
decides how children are sampled (partial re-noising in production).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional

from .backends.base import Generator
from .confidence import confidence_score
from .types import Design, DesignSpec


@dataclass
class BeamConfig:
    beam_width: int = 16
    branching: int = 4          # children per survivor each step
    steps: int = 3
    init_pool: int = 64         # initial random draw before the first prune
    seed: Optional[int] = None


def _score(d: Design, score_fn: Callable[[Design], float]) -> float:
    return score_fn(d)


def beam_search(spec: DesignSpec, generator: Generator, cfg: BeamConfig,
                score_fn: Optional[Callable[[Design], float]] = None,
                trace: Optional[list[dict]] = None) -> list[Design]:
    """Return the final beam (best ``beam_width`` designs).

    ``score_fn`` maps a Design to a scalar (higher better). Default uses the
    generation-time confidence attached to each design.
    """
    if score_fn is None:
        def score_fn(d: Design) -> float:
            return confidence_score(d.confidence) if d.confidence else 0.0

    # Initial broad draw, then prune to the beam.
    pool = generator.generate(spec, cfg.init_pool, seed=cfg.seed)
    pool.sort(key=lambda d: score_fn(d), reverse=True)
    beam = pool[: cfg.beam_width]
    if trace is not None:
        trace.append(_step_stats(0, beam, score_fn))

    for step in range(1, cfg.steps + 1):
        children = generator.refine(beam, spec, cfg.branching,
                                    seed=None if cfg.seed is None else cfg.seed + step)
        # Survivors compete with their own offspring (elitism keeps best-so-far).
        candidates = beam + children
        candidates.sort(key=lambda d: score_fn(d), reverse=True)
        beam = candidates[: cfg.beam_width]
        if trace is not None:
            trace.append(_step_stats(step, beam, score_fn))

    return beam


def _step_stats(step: int, beam: list[Design], score_fn) -> dict:
    scores = [score_fn(d) for d in beam]
    return {
        "step": step,
        "beam_size": len(beam),
        "best_score": round(max(scores), 4),
        "mean_score": round(sum(scores) / len(scores), 4),
    }
