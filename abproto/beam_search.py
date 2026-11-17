"""Test-time compute scaling: confidence-guided search over the generator.

Two strategies:
  - beam_search : keep a beam of the best candidates; each step every survivor
                  spawns refined children, then prune. With ``diverse=True`` the
                  prune is MMR over backbone RMSD so the beam keeps DISTINCT
                  binding poses instead of collapsing to one (optimisation #2).
  - smc_search  : sequential Monte Carlo / FK-steering. Particles are reweighted
                  by exp(score/temp) and resampled, concentrating compute on
                  high-confidence modes while preserving a population.

Both take a scalar ``score_fn`` (default = interface-weighted confidence).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np

from .backends.base import Generator
from .confidence import confidence_score
from .geometry import kabsch_rmsd
from .types import Design, DesignSpec


@dataclass
class BeamConfig:
    beam_width: int = 16
    branching: int = 4
    steps: int = 3
    init_pool: int = 64
    seed: Optional[int] = None
    diverse: bool = True            # MMR-style diversity-preserving prune
    diversity_weight: float = 0.5   # lambda in MMR (0 = pure score, 1 = pure diversity)
    rmsd_scale: float = 3.0         # angstrom; similarity = exp(-rmsd/scale)


# --------------------------------------------------------------------------
def _default_score(d: Design) -> float:
    return confidence_score(d.confidence) if d.confidence else 0.0


def backbone_similarity(a: Design, b: Design, scale: float = 3.0) -> float:
    """Similarity in [0,1] from antibody-CA RMSD; 0 if shapes differ (distinct)."""
    ca, cb = a.antibody_ca(), b.antibody_ca()
    if ca.shape != cb.shape:
        return 0.0
    return float(np.exp(-kabsch_rmsd(ca, cb) / scale))


def diverse_topk(cands: list[Design], score_fn: Callable[[Design], float], k: int,
                 diversity_weight: float = 0.5, rmsd_scale: float = 3.0) -> list[Design]:
    """Maximal-marginal-relevance selection: greedily pick the candidate that
    maximises score - lambda * max_similarity_to_already_selected."""
    if k >= len(cands):
        return sorted(cands, key=score_fn, reverse=True)
    remaining = sorted(cands, key=score_fn, reverse=True)
    selected = [remaining.pop(0)]
    while remaining and len(selected) < k:
        best_i, best_val = 0, -np.inf
        for i, c in enumerate(remaining):
            max_sim = max(backbone_similarity(c, s, rmsd_scale) for s in selected)
            val = score_fn(c) - diversity_weight * max_sim
            if val > best_val:
                best_val, best_i = val, i
        selected.append(remaining.pop(best_i))
    return selected


def _prune(cands: list[Design], score_fn, cfg: BeamConfig) -> list[Design]:
    if cfg.diverse:
        return diverse_topk(cands, score_fn, cfg.beam_width,
                            cfg.diversity_weight, cfg.rmsd_scale)
    return sorted(cands, key=score_fn, reverse=True)[: cfg.beam_width]


def beam_search(spec: DesignSpec, generator: Generator, cfg: BeamConfig,
                score_fn: Optional[Callable[[Design], float]] = None,
                trace: Optional[list[dict]] = None) -> list[Design]:
    score_fn = score_fn or _default_score
    pool = generator.generate(spec, cfg.init_pool, seed=cfg.seed)
    beam = _prune(pool, score_fn, cfg)
    if trace is not None:
        trace.append(_stats(0, beam, score_fn))

    for step in range(1, cfg.steps + 1):
        children = generator.refine(beam, spec, cfg.branching,
                                    seed=None if cfg.seed is None else cfg.seed + step)
        beam = _prune(beam + children, score_fn, cfg)   # elitism: parents compete
        if trace is not None:
            trace.append(_stats(step, beam, score_fn))
    return beam


# --------------------------------------------------------------------------
@dataclass
class SMCConfig:
    n_particles: int = 32
    steps: int = 3
    temperature: float = 0.5        # lower => sharper reweighting
    ess_frac: float = 0.5           # resample when ESS < ess_frac * N
    seed: Optional[int] = None


def smc_search(spec: DesignSpec, generator: Generator, cfg: SMCConfig,
               score_fn: Optional[Callable[[Design], float]] = None,
               trace: Optional[list[dict]] = None) -> list[Design]:
    """Sequential Monte Carlo over the generator (FK-steering).

    Particles are reweighted by exp(score/temp) and resampled when the effective
    sample size drops, then propagated by refinement. Returns the final particles
    sorted by score.
    """
    score_fn = score_fn or _default_score
    rng = np.random.default_rng(cfg.seed)
    particles = generator.generate(spec, cfg.n_particles, seed=cfg.seed)

    for step in range(cfg.steps + 1):
        scores = np.array([score_fn(p) for p in particles])
        w = np.exp((scores - scores.max()) / max(cfg.temperature, 1e-6))
        w = w / w.sum()
        ess = 1.0 / np.sum(w ** 2)
        if trace is not None:
            trace.append({"step": step, "best_score": round(float(scores.max()), 4),
                          "mean_score": round(float(scores.mean()), 4),
                          "ess": round(float(ess), 2)})
        if step == cfg.steps:
            break
        # Resample (systematic) when ESS is low, else keep particles.
        if ess < cfg.ess_frac * cfg.n_particles:
            idx = _systematic_resample(w, rng)
            particles = [particles[i] for i in idx]
        # Propagate: one refined child per particle.
        particles = generator.refine(particles, spec, 1,
                                     seed=None if cfg.seed is None else cfg.seed + step + 1)

    return sorted(particles, key=score_fn, reverse=True)


def _systematic_resample(w: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    n = len(w)
    positions = (rng.random() + np.arange(n)) / n
    cumsum = np.cumsum(w)
    idx = np.zeros(n, dtype=int)
    i = j = 0
    while i < n:
        if positions[i] < cumsum[j]:
            idx[i] = j; i += 1
        else:
            j += 1
    return idx


# --------------------------------------------------------------------------
def _stats(step: int, beam: list[Design], score_fn) -> dict:
    scores = [score_fn(d) for d in beam]
    return {"step": step, "beam_size": len(beam),
            "best_score": round(max(scores), 4),
            "mean_score": round(sum(scores) / len(scores), 4)}
