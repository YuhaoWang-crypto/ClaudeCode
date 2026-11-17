"""End-to-end orchestration:

    generate (Beam Search)  ->  CTMC sequence variants  ->  self-consistency
    refold + liability/physics filtering  ->  rank  ->  3A clustering  ->  <=N

Everything here is backend-agnostic: pass any Generator + Folder.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from .backends.base import Folder, Generator
from .beam_search import BeamConfig, beam_search
from .cluster import cluster, pick_representatives
from .confidence import confidence_score
from .ctmc import redesign
from .filters import cdr_liabilities, clash_count, passes
from .selfconsistency import pose_recovery
from .types import Design, DesignSpec, ScoredDesign


@dataclass
class PipelineConfig:
    beam: BeamConfig = field(default_factory=BeamConfig)
    ctmc_variants: int = 3            # extra sequences per surviving backbone
    cluster_cutoff: float = 3.0       # angstrom
    max_candidates: int = 48          # wet-lab budget
    filter_thresholds: dict = field(default_factory=dict)
    seed: Optional[int] = 0


def _final_score(conf_score: float, recovery_rmsd: float) -> float:
    """Combine confidence (higher better) and recovery RMSD (lower better)."""
    rmsd_term = 1.0 / (1.0 + max(recovery_rmsd, 0.0))
    return 0.7 * conf_score + 0.3 * rmsd_term


class DesignPipeline:
    def __init__(self, generator: Generator, folder: Folder,
                 config: Optional[PipelineConfig] = None):
        self.generator = generator
        self.folder = folder
        self.cfg = config or PipelineConfig()
        self.trace: dict = {}

    def run(self, spec: DesignSpec) -> list[ScoredDesign]:
        # ---- 1-2. generation under confidence-guided Beam Search ----------
        beam_trace: list[dict] = []
        beam = beam_search(spec, self.generator, self.cfg.beam, trace=beam_trace)
        self.trace["beam"] = beam_trace

        # ---- 3. CTMC sequence redesign on the surviving backbones ---------
        population: list[Design] = list(beam)
        rng = np.random.default_rng(self.cfg.seed)
        for d in beam:
            population.extend(redesign(d, self.cfg.ctmc_variants, rng=rng))
        self.trace["population_size"] = len(population)

        # ---- 4. self-consistency refold + 5. filtering --------------------
        scored: list[ScoredDesign] = []
        for d in population:
            recovery, refold_conf = pose_recovery(d, self.folder, spec, seed=self.cfg.seed)
            liab = cdr_liabilities(d)
            clashes = clash_count(d)
            ok = passes(recovery, liab, clashes, self.cfg.filter_thresholds)
            cs = confidence_score(refold_conf)
            scored.append(ScoredDesign(
                design=d, confidence=refold_conf, recovery_rmsd=recovery,
                liabilities=liab, clashes=clashes, passed_filters=ok,
                score=_final_score(cs, recovery),
            ))

        survivors = [s for s in scored if s.passed_filters]
        self.trace["generated"] = len(population)
        self.trace["passed_filters"] = len(survivors)
        if not survivors:
            return []

        # ---- 5. 3A clustering, keep one representative per pose ------------
        survivors.sort(key=lambda s: s.score, reverse=True)
        labels = cluster([s.design for s in survivors], cutoff=self.cfg.cluster_cutoff)
        for s, lab in zip(survivors, labels):
            s.cluster_id = int(lab)
        reps = pick_representatives(
            [s.design for s in survivors], labels, score=lambda i: survivors[i].score
        )
        for i in reps:
            survivors[i].is_representative = True
        self.trace["n_clusters"] = int(labels.max()) + 1 if len(labels) else 0

        # Final selection: representatives first (diverse poses), then best
        # remaining by score, capped at the wet-lab budget.
        chosen = [survivors[i] for i in reps]
        if len(chosen) < self.cfg.max_candidates:
            extra = [s for s in survivors if not s.is_representative]
            chosen.extend(extra[: self.cfg.max_candidates - len(chosen)])
        chosen = chosen[: self.cfg.max_candidates]
        self.trace["selected"] = len(chosen)
        return chosen
