"""
biosensor_pipeline
==================

Computational reproduction of the single-component allosteric protein-switch
workflow of:

    Guo, Smutok, Lee, ... Baker & Alexandrov,
    "Artificial allosteric protein switches with machine-learning-designed
     receptors", Nature Biotechnology (2026).
    https://doi.org/10.1038/s41587-026-03081-9

The paper builds biosensors by (1) taking a small, rigid ML-designed
ligand-binding domain (the *receptor*), (2) circularly permuting it, and
(3) inserting the permuted receptor into a loop of a reporter enzyme
(TEM-1 beta-lactamase, or a de novo luciferase). Ligand binding modulates
the conformational entropy of the chimera and switches the reporter's
catalytic activity ON.

This package reproduces the *design* half of that workflow deterministically
(pure sequence engineering, exactly as described in Methods "Design of
chimeric proteins") and provides an *in-silico validation* half that scores
candidate chimeras with structure prediction (Boltz).

HONESTY LABELING (inherited from the repo's discipline):
  ✅ rigorous   -> deterministic sequence construction; real public sequences;
                   a real structural metric computed from a predicted model.
  ⚠️ hypothesis -> anything that stands in for a wet-lab measurement
                   (e.g. an in-silico "dynamic-range proxy" score, a
                   loop-site heuristic). NEVER conflated with a measured DR.
"""

from .design import (
    circular_permute,
    build_chimera,
    Chimera,
)
from .systems import REPORTERS, RECEPTORS, SYSTEMS, get_system

__all__ = [
    "circular_permute",
    "build_chimera",
    "Chimera",
    "REPORTERS",
    "RECEPTORS",
    "SYSTEMS",
    "get_system",
]
