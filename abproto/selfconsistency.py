"""Self-consistency refolding = pose recovery.

Refold a design's sequence unconditionally with the Folder, then measure how
far the refolded antibody backbone drifts from the designed one (Kabsch CA
RMSD). Low RMSD => the sequence folds back to the intended binding geometry.
This replaces AF2-based validation, which is unreliable for antibody complexes.
"""

from __future__ import annotations

from typing import Optional

from .backends.base import Folder
from .geometry import kabsch_rmsd
from .types import Confidence, Design, DesignSpec


def pose_recovery(design: Design, folder: Folder, spec: DesignSpec,
                  seed: Optional[int] = None) -> tuple[float, Confidence]:
    """Return (recovery_rmsd in angstrom, refolded confidence)."""
    pose = folder.fold(design, spec, seed=seed)

    # Align the antibody (heavy+light) CA traces; antigen is the fixed frame.
    design_ca = design.antibody_ca()
    refold_ca = _antibody_ca(pose.backbone, design)

    if design_ca.shape != refold_ca.shape:
        # Length mismatch => folding produced a different topology: treat as a
        # maximal (failing) recovery rather than crashing the pipeline.
        return float("inf"), pose.confidence

    return kabsch_rmsd(design_ca, refold_ca), pose.confidence


def _antibody_ca(backbone: dict, design: Design):
    import numpy as np
    from .types import HEAVY, LIGHT

    parts = [backbone[c] for c in (HEAVY, LIGHT) if c in backbone]
    return np.concatenate(parts, axis=0)
