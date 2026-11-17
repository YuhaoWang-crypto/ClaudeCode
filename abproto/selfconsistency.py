"""Self-consistency refolding = pose recovery.

Refold a design's sequence unconditionally, then measure how far the refolded
antibody drifts from the designed one. Two metrics:

  - global : Kabsch RMSD over the whole antibody (heavy+light) CA trace.
  - interface (optimisation #6): align the complex on the ANTIGEN, then measure
    RMSD over the PARATOPE (antibody residues contacting the antigen). This is
    sensitive to the binding MODE — a design can pass global RMSD while the
    paratope has rotated off the epitope; interface RMSD catches that.

The pipeline uses interface recovery when the antigen is present in both the
designed and refolded structures, falling back to global otherwise.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from .backends.base import Folder
from .geometry import apply_transform, kabsch_fit, kabsch_rmsd, rmsd
from .types import ANTIGEN, HEAVY, LIGHT, Confidence, Design, DesignSpec


def pose_recovery(design: Design, folder: Folder, spec: DesignSpec,
                  seed: Optional[int] = None, contact: float = 10.0
                  ) -> tuple[float, Confidence]:
    """Return (recovery RMSD in angstrom, refolded confidence).

    Uses interface RMSD when possible (antigen present in both structures), else
    global antibody RMSD.
    """
    pose = folder.fold(design, spec, seed=seed)

    design_ab = _antibody_ca(design.backbone)
    refold_ab = _antibody_ca(pose.backbone)
    if design_ab.shape != refold_ab.shape:
        return float("inf"), pose.confidence   # topology changed => failed recovery

    if ANTIGEN in design.backbone and ANTIGEN in pose.backbone:
        r = _interface_rmsd(design.backbone, pose.backbone, contact=contact)
        if r is not None:
            return r, pose.confidence

    return kabsch_rmsd(design_ab, refold_ab), pose.confidence


def _interface_rmsd(design_bb: dict, refold_bb: dict, contact: float) -> Optional[float]:
    """Align refold onto design by the antigen, then RMSD over the paratope."""
    d_ag, r_ag = design_bb[ANTIGEN], refold_bb[ANTIGEN]
    if d_ag.shape != r_ag.shape or len(d_ag) < 3:
        return None
    d_ab = _antibody_ca(design_bb)
    r_ab = _antibody_ca(refold_bb)

    # Paratope = antibody residues within ``contact`` A of any antigen CA in the
    # DESIGNED structure.
    dmat = np.linalg.norm(d_ab[:, None, :] - d_ag[None, :, :], axis=-1)
    paratope = np.where(dmat.min(axis=1) <= contact)[0]
    if paratope.size < 3:
        paratope = np.argsort(dmat.min(axis=1))[:max(3, len(d_ab) // 5)]

    # Superpose on the antigen, apply to the refolded antibody, measure paratope.
    R, t = kabsch_fit(r_ag, d_ag)
    r_ab_aligned = apply_transform(r_ab, R, t)
    return rmsd(d_ab[paratope], r_ab_aligned[paratope])


def _antibody_ca(backbone: dict) -> np.ndarray:
    parts = [backbone[c] for c in (HEAVY, LIGHT) if c in backbone]
    return np.concatenate(parts, axis=0)
