"""Core data types passed between pipeline stages.

Backbone geometry is kept as C-alpha traces (``dict[chain_id -> (L, 3) array``)
to keep the scaffold dependency-free. A production system would carry full
all-atom coordinates; every interface here is written so that swap is local.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

# Chain id conventions used throughout.
HEAVY = "H"
LIGHT = "L"
ANTIGEN = "A"

# The six CDRs. Framework regions are everything else and are held fixed to the
# human germline (we only ever design/mutate residues inside these spans).
CDRS = ("H1", "H2", "H3", "L1", "L2", "L3")


@dataclass
class DesignSpec:
    """Everything the generator is conditioned on for one target."""

    target_name: str
    antigen_seq: str
    # C-alpha coords of the antigen, (L_ag, 3). None is allowed for mock runs.
    antigen_coords: Optional[np.ndarray] = None
    # 0-based antigen residue indices that define the epitope (hand-refined,
    # may include GPCR orthosteric anchor residues).
    epitope_resids: tuple[int, ...] = ()
    # Human germline framework the CDRs are grafted onto.
    heavy_framework: str = ""
    light_framework: str = ""
    # Sampled CDR lengths (from OAS length distributions). Maps "H3" -> int etc.
    cdr_lengths: dict[str, int] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.antigen_coords is not None:
            self.antigen_coords = np.asarray(self.antigen_coords, dtype=float)
            if self.antigen_coords.ndim != 2 or self.antigen_coords.shape[1] != 3:
                raise ValueError("antigen_coords must be (L, 3)")


@dataclass
class Confidence:
    """Model-predicted confidence for one complex (Boltz/AF3-style heads).

    Global heads plus optional INTERFACE-restricted heads (epitope x paratope).
    For antibody binders the interface metrics predict binding far better than
    global plDDT, so the scorer prefers them when present (optimisation #1).
    """

    plddt: float          # mean predicted lDDT, higher is better (0-100 or 0-1)
    pae: float            # mean predicted aligned error (A), lower is better
    pde: float            # mean predicted distance error (A), lower is better
    iptm: float = 0.0     # interface pTM proxy, higher is better
    # Optional interface-restricted heads (None => not available).
    interface_pae: Optional[float] = None    # mean PAE on epitope x paratope block
    interface_plddt: Optional[float] = None  # mean plDDT over interface residues
    contact_prob: Optional[float] = None     # mean predicted contact prob at interface

    def has_interface(self) -> bool:
        return self.interface_pae is not None

    def as_dict(self) -> dict[str, float]:
        d = {"plddt": self.plddt, "pae": self.pae, "pde": self.pde, "iptm": self.iptm}
        for k in ("interface_pae", "interface_plddt", "contact_prob"):
            v = getattr(self, k)
            if v is not None:
                d[k] = v
        return d


@dataclass
class Design:
    """A candidate antibody-antigen complex produced by the generator."""

    design_id: str
    heavy_seq: str
    light_seq: str
    # C-alpha traces per chain, e.g. {"H": (Lh,3), "L": (Ll,3), "A": (La,3)}.
    backbone: dict[str, np.ndarray]
    # 0-based residue spans of each CDR within its chain sequence.
    cdr_spans: dict[str, tuple[int, int]] = field(default_factory=dict)
    # Confidence attached at generation time (step 1/2); may be refined later.
    confidence: Optional[Confidence] = None
    parent_id: Optional[str] = None       # beam-search lineage
    provenance: str = "generator"         # "generator" | "ctmc" | ...

    def antibody_ca(self) -> np.ndarray:
        """Concatenated heavy+light C-alpha trace, (Lh+Ll, 3)."""
        parts = [self.backbone[c] for c in (HEAVY, LIGHT) if c in self.backbone]
        return np.concatenate(parts, axis=0)

    def cdr_sequence(self, cdr: str) -> str:
        chain_seq = self.heavy_seq if cdr.startswith("H") else self.light_seq
        lo, hi = self.cdr_spans[cdr]
        return chain_seq[lo:hi]


@dataclass
class Pose:
    """Result of (re)folding a design's sequence: coords + confidence."""

    backbone: dict[str, np.ndarray]
    confidence: Confidence


@dataclass
class ScoredDesign:
    """A design carried through filtering/ranking with all its metrics."""

    design: Design
    confidence: Confidence
    recovery_rmsd: float               # self-consistency, A (lower is better)
    liabilities: dict[str, int]        # motif -> count
    clashes: int                       # physics proxy
    passed_filters: bool
    score: float                       # final ranking score (higher is better)
    cluster_id: int = -1
    is_representative: bool = False

    def summary(self) -> dict:
        return {
            "design_id": self.design.design_id,
            "provenance": self.design.provenance,
            "parent_id": self.design.parent_id,
            "heavy_seq": self.design.heavy_seq,
            "light_seq": self.design.light_seq,
            "confidence": self.confidence.as_dict(),
            "recovery_rmsd": round(self.recovery_rmsd, 3),
            "liabilities": self.liabilities,
            "clashes": self.clashes,
            "passed_filters": self.passed_filters,
            "score": round(self.score, 4),
            "cluster_id": self.cluster_id,
            "is_representative": self.is_representative,
        }
