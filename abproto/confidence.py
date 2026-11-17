"""Aggregate a complex's confidence heads into one scalar used for ranking
and for Beam Search pruning.

This is the "independent confidence model" of the plan, reduced to a scalar
objective. In production the three heads come from the Boltz-2/AF3 confidence
module; here we only combine them. Higher is better.
"""

from __future__ import annotations

from .types import Confidence


# Weights are deliberately explicit so they can be tuned / learned later.
# plDDT is normalised to 0-1 before weighting (accepts either 0-1 or 0-100).
W_PLDDT = 1.0
W_PAE = 0.15      # subtracted (error)
W_PDE = 0.10      # subtracted (error)
W_IPTM = 0.6      # interface quality, the signal that matters most for binders


def _norm_plddt(plddt: float) -> float:
    return plddt / 100.0 if plddt > 1.5 else plddt


def confidence_score(conf: Confidence) -> float:
    """Scalar objective: reward plDDT and iPTM, penalise PAE and PDE."""
    return (
        W_PLDDT * _norm_plddt(conf.plddt)
        + W_IPTM * conf.iptm
        - W_PAE * conf.pae
        - W_PDE * conf.pde
    )
