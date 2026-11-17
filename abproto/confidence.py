"""Aggregate a complex's confidence heads into one scalar used for ranking
and for Beam Search / SMC pruning.

Optimisation #1 (interface-weighted): for antibody binders the global plDDT is
weakly predictive of binding; the interface-restricted heads (epitope x paratope
PAE, interface plDDT, contact probability) are what matter. When those are
present we score primarily on them and only lightly on the global terms;
otherwise we fall back to the global-only objective. Higher is better.
"""

from __future__ import annotations

from .types import Confidence

# --- global-only weights (fallback) ---
W_PLDDT = 1.0
W_PAE = 0.15
W_PDE = 0.10
W_IPTM = 0.6

# --- interface weights (preferred when available) ---
WI_IFACE_PLDDT = 1.0
WI_IFACE_PAE = 0.35        # subtracted; interface PAE is the strongest signal
WI_CONTACT = 0.8
WI_IPTM = 0.8
WI_GLOBAL_PLDDT = 0.25     # keep a little global signal for overall foldability


def _norm_plddt(plddt: float) -> float:
    return plddt / 100.0 if plddt > 1.5 else plddt


def _global_score(c: Confidence) -> float:
    return (W_PLDDT * _norm_plddt(c.plddt) + W_IPTM * c.iptm
            - W_PAE * c.pae - W_PDE * c.pde)


def _interface_score(c: Confidence) -> float:
    return (WI_IFACE_PLDDT * _norm_plddt(c.interface_plddt)
            + WI_CONTACT * (c.contact_prob or 0.0)
            + WI_IPTM * c.iptm
            + WI_GLOBAL_PLDDT * _norm_plddt(c.plddt)
            - WI_IFACE_PAE * c.interface_pae)


def confidence_score(conf: Confidence) -> float:
    """Interface-weighted objective when interface heads exist, else global."""
    if conf.has_interface():
        return _interface_score(conf)
    return _global_score(conf)
