r"""Sampling-time guidance potentials for the flow ODE (optimisation #3).

The headline potential is FR-zero-contact: penalise contacts between FRAMEWORK
residues and the antigen so that generation is actively steered toward designs
where only the six CDRs enter the pocket. CDR residues are never penalised (they
are meant to bind).

Used through ``interpolants.sample_ode(vfield, x0, cond, guidance=grad,
guidance_scale=s)``, which integrates dx/dt = v_theta - s * grad U. Descending U
pushes framework residues off the antigen.

CAVEAT (stated in the plan): "framework must not touch the antigen" is a
HYPOTHESIS. Some epitopes need framework contacts for affinity, so ``scale`` is a
knob to tune on the IL11 / beta-NGF benchmark, not a hard law. At scale=0 the
guidance is inert.
"""

from __future__ import annotations

from typing import Callable

import numpy as np


def framework_mask_from_spans(length: int, cdr_spans: dict) -> np.ndarray:
    """Boolean (L,) True at FRAMEWORK positions (everything outside the CDRs)."""
    m = np.ones(length, dtype=bool)
    for lo, hi in cdr_spans.values():
        m[lo:hi] = False
    return m


def fr_contact_potential(x: np.ndarray, antigen: np.ndarray, fr_mask: np.ndarray,
                         contact: float = 8.0) -> float:
    r"""U = 0.5 * sum_{f in FR, a} [max(0, contact - d_fa)]^2.

    Zero when no framework residue is within ``contact`` A of any antigen atom;
    grows as framework residues intrude. CDR residues (fr_mask False) are ignored.
    """
    fx = x[fr_mask]
    if fx.size == 0 or antigen.size == 0:
        return 0.0
    d = np.linalg.norm(fx[:, None, :] - antigen[None, :, :], axis=-1)
    viol = np.clip(contact - d, 0.0, None)
    return float(0.5 * (viol ** 2).sum())


def fr_contact_grad(x: np.ndarray, antigen: np.ndarray, fr_mask: np.ndarray,
                    contact: float = 8.0, eps: float = 1e-6) -> np.ndarray:
    r"""grad_x U. Shape (L, 3); nonzero only on framework residues in violation.

    dU/dx_f = sum_a -(contact - d_fa)_+ * (x_f - x_a)/d_fa  (points toward the
    antigen), so subtracting it in the ODE moves framework residues away.
    """
    grad = np.zeros_like(x)
    idx = np.where(fr_mask)[0]
    if idx.size == 0 or antigen.size == 0:
        return grad
    fx = x[idx]
    diff = fx[:, None, :] - antigen[None, :, :]          # (F, A, 3)
    d = np.linalg.norm(diff, axis=-1) + eps              # (F, A)
    viol = np.clip(contact - d, 0.0, None)               # (F, A)
    coeff = -(viol / d)                                  # (F, A)
    g = (coeff[..., None] * diff).sum(axis=1)            # (F, 3)
    grad[idx] = g
    return grad


class FRContactGuidance:
    """A GuidanceGrad closure for sample_ode.

    ``cond`` passed to the closure may override antigen/fr_mask, but defaults are
    captured here so it drops straight into ``sample_ode``.
    """

    def __init__(self, antigen_coords: np.ndarray, fr_mask: np.ndarray,
                 contact: float = 8.0):
        self.antigen = np.asarray(antigen_coords, float)
        self.fr_mask = np.asarray(fr_mask, bool)
        self.contact = contact

    def __call__(self, x: np.ndarray, t: float, cond=None) -> np.ndarray:
        return fr_contact_grad(x, self.antigen, self.fr_mask, self.contact)

    def potential(self, x: np.ndarray) -> float:
        return fr_contact_potential(x, self.antigen, self.fr_mask, self.contact)


def count_fr_contacts(x: np.ndarray, antigen: np.ndarray, fr_mask: np.ndarray,
                      contact: float = 8.0) -> int:
    """Number of framework-residue / antigen pairs within ``contact`` A."""
    fx = x[fr_mask]
    if fx.size == 0 or antigen.size == 0:
        return 0
    d = np.linalg.norm(fx[:, None, :] - antigen[None, :, :], axis=-1)
    return int((d < contact).sum())


def guidance_for_design(design, spec, contact: float = 8.0) -> Callable:
    """Build FR-contact guidance for a Design + DesignSpec (antigen from spec or
    the design backbone). Framework = positions outside the heavy/light CDRs."""
    from ..types import ANTIGEN, HEAVY

    antigen = (np.asarray(spec.antigen_coords, float)
               if spec.antigen_coords is not None
               else design.backbone.get(ANTIGEN, np.zeros((0, 3))))
    # Guidance here is illustrated on the heavy chain (extend per chain as needed).
    heavy_len = len(design.backbone[HEAVY])
    h_spans = {k: v for k, v in design.cdr_spans.items() if k.startswith("H")}
    fr_mask = framework_mask_from_spans(heavy_len, h_spans)
    return FRContactGuidance(antigen, fr_mask, contact)
