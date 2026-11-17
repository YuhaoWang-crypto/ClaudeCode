r"""Continuous Stochastic-Interpolant flow matching for all-atom coordinates.

Reference: Albergo, Boffi & Vanden-Eijnden, "Stochastic Interpolants: A Unifying
Framework for Flows and Diffusions" (2023). We use the linear interpolant with an
optional latent, which subsumes flow matching.

Setup
-----
Let x_1 ~ p_data be a (clean) coordinate set for the designed region (CDR
backbone + side-chain atoms), conditioned on c = (antigen coords, epitope,
germline framework). Let x_0 ~ N(0, I) be the base. Define the interpolant

    x_t = alpha(t) x_1 + beta(t) x_0 + gamma(t) z ,    z ~ N(0, I),  t in [0,1]

with boundary conditions alpha(0)=0, beta(0)=1, alpha(1)=1, beta(1)=0 and
gamma(0)=gamma(1)=0 (the latent vanishes at the endpoints). Linear choice:

    alpha(t) = t ,  beta(t) = 1 - t ,  gamma(t) = c * sqrt(t (1-t)) .

The time derivative of the interpolant (the conditional velocity) is

    d/dt x_t = alpha'(t) x_1 + beta'(t) x_0 + gamma'(t) z
             = x_1 - x_0 + gamma'(t) z        (linear, c>0 adds the latent term).

Learning
--------
A network v_theta(x_t, t, c) is trained to match the *marginal* velocity
b(x,t,c) = E[ d/dt x_t | x_t = x, c ] by regressing on the conditional target
(this expectation is achieved by the per-sample MSE below):

    L_flow = E_{t, x_0, x_1, z} || v_theta(x_t, t, c) - (x_1 - x_0 + gamma'(t) z) ||^2 .

Sampling
--------
Integrate the probability-flow ODE from t=0 (noise) to t=1 (data):

    dx/dt = v_theta(x, t, c) ,     x(0) = x_0 ~ N(0, I).

A conditioning potential U(x, c) (e.g. "framework atoms must not contact the
antigen") can be added as guidance: dx/dt = v_theta - s * grad_x U (see
``sample_ode(..., guidance=...)``).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Optional

import numpy as np

# A velocity field: (x_t, t, cond) -> velocity of the same shape as x_t.
VelocityField = Callable[[np.ndarray, float, Any], np.ndarray]
# A guidance potential's gradient: (x_t, t, cond) -> grad wrt x_t.
GuidanceGrad = Callable[[np.ndarray, float, Any], np.ndarray]


@dataclass
class LinearInterpolant:
    r"""Linear stochastic interpolant x_t = t x_1 + (1-t) x_0 + gamma(t) z."""

    latent_scale: float = 0.0        # c; 0 => deterministic flow-matching interpolant

    # --- schedule coefficients ---
    def alpha(self, t: float) -> float:
        return t

    def beta(self, t: float) -> float:
        return 1.0 - t

    def gamma(self, t: float) -> float:
        if self.latent_scale == 0.0:
            return 0.0
        return self.latent_scale * np.sqrt(max(t * (1.0 - t), 0.0))

    def gamma_dot(self, t: float) -> float:
        if self.latent_scale == 0.0:
            return 0.0
        # d/dt sqrt(t(1-t)) = (1 - 2t) / (2 sqrt(t(1-t)))
        denom = 2.0 * np.sqrt(max(t * (1.0 - t), 1e-12))
        return self.latent_scale * (1.0 - 2.0 * t) / denom

    # --- interpolant path and its conditional velocity ---
    def interpolate(self, x0: np.ndarray, x1: np.ndarray, t: float,
                    z: Optional[np.ndarray] = None) -> np.ndarray:
        xt = self.alpha(t) * x1 + self.beta(t) * x0
        if self.latent_scale != 0.0:
            if z is None:
                z = np.random.standard_normal(x0.shape)
            xt = xt + self.gamma(t) * z
        return xt

    def target_velocity(self, x0: np.ndarray, x1: np.ndarray, t: float,
                        z: Optional[np.ndarray] = None) -> np.ndarray:
        v = x1 - x0
        if self.latent_scale != 0.0:
            if z is None:
                z = np.zeros_like(x0)
            v = v + self.gamma_dot(t) * z
        return v


def flow_matching_loss(vfield: VelocityField, interp: LinearInterpolant,
                       x0: np.ndarray, x1: np.ndarray, t: float, cond: Any = None,
                       z: Optional[np.ndarray] = None) -> float:
    r"""Per-sample flow-matching loss ||v_theta(x_t,t,c) - target||^2 (mean over
    coordinates). Averaging this over (t, x0, x1, z) is the Monte-Carlo estimate
    of L_flow above."""
    if z is None and interp.latent_scale != 0.0:
        z = np.random.standard_normal(x0.shape)
    xt = interp.interpolate(x0, x1, t, z)
    target = interp.target_velocity(x0, x1, t, z)
    pred = vfield(xt, t, cond)
    diff = pred - target
    return float((diff * diff).mean())


def sample_ode(vfield: VelocityField, x0: np.ndarray, cond: Any = None,
               steps: int = 100, guidance: Optional[GuidanceGrad] = None,
               guidance_scale: float = 1.0,
               trace: Optional[list] = None) -> np.ndarray:
    r"""Integrate dx/dt = v_theta(x,t,c) - s*grad U from t=0 to t=1 (Euler).

    ``guidance`` lets you steer sampling (e.g. the FR-zero-contact potential):
    its gradient is subtracted so the flow descends U while following the model.
    """
    x = np.array(x0, dtype=float)
    dt = 1.0 / steps
    for k in range(steps):
        t = k * dt
        v = vfield(x, t, cond)
        if guidance is not None:
            v = v - guidance_scale * guidance(x, t, cond)
        x = x + dt * v
        if trace is not None:
            trace.append(x.copy())
    return x


# --------------------------------------------------------------------------
# Coupling continuous structure with the discrete sequence CTMC.
# --------------------------------------------------------------------------
def structure_first_time(t: float, delta: float = 0.2) -> float:
    r"""Map global time ``t`` to the sequence's effective time.

        t_seq(t) = clip((t - delta) / (1 - delta), 0, 1)

    With delta>0 the sequence is *more noised* than the structure at any global
    t (its unmasking lags), so the backbone resolves first and the sequence is
    predicted against an increasingly definite structure. Empirically this
    "structure-first" schedule improves co-design self-consistency.
    """
    if not 0.0 <= delta < 1.0:
        raise ValueError("delta must be in [0, 1)")
    return float(np.clip((t - delta) / (1.0 - delta), 0.0, 1.0))


def coupled_loss(flow_term: float, ctmc_term: float, lam: float = 1.0) -> float:
    r"""Total training objective L = L_flow + lam * L_ctmc.

    Both terms share the sampled time t (the structure uses t, the sequence uses
    structure_first_time(t)) and, in production, a shared backbone network that
    emits both the velocity and the sequence posterior. Keeping them in one loss
    is what makes generation *jointly* self-consistent rather than staged.
    """
    return float(flow_term + lam * ctmc_term)
