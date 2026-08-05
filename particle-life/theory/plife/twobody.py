r"""
Exact small-cluster theory: where self-propulsion comes from.
============================================================

The two-particle problem of Particle Life is exactly solvable and it already
contains the mechanism behind every moving structure in the simulator.

Take one particle of type *a* at :math:`x_a` and one of type *b* at :math:`x_b`,
:math:`r=|x_b-x_a|`, :math:`\hat r=(x_b-x_a)/r`.  The equations of motion are

.. math::
    \ddot x_a = 60\kappa\,f_{ab}(r)\,\hat r - \gamma \dot x_a ,\qquad
    \ddot x_b = -60\kappa\,f_{ba}(r)\,\hat r - \gamma \dot x_b .

Split into centre of mass :math:`X=(x_a+x_b)/2` and separation:

.. math::
    \ddot{\vec r} &= -120\kappa\,f_S(r)\,\hat r - \gamma\dot{\vec r},
        &f_S=\tfrac12\!\left(f_{ab}+f_{ba}\right)\quad\text{(reciprocal part)}\\
    \ddot X &= 60\kappa\,f_A(r)\,\hat r - \gamma\dot X,
        &f_A=\tfrac12\!\left(f_{ab}-f_{ba}\right)\quad\text{(non-reciprocal part)}

Two exact statements follow immediately.

**Bond length is set by the reciprocal part.**  The separation relaxes to
:math:`r^\*` with :math:`f_S(r^\*)=0` and :math:`f_S'(r^\*)>0`.  When both
:math:`A_{ab},A_{ba}>0` and the radii agree this is just :math:`r^\*=\alpha`:
the min-radius matrix *is* the bond-length matrix.

**Speed is set by the non-reciprocal part.**  Once the bond has relaxed, the
centre of mass drifts at the *constant* terminal velocity

.. math::  \boxed{\;\vec V = \frac{60\kappa}{\gamma}\,f_A(r^\*)\;\hat r\;}

along its own (spontaneously chosen) axis.  A bound pair with
:math:`A_{ab}\neq A_{ba}` is a **self-propelled dimer** -- an internal force that
does not cancel, because Newton's third law does not hold.  Every crawling snake,
every chaser, every migrating cell in Particle Life is this mechanism repeated
inside a bigger cluster:

.. math::
    \vec F^{\rm int}_{\rm cluster}=\sum_{i<j\in C} 2 f_A(r_{ij})\,\hat r_{ij},
    \qquad
    \tau^{\rm int}_{\rm cluster}=\sum_{i<j\in C}
        \bigl[(x_i-X)f_{ij}-(x_j-X)f_{ji}\bigr]\times\hat r_{ij},

both of which vanish identically for a reciprocal matrix and are generically
non-zero otherwise.  Non-zero net internal force = translation; non-zero net
internal torque = rotation (the spinning cells produced by cyclic
"rock-paper-scissors" matrices).
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from .kernels import pair_force

__all__ = [
    "f_sym", "f_asym", "bond_length_pair", "dimer",
    "cluster_internal_force", "cluster_internal_torque",
]


# --------------------------------------------------------------------------- #
def _pf(r, params, a, b):
    return pair_force(r, params.A[a, b], params.alpha[a, b], params.beta[a, b],
                      params.repel)


def f_sym(r, params, a, b):
    r"""Reciprocal part :math:`f_S=\tfrac12(f_{ab}+f_{ba})` -- sets the bond length."""
    return 0.5 * (_pf(r, params, a, b) + _pf(r, params, b, a))


def f_asym(r, params, a, b):
    r"""Non-reciprocal part :math:`f_A=\tfrac12(f_{ab}-f_{ba})` -- sets the speed."""
    return 0.5 * (_pf(r, params, a, b) - _pf(r, params, b, a))


# --------------------------------------------------------------------------- #
def bond_length_pair(params, a, b, n=4001):
    r"""Stable separation :math:`r^\*`: the largest root of :math:`f_S` with
    :math:`f_S'>0` inside the interaction range.  ``None`` if the pair is unbound.
    """
    rmax = max(params.beta[a, b], params.beta[b, a])
    r = np.linspace(1e-6, rmax * 0.999, n)
    f = f_sym(r, params, a, b)
    roots = []
    for i in range(n - 1):
        if f[i] == 0.0 and 0 < i:
            roots.append(r[i])
        elif f[i] * f[i + 1] < 0:
            try:
                roots.append(brentq(lambda x: f_sym(np.array([x]), params, a, b)[0],
                                    r[i], r[i + 1]))
            except ValueError:
                pass
    stable = []
    for rr in roots:
        h = max(1e-4, rr * 1e-4)
        d = (f_sym(np.array([rr + h]), params, a, b)[0]
             - f_sym(np.array([rr - h]), params, a, b)[0]) / (2 * h)
        if d > 0:                       # restoring
            stable.append(rr)
    return max(stable) if stable else None


def dimer(params, a, b, discrete=True):
    r"""Full prediction for the (a, b) dimer.

    Returns the bond length ``r_star``, the terminal self-propulsion speed ``V``,
    the bond stiffness, the internal oscillation frequency and whether the bond is
    under- or over-damped.  ``discrete=True`` uses the exact fixed point of the
    simulator's finite-``dt`` update (``Params.mobility_discrete``); set it False
    for the continuum value ``60 kappa / gamma``.

    Degenerate case worth knowing: if the two radii matrices are symmetric
    (:math:`\alpha_{ab}=\alpha_{ba}`, :math:`\beta_{ab}=\beta_{ba}`) then
    :math:`f_{ab}` and :math:`f_{ba}` share their zeros, so the bond length
    :math:`r^\*=\alpha` is *also* a zero of :math:`f_A` and ``V = 0``: a two-body
    chase always stalls.  Sustained self-propulsion therefore needs either
    asymmetric radii or >= 3 particles, where geometric frustration prevents every
    bond from sitting at its own zero (see ``cluster_internal_force``).
    """
    rs = bond_length_pair(params, a, b)
    out = dict(a=a, b=b, r_star=rs, bound=rs is not None)
    if rs is None:
        out.update(V=0.0, stiffness=np.nan, omega0=np.nan, underdamped=False)
        return out

    kappa, gamma = params.force_factor, params.gamma
    mob = params.mobility_discrete if discrete else params.mobility
    V = mob * float(f_asym(np.array([rs]), params, a, b)[0])

    h = max(1e-4, rs * 1e-4)
    dfS = float(
        (f_sym(np.array([rs + h]), params, a, b)[0]
         - f_sym(np.array([rs - h]), params, a, b)[0]) / (2 * h)
    )
    stiff = 120.0 * kappa * dfS                      # d(-r_ddot)/dr for the separation
    omega0 = float(np.sqrt(max(stiff, 0.0)))
    out.update(V=V, stiffness=stiff, omega0=omega0,
               underdamped=bool(omega0 > gamma / 2))
    return out


# --------------------------------------------------------------------------- #
def cluster_internal_force(pos, types, params, box=None):
    r"""Net **internal** force of a set of particles, :math:`\sum_{i<j}2f_A\hat r_{ij}`.

    Exactly zero for a reciprocal matrix (Newton's third law).  Otherwise it is
    the cluster's self-propulsion drive: terminal velocity ``F_int * 60 kappa /
    (gamma * N)``.
    """
    pos = np.asarray(pos, float)
    types = np.asarray(types, int)
    n = len(pos)
    F = np.zeros(2)
    for i in range(n):
        for j in range(i + 1, n):
            d = pos[j] - pos[i]
            if box is not None:
                d -= np.asarray(box) * np.round(d / np.asarray(box))
            r = np.hypot(*d)
            if r < 1e-12:
                continue
            u = d / r
            fa = float(f_asym(np.array([r]), params, types[i], types[j])[0])
            F += 2.0 * fa * u
    return F


def largest_cluster(sim, cutoff_scale=1.35, max_size=400, min_size=4):
    """Indices of the biggest connected cluster that still fits in ``max_size``.

    Truncating a bigger cluster would corrupt the internal-force sum (the
    prediction divides by the cluster size), so an oversized cluster is skipped
    in favour of the largest one that can be summed exactly.  Returns an empty
    array if nothing qualifies.
    """
    from .observables import clusters

    box = np.asarray(sim.p.L, float)
    cutoff = cutoff_scale * float(np.median(sim.p.alpha))
    labels, _ = clusters(sim.pos, box, cutoff)
    counts = np.bincount(labels)
    ok = np.where((counts >= min_size) & (counts <= max_size))[0]
    if ok.size == 0:
        return np.array([], dtype=int)
    best = ok[int(np.argmax(counts[ok]))]
    return np.where(labels == best)[0]


def cluster_prediction(sim, cutoff_scale=1.35, max_size=400):
    r"""Predict a real cluster's drift and spin from its **internal** force/torque.

    For the biggest cluster in the snapshot, unwrapped across the periodic seam:

    .. math::
        V_{\rm pred} = \mu_{\rm disc}\,\frac{|F^{\rm int}|}{N_C},\qquad
        \omega_{\rm pred} = \mu_{\rm disc}\,\frac{\tau^{\rm int}}{I},\quad
        I=\sum_{i\in C}|x_i-X|^2 ,

    and compares against the measured centre-of-mass speed and angular velocity.
    Both predictions are identically zero for a reciprocal interaction.
    """
    p = sim.p
    box = np.asarray(p.L, float)
    idx = largest_cluster(sim, cutoff_scale, max_size=max_size)
    if len(idx) < 4:
        return dict(size=len(idx), V_pred=0.0, V_meas=0.0, omega_pred=0.0,
                    omega_meas=0.0, F_int=0.0, tau_int=0.0)

    pos = sim.pos[idx].copy()
    d = pos - pos[0]
    d -= box * np.round(d / box)                  # unwrap around the first member
    pos = pos[0] + d
    vel = sim.vel[idx]
    types = sim.types[idx]

    X, V = pos.mean(axis=0), vel.mean(axis=0)
    rel, vrel = pos - X, vel - V
    I = float((rel ** 2).sum())
    Lz = float((rel[:, 0] * vrel[:, 1] - rel[:, 1] * vrel[:, 0]).sum())

    F = cluster_internal_force(pos, types, p)
    tau = cluster_internal_torque(pos, types, p)
    mob = p.mobility_discrete
    return dict(
        size=int(len(idx)),
        F_int=float(np.hypot(*F)), tau_int=float(tau),
        V_pred=float(mob * np.hypot(*F) / len(idx)), V_meas=float(np.hypot(*V)),
        omega_pred=float(mob * tau / I) if I > 0 else 0.0,
        omega_meas=float(Lz / I) if I > 0 else 0.0,
    )


def cluster_internal_torque(pos, types, params, box=None):
    r"""Net **internal** torque about the centre of mass (scalar, 2-D).

    Vanishes identically for a reciprocal matrix; non-zero torque is the origin
    of spinning clusters and vortex lattices.
    """
    pos = np.asarray(pos, float)
    types = np.asarray(types, int)
    X = pos.mean(axis=0)
    n = len(pos)
    tau = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            d = pos[j] - pos[i]
            if box is not None:
                d -= np.asarray(box) * np.round(d / np.asarray(box))
            r = np.hypot(*d)
            if r < 1e-12:
                continue
            u = d / r
            fij = float(pair_force(r, params.A[types[i], types[j]],
                                   params.alpha[types[i], types[j]],
                                   params.beta[types[i], types[j]], params.repel))
            fji = float(pair_force(r, params.A[types[j], types[i]],
                                   params.alpha[types[j], types[i]],
                                   params.beta[types[j], types[i]], params.repel))
            ri, rj = pos[i] - X, pos[j] - X
            vec = ri * fij - rj * fji
            tau += vec[0] * u[1] - vec[1] * u[0]
    return tau
