"""
The Particle Life integrator - a faithful CPU port of the platform's GPU kernel.
===============================================================================

Discrete update (exactly the order of operations in the site's WGSL shaders
``forceCompute.main`` followed by ``particleAdvance``)::

    dv_i   = (kappa * 60 * dt) * sum_j f_{s_i s_j}(r_ij) * rhat_ij     # force pass
    v_i   <- (v_i + dv_i) * (1 - lambda)^(60 dt)                       # drag
    x_i   <- x_i + v_i * dt                                            # advect

Continuum-time form
-------------------
Letting ``dt -> 0`` with ``(1-lambda)^(60 dt) = exp(-gamma dt)``:

.. math::
    \\ddot x_i = 60\\kappa \\sum_j f_{s_i s_j}(r_{ij})\\,\\hat r_{ij} - \\gamma\\,\\dot x_i,
    \\qquad \\gamma \\equiv -60\\ln(1-\\lambda).

So Particle Life is a **damped, non-reciprocal, finite-range Newtonian N-body
system** with unit masses.  With the platform default ``lambda = 0.3`` one gets
``gamma = 21.4 s^-1``: the velocity memory time is 47 ms while structures live
for seconds, so the system is *close to* - but not exactly - overdamped.  That
small residual inertia is what lets chasers overshoot and orbit; set
``inertialess=True`` to integrate the strictly overdamped limit

.. math::  \\dot x_i = (60\\kappa/\\gamma)\\sum_j f\\,\\hat r_{ij}

and see which phenomena survive.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import cKDTree

from .kernels import pair_force, pair_potential

__all__ = ["Params", "ParticleLife"]


# --------------------------------------------------------------------------- #
@dataclass
class Params:
    """Every knob of the platform that has physical meaning."""

    A: np.ndarray                      # (S,S) force matrix, entries in [-1,1]
    alpha: np.ndarray                  # (S,S) min radius  (bond length)
    beta: np.ndarray                   # (S,S) max radius  (interaction cutoff)
    repel: float = 1.0                 # global core-repulsion strength R
    force_factor: float = 1.0          # kappa
    friction: float = 0.3              # lambda, per 1/60 s frame
    dt: float = 1.0 / 60.0
    L: tuple = (600.0, 600.0)          # periodic box
    wrap: bool = True                  # True: torus.  False: reflecting walls
    inertialess: bool = False          # integrate the overdamped limit instead

    def __post_init__(self):
        self.A = np.asarray(self.A, dtype=float)
        self.alpha = np.asarray(self.alpha, dtype=float)
        self.beta = np.asarray(self.beta, dtype=float)
        assert self.A.shape == self.alpha.shape == self.beta.shape
        assert self.A.shape[0] == self.A.shape[1]

    @property
    def S(self) -> int:
        return self.A.shape[0]

    @property
    def gamma(self) -> float:
        """Continuum drag rate ``gamma = -60 ln(1 - lambda)`` [1/s]."""
        return -60.0 * np.log(max(1e-12, 1.0 - self.friction))

    @property
    def mobility(self) -> float:
        """Continuum overdamped mobility ``mu = 60 kappa / gamma`` (velocity per force)."""
        return 60.0 * self.force_factor / self.gamma

    @property
    def mobility_discrete(self) -> float:
        r"""Terminal velocity per unit force of the **discrete** map.

        Iterating ``v <- (v + ff a) d`` with ``ff = 60 kappa dt`` and
        ``d = (1-lambda)^{60 dt}`` has the exact fixed point

        .. math::  v^\* = 60\kappa\,\Delta t\,\frac{d}{1-d}\,a .

        As ``dt -> 0`` this tends to ``60 kappa / gamma``; at the platform default
        (``dt = 1/60``, ``lambda = 0.3``) it is 2.3333 against a continuum value of
        2.8037, i.e. the simulator runs 16.8 % "slower" than the continuum theory.
        Use this when comparing predictions against the actual simulator.
        """
        d = (1.0 - self.friction) ** (60.0 * self.dt)
        return 60.0 * self.force_factor * self.dt * d / (1.0 - d)

    @property
    def rmax(self) -> float:
        return float(self.beta.max())

    def copy(self, **kw):
        d = dict(A=self.A.copy(), alpha=self.alpha.copy(), beta=self.beta.copy(),
                 repel=self.repel, force_factor=self.force_factor,
                 friction=self.friction, dt=self.dt, L=self.L, wrap=self.wrap,
                 inertialess=self.inertialess)
        d.update(kw)
        return Params(**d)


# --------------------------------------------------------------------------- #
class ParticleLife:
    """N-body Particle Life simulator (vectorised, cKDTree neighbour lists)."""

    def __init__(self, params: Params, N: int = 4000, seed: int = 0,
                 init: str = "random", types: np.ndarray | None = None):
        self.p = params
        self.N = int(N)
        self.rng = np.random.default_rng(seed)
        Lx, Ly = params.L

        if types is None:
            self.types = self.rng.integers(0, params.S, size=self.N)
        else:
            self.types = np.asarray(types, dtype=int)

        self.pos = self._init_positions(init, Lx, Ly)
        self.vel = np.zeros((self.N, 2))
        self.t = 0.0
        self.step_count = 0

        # flattened per-ordered-pair parameters, indexed by s_i * S + s_j
        S = params.S
        self._Af = params.A.reshape(-1)
        self._alf = params.alpha.reshape(-1)
        self._bef = params.beta.reshape(-1)
        self._S = S

    # ------------------------------------------------------------------ #
    def _init_positions(self, init, Lx, Ly):
        r = self.rng
        if init == "random":
            return np.column_stack([r.uniform(0, Lx, self.N), r.uniform(0, Ly, self.N)])
        if init == "disk":
            th = r.uniform(0, 2 * np.pi, self.N)
            rad = 0.46 * min(Lx, Ly) * np.sqrt(r.uniform(0, 1, self.N))
            return np.column_stack([Lx / 2 + rad * np.cos(th), Ly / 2 + rad * np.sin(th)])
        if init == "ring":
            th = r.uniform(0, 2 * np.pi, self.N)
            rad = 0.42 * min(Lx, Ly) * (1 - 0.06 * r.uniform(0, 1, self.N))
            return np.column_stack([Lx / 2 + rad * np.cos(th), Ly / 2 + rad * np.sin(th)])
        if init == "grid":
            n = int(np.ceil(np.sqrt(self.N)))
            gx, gy = np.meshgrid(np.linspace(0, Lx, n, endpoint=False),
                                 np.linspace(0, Ly, n, endpoint=False))
            pts = np.column_stack([gx.ravel(), gy.ravel()])[: self.N]
            return pts + r.normal(0, 0.5, pts.shape)
        raise ValueError(f"unknown init {init!r}")

    # ------------------------------------------------------------------ #
    def pairs(self):
        """Neighbour pairs ``(i, j)`` with ``i < j`` and ``|x_i - x_j| < rmax``."""
        p = self.p
        if p.wrap:
            box = np.array(p.L, dtype=float)
            pos = np.mod(self.pos, box)
            tree = cKDTree(pos, boxsize=box)
        else:
            tree = cKDTree(self.pos)
        return tree.query_pairs(p.rmax, output_type="ndarray")

    def accelerations(self):
        r"""Net interaction force per particle, :math:`\sum_j f\,\hat r_{ij}` (mass = 1).

        Both ordered directions of every pair are evaluated separately - this is
        precisely where the third law is broken.
        """
        p, N, S = self.p, self.N, self._S
        acc = np.zeros((N, 2))
        ij = self.pairs()
        if ij.size == 0:
            return acc
        i, j = ij[:, 0], ij[:, 1]

        d = self.pos[j] - self.pos[i]
        if p.wrap:
            box = np.array(p.L, dtype=float)
            d -= box * np.round(d / box)
        r = np.hypot(d[:, 0], d[:, 1])
        ok = r > 1e-9
        i, j, d, r = i[ok], j[ok], d[ok], r[ok]
        if r.size == 0:
            return acc
        u = d / r[:, None]                      # unit vector i -> j

        si, sj = self.types[i], self.types[j]
        k_ij = si * S + sj                      # force felt by i due to j
        k_ji = sj * S + si                      # force felt by j due to i

        f_i = pair_force(r, self._Af[k_ij], self._alf[k_ij], self._bef[k_ij], p.repel)
        f_j = pair_force(r, self._Af[k_ji], self._alf[k_ji], self._bef[k_ji], p.repel)

        fx_i, fy_i = f_i * u[:, 0], f_i * u[:, 1]
        fx_j, fy_j = -f_j * u[:, 0], -f_j * u[:, 1]

        acc[:, 0] = (np.bincount(i, fx_i, minlength=N) + np.bincount(j, fx_j, minlength=N))
        acc[:, 1] = (np.bincount(i, fy_i, minlength=N) + np.bincount(j, fy_j, minlength=N))
        return acc

    # ------------------------------------------------------------------ #
    def step(self, n: int = 1):
        """Advance ``n`` frames using the platform's exact update order."""
        p = self.p
        box = np.array(p.L, dtype=float)
        drag = (1.0 - p.friction) ** (60.0 * p.dt)
        ff = p.force_factor * 60.0 * p.dt

        for _ in range(n):
            acc = self.accelerations()
            if p.inertialess:
                self.vel = p.mobility * acc
            else:
                self.vel = (self.vel + ff * acc) * drag
            self.pos = self.pos + self.vel * p.dt

            if p.wrap:
                self.pos = np.mod(self.pos, box)
            else:                                      # reflecting walls
                for ax in (0, 1):
                    lo = self.pos[:, ax] < 0
                    hi = self.pos[:, ax] > box[ax]
                    self.vel[lo | hi, ax] *= -1.0
                    self.pos[:, ax] = np.clip(self.pos[:, ax], 0, box[ax])

            self.t += p.dt
            self.step_count += 1
        return self

    def run(self, seconds: float):
        return self.step(int(round(seconds / self.p.dt)))

    # ------------------------------------------------------------------ #
    #  diagnostics
    # ------------------------------------------------------------------ #
    def potential_energy(self):
        r"""Total pair potential :math:`V=\sum_{i<j} \tfrac12 [U_{ab}+U_{ba}]`.

        For a **reciprocal** matrix this is an exact Lyapunov functional of the
        overdamped dynamics (:math:`\dot V \le 0`).  For a non-reciprocal matrix
        the same expression is *not* conserved or monotone - measuring its
        non-monotonicity is a direct test of non-variational dynamics.
        """
        p, S = self.p, self._S
        ij = self.pairs()
        if ij.size == 0:
            return 0.0
        i, j = ij[:, 0], ij[:, 1]
        d = self.pos[j] - self.pos[i]
        if p.wrap:
            box = np.array(p.L, dtype=float)
            d -= box * np.round(d / box)
        r = np.hypot(d[:, 0], d[:, 1])
        si, sj = self.types[i], self.types[j]
        k_ij, k_ji = si * S + sj, sj * S + si
        u_ij = pair_potential(r, self._Af[k_ij], self._alf[k_ij], self._bef[k_ij], p.repel)
        u_ji = pair_potential(r, self._Af[k_ji], self._alf[k_ji], self._bef[k_ji], p.repel)
        return float(0.5 * (u_ij + u_ji).sum())

    def kinetic_energy(self):
        return float(0.5 * (self.vel ** 2).sum())

    def injected_power(self):
        r"""Power pumped in by the **antisymmetric** part of the interaction.

        Split each ordered force into reciprocal and non-reciprocal halves,
        :math:`f_{ab}=f_S\pm f_A` with :math:`f_S=\tfrac12(f_{ab}+f_{ba})`,
        :math:`f_A=\tfrac12(f_{ab}-f_{ba})`.  The reciprocal half is exactly
        :math:`-\nabla V`, so the only non-variational contribution to
        :math:`\dot E` is

        .. math::
            P_{\rm inj} \;=\; 60\kappa \sum_{i<j} f_A(r_{ij})\,
                              \hat r_{ij}\cdot(v_i+v_j) ,

        which vanishes **identically** for a reciprocal matrix.  In a
        non-reciprocal steady state :math:`\langle P_{\rm inj}\rangle` must equal
        the frictional dissipation :math:`\gamma\sum_i|v_i|^2`: the broken third
        law is literally the engine that pays for the motion.
        """
        p, S = self.p, self._S
        ij = self.pairs()
        if ij.size == 0:
            return 0.0
        i, j = ij[:, 0], ij[:, 1]
        d = self.pos[j] - self.pos[i]
        if p.wrap:
            box = np.array(p.L, dtype=float)
            d -= box * np.round(d / box)
        r = np.hypot(d[:, 0], d[:, 1])
        ok = r > 1e-9
        i, j, d, r = i[ok], j[ok], d[ok], r[ok]
        if r.size == 0:
            return 0.0
        u = d / r[:, None]
        si, sj = self.types[i], self.types[j]
        k_ij, k_ji = si * S + sj, sj * S + si
        f_ij = pair_force(r, self._Af[k_ij], self._alf[k_ij], self._bef[k_ij], p.repel)
        f_ji = pair_force(r, self._Af[k_ji], self._alf[k_ji], self._bef[k_ji], p.repel)
        f_a = 0.5 * (f_ij - f_ji)
        vsum = self.vel[i] + self.vel[j]
        return float(60.0 * p.force_factor * (f_a * (u * vsum).sum(axis=1)).sum())

    def dissipated_power(self):
        r"""Frictional dissipation :math:`\gamma\sum_i|v_i|^2`."""
        return float(self.p.gamma * (self.vel ** 2).sum())

    def mean_speed(self):
        return float(np.hypot(self.vel[:, 0], self.vel[:, 1]).mean())

    def snapshot(self):
        return dict(pos=self.pos.copy(), vel=self.vel.copy(),
                    types=self.types.copy(), t=self.t)


# --------------------------------------------------------------------------- #
def default_params(S=6, seed=0, L=600.0, **kw):
    """Platform-default universe: 7 species, alpha in [12,24], beta in [32,64]."""
    from .matrices import generate, random_radii

    rng = np.random.default_rng(seed)
    A = generate(kw.pop("matrix", "random"), S, rng)
    alpha, beta = random_radii(S, rng=rng)
    return Params(A=A, alpha=alpha, beta=beta, L=(L, L), **kw)
