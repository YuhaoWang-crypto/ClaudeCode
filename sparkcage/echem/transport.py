"""
Analyte mass transport to the SparkCage electrode + surface capture kinetics.

Open-source stand-in for the COMSOL "Transport of Diluted Species" +
surface-reaction physics named in the quote.  One spatial dimension normal to a
planar electrode is the honest geometry for a flat working electrode in a
quiescent or well-stirred cell; 2D/3D FEA only becomes necessary for recessed
electrodes, microelectrode arrays, or a flowing microfluidic channel (see
FEASIBILITY.md, section on when COMSOL is actually required).

Model
-----
    dc/dt = D d2c/dx2                              0 < x < L
    c(L, t) = c_bulk                               (bulk reservoir / stirred layer)
    D dc/dx |_{x=0} = Gamma_max dphi/dt            (flux consumed by capture)
    dphi/dt = k_on c(0,t) (1 - phi) - k_off phi    (Langmuir surface binding)

phi is the fraction of immobilised SparkCage molecules with analyte bound; it
feeds straight into the thermodynamic switch model, whose equilibrium answer
this reproduces as t -> infinity.

Also provides the Squires-Manalis style diagnostic that tells you whether a
given sensor is reaction-limited (kinetics set the response time, geometry does
not matter, 1D is plenty) or diffusion-limited (transport sets it, and electrode
geometry / flow starts to matter).
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

N_A = 6.02214076e23


def geometric_grid(length_cm=0.05, n=160, growth=1.06):
    """Non-uniform grid, fine at the electrode where the depletion layer forms."""
    w = growth ** np.arange(n)
    x = np.concatenate([[0.0], np.cumsum(w)])
    return x / x[-1] * length_cm


def simulate(c_bulk_m, k_on=1.0e6, k_off=1.0e-3, gamma_max_mol_cm2=2.0e-11,
             d_cm2_s=1.0e-6, length_cm=0.05, t_end_s=600.0, n_grid=160,
             n_out=400):
    """Integrate transport + capture. Concentrations in mol/L, k_on in 1/(M s).

    Returns dict with t (s), phi(t), surface concentration c0(t), and the
    final profile.
    """
    x = geometric_grid(length_cm, n_grid)
    n = len(x)
    h = np.diff(x)
    # litre <-> cm^3 : 1 L = 1000 cm^3, so mol/L = 1e-3 mol/cm^3
    c_bulk = c_bulk_m * 1.0e-3

    def rhs(_t, y):
        c = y[:n]
        phi = y[n]
        dphi = k_on * (c[0] * 1.0e3) * (1.0 - phi) - k_off * phi  # c[0] back to mol/L
        dc = np.zeros(n)
        # surface control volume: half of first cell
        j_in = d_cm2_s * (c[1] - c[0]) / h[0]
        dc[0] = (j_in - gamma_max_mol_cm2 * dphi) / (0.5 * h[0])
        # interior nodes, non-uniform second derivative (finite volume)
        flux = d_cm2_s * (c[1:] - c[:-1]) / h
        dc[1:-1] = (flux[1:] - flux[:-1]) / (0.5 * (h[1:] + h[:-1]))
        dc[-1] = 0.0  # Dirichlet bulk
        return np.concatenate([dc, [dphi]])

    y0 = np.concatenate([np.full(n, c_bulk), [0.0]])
    t_eval = np.linspace(0.0, t_end_s, n_out)
    sol = solve_ivp(rhs, (0.0, t_end_s), y0, t_eval=t_eval, method="BDF",
                    rtol=1e-8, atol=1e-16)
    phi = sol.y[n]
    return {
        "t": sol.t,
        "phi": phi,
        "c_surface_M": sol.y[0] * 1.0e3,
        "x_cm": x,
        "c_final_M": sol.y[:n, -1] * 1.0e3,
        "phi_eq": float(phi[-1]),
        "success": bool(sol.success),
    }


def time_to_fraction(res, frac=0.9):
    """Time to reach `frac` of the final surface occupancy (seconds)."""
    phi = res["phi"]
    target = frac * phi[-1]
    idx = np.searchsorted(phi, target)
    if idx == 0 or idx >= len(phi):
        return float("nan")
    t0, t1 = res["t"][idx - 1], res["t"][idx]
    p0, p1 = phi[idx - 1], phi[idx]
    return float(t0 + (target - p0) * (t1 - t0) / (p1 - p0))


def damkohler(k_on=1.0e6, gamma_max_mol_cm2=5.0e-12, d_cm2_s=1.0e-6,
              length_cm=0.05):
    """Ratio of intrinsic capture rate to diffusive supply rate.

    Both rates are first order in analyte concentration, so Da is CONCENTRATION
    INDEPENDENT -- it is a property of the sensor and the cell geometry, not of
    the sample.

    Da << 1 : reaction-limited. Transport is irrelevant, a well-mixed ODE is
              adequate, and electrode geometry does not matter.
    Da >> 1 : diffusion-limited. A depletion layer forms in front of the
              electrode, the response time is set by transport, and stirring,
              flow or electrode geometry become first-order design variables.
              This is the regime where 2D/3D FEA earns its cost.
    """
    capture = k_on * gamma_max_mol_cm2                     # cm/s, per unit conc
    supply = d_cm2_s * 1.0e-3 / length_cm                  # cm/s, per unit conc
    return float(capture / supply)


def diffusion_limited_arrival_time(c_bulk_m, gamma_max_mol_cm2=5.0e-12,
                                   d_cm2_s=1.0e-6, occupancy=0.1):
    """Delivery time in a QUIESCENT, SEMI-INFINITE sample (seconds).

    Even an infinitely fast binder cannot report until diffusion has delivered
    enough analyte to occupy the fraction `occupancy` of the immobilised sensor
    monolayer.  For semi-infinite planar diffusion the cumulative delivery is
    2 c sqrt(D t / pi), so

        t_floor = pi ( occupancy * Gamma_max / (2 c) )^2 / D

    Note this is quadratic in 1/c: dropping the analyte concentration by 10x
    costs 100x in response time.  That scaling, not the electrochemistry, is
    what limits a planar affinity sensor at low picomolar concentrations, and
    it is the single most useful number the transport model produces.

    IMPORTANT: this is the answer for a truly quiescent, unbounded sample.  It
    is NOT a lower bound on `simulate`, which places a well-mixed reservoir at
    x = length_cm and therefore delivers analyte at a constant rate D c / L once
    the profile is established, instead of the decaying t^-1/2 flux assumed here.
    The two idealisations bracket the real answer, and the gap between them is
    exactly the hydrodynamic question a device-level simulation has to settle.
    """
    c = c_bulk_m * 1.0e-3  # mol/cm^3
    return float(np.pi * (occupancy * gamma_max_mol_cm2 / (2.0 * c)) ** 2 / d_cm2_s)
