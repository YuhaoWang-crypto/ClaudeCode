"""
Square-wave voltammetry of a SURFACE-CONFINED methylene-blue reporter.

This is the open-source stand-in for quote item Phase I / Step 2:
"Computational modeling and simulation of SparkCage biosensors ... finite
element analysis (FEA) software (COMSOL Multiphysics 5.6)".

Physics implemented
-------------------
The methylene blue (MB) reporter is covalently tethered to the SparkCage latch,
so it is an ADSORBED (surface-confined) redox couple, not a diffusing one.  Its
voltammetry is therefore governed by a surface coverage ODE with Butler-Volmer
interfacial kinetics -- no bulk transport term at all:

    theta = fraction of MB in the reduced (leuco-MB) form
    d(theta)/dt = k_c (1 - theta) - k_a theta

    k_c = k0 exp( -alpha  * n F/RT * (E - E0) )      (reduction)
    k_a = k0 exp( (1-alpha) * n F/RT * (E - E0) )    (oxidation)

The electron count n belongs in the exponent, not only in the charge prefactor.
This is Laviron's surface-confined form, and it is why a two-electron couple
gives a peak about half as wide as a one-electron one (roughly 90/n mV at half
height in the reversible limit).  Leaving n out of the exponent -- which this
module did until the Biomni report was reproduced -- makes n a pure scale factor
and silently predicts a one-electron peak shape for methylene blue.

    i = n F A Gamma d(theta)/dt

The SENSING mechanism is that k0 depends exponentially on the electron-tunnelling
distance d between MB and the electrode surface (Marcus / superexchange):

    k0(d) = k0_contact * exp( -beta * (d - d_contact) ),   beta ~ 1.0-1.4 / Angstrom

When the SparkCage latch is docked, MB is held near the carbon surface (fast
electron transfer, large square-wave peak).  When analyte binding releases the
latch, MB is displaced (slow electron transfer, small peak).  The measured peak
current is therefore a direct readout of the fraction of switched sensors that
the thermodynamic model in ../thermo/switch_model.py predicts.

Numerics
--------
During each half-cycle of the square wave the potential is constant, so the
coverage ODE has a closed-form solution:

    theta_ss = k_c / (k_c + k_a)
    theta(t) = theta_ss + (theta_0 - theta_ss) exp( -(k_c + k_a) t )
    d(theta)/dt = (k_c + k_a) (theta_ss - theta(t))

The waveform is integrated exactly, half-cycle by half-cycle.  This is the same
answer a finite-element solver returns for a surface species, computed without
one, because there is no spatial dimension in the problem.

Sign convention: reduction current is reported positive.

Validation anchors
------------------
Dauphin-Ducharme, Arroyo-Curras, Kurnik, Ortega, Li, Plaxco, "Simulation-Based
Approach to Determining Electron Transfer Rates Using Square-Wave Voltammetry",
Langmuir 2017, 33, 4407-4413, doi:10.1021/acs.langmuir.7b00359.  That paper is a
COMSOL 5.2 model of this exact problem.  Notably, it has no way to represent a
surface-confined species directly, so it approximates the monolayer with a
150 um thin-layer cell chosen so that the diffusion layer never reaches the far
wall.  The closed-form treatment used here is not an approximation to that
model; it is the exact answer that model is approximating.

Theory: Mirceski, Komorsky-Lovric, Lovric, "Square-Wave Voltammetry: Theory and
Application", Springer 2007; Lovric & Komorsky-Lovric, J. Electroanal. Chem.
1988, 248, 239 (surface-confined methylene blue specifically).
"""

from __future__ import annotations

import numpy as np

F = 96485.332          # C / mol
R = 8.314462           # J / (mol K)
N_A = 6.02214076e23


def k0_from_distance(distance_ang, k0_contact=5.0e3, beta_per_ang=1.0,
                     d_contact_ang=5.0):
    """Standard heterogeneous rate constant vs tunnelling distance (1/s).

    beta = 1.0 /Angstrom is the value reported for alkanethiol self-assembled
    monolayers and used by Dauphin-Ducharme et al. (Langmuir 2017, 33, 4407) to
    convert their measured MB rate constants between tether lengths.

    k0_contact is a CALIBRATION parameter, not a prediction.  Anchor it to the
    measured values in that paper before quoting any absolute current:
        MB-C6 thiol SAM on gold        k0 = 1.5 +/- 0.5 cm/s
        MB on 20-nt ssDNA, 5'-hexanethiol   k0 = 0.4 +/- 0.1 cm/s
    (converted to 1/s using an 8.0 Angstrom monolayer thickness).  No equivalent
    measurement exists for MB on a CARBON electrode with a protein tether, which
    is a genuine open parameter of this project, not a modelling shortcut.
    """
    d = np.asarray(distance_ang, dtype=float)
    return k0_contact * np.exp(-beta_per_ang * (d - d_contact_ang))


def square_wave_potentials(e_start=0.0, e_end=-0.5, e_step=0.004,
                           e_sw=0.025):
    """Staircase + square-wave waveform.

    Returns (e_stair, e_forward, e_reverse) arrays, volts vs reference.
    Defaults are a typical MB scan window (0 to -0.5 V), 4 mV step,
    25 mV square-wave amplitude.
    """
    n = int(round(abs(e_end - e_start) / e_step)) + 1
    e_stair = np.linspace(e_start, e_end, n)
    sign = -1.0 if e_end < e_start else 1.0
    # forward pulse drives the scan direction, reverse pulse opposes it
    return e_stair, e_stair + sign * e_sw, e_stair - sign * e_sw


def _half_cycle(theta0, e, tau, k0, e0, alpha, temperature_k,
                sample_window=0.25, n_electrons=1):
    """Exact integration of one constant-potential half cycle.

    Returns (theta_end, mean_current_density_factor).

    A real potentiostat does not read an instantaneous current: it integrates
    over a sampling window at the END of each pulse.  Sampling the average
    rather than the instantaneous value matters here, because in the reversible
    limit the faradaic transient has fully decayed by the end of the pulse and
    the instantaneous current is identically zero.  The averaged (charge-based)
    form is both what the instrument measures and numerically well behaved:

        <d(theta)/dt> = [ theta(t2) - theta(t1) ] / (t2 - t1)

    with t2 = tau and t1 = (1 - sample_window) tau.
    """
    f_rt = n_electrons * F / (R * temperature_k)
    eta = e - e0
    k_c = k0 * np.exp(-alpha * f_rt * eta)
    k_a = k0 * np.exp((1.0 - alpha) * f_rt * eta)
    ksum = k_c + k_a
    theta_ss = k_c / ksum
    t1 = (1.0 - sample_window) * tau
    th1 = theta_ss + (theta0 - theta_ss) * np.exp(-ksum * t1)
    th2 = theta_ss + (theta0 - theta_ss) * np.exp(-ksum * tau)
    return th2, (th2 - th1) / (tau - t1)


def _refine_peak(e, i, j):
    """Sub-grid peak height and position by parabolic interpolation.

    Without this the reported peak height jumps between staircase grid points as
    parameters vary, which shows up as spurious non-monotonicity when sweeping
    the reporter distance.  The peak is a smooth function of the parameters; the
    grid is not.
    """
    if j == 0 or j == len(i) - 1:
        return float(i[j]), float(e[j])
    y0, y1, y2 = float(i[j - 1]), float(i[j]), float(i[j + 1])
    denom = y0 - 2.0 * y1 + y2
    if denom == 0.0:
        return y1, float(e[j])
    delta = 0.5 * (y0 - y2) / denom
    delta = max(-1.0, min(1.0, delta))
    step = float(e[1] - e[0])
    return float(y1 - 0.25 * (y0 - y2) * delta), float(e[j] + delta * step)


def swv_scan(k0, gamma_mol_cm2=2.0e-11, area_cm2=0.07, e0=-0.27,
             alpha=0.37, freq_hz=100.0, e_start=0.0, e_end=-0.5,
             e_step=0.004, e_sw=0.025, temperature_k=298.15,
             n_electrons=2, theta_init=None, sample_window=0.25):
    """Simulate one square-wave voltammogram of a surface-confined couple.

    Parameters
    ----------
    k0              : heterogeneous rate constant, 1/s (use k0_from_distance)
    gamma_mol_cm2   : MB surface coverage, mol/cm^2 (1e-11 - 1e-10 typical for
                      a packed protein monolayer with one label per protein)
    area_cm2        : electrode area
    e0              : MB formal potential vs Ag/AgCl (about -0.27 V)
    alpha           : transfer coefficient; 0.37 +/- 0.02 measured for MB
                      (Dauphin-Ducharme 2017), not the textbook 0.5
    freq_hz         : square-wave frequency; E-AB sensors are interrogated at
                      20-1000 Hz, with E_step 1-4 mV and E_sw 25-75 mV
    n_electrons     : MB/leucoMB is a 2-electron, 1-proton couple

    Returns
    -------
    dict with keys e (staircase potential), i_forward, i_reverse, i_net (A),
    peak_current (A), peak_potential (V).
    """
    e_stair, e_fwd, e_rev = square_wave_potentials(e_start, e_end, e_step, e_sw)
    tau = 1.0 / (2.0 * freq_hz)
    scale = n_electrons * F * area_cm2 * gamma_mol_cm2

    theta = 1.0 if theta_init is None else float(theta_init)
    if theta_init is None:
        # equilibrate at the initial potential (fully oxidised at E >> E0)
        f_rt = n_electrons * F / (R * temperature_k)
        eta0 = e_start - e0
        kc = np.exp(-alpha * f_rt * eta0)
        ka = np.exp((1.0 - alpha) * f_rt * eta0)
        theta = kc / (kc + ka)

    theta_start = theta
    i_f = np.empty_like(e_stair)
    i_r = np.empty_like(e_stair)
    for j in range(len(e_stair)):
        theta, d1 = _half_cycle(theta, e_fwd[j], tau, k0, e0, alpha,
                                temperature_k, sample_window, n_electrons)
        i_f[j] = scale * d1
        theta, d2 = _half_cycle(theta, e_rev[j], tau, k0, e0, alpha,
                                temperature_k, sample_window, n_electrons)
        i_r[j] = scale * d2

    i_net = i_f - i_r
    j = int(np.argmax(np.abs(i_net)))
    pk_i, pk_e = _refine_peak(e_stair, i_net, j)
    return {
        "theta_initial": float(theta_start),
        "theta_final": float(theta),
        "e": e_stair,
        "i_forward": i_f,
        "i_reverse": i_r,
        "i_net": i_net,
        "peak_current": pk_i,
        "peak_potential": pk_e,
        "k0": float(k0),
    }


def sensor_response(fraction_switched, d_closed_ang=6.0, d_open_ang=22.0,
                    gamma_total_mol_cm2=2.0e-11, **kw):
    """Net SWV peak for a monolayer that is a mixture of docked and released sensors.

    The two sub-populations are independent surface-confined couples with
    different k0, so their currents add.  Returns the peak of the summed
    voltammogram and the signal change relative to zero analyte.
    """
    f = float(fraction_switched)
    k0_closed = k0_from_distance(d_closed_ang)
    k0_open = k0_from_distance(d_open_ang)
    a = swv_scan(k0_closed, gamma_mol_cm2=gamma_total_mol_cm2 * (1.0 - f), **kw)
    b = swv_scan(k0_open, gamma_mol_cm2=gamma_total_mol_cm2 * f, **kw)
    i_net = a["i_net"] + b["i_net"]
    j = int(np.argmax(np.abs(i_net)))
    pk_i, pk_e = _refine_peak(a["e"], i_net, j)
    return {
        "e": a["e"],
        "i_net": i_net,
        "peak_current": pk_i,
        "peak_potential": pk_e,
        "peak_docked_population": float(a["peak_current"]),
        "peak_released_population": float(b["peak_current"]),
        "k0_closed": float(k0_closed),
        "k0_open": float(k0_open),
    }


def capacitive_background(freq_hz=100.0, area_cm2=0.07, c_dl_uf_cm2=20.0,
                          r_solution_ohm=100.0, e_sw=0.025, e_step=0.004,
                          sample_window=0.25):
    """Double-layer charging current surviving into the SWV sampling window (A).

    Square-wave voltammetry suppresses, but does not remove, the capacitive
    background: after each potential step the charging current decays with the
    cell time constant tau_RC = R_s C_dl A, and only what is left at the end of
    the pulse is measured.  Forward and reverse steps have opposite sign, so in
    the net (forward minus reverse) current their contributions ADD rather than
    cancel.

    This term is what sets the usable frequency window.  Without it the faradaic
    peak grows monotonically with frequency and a frequency optimisation is
    meaningless; with it there is a genuine optimum, which is the experimentally
    familiar behaviour of E-AB sensors.
    """
    c_dl = c_dl_uf_cm2 * 1e-6 * area_cm2          # farad
    tau_rc = r_solution_ohm * c_dl                # seconds
    tau = 1.0 / (2.0 * freq_hz)
    t1 = (1.0 - sample_window) * tau
    delta_e = 4.0 * e_sw + e_step                 # net step seen by i_forward - i_reverse
    charge = c_dl * delta_e * (np.exp(-t1 / tau_rc) - np.exp(-tau / tau_rc))
    return float(charge / (tau - t1))


def sensor_response_distributed(fraction_switched, d_closed_ang=6.0,
                                d_open_ang=14.0, sigma_ang=2.5,
                                gamma_total_mol_cm2=2.0e-12, n_quad=61, **kw):
    """Same as sensor_response, but each state samples a DISTRIBUTION of distances.

    A methylene blue tethered through a flexible linker to a protein does not sit
    at one distance from the electrode; it samples a conformational ensemble.
    Because k0 depends exponentially on distance, averaging over that ensemble is
    not the same as evaluating at the mean distance -- the fast (close) tail
    dominates the current, which blunts the predicted signal change considerably.

    This matters quantitatively.  A single-distance model predicts that
    sub-Angstrom displacements give tens of percent signal change, which is far
    more sensitive than the 4-30% that Kang/Plaxco (JACS 2017, 139, 12113)
    measure for real protein-based E-AB sensors across genuinely large
    conformational changes.  Distance dispersion is the leading candidate for
    that gap, and it is also what RedoxPySolid models as a k0 distribution.

    sigma_ang is the width of the distance distribution, and it is a fitted
    quantity, not a predicted one.

    n_quad matters more than it looks.  Individual quadrature nodes cross the
    quasi-reversible maximum as the mean distance is swept, so a coarse rule
    produces a spuriously non-monotonic response.  15 nodes is visibly
    unconverged; 41 and 81 agree to about 1%, so 61 is the default.
    """
    # Gauss-Hermite quadrature over a Gaussian distance distribution
    nodes, weights = np.polynomial.hermite_e.hermegauss(n_quad)
    weights = weights / weights.sum()
    f = float(fraction_switched)
    i_net = None
    e_axis = None
    for mean_d, pop in [(d_closed_ang, 1.0 - f), (d_open_ang, f)]:
        if pop <= 0.0:
            continue
        for node, w in zip(nodes, weights):
            d = mean_d + sigma_ang * node
            r = swv_scan(k0_from_distance(d),
                         gamma_mol_cm2=gamma_total_mol_cm2 * pop * w, **kw)
            i_net = r["i_net"] if i_net is None else i_net + r["i_net"]
            e_axis = r["e"]
    j = int(np.argmax(np.abs(i_net)))
    pk_i, pk_e = _refine_peak(e_axis, i_net, j)
    return {"e": e_axis, "i_net": i_net,
            "peak_current": pk_i, "peak_potential": pk_e}


def signal_gain(fraction_switched, **kw):
    """Relative peak-current change, the quantity an E-AB sensor actually reports.

    Returns NaN when the baseline peak is numerically zero, which happens when
    the docked reporter is fully reversible at the chosen frequency (no faradaic
    current survives into the sampling window).  That is a real operating regime,
    not a numerical artefact: it is why the interrogation frequency has to be
    tuned, and it is the reason the frequency sweep below is worth running.
    """
    i0 = sensor_response(0.0, **kw)["peak_current"]
    i = sensor_response(fraction_switched, **kw)["peak_current"]
    if abs(i0) < 1e-18:
        return float("nan")
    return (i - i0) / i0
