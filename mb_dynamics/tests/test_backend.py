"""Numerical-correctness tests for the SWV backend (run: python -m pytest mb_dynamics/tests -q)."""
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import swv_backend as sb  # noqa: E402


def random_reversible(M, rng, scale=1e4):
    pi = rng.dirichlet(np.ones(M))
    F = rng.uniform(0, 1, (M, M))
    F = np.triu(F, 1)
    F = F + F.T                         # symmetric equilibrium fluxes pi_j L_ij
    L = F / pi[None, :] * scale          # L_ij = F_ij / pi_j (rate j -> i)
    L -= np.diag(L.sum(0))
    return pi, L


def test_half_pulse_matches_full_expm():
    rng = np.random.default_rng(1)
    pi, L = random_reversible(6, rng)
    w = np.exp(-5 * rng.uniform(0, 2, 6))
    kin = sb.Kinetics(pi, w, L)
    wf = sb.Waveform(100.0, 0.05, e_start=-0.20, e_end=-0.20)   # one step
    k_c, a, e0 = 3e3, 0.37, -0.26
    out = sb.simulate(kin, [wf], k_c, a, e0)[0]
    fn = sb.N_ELEC * sb.F_CONST / (sb.R_GAS * sb.T_DEFAULT)
    x = kin.pi / (1 + np.exp(-fn * (wf.e_start - e0)))
    cur = []
    for mv in (wf.mv_fwd[0], wf.mv_bwd[0]):
        x = sb.full_system_expm_check(kin, k_c, a, e0, mv / 1000, wf.h, x)
        eta = mv / 1000 - e0
        cur.append(np.sum(k_c * w * (np.exp(-a * fn * eta) * x - np.exp((1 - a) * fn * eta) * (kin.pi - x))))
    # reference expm on the stiff full system cancels near zero current: absolute tolerance
    assert np.allclose([out["fwd"][0], out["bwd"][0]], cur, rtol=1e-7, atol=1e-8 * np.abs(cur).max())


def test_full_scan_matches_stiff_ode():
    rng = np.random.default_rng(2)
    pi, L = random_reversible(4, rng, scale=300.0)
    w = np.array([1.0, 0.3, 0.05, 0.01])
    kin = sb.Kinetics(pi, w, L)
    k_c, a, e0, f, esw = 2e3, 0.4, -0.26, 100.0, 0.05
    wf = sb.Waveform(f, esw, e_start=-0.15, e_end=-0.35)
    out = sb.simulate(kin, [wf], k_c, a, e0)[0]
    fn = sb.N_ELEC * sb.F_CONST / (sb.R_GAS * sb.T_DEFAULT)
    p = np.r_[kin.pi / (1 + np.exp(-fn * (wf.e_start - e0))), 0 * kin.pi]
    p[kin.M:] = kin.pi - p[:kin.M]
    i_f, i_b = [], []
    kl = k_c * w
    for k in range(len(wf.e_stair)):
        for half, mv in ((0, wf.mv_fwd[k]), (1, wf.mv_bwd[k])):
            eta = mv / 1000 - e0
            Kr, Ko = kl * np.exp(-a * fn * eta), kl * np.exp((1 - a) * fn * eta)

            def rhs(t, q):
                pO, pR = q[:4], q[4:]
                return np.r_[L @ pO - Kr * pO + Ko * pR, L @ pR + Kr * pO - Ko * pR]
            sol = solve_ivp(rhs, (0, wf.h), p, method="Radau", rtol=1e-10, atol=1e-13)
            p = sol.y[:, -1]
            cur = np.sum(Kr * p[:4] - Ko * p[4:])
            (i_f if half == 0 else i_b).append(cur)
    assert np.allclose(out["fwd"], i_f, rtol=1e-6, atol=1e-9 * np.abs(i_f).max())
    assert np.allclose(out["bwd"], i_b, rtol=1e-6, atol=1e-9 * np.abs(i_b).max())


def test_fast_exchange_limit_equals_R0_and_slow_limit_equals_static():
    rng = np.random.default_rng(3)
    z = 0.35 + rng.gamma(3.0, 0.6, 20000)
    kin1 = sb.smoluchowski_kinetics(z, 40, D_nm2_per_s=7e11, beta_nm=5.0, z_min=0.35)
    kin0 = sb.fast_eq_kinetics(kin1.pi, kin1.w)
    wf = [sb.Waveform(100.0, 0.05)]
    k_c = 1e5
    r0 = sb.simulate(kin0, wf, k_c, 0.37, -0.26)[0]["net"]
    rfast = sb.simulate(kin1.scaled(1e3), wf, k_c, 0.37, -0.26)[0]["net"]
    assert np.abs(rfast - r0).max() < 1e-3 * np.abs(r0).max()
    static = sb.Kinetics(kin1.pi, kin1.w, None)
    rslow = sb.simulate(kin1.scaled(1e-20), wf, k_c, 0.37, -0.26)[0]["net"]
    rstat = sb.simulate(static, wf, k_c, 0.37, -0.26)[0]["net"]
    assert np.abs(rslow - rstat).max() < 1e-6 * np.abs(rstat).max()


def test_reversible_and_irreversible_limits_vanish():
    kin = sb.fast_eq_kinetics([1.0], [1.0])
    wf = [sb.Waveform(100.0, 0.025)]
    peak = lambda k: sb.simulate(kin, wf, k, 0.5, -0.26)[0]["net"].max() / (100.0)
    ks = np.logspace(-1, 6, 29)
    psi = np.array([peak(k) for k in ks])
    imax = psi.argmax()
    # quasi-reversible maximum: interior, near k/f ~ 1 for n*Esw = 50 mV
    assert 0 < imax < len(ks) - 1
    assert 0.2 < ks[imax] / 100.0 < 5.0
    assert psi[-1] < 0.05 * psi[imax] and psi[0] < 0.1 * psi[imax]


def test_smoluchowski_generator_properties():
    rng = np.random.default_rng(4)
    z = 0.35 + rng.gamma(2.0, 1.0, 50000)
    for n in (40, 80):
        kin = sb.smoluchowski_kinetics(z, n, 7e11, 5.0, 0.35)
        sq = np.sqrt(kin.pi)
        L = (kin.S * sq[:, None]) / sq[None, :]
        assert np.allclose(L.sum(0), 0, atol=1e-6 * np.abs(L).max())
        off = L - np.diag(np.diag(L))
        assert (off >= 0).all()
        assert np.allclose(L @ kin.pi, 0, atol=1e-8 * np.abs(L).max())
        # detailed balance
        J = L * kin.pi[None, :]
        assert np.allclose(J, J.T, rtol=1e-10, atol=1e-12 * np.abs(J).max())


def test_long_hold_reaches_nernst():
    rng = np.random.default_rng(5)
    pi, L = random_reversible(5, rng, scale=1e3)
    kin = sb.Kinetics(pi, rng.uniform(0.01, 1, 5), L)
    wf = sb.Waveform(0.01, 0.0, e_start=-0.25, e_end=-0.25)   # 50 s hold
    out = sb.simulate(kin, [wf], 1e3, 0.4, -0.26)[0]
    assert abs(out["fwd"][0]) < 1e-10
