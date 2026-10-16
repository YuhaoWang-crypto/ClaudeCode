"""Multi-species (stepwise) chemistry backend tests (run from scripts/: python -m pytest tests -q)."""
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import chem_backend as cb  # noqa: E402
import swv_backend as sb  # noqa: E402
from test_backends import random_reversible  # noqa: E402

P1 = dict(k1=2e3, k2=5e2, alpha=0.4, Em=-0.26, dE=-0.05)
P2 = dict(k1=2e3, k2=5e2, alpha=0.4, Em=-0.26, dE=0.04, kp=3e3)


def _kin(M, seed, scale=300.0):
    rng = np.random.default_rng(seed)
    pi, L = random_reversible(M, rng, scale=scale)
    return sb.Kinetics(pi, np.exp(-rng.uniform(0, 3, M)), L)


def test_C0_identical_to_round1_backend():
    rng = np.random.default_rng(0)
    z = 0.35 + rng.gamma(3.0, 0.6, 20000)
    kin = sb.smoluchowski_kinetics(z, 12, 7e8, 5.0, 0.35)
    wfs = [sb.Waveform(100.0, 0.05), sb.Waveform(500.0, 0.025)]
    for k in (sb.fast_eq_kinetics(kin.pi, kin.w), kin, kin.scaled(1e-3)):
        a = sb.simulate(k, wfs, 3e3, 0.37, -0.26)
        b = cb.simulate(k, cb.Chemistry("C0"), dict(k=3e3, alpha=0.37, E0=-0.26), wfs)
        for x, y in zip(a, b):
            for c in ("fwd", "bwd", "net"):
                assert np.allclose(x[c], y[c], rtol=1e-8, atol=1e-10 * np.abs(x[c]).max())


def _reference_scan(kin, chem, p, wf):
    """Stiff Radau integration of the dense joint generator over a whole scan."""
    pvec = cb.equilibrium(kin, chem, p, wf.e_start)
    i_f, i_b = [], []
    for k in range(len(wf.e_stair)):
        for half, mv in ((0, wf.mv_fwd[k]), (1, wf.mv_bwd[k])):
            Gm = cb.full_generator(kin, chem, p, mv / 1000)
            sol = solve_ivp(lambda t, y: Gm @ y, (0, wf.h), pvec, method="Radau", jac=Gm, rtol=1e-10, atol=1e-14)
            pvec = sol.y[:, -1]
            (i_f if half == 0 else i_b).append(cb.current_from_p(kin, chem, p, mv / 1000, pvec))
    return np.array(i_f), np.array(i_b), pvec


def test_C1_C2_match_stiff_ode_and_expm():
    for chem, p in ((cb.Chemistry("C1"), P1), (cb.Chemistry("C2"), P2),
                    (cb.Chemistry("E1"), dict(k=3e3, alpha=0.5, E0=-0.26))):
        kin = _kin(3, 11)
        wf = sb.Waveform(100.0, 0.05, e_start=-0.17, e_end=-0.33)
        out = cb.simulate(kin, chem, p, [wf])[0]
        rf, rb, pend = _reference_scan(kin, chem, p, wf)
        tol = 1e-6 * np.abs(rf).max()
        assert np.allclose(out["fwd"], rf, rtol=1e-6, atol=tol)
        assert np.allclose(out["bwd"], rb, rtol=1e-6, atol=tol)
        # per-conformation population conserved by the dense reference
        assert np.allclose(pend.reshape(kin.M, chem.S).sum(1), kin.pi, atol=1e-10)
        # one half-pulse by expm from equilibrium at another potential
        p0 = cb.equilibrium(kin, chem, p, -0.20)
        Gm = cb.full_generator(kin, chem, p, -0.30)
        p1 = expm(Gm * 0.005) @ p0
        w1 = sb.Waveform(100.0, 0.10, e_start=-0.20, e_end=-0.20)   # forward half at -0.30 V
        o1 = cb.simulate(kin, chem, p, [w1])[0]
        assert np.isclose(o1["fwd"][0], cb.current_from_p(kin, chem, p, -0.30, p1), rtol=1e-7)


def test_exchange_limits_C1():
    rng = np.random.default_rng(3)
    z = 0.35 + rng.gamma(3.0, 0.6, 20000)
    kin1 = sb.smoluchowski_kinetics(z, 20, 7e8, 5.0, 0.35)
    kin0 = sb.fast_eq_kinetics(kin1.pi, kin1.w)
    chem = cb.Chemistry("C1")
    p = dict(k1=1e4, k2=3e3, alpha=0.4, Em=-0.26, dE=-0.05)
    wf = [sb.Waveform(100.0, 0.05)]
    r0 = cb.simulate(kin0, chem, p, wf)[0]["net"]
    rf = cb.simulate(kin1.scaled(1e3), chem, p, wf)[0]["net"]
    assert np.abs(rf - r0).max() < 1e-4 * np.abs(r0).max()
    st = cb.simulate(sb.Kinetics(kin1.pi, kin1.w, None), chem, p, wf)[0]["net"]
    rs = cb.simulate(kin1.scaled(1e-20), chem, p, wf)[0]["net"]
    assert np.abs(rs - st).max() < 1e-6 * np.abs(st).max()


def test_reversible_limit_and_large_dE_C1():
    chem = cb.Chemistry("C1")
    kin = sb.fast_eq_kinetics([1.0], [1.0])
    wf = [sb.Waveform(100.0, 0.025)]
    fast = cb.simulate(kin, chem, dict(k1=1e9, k2=1e9, alpha=0.5, Em=-0.26, dE=0.0), wf)[0]["net"]
    mid = cb.simulate(kin, chem, dict(k1=1e2, k2=1e2, alpha=0.5, Em=-0.26, dE=0.0), wf)[0]["net"]
    assert np.abs(fast).max() < 1e-3 * np.abs(mid).max()
    # widely separated potentials -> two separate one-electron waves near E1 and E2
    sep = cb.simulate(kin, chem, dict(k1=1e2, k2=1e2, alpha=0.5, Em=-0.26, dE=0.20), wf)[0]
    E, n = sep["E"], sep["net"]
    peaks = [E[i] for i in range(1, len(n) - 1) if n[i] > n[i - 1] and n[i] >= n[i + 1] and n[i] > 0.3 * n.max()]
    assert len(peaks) == 2 and abs(peaks[0] + 0.16) < 0.02 and abs(peaks[1] + 0.36) < 0.02, peaks
