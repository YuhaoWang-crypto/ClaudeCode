"""
Fast self-checks for the plife package.  Run with `python3 test_plife.py`.

These are the invariants the whole analysis rests on; they take a few seconds
and catch any regression in the force law, the integrator or the theory.
"""

import sys

import numpy as np

from plife import kernels as K, matrices as M, model, observables as OB
from plife import stability as ST, twobody as TB

fails = []


def check(name, cond, detail=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"   {detail}" if detail else ""))
    if not cond:
        fails.append(name)


print("kernel")
r = np.linspace(1e-4, 60, 200001)
A, al, be, R = 0.7, 15.0, 45.0, 1.0
U, f = K.pair_potential(r, A, al, be, R), K.pair_force(r, A, al, be, R)
m = (r > 0.4) & (r < be - 0.4) & (abs(r - al) > 0.3) & (abs(r - (al + be) / 2) > 0.3)
check("U' == f", np.abs(np.gradient(U, r)[m] - f[m]).max() < 1e-6,
      f"max err {np.abs(np.gradient(U, r)[m] - f[m]).max():.2e}")
check("well depth = -A(beta-alpha)/2",
      abs(K.pair_potential(al, A, al, be, R) - (-A * (be - al) / 2)) < 1e-9)
check("f(alpha) == 0 (bond length)", abs(K.pair_force(al, A, al, be, R)) < 1e-12)
check("Uhat_core(0) == pi R alpha^3/12",
      abs(K.integrated_interaction(0.0, al, be, R) - np.pi * R * al**3 / 12) < 1e-6)

print("\nreciprocity algebra")
Arand = M.g_random(5, np.random.default_rng(0))
check("nu(symmetric) == 0", M.nonreciprocity(M.symmetrise(Arand)) < 1e-12)
check("nu(antisymmetric) == 1", M.nonreciprocity(M.antisymmetrise(Arand)) > 1 - 1e-12)
check("A == A_S + A_N",
      np.allclose(Arand, M.symmetrise(Arand) + M.antisymmetrise(Arand)))

print("\ntwo-body theory vs direct integration")
alm = np.array([[15.0, 15.0], [25.0, 15.0]])
bem = np.full((2, 2), 45.0)
p = model.Params(A=M.chase_pair(0.9, -0.4), alpha=alm, beta=bem, L=(9000.0, 9000.0))
pred = TB.dimer(p, 0, 1)
sim = model.ParticleLife(p, N=2, seed=1)
sim.types = np.array([0, 1])
sim.pos = np.array([[4500.0, 4500.0], [4520.0, 4500.0]])
sim.vel = np.zeros((2, 2))
sim.run(80.0)
X0, t0 = sim.pos.mean(0).copy(), sim.t
sim.run(40.0)
V = np.hypot(*(sim.pos.mean(0) - X0)) / (sim.t - t0)
d = np.hypot(*(sim.pos[1] - sim.pos[0]))
check("bond length r*", abs(d - pred["r_star"]) < 1e-3, f"{pred['r_star']:.4f} vs {d:.4f}")
check("drift speed V = mu_disc f_A(r*)", abs(V - abs(pred["V"])) / abs(pred["V"]) < 5e-3,
      f"{abs(pred['V']):.5f} vs {V:.5f}")

print("\nreciprocity theorem: injected power vanishes identically")
alm2, bem2 = M.random_radii(4, rng=np.random.default_rng(3), symmetric=True)
A_sym = M.symmetrise(M.g_random(4, np.random.default_rng(4)))
ps = model.Params(A=A_sym, alpha=alm2, beta=bem2, L=(500.0, 500.0))
s = model.ParticleLife(ps, N=600, seed=0)
s.run(3.0)
check("P_inj == 0 for reciprocal", abs(s.injected_power()) < 1e-9,
      f"{s.injected_power():.3e}")
pn = ps.copy(A=M.g_random(4, np.random.default_rng(4)))
s2 = model.ParticleLife(pn, N=600, seed=0)
s2.run(3.0)
check("P_inj != 0 for non-reciprocal", abs(s2.injected_power()) > 1e-6,
      f"{s2.injected_power():.3e}")

print("\nlinear stability")
rho = np.full(4, 600 / 4 / 500.0**2)
res = ST.scan(ps, rho, nk=200)
check("reciprocal matrix has no growing oscillatory mode", res.osc_growth <= 0,
      f"osc growth {res.osc_growth:+.4f}")
lt = ST.LinearTheory.from_params(ps, nk=200)
check("LinearTheory matches scan",
      abs(lt.analyse(ps.A, rho).growth - res.growth) < 1e-6)

print("\nobservables")
o = OB.observables(s)
check("observables complete", all(k in o for k in
      ("speed", "motility", "dim2", "psi6", "k_peak", "n_clusters")))
check("classify returns a label", isinstance(OB.classify(o), str), OB.classify(o))

print()
if fails:
    print(f"{len(fails)} FAILED: {fails}")
    sys.exit(1)
print("all checks passed")
