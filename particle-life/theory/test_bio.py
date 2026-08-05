"""
Self-checks for plife.bio.  Run with `python3 test_bio.py`.

These are limit tests, not regression tests: each one pins the implementation to
a result that is known independently of the code (law of mass action, the
valence gate, the ideal-gas limit, Newton's third law).
"""

import sys

import numpy as np

from plife.bio import inverse as I, phases as P, thermo as T, units as U

fails = []


def check(name, cond, detail=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"   {detail}" if detail else ""))
    if not cond:
        fails.append(name)


print("units")
check("1 M = 0.60221 nm^-3", abs(U.M_TO_NM3 - 0.602214076) < 1e-9)
check("R(50 kDa) ~ 2.44 nm", abs(float(U.radius_from_mw(50)) - 2.437) < 0.01)
check("uM <-> density round trip",
      abs(float(U.density_to_conc(U.conc_to_density(12.3))) - 12.3) < 1e-9)
phi = float(U.volume_fraction(U.conc_to_density(10), 2.5))
check("10 uM of a 2.5 nm protein is dilute", phi < 1e-3, f"phi = {phi:.2e}")

print("\nthird law is enforced")
try:
    T.Mixture(names=["A", "B"], diameters=[5, 5], valence=[2, 2],
              kd=np.array([[1e-6, 1e-6], [1e-5, 1e-6]]))
    check("asymmetric Kd rejected", False)
except ValueError as e:
    check("asymmetric Kd rejected", "third law" in str(e))

print("\nlaw of mass action, self-association  (Kd = [A]^2/[A2])")
worst = 0.0
for kd in (1e-7, 1e-6, 1e-5):
    mix = T.Mixture(names=["A"], diameters=[5.0], valence=[1.0], kd=np.array([[kd]]))
    for c in (0.1, 1.0, 10.0):
        rho = U.conc_to_density(np.array([[c]]))
        X = T.unbonded_fractions(rho, mix)[0, 0]
        mono, dim = rho[0, 0] * X, rho[0, 0] * (1 - X) / 2
        worst = max(worst, abs(dim / (mono**2 / (kd * U.M_TO_NM3)) - 1))
check("dimer count matches mass action", worst < 5e-3, f"max rel err {worst:.2e}")

print("\nlaw of mass action, hetero-association  (Kd = [A][B]/[AB])")
worst = 0.0
kd = 1e-6
mix = T.Mixture(names=["A", "B"], diameters=[5.0, 5.0], valence=[1.0, 1.0],
                kd=np.array([[np.inf, kd], [kd, np.inf]]))
for cA, cB in ((1.0, 1.0), (10.0, 1.0), (5.0, 20.0)):
    rho = U.conc_to_density(np.array([[cA, cB]]))
    X = T.unbonded_fractions(rho, mix)[0]
    mA, mB = rho[0, 0] * X[0], rho[0, 1] * X[1]
    bonds = rho[0, 0] * (1 - X[0])
    worst = max(worst, abs(bonds / (mA * mB / (kd * U.M_TO_NM3)) - 1))
check("complex count matches mass action", worst < 5e-3, f"max rel err {worst:.2e}")

print("\nvalence gate")
mono_never = True
for kd in np.geomspace(1e-9, 1e-3, 9):
    m = T.Mixture(names=["A"], diameters=[5.0], valence=[1.0], kd=np.array([[kd]]))
    if np.isfinite(T.csat_scan(m, 0, lo_uM=1e-3, hi_uM=1e4, n=140)["csat_uM"]):
        mono_never = False
check("monovalent never phase separates (9 Kd from 1 nM to 1 mM)", mono_never)
m4 = T.Mixture(names=["A"], diameters=[5.0], valence=[4.0], kd=np.array([[3e-4]]))
cs4 = T.csat_scan(m4, 0, lo_uM=1e-3, hi_uM=1e4, n=200)["csat_uM"]
check("tetravalent does phase separate", np.isfinite(cs4), f"c_sat = {cs4:.3f} uM")

print("\nstability of the dilute limit")
mix = T.Mixture(names=["A", "B"], diameters=[5.2, 4.3], valence=[5.0, 3.0],
                kd=np.array([[3e-4, 3e-5], [3e-5, np.inf]]))
rho_dil = U.conc_to_density(np.array([[1e-4, 1e-4]]))
check("very dilute mixture is stable", float(T.spinodal_eig(rho_dil, mix)[0]) > 0,
      f"min eig = {float(T.spinodal_eig(rho_dil, mix)[0]):.3e}")

print("\nphase behaviour")
# The hull is built over the sampled columns, so c_sat is only resolved to the
# grid: a 3-point partner axis is not enough to see the shallow minimum.
cB_probe = np.concatenate([[0.0], np.geomspace(0.1, 3000, 39)])
cs = I.csat_curve(mix, cB_probe)[0]
check("c_sat is finite without partner", np.isfinite(cs[0]), f"{cs[0]:.3f} uM")
jm = int(np.nanargmin(cs))
check("a little partner lowers c_sat (cross-linking)", cs[jm] < cs[0],
      f"min {cs[jm]:.3f} at [B]={cB_probe[jm]:.2f} uM  <  {cs[0]:.3f}")
check("excess partner raises it again (re-entrance)", cs[-1] > 10 * cs[0],
      f"{cs[-1]:.1f} >> {cs[0]:.3f}")

pd = P.phase_diagram(mix, np.geomspace(0.1, 300, 24), np.geomspace(0.1, 300, 24))
check("phase diagram has homogeneous and condensed regions",
      (pd["label"] == 0).any() and (pd["label"] > 0).any())
w = P.partner_window(pd, 10.0)
check("partner_window returns an interval", len(w) >= 1 and w[0][1] > w[0][0],
      f"{w}")

print()
if fails:
    print(f"{len(fails)} FAILED: {fails}")
    sys.exit(1)
print("all checks passed")
