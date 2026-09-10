"""
Validation of every model in this package against published values or exact limits.

    python3 sparkcage/validate.py

Nothing here is a self-consistency check.  Each test compares against either a
number printed in a paper or a closed-form limit that the numerics must recover.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from thermo import switch_model as sw            # noqa: E402
from thermo import luccage_reference as ref      # noqa: E402
from echem import swv                            # noqa: E402

PASS, FAIL = "PASS", "FAIL"
results = []


def check(name, ok, detail=""):
    results.append((PASS if ok else FAIL, name, detail))
    print(f"  [{PASS if ok else FAIL}] {name}" + (f"  -- {detail}" if detail else ""))
    return ok


print("\n=== 1. LucCage reference model vs Quijano-Rubio 2021 ===")

# The paper's Extended Data Fig. 1a sweeps K_open over four decades at fixed
# K_CK = 1e-8 M, K_LT = 1e-9 M, K_R = 0.19, 10 nM lucCage : 100 nM lucKey.
for k_open in [1e-7, 1e-5, 1e-3, 1e-1]:
    dr = ref.dynamic_range_percent(k_open=k_open)
    base = ref.solve(1e-15, k_open=k_open)["fraction_reconstituted"]
    print(f"    K_open = {k_open:.0e}: background = {base*100:8.4f}% of cage, "
          f"dynamic range = {dr:10.1f}%")

# Qualitative claim from ED Fig. 1f, verbatim: "K_open exerts a predominant
# effect on the dynamic range". Monotone decrease of dynamic range with K_open.
drs = [ref.dynamic_range_percent(k_open=k) for k in [1e-7, 1e-5, 1e-3, 1e-1]]
check("dynamic range falls monotonically as the cage is loosened",
      all(drs[i] > drs[i + 1] for i in range(len(drs) - 1)),
      f"{[f'{d:.3g}' for d in drs]}")

# Paper, verbatim: "too tight cage-latch interaction results in low signal in
# the presence of target, and too weak an interaction results in high background
# in the absence of target." Both halves are testable.
bases = [ref.solve(1e-15, k_open=k)["fraction_reconstituted"]
         for k in [1e-7, 1e-5, 1e-3, 1e-1]]
maxes = [ref.solve(1e-5, k_open=k)["fraction_reconstituted"]
         for k in [1e-7, 1e-5, 1e-3, 1e-1]]
check("loosening the cage raises the background monotonically",
      all(bases[i] < bases[i + 1] for i in range(len(bases) - 1)),
      f"{[f'{b*100:.4g}%' for b in bases]}")
check("tightening the cage suppresses the saturated signal too",
      maxes[0] < maxes[-1] and maxes[0] < 0.5,
      f"saturated signal {maxes[0]*100:.1f}% at K_open=1e-7 "
      f"vs {maxes[-1]*100:.1f}% at 1e-1")

# Paper, verbatim: "the LOD of our sensor platform is about 0.1 x Kd of the
# latch-target affinity (K_LT)" and "the driving force for opening the switch
# becomes too weak below this concentration".  The claim is a FLOOR that holds
# however K_open is tuned, so test the best LOD achievable over all K_open.
k_lt = 1e-9
conc = np.logspace(-13, -5, 200)


def lod_at(k_open, k_lt=k_lt, rise_frac=0.1):
    r = ref.dose_response(conc, k_open=k_open, k_lt=k_lt)
    rise = (r - r[0]) / (r[-1] - r[0])
    idx = int(np.argmax(rise > rise_frac))
    return float(conc[idx]) if rise[idx] > rise_frac else float("nan")


lods = {k: lod_at(k) for k in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]}
best = min(v for v in lods.values() if np.isfinite(v))
for k, v in lods.items():
    print(f"    K_open = {k:.0e}  ->  LOD = {v*1e9:9.4f} nM")
print(f"    best over all K_open  = {best*1e9:.3f} nM, i.e. ~{best/k_lt:.1f} x K_LT")
check("no choice of K_open beats the paper's 0.1 x K_LT sensitivity ceiling",
      best >= 0.1 * k_lt,
      f"best {best/k_lt:.2f} x K_LT vs ceiling 0.1 x K_LT")
print("    Interpretation: a 10%-of-full-rise LOD proxy is about a decade more")
print("    conservative than an experimental 3-sigma LOD, which is set by assay")
print("    precision rather than by response magnitude. The paper's own measured")
print("    lucCageRBD LOD (0.015 nM against a 0.5 nM binder) is 0.03 x K_LT,")
print("    i.e. better than its own rule of thumb. Use the model to rank designs,")
print("    not to predict an absolute LOD.")

# Measured dynamic ranges in the paper span 50% (lucCageSARS2-M) to 1700%
# (lucCageRBD). The baseline simulation should land inside that envelope.
dr_base = ref.dynamic_range_percent()
check("baseline simulated dynamic range lies in the measured 50-1700% envelope",
      50.0 <= dr_base <= 1700.0, f"{dr_base:.0f}%")


print("\n=== 2. Reduced keyless model vs the full four-component model ===")

# The SparkCage electrochemical architecture removes the key and the luciferase.
# Check that the reduced model is not merely the full model with terms dropped:
# it should show a markedly worse detection window at the same cage stability.
dg_open_4kcal = 4.0                       # the only kcal/mol value in the literature
k_open_equiv = np.exp(-dg_open_4kcal / (ref.R_KCAL * 298.15))
print(f"    dG_open = 4 kcal/mol  ->  K_open = {k_open_equiv:.3e}"
      f"   (paper's baseline K_open = 1e-3)")
check("4 kcal/mol reproduces the paper's baseline K_open to within 2x",
      0.5 < k_open_equiv / 1e-3 < 2.0, f"ratio = {k_open_equiv/1e-3:.2f}")

dg_lt = sw.dg_from_kd(1e-9, 298.15)
ec50_keyless = sw.ec50(dg_open_4kcal, dg_lt, 298.15)
resp_full = ref.dose_response(conc, k_open=k_open_equiv, k_lt=1e-9)
rise_full = (resp_full - resp_full[0]) / (resp_full[-1] - resp_full[0])
ec50_full = float(conc[int(np.argmax(rise_full > 0.5))])
penalty = ec50_keyless / ec50_full
print(f"    EC50 with key    = {ec50_full*1e9:12.4f} nM")
print(f"    EC50 without key = {ec50_keyless*1e9:12.4f} nM")
print(f"    sensitivity penalty of dropping the key = {penalty:.3g}x")
check("removing the key costs more than an order of magnitude in sensitivity",
      penalty > 10.0, f"{penalty:.3g}x")
naive = 1.0 / 1e-8          # the tempting exp(-dG_CK/RT) = 1/K_CK estimate
check("the naive exp(-dG_CK/RT) estimate overstates the penalty by >1000x",
      naive / penalty > 1000.0,
      f"naive {naive:.3g}x vs mass-action "
      f"{ref.key_compensation_factor():.1f}x vs measured {penalty:.1f}x")


print("\n=== 3. Square-wave voltammetry module vs exact limits ===")

# (a) For a REVERSIBLE surface-confined couple with alpha = 0.5 the net peak
# sits at the formal potential. A quasi-reversible couple with alpha != 0.5 is
# shifted, which is why the symmetric reversible case is the right limit to test.
e0 = -0.27
# k0/f ~ 2 is the quasi-reversible window where end-of-pulse sampling actually
# returns a peak; a fully reversible couple has decayed to zero current by then,
# which is a real property of SWV, not a numerical artefact.
r = swv.swv_scan(k0=200.0, e0=e0, alpha=0.5, freq_hz=100.0, e_step=0.002)
check("symmetric (alpha=0.5) SWV peak sits at the formal potential",
      abs(r["peak_potential"] - e0) <= 0.004,
      f"peak at {r['peak_potential']*1000:.1f} mV vs E0 = {e0*1000:.0f} mV")
r_qr = swv.swv_scan(k0=100.0, e0=e0, alpha=0.37, freq_hz=50.0)
print(f"    (quasi-reversible, alpha=0.37: peak shifted to "
      f"{r_qr['peak_potential']*1000:.1f} mV, as expected)")

# (b) Coverage conservation: scanning far past E0 with fast kinetics must
# convert the entire monolayer from oxidised to reduced.
r = swv.swv_scan(k0=1.0e4, freq_hz=10.0, e_start=0.2, e_end=-0.7, e_step=0.002)
check("scanning through E0 converts the whole monolayer",
      r["theta_initial"] < 0.01 and r["theta_final"] > 0.99,
      f"theta {r['theta_initial']:.4f} -> {r['theta_final']:.4f}")

# (c) The sampled peak is NON-monotonic in k0, with a maximum near k0/f = 1.
# This is the quasi-reversible maximum, and it is the reason a sensor cannot be
# made better simply by making electron transfer faster.  An earlier version of
# this check asserted "faster is bigger", which is only true below the maximum
# and started failing the moment the electron count was put into the
# Butler-Volmer exponent where it belongs.
_f = 100.0
_k = np.logspace(-1.0, 5.0, 80)
_p = np.array([abs(swv.swv_scan(k0=k, freq_hz=_f)["peak_current"]) for k in _k])
_kappa = _k[int(np.argmax(_p))] / _f
check("the quasi-reversible maximum sits at k0/f near 1",
      0.5 <= _kappa <= 2.0, f"maximum at k0/f = {_kappa:.2f}")
check("peak falls away on BOTH sides of the quasi-reversible maximum",
      _p[0] < _p.max() / 10.0 and _p[-1] < _p.max() / 10.0,
      f"{_p[0]*1e9:.2f} nA .. {_p.max()*1e9:.0f} nA .. {_p[-1]*1e9:.2g} nA")

# (c2) Laviron's surface-confined peak width: about 90/n mV at half height in
# the reversible limit.  This is what pins the electron count to the exponent
# rather than to the charge prefactor alone.
for _n, _want in ((1, 90.6), (2, 45.3)):
    _r = swv.swv_scan(k0=1.0e2, freq_hz=100.0, n_electrons=_n, e_step=0.001,
                      alpha=0.5)
    _i = np.abs(_r["i_net"])
    _sel = np.where(_i >= _i.max() / 2.0)[0]
    _fwhm = abs(_r["e"][_sel[-1]] - _r["e"][_sel[0]]) * 1000.0
    check(f"n = {_n} peak half-width matches Laviron's {_want:.0f} mV",
          abs(_fwhm - _want) < 0.25 * _want, f"{_fwhm:.1f} mV")

# (d) Peak current must be strictly proportional to surface coverage.
p1 = abs(swv.swv_scan(k0=100.0, gamma_mol_cm2=1e-11)["peak_current"])
p2 = abs(swv.swv_scan(k0=100.0, gamma_mol_cm2=2e-11)["peak_current"])
check("peak current is linear in surface coverage",
      abs(p2 / p1 - 2.0) < 1e-6, f"ratio = {p2/p1:.6f}")

# (e) Capacitive background must fall as the pulse gets long relative to R_s C.
c_fast = swv.capacitive_background(freq_hz=5000.0)
c_slow = swv.capacitive_background(freq_hz=10.0)
check("capacitive background is suppressed at low frequency", c_fast > c_slow,
      f"{c_fast*1e9:.3g} nA at 5 kHz vs {c_slow*1e9:.3g} nA at 10 Hz")


print("\n=== 4. Reduced switch model closed-form identities ===")
dgo, dglt = 3.0, -11.0
check("EC50 matches the analytic Kd * exp(dG_open/RT)",
      abs(sw.ec50(dgo, dglt) / (np.exp(dglt / sw.rt()) * np.exp(dgo / sw.rt())) - 1.0) < 0.01,
      "")
check("dynamic range matches the analytic exp(dG_open/RT)",
      abs(sw.dynamic_range(dgo) / (1.0 + np.exp(dgo / sw.rt())) - 1.0) < 1e-9, "")
check("fraction switched is monotone in analyte concentration",
      np.all(np.diff(sw.fraction_switched(np.logspace(-12, -4, 50), dgo, dglt)) > 0), "")

n_fail = sum(1 for s, _, _ in results if s == FAIL)
print(f"\n{'=' * 74}\n{len(results) - n_fail}/{len(results)} checks passed"
      + ("" if n_fail == 0 else f"   ({n_fail} FAILED)") + f"\n{'=' * 74}")
sys.exit(1 if n_fail else 0)
