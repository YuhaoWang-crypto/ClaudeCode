"""
Module template — copy to grn_pipeline/mXX_<name>.py and fill in.

Pick ONE dynamical class below, delete the others, and reuse the matching
generic engine (m19 for switches, m21 for oscillators, m22 for SNIC). Do NOT
re-implement Lyapunov / early-warning statistics — import the validated ones.

State honestly in this docstring:
  ✅ rigorous  = deterministic computation / decimal-exact reproduction / real stat
  ⚠️ hypothesis = illustrative parameterization or an untested claim
"""
import numpy as np


# ── CLASS A: bistable switch (saddle-node) ──────────────────────────────────
# Reuse m19_switch_library. Write a 1-D rate law with a control parameter p.
def f_switch(x, p):
    # example: positive feedback via a Hill term minus linear removal + drive p
    def hill(x, K, n): return x**n / (K**n + x**n)
    return 1.5 * hill(x, K=1.0, n=4) - x + p


def report_switch(fig_path="figures/mXX_myswitch.png"):
    from grn_pipeline import m19_switch_library as sw
    prange = (0.0, 0.9)
    window = sw.bistable_window(f_switch, prange, xmax=3.0)   # fold edges
    model = {"id": "my switch", "f": f_switch, "prange": prange, "xmax": 3.0}
    # sw.titration(model, ...) → Langevin traces + SD/AR1 across the window
    # (m19 model dict keys: "id", "f", "xmax", "prange")
    print("bistable window:", window)
    return {"window": window}


# ── CLASS B: oscillator (Hopf) ──────────────────────────────────────────────
# Reuse m21_oscillators. Write an ODE system f(state, p); control parameter p.
def f_osc(state, p):
    x, y = state
    # replace with your negative-feedback oscillator; p is the bifurcation knob
    return np.array([p - (1 + 1.0) * x + x * x * y,   # Brusselator-like skeleton
                     1.0 * x - x * x * y])


def report_osc(fig_path="figures/mXX_myosc.png"):
    from grn_pipeline import m21_oscillators as osc
    # m21 model dict keys: "f", "x0", "prange"; power_spectrum also needs "var"
    # (index of the state variable whose spectrum to take)
    model = {"f": f_osc, "x0": np.array([1.0, 1.0]), "prange": (1.5, 2.5),
             "var": 0}
    hopf = osc.find_hopf(model)                 # where Re(λ)→0, Im(λ)≠0
    # osc.power_spectrum(model, p, ...) → spectral peak at the intrinsic freq
    print("Hopf at:", hopf)
    return {"hopf": hopf}


# ── CLASS C: exact literature model (decimal-exact) ─────────────────────────
# Reuse m20b. Needs libroadrunner + network (guard in run_all.py).
def report_exact(biomd_id="BIOMD0000000027"):
    from grn_pipeline import m20b_biomodels_exact as ex
    rr = ex._rr(ex.fetch_sbml(biomd_id))
    lam = float(np.linalg.eigvals(np.array(rr.getFullJacobian())).real.max())
    print(f"{biomd_id}: λ_max = {lam:+.5f}")
    return {"biomd": biomd_id, "lambda_max": lam}


def report(fig_path="figures/mXX_mypathway.png"):
    """Entry point aggregated by run_all.py. Delete unused classes above."""
    print("=" * 70)
    print("MODULE XX — <pathway name> (<bifurcation class>)")
    print("=" * 70)
    r = report_switch(fig_path)      # or report_osc / report_exact
    return r


if __name__ == "__main__":
    report()
