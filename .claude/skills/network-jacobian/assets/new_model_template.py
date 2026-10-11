"""
Template: onboard a new kinetic model into netjac and run the three
diagnostics on it. Copy to netjac/demos/dX_<name>.py and fill in the TODOs.
"""
import numpy as np

from netjac.core import KineticModel
from netjac import spectrum as sp, relax as rx

MODEL = "BIOMD0000000064"        # TODO: BIOMD id or path to an SBML file


def onboard():
    m = KineticModel(MODEL)
    print("state:", m.state_ids, "| conservation laws:", m.n_cons)
    print("enzyme knobs:", m.enzyme_params())
    try:
        x = m.steady_state()
    except RuntimeError:                 # oscillating: unstable fixed point
        m.rr.simulate(0, 50, 10)
        x = m.fixed_point(m.x())
    cl = sp.classify(m.reduced_jacobian(x))
    print("state class:", cl["kind"], "λ_lead =", cl["lambda_lead"])
    if cl["lambda_lead"].real < 0:
        c = m.control(x)
        print("MCA summation / connectivity residuals:",
              c["summation_residual"], c["connectivity_residual"])  # ✅ < 1e-6
    return m, x, cl


def three_diagnostics(m, x, cl):
    # 1) edges: elasticities (rows = reactions, cols = species)
    eps = m.elasticities(x)
    # 2) relaxation: pulse the first species, observe two species
    t, Y = rx.pulse_response(m, x, m.state_ids[0], 0.1 * x[0],
                             t_end=10 * cl["tau_relax"], n=500,
                             observe=m.state_ids[:2])
    lam_dmd = rx.hankel_dmd(Y, t[1] - t[0], rank=4)
    # 3) fragility: total, per enzyme (choose solver to match stability)
    solver = "integrate" if cl["lambda_lead"].real < 0 else "fixed_point"
    F = {r: sp.total_fragility(m, r, x, solver=solver)
         for r in m.enzyme_params()}
    _, _, part = sp.mode_shape(m, x)
    return eps, lam_dmd, F, dict(zip(m.state_ids, part))


if __name__ == "__main__":
    m, x, cl = onboard()
    eps, lam, F, part = three_diagnostics(m, x, cl)
    print("DMD eigenvalues:", np.round(lam[:3], 3))
    print("top fragile:", sorted(F.items(), key=lambda kv: -abs(kv[1]))[:5])
    print("best probes:", sorted(part.items(), key=lambda kv: -kv[1])[:5])
