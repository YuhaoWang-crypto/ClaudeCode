"""
spectrum — what the eigenvalues of the dynamical Jacobian Jz do and do NOT say,
plus eigenvalue sensitivities ("fragility") with finite-difference checks.

Rigour notes (keep these in every write-up):
  * Eigenvalues of Jz at a fixed point classify LOCAL stability of that fixed
    point (Hartman-Grobman, hyperbolic case). They do not by themselves certify
    a limit cycle (needs Floquet multipliers) or chaos (needs a positive
    Lyapunov exponent along the attractor). classify() says which case applies.
  * The "dFlux/dEnzyme" control matrix is NOT the stability Jacobian; it is
    built FROM its inverse (core.control). Its eigenvalues have no stability
    meaning; use Jz.

Fragility (first-order perturbation theory, simple eigenvalue lambda with right
eigenvector u and left eigenvector w of Jz):
    d lambda / d Jz_ab        = conj(w_a) u_b / (w^H u)
  An EDGE is an elasticity eps[r,k] (reaction r senses species k). Since
  Jz = B^T N eps B:
    d lambda / d eps[r,k]     = (w^H B^T N)_r (B u)_k / (w^H u)
  scaled edge fragility       F_rk = eps[r,k] * d lambda / d eps[r,k]
                              (= change of lambda per unit RELATIVE change of
                                 that interaction's strength)
  enzyme fragility            F_r  = sum_k F_rk = d lambda / d ln E_r
"""
import numpy as np
from scipy.linalg import eig


def leading(Jz, target=None):
    """Eigen-triple of the leading eigenvalue, or of the one closest to
    `target` (e.g. to follow the critical complex pair of a Hopf even when a
    non-critical real mode is slower)."""
    vals, wl, vr = eig(Jz, left=True, right=True)
    if target is None:
        i = int(np.argmax(vals.real))
    else:
        i = int(np.argmin(np.abs(vals - target)))
    return vals, i, wl[:, i], vr[:, i]


def critical_pair(Jz):
    """Complex eigenvalue (Im>0) with the largest real part, or None."""
    vals = np.linalg.eigvals(Jz)
    cp = vals[vals.imag > 1e-9 * (np.abs(vals).max() + 1)]
    return cp[np.argmax(cp.real)] if len(cp) else None


def classify(Jz, tol=1e-9):
    vals = np.linalg.eigvals(Jz)
    lead = vals[np.argmax(vals.real)]
    re, im = lead.real, abs(lead.imag)
    # median, not max: stiff models (|λ|max/|λ|min ~ 1e10) would otherwise
    # mislabel a clearly stable slow mode as "critical"
    scale = np.median(np.abs(vals))
    if abs(re) < tol * scale:
        kind = "critical (non-hyperbolic: linearisation inconclusive)"
    elif re < 0:
        kind = "stable focus (damped oscillatory return)" if im > 0 \
            else "stable node (monotone return)"
    else:
        kind = ("unstable focus -> oscillation (confirm limit cycle by "
                "simulation / Floquet)" if im > 0
                else "unstable saddle/node -> switch to another state")
    return {"kind": kind, "lambda_lead": lead,
            "tau_relax": (1.0 / abs(re)) if re != 0 else np.inf,
            "period": (2 * np.pi / im) if im > 0 else None,
            "damping_ratio": (-re / abs(lead)) if abs(lead) > 0 else None,
            "eigenvalues": vals}


def edge_fragility(model, x, eps=None, target=None):
    """Return (F, dlam_deps, lam): F[r,k] scaled fragility of edge r<-k on the
    real part of the leading eigenvalue (or the eigenvalue nearest `target`)."""
    eps = model.elasticities(x) if eps is None else eps
    Jz = model.B.T @ model.N @ eps @ model.B
    vals, i, w, u = leading(Jz, target)
    denom = np.conj(w) @ u
    left = (np.conj(w) @ model.B.T @ model.N)          # over reactions
    right = model.B @ u                                 # over species
    dlam = np.outer(left, right) / denom
    F = (eps * dlam).real
    return F, dlam, vals[i]


def enzyme_fragility(model, x, eps=None, target=None):
    F, _, lam = edge_fragility(model, x, eps, target)
    return F.sum(1), lam


def mode_shape(model, x, target=None):
    """Right eigenvector of the chosen eigenvalue in species coordinates,
    |u_k| / x_k = relative participation of species k (probe placement)."""
    Jz = model.reduced_jacobian(x)
    vals, i, w, u = leading(Jz, target)
    ux = model.B @ u
    rel = np.abs(ux) / np.maximum(np.abs(x), 1e-12)
    return vals[i], ux, rel / rel.max()


def enzyme_fragility_fd(model, x, rel=1e-4):
    """Finite-difference check of d Re(lambda_lead)/d ln E_r: scale enzyme r by
    exp(+-rel) at FIXED state (local effect on the Jacobian at x)."""
    out = {}
    ez = model.enzyme_params()
    for rid in ez:
        lams = []
        for s in (+1, -1):
            model.scale_enzyme(rid, np.exp(s * rel))
            Jz = model.reduced_jacobian(x)
            lams.append(np.linalg.eigvals(Jz).real.max())
            model.scale_enzyme(rid, np.exp(-s * rel))
        out[rid] = (lams[0] - lams[1]) / (2 * rel)
    model.set_x(x)
    return out


def total_parameter_fragility(model, p, x, rel=1e-3, t_end=500.0):
    """GLOBAL sensitivity d Re(lambda_lead)/d ln p including the shift of the
    steady state itself (re-solve). This is what a real perturbation does."""
    p0 = model.get(p)
    lams = []
    for s in (+1, -1):
        model.set_params({p: p0 * np.exp(s * rel)})
        xs = model.steady_state(t_end=t_end, x0=x)
        lams.append(np.linalg.eigvals(model.reduced_jacobian(xs)).real.max())
    model.set_params({p: p0})
    model.set_x(x)
    return (lams[0] - lams[1]) / (2 * rel)


def total_fragility(model, rid, x, target=None, rel=1e-3, solver="fixed_point"):
    """TOTAL d Re(λ)/d ln E_rid: local Jacobian change PLUS the shift of the
    steady state the enzyme change causes (the local fragility alone can even
    have the wrong sign when the state shift dominates, e.g. an influx).
    The eigenvalue followed is the leading one, or the one nearest `target`."""
    lam0 = (leading(model.reduced_jacobian(x), target)[0]
            [leading(model.reduced_jacobian(x), target)[1]])
    out = []
    for s in (+1, -1):
        model.scale_enzyme(rid, np.exp(s * rel))
        xs = (model.fixed_point(x) if solver == "fixed_point"
              else model.steady_state(x0=x, max_rounds=3))
        vals = np.linalg.eigvals(model.reduced_jacobian(xs))
        out.append(vals[np.argmin(np.abs(vals - lam0))].real)
        model.scale_enzyme(rid, np.exp(-s * rel))
    model.set_x(x)
    return (out[0] - out[1]) / (2 * rel)
