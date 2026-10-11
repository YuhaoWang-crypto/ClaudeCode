"""
relax — perturbation-relaxation diagnostics: push the system, read the
eigenvalues off the way it comes back.

pulse_response()  exact model: bolus on one species, record a few observables
hankel_dmd()      eigenvalues from a (possibly partial, noisy) relaxation
                  record via time-delay embedding + DMD (= multichannel Prony /
                  ERA). With only 1-2 measured analytes it still recovers the
                  slow, dominant eigenvalues because every mode leaks into
                  every observable (unless the observable is orthogonal to the
                  mode's eigenvector -> see D3 probe placement).
fit_decay()       single-exponential relaxation rate for short clinical traces
                  (CGM after a meal), with a sanity flag.
"""
import numpy as np


def pulse_response(model, x_ss, species, amount, t_end, n, observe):
    """Add `amount` (concentration units) to `species` at t=0 and simulate.
    Returns (t, Y) with Y[:, j] = deviation of observable j from steady state."""
    x = x_ss.copy()
    x[model.state_ids.index(species)] += amount
    model.set_x(x)
    res = model.rr.simulate(0, t_end, n, selections=["time"] + list(observe))
    t = res[:, 0]
    base = np.array([_value_at_ss(model, x_ss, o) for o in observe])
    model.set_x(x_ss)
    return t, res[:, 1:] - base[None, :]


def _value_at_ss(model, x_ss, obs):
    model.set_x(x_ss)
    return model.rr[obs]


def hankel_dmd(Y, dt, rank, delays=None):
    """Y: (T, p) samples of p observables at spacing dt. Returns continuous-
    time eigenvalues sorted by real part (slowest first)."""
    Y = np.atleast_2d(np.asarray(Y, float))
    if Y.shape[0] < Y.shape[1]:
        Y = Y.T
    T, p = Y.shape
    d = delays or max(2, int(np.ceil(2 * rank / p)) + 2)
    H = np.hstack([Y[i:T - d + i + 1] for i in range(d)]).T   # (d*p, T-d+1)
    X0, X1 = H[:, :-1], H[:, 1:]
    U, s, Vh = np.linalg.svd(X0, full_matrices=False)
    r = min(rank, (s > 1e-12 * s[0]).sum())
    U, s, Vh = U[:, :r], s[:r], Vh[:r]
    A = U.T @ X1 @ Vh.T / s
    mu = np.linalg.eigvals(A)
    lam = np.log(mu.astype(complex)) / dt
    return lam[np.argsort(-lam.real)]


def fit_decay(t, y):
    """Fit y(t) = c + a*exp(-k t) by profile least squares on k.
    Returns dict(k, tau, a, c, r2)."""
    t, y = np.asarray(t, float), np.asarray(y, float)
    best = None
    for k in np.exp(np.linspace(np.log(1e-3), np.log(1.0), 300)):
        e = np.exp(-k * (t - t[0]))
        A = np.vstack([np.ones_like(e), e]).T
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        sse = ((A @ coef - y) ** 2).sum()
        if best is None or sse < best[0]:
            best = (sse, k, coef)
    sse, k, (c, a) = best
    sst = ((y - y.mean()) ** 2).sum() + 1e-12
    return {"k": k, "tau": 1.0 / k, "a": a, "c": c, "r2": 1 - sse / sst}


def fit_modes(t, Y, rates0, osc0=None):
    """Global nonlinear least-squares fit of a relaxation record
        y_j(t) = sum_i a_ij exp(r_i t) + sum_q exp(s_q t)(b_qj cos w_q t + c_qj sin w_q t)
    with rates SHARED across observables (they are eigenvalues of one
    Jacobian; amplitudes are observable-specific). Unbiased under white noise,
    unlike plain DMD. rates0: initial real rates (<0); osc0: list of (s, w).
    Amplitudes are solved linearly (variable projection). Returns dict with
    'rates', 'osc', 'rss'."""
    from scipy.optimize import least_squares
    t = np.asarray(t, float) - t[0]
    Y = np.atleast_2d(np.asarray(Y, float))
    if Y.shape[0] != len(t):
        Y = Y.T
    osc0 = osc0 or []
    nr, no = len(rates0), len(osc0)

    def basis(th):
        cols = [np.exp(-np.exp(th[i]) * t) for i in range(nr)]
        for q in range(no):
            s, w = -np.exp(th[nr + 2 * q]), th[nr + 2 * q + 1]
            cols += [np.exp(s * t) * np.cos(w * t), np.exp(s * t) * np.sin(w * t)]
        return np.array(cols).T

    def resid(th):
        Phi = basis(th)
        coef, *_ = np.linalg.lstsq(Phi, Y, rcond=None)
        return (Phi @ coef - Y).ravel()
    th0 = [np.log(-r) for r in rates0]
    for s, w in osc0:
        th0 += [np.log(-s), w]
    sol = least_squares(resid, th0, method="lm", max_nfev=4000)
    th = sol.x
    rates = [-np.exp(th[i]) for i in range(nr)]
    osc = [(-np.exp(th[nr + 2 * q]), abs(th[nr + 2 * q + 1])) for q in range(no)]
    return {"rates": sorted(rates, reverse=True), "osc": osc,
            "rss": (sol.fun ** 2).sum()}
