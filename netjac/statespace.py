"""
statespace — state judgments BEYOND "eigenvalues -> stable/unstable".

For a linearised system dx/dt = J x (+ B u + noise), each function below
answers a different, testable question. Eigenvalues only tell the long-run
fate; these tell how the system gets there, how it can be seen, and how it
can be steered.

  phase_map            where each parameter pair sits in the trace-determinant
                       plane (node / focus / centre / saddle / source)
  reactivity           can a STABLE system still amplify a shock at first?
                       (max eigenvalue of (J+J^T)/2 > 0; Neubert & Caswell 1997)
  max_amplification    how much, and when (max_t ||e^{Jt}||, worst shock)
  quadratic_energy     a Lyapunov "energy" V = x^T P x with dV/dt = -x^T Q x;
                       diagonal P exists -> each variable's self-damping is a
                       dissipation channel (the "buffer" gamma/delta, exactly)
  noise_response       stationary covariance under random shocks, spectral
                       peak (noise-sustained quasi-cycles from a STABLE focus)
                       and the lead-lag between two variables
  gramians             controllability / observability: which input steers,
                       which single readout reveals the whole state
  pole_placement, lqr  the feedback gains that give a chosen damping
  critical_delay       the delay on one coupling at which a stable loop turns
                       into a sustained oscillation (delay-induced Hopf)
  coupling_sensitivity d lambda / d J_ij for every entry (which interaction
                       sets damping vs frequency)

All are exact linear algebra on J (✅) except critical_delay (numerical root +
DDE simulation check) and anything about what J means in a real system (⚠️).
"""
import numpy as np
from scipy.linalg import expm, solve_continuous_lyapunov, eig


# --------------------------------------------------------- phase portrait
def classify_2x2(J):
    tr, det = np.trace(J), np.linalg.det(J)
    disc = tr ** 2 - 4 * det
    if det < 0:
        return "saddle"
    if abs(tr) < 1e-12:
        return "centre"
    kind = "focus" if disc < 0 else "node"
    return ("stable " if tr < 0 else "unstable ") + kind


def phase_map(Jfun, a_vals, b_vals):
    """Grid of classify_2x2(Jfun(a, b)) -> labels array [len(b), len(a)]."""
    return np.array([[classify_2x2(Jfun(a, b)) for a in a_vals]
                     for b in b_vals])


# -------------------------------------------------- transient behaviour
def reactivity(J):
    """Initial maximal growth rate of ||x||: lambda_max((J+J^T)/2)."""
    return float(np.linalg.eigvalsh((J + J.T) / 2).max())


def max_amplification(J, t_max=None, n=2000):
    """max_t ||e^{Jt}||_2, the time it happens, and the worst initial shock
    direction (top right singular vector at that time)."""
    lam = np.linalg.eigvals(J)
    if t_max is None:
        slow = np.abs(lam.real).min()
        t_max = 10 / slow if slow > 0 else 100.0
    ts = np.linspace(0, t_max, n)
    norms = [np.linalg.norm(expm(J * t), 2) for t in ts]
    k = int(np.argmax(norms))
    U, s, Vt = np.linalg.svd(expm(J * ts[k]))
    return {"amp": norms[k], "t": ts[k], "worst_shock": Vt[0],
            "ts": ts, "norms": np.array(norms)}


# ----------------------------------------------------- energy function
def quadratic_energy(J, Q=None):
    """Solve J^T P + P J = -Q (Q = I by default). P > 0 <=> stable; V = x^T P x
    is a strict Lyapunov function. Also search for a DIAGONAL P (diagonal
    stability): then V = sum p_i x_i^2 and dV/dt splits into per-variable
    dissipation terms set by the diagonal of J."""
    n = len(J)
    Q = np.eye(n) if Q is None else Q
    P = solve_continuous_lyapunov(J.T, -Q)
    pos = bool(np.linalg.eigvalsh((P + P.T) / 2).min() > 0)
    diag = None
    # 2x2 closed form: off-diagonal cross terms cancel when
    # p1*J12 + p2*J21 = 0  -> p2/p1 = -J12/J21 (needs opposite signs)
    if n == 2 and J[0, 1] * J[1, 0] < 0:
        p = np.array([1.0, -J[0, 1] / J[1, 0]])
        Qd = -(J.T @ np.diag(p) + np.diag(p) @ J)
        diag = {"p": p, "dissipation": np.diag(Qd),
                "cross_terms_cancel": bool(abs(Qd[0, 1]) < 1e-12)}
    return {"P": P, "positive": pos, "diag": diag}


# ------------------------------------------------- stochastic forcing
def noise_response(J, D=None, freqs=None, i=0, j=1):
    """Random shocks with covariance D: stationary covariance C (Lyapunov
    eq.), power spectrum of each variable, spectral peak (quasi-cycle) and
    the lag at which corr(x_i(t), x_j(t+tau)) peaks."""
    n = len(J)
    D = np.eye(n) if D is None else D
    C = solve_continuous_lyapunov(J, -D)
    lam = np.linalg.eigvals(J)
    if freqs is None:
        w0 = max(np.abs(lam.imag).max(), np.abs(lam.real).min(), 1e-3)
        freqs = np.linspace(1e-4, 4 * w0, 2000)
    S = []
    for w in freqs:
        R = np.linalg.inv(1j * w * np.eye(n) - J)
        S.append(np.real(np.diag(R @ D @ R.conj().T)))
    S = np.array(S)
    peaks = [freqs[np.argmax(S[:, k])] for k in range(n)]
    # cross-covariance C(tau) = e^{J tau} C for tau >= 0
    w_slow = np.abs(lam.real).min()
    taus = np.linspace(-6 / max(w_slow, 1e-3), 6 / max(w_slow, 1e-3), 2401)
    corr = []
    for tau in taus:
        Ct = expm(J * tau) @ C if tau >= 0 else (expm(J * -tau) @ C).T
        corr.append(Ct[j, i] / np.sqrt(C[i, i] * C[j, j]))
    corr = np.array(corr)
    k = int(np.argmax(np.abs(corr)))
    return {"C": C, "freqs": freqs, "S": S, "peak_freq": peaks,
            "taus": taus, "xcorr": corr, "lag": taus[k], "lag_corr": corr[k],
            "contemp_corr": C[i, j] / np.sqrt(C[i, i] * C[j, j])}


# ------------------------------------------- what can be seen / steered
def gramians(J, B=None, Cm=None):
    """Controllability (from inputs B) and observability (from readouts Cm):
    rank tests + Gramian eigenvalues (how much 'effort' / 'signal')."""
    n = len(J)
    out = {}
    if B is not None:
        B = np.atleast_2d(B).reshape(n, -1)
        K = np.hstack([np.linalg.matrix_power(J, k) @ B for k in range(n)])
        Wc = solve_continuous_lyapunov(J, -B @ B.T)
        out["controllable"] = bool(np.linalg.matrix_rank(K) == n)
        out["Wc_eigs"] = np.linalg.eigvalsh(Wc)
    if Cm is not None:
        Cm = np.atleast_2d(Cm).reshape(-1, n)
        O = np.vstack([Cm @ np.linalg.matrix_power(J, k) for k in range(n)])
        Wo = solve_continuous_lyapunov(J.T, -Cm.T @ Cm)
        out["observable"] = bool(np.linalg.matrix_rank(O) == n)
        out["Wo_eigs"] = np.linalg.eigvalsh(Wo)
    return out


def pole_placement(J, B, poles):
    """Feedback u = -K x giving closed-loop eigenvalues `poles`."""
    from scipy.signal import place_poles
    B = np.atleast_2d(B).reshape(len(J), -1)
    K = place_poles(J, B, poles).gain_matrix
    return K, np.linalg.eigvals(J - B @ K)


def lqr(J, B, Q, R):
    from scipy.linalg import solve_continuous_are
    B = np.atleast_2d(B).reshape(len(J), -1)
    P = solve_continuous_are(J, B, Q, R)
    K = np.linalg.solve(np.atleast_2d(R), B.T @ P)
    return K, np.linalg.eigvals(J - B @ K)


# ---------------------------------------------------------- delays
def critical_delay(J, i, j, w_max=None):
    """Delay tau on the single coupling J[i,j] (x_i' depends on
    x_j(t - tau)). Characteristic equation det(lambda I - J0 - J1 e^{-lambda
    tau}) = 0 with J1 holding only that entry. Returns the smallest tau at
    which a root crosses the imaginary axis (lambda = i w), and w."""
    from scipy.optimize import brentq
    J0 = J.copy()
    J0[i, j] = 0.0
    J1 = np.zeros_like(J)
    J1[i, j] = J[i, j]
    n = len(J)
    # det(iw I - J0 - J1 z) is affine in z = e^{-i w tau}: a(w) + b(w) z
    def ab(w):
        M = 1j * w * np.eye(n) - J0
        a = np.linalg.det(M)
        b = np.linalg.det(M - J1) - a
        return a, b
    w_max = w_max or 10 * (np.abs(np.linalg.eigvals(J)).max() + 1)
    ws = np.linspace(1e-6, w_max, 20000)
    g = np.array([abs(ab(w)[0]) - abs(ab(w)[1]) for w in ws])   # |z| = 1
    sols = []
    for k in np.where(np.sign(g[:-1]) != np.sign(g[1:]))[0]:
        w = brentq(lambda q: abs(ab(q)[0]) - abs(ab(q)[1]), ws[k], ws[k + 1])
        a, b = ab(w)
        z = -a / b                                   # = e^{-i w tau}
        tau = (-np.angle(z)) % (2 * np.pi) / w
        sols.append((tau, w))
    if not sols:
        return None
    return min(sols)


def simulate_dde(J, i, j, tau, x0, t_end, dt=0.01):
    """Euler integration with the J[i,j] coupling delayed by tau."""
    n = len(J)
    J0 = J.copy(); J0[i, j] = 0.0
    lag = int(round(tau / dt))
    steps = int(t_end / dt)
    X = np.zeros((steps + 1, n))
    X[0] = x0
    for k in range(steps):
        xd = X[k - lag] if k >= lag else X[0]
        dx = J0 @ X[k]
        dx[i] += J[i, j] * xd[j]
        X[k + 1] = X[k] + dt * dx
    return np.arange(steps + 1) * dt, X


# --------------------------------------------------- sensitivities
def coupling_sensitivity(J):
    """d lambda_lead / d J_ij = conj(w_i) v_j / (w^H v) for every entry:
    real part -> effect on damping, imaginary part -> effect on frequency."""
    vals, wl, vr = eig(J, left=True, right=True)
    k = int(np.argmax(vals.real + 1e-12 * vals.imag))
    w, v = wl[:, k], vr[:, k]
    return np.outer(np.conj(w), v) / (np.conj(w) @ v), vals[k]


def var_jacobian(Z, dt, p=1):
    """Fit Z_t = c + A Z_{t-1} (+...) by OLS; J = logm(A)/dt for p = 1.
    Returns J, A, residual covariance (-> noise matrix D), equilibrium."""
    from scipy.linalg import logm
    Y = Z[p:]
    X = np.hstack([np.ones((len(Y), 1))] + [Z[p - k:len(Z) - k]
                                            for k in range(1, p + 1)])
    Bc, *_ = np.linalg.lstsq(X, Y, rcond=None)
    E = Y - X @ Bc
    A = Bc[1:1 + Z.shape[1]].T
    J = np.real(logm(A)) / dt
    c = Bc[0]
    eq = np.linalg.solve(np.eye(len(A)) - A, c) if p == 1 else None
    return {"J": J, "A": A, "Sigma": E.T @ E / len(E), "eq": eq,
            "resid": E}
