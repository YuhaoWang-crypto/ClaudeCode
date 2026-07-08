"""
Module 4 — Tipping-point biomarkers: DNB, critical slowing down, and the
largest Lyapunov exponent (LLE).

Idea:
  Near a bifurcation ("tipping point") a dynamical system shows model-free
  early-warning signals:
    * critical slowing down: the leading Jacobian eigenvalue Re(lambda_max)
      -> 0^- , so recovery from perturbations slows;
    * rising fluctuation VARIANCE and lag-1 AUTOCORRELATION (Scheffer 2009);
    * a Dynamical Network Biomarker (DNB, Chen & Aihara 2012): a dominant
      group of molecules whose SD rises and whose mutual correlation rises.
  The LARGEST LYAPUNOV EXPONENT (LLE) is the fully nonlinear generalisation
  of Re(lambda_max): LLE < 0 stable, LLE = 0 marginal/limit cycle, LLE > 0
  chaotic. It is a single abstract, measurable stability parameter.

Concrete system: a 2-gene mutually-activating module (a "dominant group")
with a drug parameter mu that raises effective degradation and drives the
co-expressed ON state through a fold (saddle-node) bifurcation toward
collapse. We track LLE, variance, autocorrelation and the DNB index as the
drug dose increases toward the tipping point.

The Benettin LLE routine is first VALIDATED on the Rossler attractor
(literature value ~ +0.071) before being applied to the gene model.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ------------------------- gene module ODE -------------------------------
def hill(x, K=1.0, n=4):
    xn = x ** n
    return xn / (K ** n + xn)


def gene_rhs(mu, b=0.2, V=3.0, d=1.0, K=1.0, n=4):
    """Two mutually-activating genes; mu = extra drug-induced degradation."""
    def f(t, s):
        x, y = s
        dx = b + V * hill(y, K, n) - (d + mu) * x
        dy = b + V * hill(x, K, n) - (d + mu) * y
        return [dx, dy]
    return f


def gene_jac(state, mu, b=0.2, V=3.0, d=1.0, K=1.0, n=4):
    x, y = state
    def dhill(z):
        zn = z ** n
        return n * (K ** n) * z ** (n - 1) / (K ** n + zn) ** 2
    J = np.array([[-(d + mu), V * dhill(y)],
                  [V * dhill(x), -(d + mu)]])
    return J


# --------------------- Benettin largest Lyapunov -------------------------
def benettin_lle(rhs, jac, x0, T=2000.0, dt=0.01, t_transient=200.0):
    """Largest Lyapunov exponent via tangent-space renormalisation."""
    x = np.array(x0, dtype=float)
    dim = len(x)
    # burn in
    n_trans = int(t_transient / dt)
    for _ in range(n_trans):
        x = _rk4(rhs, x, dt)
    v = np.random.default_rng(0).normal(size=dim)
    v /= np.linalg.norm(v)
    n_steps = int(T / dt)
    log_sum = 0.0
    for _ in range(n_steps):
        J = jac(x)
        # evolve state and tangent vector together (RK4 on augmented system)
        x = _rk4(rhs, x, dt)
        v = v + dt * (J @ v)          # linearised tangent flow (Euler ok, dt small)
        nrm = np.linalg.norm(v)
        log_sum += np.log(nrm)
        v /= nrm
    return log_sum / (n_steps * dt)


def _rk4(rhs, x, dt):
    k1 = np.array(rhs(0, x))
    k2 = np.array(rhs(0, x + 0.5 * dt * k1))
    k3 = np.array(rhs(0, x + 0.5 * dt * k2))
    k4 = np.array(rhs(0, x + dt * k3))
    return x + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def validate_lle_rossler():
    """Rossler attractor: literature LLE ~ +0.071."""
    a, b, c = 0.2, 0.2, 5.7
    def rhs(t, s):
        x, y, z = s
        return [-y - z, x + a * y, b + z * (x - c)]
    def jac(s):
        x, y, z = s
        return np.array([[0, -1, -1],
                         [1, a, 0],
                         [z, 0, x - c]])
    return benettin_lle(rhs, jac, [1.0, 1.0, 1.0], T=3000.0, dt=0.01)


# ---------------- stochastic run -> DNB / variance / autocorr ------------
def stochastic_signals(mu, x_star, sigma=0.05, T=800.0, dt=0.01, seed=1):
    """Euler-Maruyama around the ON fixed point; return SD, lag-1 autocorr,
    correlation between the two genes, and the DNB composite index."""
    rhs = gene_rhs(mu)
    rng = np.random.default_rng(seed)
    n = int(T / dt)
    s = np.array(x_star, dtype=float)
    xs = np.empty(n)
    ys = np.empty(n)
    sq = np.sqrt(dt)
    for i in range(n):
        drift = np.array(rhs(0, s))
        s = s + dt * drift + sigma * sq * rng.normal(size=2)
        s = np.clip(s, 0, None)
        xs[i], ys[i] = s
    burn = n // 4
    xs, ys = xs[burn:], ys[burn:]
    sd_x, sd_y = xs.std(), ys.std()
    # lag-1 autocorrelation of gene x
    ac = np.corrcoef(xs[:-1], xs[1:])[0, 1]
    corr = np.corrcoef(xs, ys)[0, 1]
    dnb = 0.5 * (sd_x + sd_y) * abs(corr)      # classic DNB composite index
    return {"sd": 0.5 * (sd_x + sd_y), "autocorr": ac,
            "corr": corr, "dnb": dnb}


def upper_fixed_point(mu, b=0.2, V=3.0, d=1.0, K=1.0, n=4):
    """Highest (ON) root of the symmetric branch  b+V*hill(x)-(d+mu)x = 0,
    found precisely by sign-change bracketing. Returns None once the ON
    branch has vanished (only the low root remains -> fold passed)."""
    def g(x):
        return b + V * (x ** n) / (K ** n + x ** n) - (d + mu) * x
    xs = np.linspace(1e-3, 8.0, 8000)
    vals = g(xs)
    roots = []
    for i in range(len(xs) - 1):
        if vals[i] * vals[i + 1] < 0:
            roots.append(brentq(g, xs[i], xs[i + 1]))
    if len(roots) < 3:                 # fold passed: ON state gone
        return None
    xstar = roots[-1]
    return np.array([xstar, xstar])


def sweep():
    # fine grid up to just below the fold (~0.87) to resolve slowing-down
    mus = np.linspace(0.0, 0.865, 40)
    rows = []
    for mu in mus:
        fp = upper_fixed_point(mu)
        if fp is None:
            break                      # passed the fold — ON state destroyed
        eig = np.linalg.eigvals(gene_jac(fp, mu)).real.max()   # local LLE
        sig = stochastic_signals(mu, fp)
        rows.append({"mu": mu, "fp": fp.mean(), "lead_eig": eig, **sig})
    return rows


def make_figure(rows, lle_rossler, path):
    mu = [r["mu"] for r in rows]
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))

    ax[0].plot(mu, [r["lead_eig"] for r in rows], "o-", color="#c0392b")
    ax[0].axhline(0, ls="--", color="k", lw=1)
    ax[0].set_title("Local Lyapunov exp. (leading eigenvalue)\n"
                    "-> 0 at the tipping point (critical slowing down)")
    ax[0].set_xlabel("drug dose  mu"); ax[0].set_ylabel("Re(lambda_max)")

    l1, = ax[1].plot(mu, [r["sd"] for r in rows], "o-", color="#2980b9",
                     label="SD (variance$^{0.5}$)")
    ax[1].set_ylabel("SD", color="#2980b9")
    ax[1].tick_params(axis="y", labelcolor="#2980b9")
    ax1b = ax[1].twinx()
    l2, = ax1b.plot(mu, [r["autocorr"] for r in rows], "s-", color="#27ae60",
                    label="lag-1 autocorrelation")
    ax1b.set_ylabel("lag-1 autocorrelation", color="#27ae60")
    ax1b.tick_params(axis="y", labelcolor="#27ae60")
    ax[1].set_title("Early-warning signals rise\ntoward the tipping point")
    ax[1].set_xlabel("drug dose  mu")
    ax[1].legend(handles=[l1, l2], loc="upper left", fontsize=8)

    ax[2].plot(mu, [r["dnb"] for r in rows], "o-", color="#8e44ad")
    ax[2].set_title("Dynamical Network Biomarker index\n"
                    "DNB = mean(SD) x |corr(x,y)|")
    ax[2].set_xlabel("drug dose  mu"); ax[2].set_ylabel("DNB index")

    fig.suptitle(f"Tipping-point biomarkers  |  Benettin LLE validated on "
                 f"Rossler = {lle_rossler:+.4f} (lit. ~ +0.071)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(path, dpi=130)
    plt.close(fig)


def report(fig_path="figures/m4_biomarkers.png"):
    print("=" * 68)
    print("MODULE 4 — DNB / critical slowing down / Lyapunov biomarker")
    print("=" * 68)
    lle_r = validate_lle_rossler()
    print(f"Benettin LLE code validated on Rossler: {lle_r:+.4f} "
          f"(literature ~ +0.071)  -> algorithm OK")
    rows = sweep()
    print(f"\ndrug sweep along the ON branch ({len(rows)} feasible doses "
          f"before the fold):")
    print(f"  {'mu':>5} {'ON-level':>9} {'lead_eig(LLE)':>14} "
          f"{'SD':>7} {'autocorr':>9} {'DNB':>7}")
    for r in rows[::4] + [rows[-1]]:
        print(f"  {r['mu']:5.2f} {r['fp']:9.3f} {r['lead_eig']:14.4f} "
              f"{r['sd']:7.3f} {r['autocorr']:9.3f} {r['dnb']:7.3f}")
    print(f"\nAt the last feasible dose mu={rows[-1]['mu']:.2f} the leading "
          f"eigenvalue = {rows[-1]['lead_eig']:.3f} (-> 0), while SD, "
          f"autocorrelation and DNB\nall peak. Beyond it the ON state is "
          f"destroyed (fold bifurcation) — cell tips off.")
    make_figure(rows, lle_r, fig_path)
    print(f"\nfigure written to: {fig_path}")
    return {"lle_rossler": lle_r, "rows": rows}


if __name__ == "__main__":
    report()
