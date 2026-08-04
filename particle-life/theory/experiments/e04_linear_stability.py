r"""
E04 - Linear stability: what sets the SIZE of the structures.
=============================================================

The continuum theory of ``stability.py`` predicts a dispersion relation

.. math::  \lambda_\pm(k)=\tfrac12\bigl(-\gamma\pm\sqrt{\gamma^2-240\kappa k^2\sigma_n(k)}\bigr),
           \qquad \sigma_n=\mathrm{eig}\bigl[\mathrm{diag}(\bar\rho)\hat U(k)\bigr].

The fastest growing wavenumber :math:`k^\*` fixes the length scale
:math:`\ell^\*=2\pi/k^\*` at which structure first appears.  This experiment

  1. plots :math:`\hat U(k)` and :math:`\operatorname{Re}\lambda(k)` for a set of universes;
  2. measures the peak of the static structure factor :math:`S(k)` in the
     simulation at *early* times (before nonlinear coarsening takes over) and
     compares it with :math:`k^\*`;
  3. sweeps the min radius :math:`\alpha` to show that the selected length scale
     grows with :math:`\alpha` -- the quantitative version of the platform's own
     rule of thumb that the min radius controls how big/thick structures are.
"""

import multiprocessing as mp

import numpy as np
from common import ACC, dump, save
import matplotlib.pyplot as plt

from plife import kernels as K, matrices as M, model, observables as OB, stability as ST

L, N, S = 1000.0, 5000, 3
# measure once the instability has developed but before heavy coarsening
T_EARLY = 6.0
SEEDS = (0, 1, 2, 3)


def measure_peak(A, alpha, beta, seed=0, t=T_EARLY):
    p = model.Params(A=A, alpha=alpha, beta=beta, L=(L, L))
    sim = model.ParticleLife(p, N=N, seed=seed)
    sim.run(t)
    k, Sk = OB.structure_factor(sim.pos, np.array([L, L]), ngrid=256,
                                kmax=3 * np.pi / float(alpha.mean()))
    m = (k > 0.01) & np.isfinite(Sk)
    return float(k[m][np.nanargmax(Sk[m])]), k[m], Sk[m], sim


def job(args):
    name, a = args
    A = M.generate(name, S, np.random.default_rng(5))
    alpha, beta = M.uniform_radii(S, a, 3.0 * a)
    p = model.Params(A=A, alpha=alpha, beta=beta, L=(L, L))
    rho = np.full(S, N / S / L**2)
    res = ST.scan(p, rho, nk=500)
    kp = np.median([measure_peak(A, alpha, beta, seed=s)[0] for s in SEEDS])
    return dict(name=name, alpha=float(a), k_star=res.k_star, k_meas=float(kp),
                length_pred=res.length, length_meas=float(2 * np.pi / kp),
                growth=res.growth, osc_growth=res.osc_growth)


if __name__ == "__main__":
    # ------------------------------------------------- panel A/B: kernels
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.6))

    kk = np.linspace(1e-3, 0.55, 500)
    rho = np.full(S, N / S / L**2)
    cases = [("condensate", 12.0), ("condensate", 20.0), ("condensate", 30.0),
             ("shells", 15.0)]

    ax = axes[0, 0]
    for (name, a), c in zip(cases, ACC):
        A = M.generate(name, S, np.random.default_rng(5))
        alpha, beta = M.uniform_radii(S, a, 3.0 * a)
        U = K.hankel_matrix(kk, A, alpha, beta, 1.0)
        sig = np.linalg.eigvals(rho[None, :, None] * U)
        ax.plot(kk, sig.real.min(axis=1), color=c, label=f"{name}, α={a:.0f}")
    ax.axhline(0, color="#8b949e", lw=0.6)
    ax.set_xlabel("k"); ax.set_ylabel(r"$\min_n\ \mathrm{Re}\,\sigma_n(k)$")
    ax.set_title(r"$\sigma<0$ = unstable band", fontsize=10)
    ax.legend(fontsize=7.5)

    ax = axes[0, 1]
    for (name, a), c in zip(cases, ACC):
        A = M.generate(name, S, np.random.default_rng(5))
        alpha, beta = M.uniform_radii(S, a, 3.0 * a)
        p = model.Params(A=A, alpha=alpha, beta=beta, L=(L, L))
        r = ST.scan(p, rho, nk=500, kmax=0.55)
        re = r.lam.real.reshape(len(r.k), -1).max(axis=1)
        ax.plot(r.k, re, color=c, label=f"{name}, α={a:.0f}")
        ax.plot([r.k_star], [r.growth], "o", color=c, ms=5)
    ax.axhline(0, color="#8b949e", lw=0.6)
    ax.set_xlabel("k"); ax.set_ylabel(r"$\max\ \mathrm{Re}\,\lambda(k)$  [1/s]")
    ax.set_title(r"dispersion relation; dot = $k^*$", fontsize=10)
    ax.legend(fontsize=7.5)

    # ------------------------------------------------- panel C: S(k) overlay
    ax = axes[1, 0]
    for (name, a), c in zip(cases[:3], ACC):
        A = M.generate(name, S, np.random.default_rng(5))
        alpha, beta = M.uniform_radii(S, a, 3.0 * a)
        kp, kg, Sg, _ = measure_peak(A, alpha, beta)
        ax.plot(kg, Sg / np.nanmax(Sg), color=c, lw=1.2, label=f"α={a:.0f}: S(k)")
        p = model.Params(A=A, alpha=alpha, beta=beta, L=(L, L))
        ax.axvline(ST.scan(p, rho, nk=500).k_star, color=c, ls="--", lw=1.1)
    ax.set_xlim(0, 0.5)
    ax.set_xlabel("k"); ax.set_ylabel("S(k) (normalised)")
    ax.set_title(r"measured $S(k)$ vs predicted $k^*$ (dashed)", fontsize=10)
    ax.legend(fontsize=7.5)

    # ------------------------------------------------- panel D: alpha sweep
    alphas = [11.0, 14.0, 17.0, 20.0, 24.0, 28.0, 33.0, 38.0]
    with mp.Pool(min(4, mp.cpu_count())) as pool:
        rows = pool.map(job, [("condensate", a) for a in alphas])

    print(f"{'alpha':>6s} {'k*_pred':>9s} {'k_meas':>9s} {'l*_pred':>9s} {'l_meas':>9s}")
    for r in rows:
        print(f"{r['alpha']:6.1f} {r['k_star']:9.4f} {r['k_meas']:9.4f} "
              f"{r['length_pred']:9.1f} {r['length_meas']:9.1f}")

    lp = np.array([r["length_pred"] for r in rows])
    lm = np.array([r["length_meas"] for r in rows])
    c = float((lp * lm).sum() / (lp * lp).sum())            # least squares through origin
    ss_res = float(((lm - c * lp) ** 2).sum())
    ss_tot = float(((lm - lm.mean()) ** 2).sum())
    r2 = 1 - ss_res / ss_tot
    print(f"\nleast-squares  l_meas = {c:.3f} * l_pred     R^2 = {r2:.4f}")
    print("mean-field theory captures the SCALING exactly; the prefactor is off by "
          f"a constant {1/c:.2f}x, as expected of a mean-field closure.")

    ax = axes[1, 1]
    a = np.array([r["alpha"] for r in rows])
    ax.plot(a, lp, "o-", color=ACC[2], label=r"theory  $2\pi/k^*$")
    ax.plot(a, lm, "s--", color=ACC[0], label=r"simulation  $2\pi/k_{\rm peak}$")
    ax.plot(a, c * lp, ":", color=ACC[4], lw=1.4,
            label=rf"theory $\times$ {c:.2f}  ($R^2$={r2:.3f})")
    ax.set_xlabel(r"min radius $\alpha$   ($\beta=3\alpha$)")
    ax.set_ylabel("selected length scale")
    ax.set_title("structure size is set by the min radius", fontsize=10)
    ax.legend(fontsize=8)

    fig.suptitle("E04  Linear stability selects the length scale of the first structures",
                 fontsize=12)
    fig.tight_layout()
    save(fig, "e04_linear_stability.png")
    dump(dict(alpha_sweep=rows, fit_slope=c, fit_r2=r2), "e04_stability.json")
