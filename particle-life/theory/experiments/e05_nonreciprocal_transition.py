r"""
E05 - The non-reciprocal transition: when structures start to move.
===================================================================

Linear theory says a growing mode can only *oscillate* -- i.e. the pattern can only
travel -- if the eigenvalue :math:`\sigma_n(k)` of :math:`\mathrm{diag}(\bar\rho)\hat U(k)`
is complex **and**

.. math::  \sigma_i^2 \;>\; \frac{\gamma^2}{60\kappa k^2}\,\sigma_r .

Complex :math:`\sigma` requires :math:`\hat U\neq\hat U^{T}`; for two species the
condition is exactly

.. math::  D(k)=\Bigl[\tfrac12(\rho_1\hat u_{11}-\rho_2\hat u_{22})\Bigr]^2
                 +\rho_1\rho_2\,\hat u_{12}\hat u_{21}\;<\;0 ,

which needs :math:`\hat u_{12}\hat u_{21}<0`: **one species must chase while the
other flees**.  ``D=0`` is an exceptional point -- two static demixing modes merge
and re-emerge as one travelling mode (the Particle-Life instance of the
non-reciprocal phase transition of Fruchart et al., Nature 592, 363 (2021)).

The inequality also contains an explicit :math:`\gamma^2` -- **friction suppresses
the travelling instability**.  That gives a sharp, falsifiable prediction: hold
the matrix fixed, lower ``frictionFactor``, and the structures must switch from
standing to travelling at a computable threshold.  Panels 3-4 test exactly that.
"""

import multiprocessing as mp

import numpy as np
from common import ACC, dump, save, scatter_snapshot
import matplotlib.pyplot as plt

from plife import matrices as M, model, observables as OB, stability as ST

L, N, S = 1200.0, 6000, 2
ALPHA, BETA = 15.0, 45.0
A_CHASE = np.array([[-0.1, 0.9], [-0.9, -0.1]])
LAMBDAS = [0.30, 0.20, 0.15, 0.12, 0.10, 0.08, 0.06, 0.04, 0.02]


def sim_motility(lam, seed=0, t_eq=25.0, t_meas=10.0):
    alpha, beta = M.uniform_radii(S, ALPHA, BETA)
    p = model.Params(A=A_CHASE, alpha=alpha, beta=beta, L=(L, L), friction=lam)
    sim = model.ParticleLife(p, N=N, seed=seed)
    sim.run(t_eq)
    o1 = OB.observables(sim)
    sim.run(t_meas)
    o2 = OB.observables(sim)
    return dict(lam=lam, seed=seed,
                motility=0.5 * (o1["motility"] + o2["motility"]),
                cluster_speed=0.5 * (o1["cluster_speed"] + o2["cluster_speed"]),
                speed=0.5 * (o1["speed"] + o2["speed"]))


def job(args):
    return sim_motility(*args)


if __name__ == "__main__":
    alpha, beta = M.uniform_radii(S, ALPHA, BETA)
    rho = np.full(S, N / S / L**2)

    # ---------------------------------------------- theory: friction threshold
    theory = []
    for lam in np.linspace(0.005, 0.35, 120):
        p = model.Params(A=A_CHASE, alpha=alpha, beta=beta, L=(L, L), friction=lam)
        r = ST.scan(p, rho, nk=500)
        theory.append(dict(lam=float(lam), gamma=p.gamma, osc=r.osc_growth,
                           osc_k=r.osc_k, c=r.osc_speed))
    lam_c = None
    for a, b in zip(theory[:-1], theory[1:]):
        if (a["osc"] > 0) != (b["osc"] > 0):
            lam_c = 0.5 * (a["lam"] + b["lam"])
    print(f"predicted travelling-instability threshold:  lambda_c ~ {lam_c}")

    # ---------------------------------------------- simulation sweep
    jobs = [(l, s) for l in LAMBDAS for s in (0, 1)]
    with mp.Pool(min(4, mp.cpu_count())) as pool:
        runs = pool.map(job, jobs)
    agg = {}
    for r in runs:
        agg.setdefault(r["lam"], []).append(r["motility"])
    sim_rows = [dict(lam=l, motility=float(np.mean(v)), sd=float(np.std(v)))
                for l, v in sorted(agg.items())]
    print(f"\n{'lambda':>7s} {'motility':>10s}")
    for r in sim_rows:
        print(f"{r['lam']:7.2f} {r['motility']:10.5f}")

    # ---------------------------------------------- (A12, A21) maps
    g = np.linspace(-1, 1, 41)
    p0 = model.Params(A=np.array([[-0.1, 0.0], [0.0, -0.1]]),
                      alpha=alpha, beta=beta, L=(L, L), friction=0.06)
    growth, osc, speed = ST.two_species_phase_boundary(g, g, p0, rho, nk=300)

    # ---------------------------------------------- figure
    fig = plt.figure(figsize=(13.4, 7.6))
    gs = fig.add_gridspec(2, 3)

    ax = fig.add_subplot(gs[0, 0])
    im = ax.pcolormesh(g, g, np.where(np.isfinite(osc), osc, np.nan),
                       cmap="magma", shading="auto")
    ax.contour(g, g, np.where(np.isfinite(osc), osc, -1), levels=[0.0],
               colors="#39d353", linewidths=1.6)
    ax.set_xlabel(r"$A_{12}$"); ax.set_ylabel(r"$A_{21}$")
    ax.set_title(r"growth rate of the oscillatory branch" + "\n"
                 r"green = exceptional/onset line, $\lambda_{\rm fric}=0.06$", fontsize=9)
    fig.colorbar(im, ax=ax, fraction=0.046)

    ax = fig.add_subplot(gs[0, 1])
    im = ax.pcolormesh(g, g, np.where(osc > 0, speed, np.nan), cmap="viridis",
                       shading="auto")
    ax.plot([-1, 1], [1, -1], color="#8b949e", ls="--", lw=0.8)
    ax.set_xlabel(r"$A_{12}$"); ax.set_ylabel(r"$A_{21}$")
    ax.set_title("phase velocity  c = Im λ / k\nof the travelling mode", fontsize=9)
    fig.colorbar(im, ax=ax, fraction=0.046)

    ax = fig.add_subplot(gs[0, 2])
    ax.plot([t["lam"] for t in theory], [max(t["osc"], -0.05) for t in theory],
            color=ACC[2], lw=1.8, label=r"theory: Re$\lambda_{\rm osc}$")
    ax.axhline(0, color="#8b949e", lw=0.6)
    if lam_c:
        ax.axvline(lam_c, color=ACC[0], ls="--", lw=1.2,
                   label=rf"$\lambda_c\approx{lam_c:.3f}$")
    ax2 = ax.twinx()
    ax2.errorbar([r["lam"] for r in sim_rows], [r["motility"] for r in sim_rows],
                 yerr=[r["sd"] for r in sim_rows], color=ACC[1], marker="o", ms=4,
                 lw=1.4, capsize=2.5, label="simulation motility")
    ax2.set_ylabel("measured motility", color=ACC[1])
    ax2.grid(False)
    ax.set_xlabel(r"friction $\lambda$")
    ax.set_ylabel(r"Re$\lambda_{\rm osc}$  [1/s]", color=ACC[2])
    ax.set_title("friction kills the travelling instability", fontsize=9)
    ax.legend(fontsize=7.5, loc="upper right")

    for col, lam in enumerate((0.30, 0.10, 0.03)):
        ax = fig.add_subplot(gs[1, col])
        p = model.Params(A=A_CHASE, alpha=alpha, beta=beta, L=(L, L), friction=lam)
        sim = model.ParticleLife(p, N=N, seed=0)
        sim.run(35.0)
        scatter_snapshot(ax, sim, size=1.2)
        o = OB.observables(sim)
        r = ST.scan(p, rho, nk=400)
        ax.set_title(rf"$\lambda$={lam:.2f}  ($\gamma$={p.gamma:.1f})   "
                     rf"motility={o['motility']:.4f}" + "\n"
                     rf"theory: Re$\lambda_{{\rm osc}}$={r.osc_growth:+.3f}, "
                     rf"c={r.osc_speed:.1f}", fontsize=8.5)

    fig.suptitle("E05  Non-reciprocity + low friction = travelling structures", fontsize=12)
    fig.tight_layout()
    save(fig, "e05_nonreciprocal_transition.png")
    dump(dict(lambda_c=lam_c, theory=theory, simulation=sim_rows), "e05_transition.json")
