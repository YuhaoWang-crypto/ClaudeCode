r"""
E03 - The reciprocity theorem: why Particle Life needs a broken third law.
=========================================================================

**Theorem.**  If the interaction is reciprocal, i.e. ``f_ab(r) = f_ba(r)`` for all
species pairs (equivalently ``A``, ``alpha`` and ``beta`` are all symmetric), then

.. math::
    E \;=\; \tfrac12\sum_i |v_i|^2 \;+\; \sum_{i<j} U_{s_is_j}(r_{ij}),
    \qquad
    \dot E \;=\; -\gamma\sum_i |v_i|^2 \;\le\; 0 .

``E`` is bounded below (``U`` is bounded and ``N`` is finite), so by LaSalle's
invariance principle every trajectory converges to the largest invariant set
inside ``{v = 0}``: the mechanical equilibria of ``V``.  **No persistent motion of
any kind is possible** -- no crawling, no chasing, no rotation, no oscillation.

Everything alive in Particle Life therefore comes from the antisymmetric part of
the interaction.  This experiment measures that statement.

Protocol
--------
Start from one random universe, split ``A = A_S + A_N`` and interpolate
``A(theta) = A_S + theta A_N`` with **symmetric radii**, so that ``theta = 0`` is
exactly reciprocal.  Sweep ``theta`` and measure (i) whether ``E(t)`` is monotone
and (ii) the steady-state cluster motility.
"""

import multiprocessing as mp

import numpy as np
from common import ACC, dump, save, scatter_snapshot
import matplotlib.pyplot as plt

from plife import matrices as M, model, observables as OB

N, L, T_EQ, T_MEAS = 2000, 800.0, 30.0, 12.0
S = 5
SEEDS = (0, 1, 2)
THETAS = np.linspace(0.0, 1.0, 9)
# The theorem holds at any friction, but lambda = 0.3 (the platform default) is
# deep in the over-damped regime where even a strongly non-reciprocal drive only
# produces a crawl.  lambda = 0.10 keeps everything else identical while making
# the effect large enough to see; E05 maps the friction axis itself.
FRICTION = 0.10
T_LONG = 150.0


def build(theta, seed):
    rng = np.random.default_rng(11)
    A0 = M.g_random(S, rng)
    A = M.symmetrise(A0) + theta * M.antisymmetrise(A0)
    alpha, beta = M.random_radii(S, rng=np.random.default_rng(12), symmetric=True)
    p = model.Params(A=A, alpha=alpha, beta=beta, L=(L, L), friction=FRICTION)
    return p, model.ParticleLife(p, N=N, seed=seed)


def long_trace(theta):
    """A single long run, recording <|v|> - shows decay-to-rest vs saturation."""
    p, sim = build(theta, 0)
    ts, spd = [], []
    for i in range(int(T_LONG / p.dt)):
        sim.step(1)
        if i % 30 == 0:
            ts.append(sim.t)
            spd.append(sim.mean_speed())
    return dict(theta=float(theta), ts=ts, spd=spd)


def one(args):
    theta, seed = args
    p, sim = build(theta, seed)
    ts, spd, pin, pdis = [], [], [], []
    nsteps = int((T_EQ + T_MEAS) / p.dt)
    for i in range(nsteps):
        sim.step(1)
        if i % 20 == 0:
            ts.append(sim.t)
            spd.append(sim.mean_speed())
            pin.append(sim.injected_power() / sim.N)
            pdis.append(sim.dissipated_power() / sim.N)
    obs = OB.observables(sim)
    tail = slice(int(0.65 * len(ts)), None)
    return dict(theta=float(theta), seed=seed, nu=M.nonreciprocity(p.A),
                motility=obs["motility"], speed=obs["speed"],
                p_inj=float(np.mean(np.array(pin)[tail])),
                p_diss=float(np.mean(np.array(pdis)[tail])),
                ts=ts, spd=spd, pin=pin, pdis=pdis)


def dispatch(args):
    return long_trace(args[1]) if args[0] == "trace" else one(args[1])


TRACE_THETAS = (0.0, 0.5, 1.0)

if __name__ == "__main__":
    jobs = [("sweep", (t, s)) for t in THETAS for s in SEEDS]
    jobs += [("trace", t) for t in TRACE_THETAS]
    with mp.Pool(min(4, mp.cpu_count())) as pool:
        res_all = pool.map(dispatch, jobs)
    out = [r for r in res_all if "seed" in r]
    traces = {round(r["theta"], 4): r for r in res_all if "seed" not in r}

    by_theta = {}
    for r in out:
        by_theta.setdefault(round(r["theta"], 4), []).append(r)

    print(f"{'theta':>6s} {'nu':>6s} {'motility':>10s} {'+/-':>8s} "
          f"{'P_inj/N':>10s} {'P_diss/N':>10s}")
    summary = []
    for th in sorted(by_theta):
        rs = by_theta[th]
        mo = np.array([r["motility"] for r in rs])
        pi = np.array([r["p_inj"] for r in rs])
        pd = np.array([r["p_diss"] for r in rs])
        summary.append(dict(theta=th, nu=rs[0]["nu"], motility=mo.mean(),
                            motility_sd=mo.std(), p_inj=pi.mean(), p_diss=pd.mean()))
        print(f"{th:6.2f} {rs[0]['nu']:6.3f} {mo.mean():10.5f} {mo.std():8.5f} "
              f"{pi.mean():10.6f} {pd.mean():10.6f}")

    # ------------------------------------------------------------- figures
    fig = plt.figure(figsize=(13.2, 7.4))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1])

    ax = fig.add_subplot(gs[0, 0])
    for th, c, lab in ((0.0, ACC[2], r"$\theta=0$  reciprocal"),
                       (0.5, ACC[3], r"$\theta=0.5$"),
                       (1.0, ACC[1], r"$\theta=1$  full non-reciprocity")):
        tr = traces[round(th, 4)]
        ax.loglog(tr["ts"], np.maximum(tr["spd"], 1e-12), color=c, label=lab)
    ax.set_xlabel("t [s]"); ax.set_ylabel(r"$\langle|v|\rangle$")
    ax.set_title(r"$\theta=0$ decays without bound;" + "\n"
                 r"$\theta>0$ levels off at a finite speed", fontsize=10)
    ax.legend(fontsize=8)

    ax = fig.add_subplot(gs[0, 1])
    th = [s["theta"] for s in summary]
    ax.errorbar(th, [s["motility"] for s in summary],
                yerr=[s["motility_sd"] for s in summary],
                color=ACC[0], marker="o", ms=5, lw=1.6, capsize=3)
    ax.axhline(0, color="#8b949e", lw=0.6)
    ax.set_xlabel(r"non-reciprocity mixing $\theta$")
    ax.set_ylabel(r"motility  $\langle |V_{\rm cluster}|\rangle/\alpha$  [1/s]")
    ax.set_title("motility order parameter", fontsize=10)

    ax = fig.add_subplot(gs[0, 2])
    ax.plot([s["nu"] for s in summary], [s["p_inj"] for s in summary], "o-",
            color=ACC[3], ms=5, label=r"$P_{\rm inj}/N$ (non-reciprocal drive)")
    ax.plot([s["nu"] for s in summary], [s["p_diss"] for s in summary], "s--",
            color=ACC[2], ms=5, label=r"$\gamma\sum|v|^2/N$ (dissipation)")
    ax.axhline(0, color="#8b949e", lw=0.6)
    ax.set_xlabel(r"non-reciprocity index $\nu$")
    ax.set_ylabel("power per particle")
    ax.set_title("the broken third law pays for the motion", fontsize=10)
    ax.legend(fontsize=7.5)

    for col, th in enumerate((0.0, 0.5, 1.0)):
        ax = fig.add_subplot(gs[1, col])
        p, sim = build(th, 0)
        sim.run(T_EQ + T_MEAS)
        scatter_snapshot(ax, sim)
        o = OB.observables(sim)
        ax.set_title(rf"$\theta$={th:.1f}   $\nu$={M.nonreciprocity(p.A):.2f}   "
                     rf"motility={o['motility']:.4f}" + "\n" + OB.classify(o),
                     fontsize=9)

    fig.suptitle("E03  Reciprocal Particle Life always freezes; life needs a broken third law",
                 fontsize=12)
    fig.tight_layout()
    save(fig, "e03_reciprocity_theorem.png")
    dump(dict(summary=summary, friction=FRICTION,
              traces={str(k): v for k, v in traces.items()}), "e03_reciprocity.json")
