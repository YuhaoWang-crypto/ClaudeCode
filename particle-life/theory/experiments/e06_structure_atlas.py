r"""
E06 - Structure atlas: which motif in A makes which morphology.
===============================================================

Runs every named interaction-matrix generator shipped by the platform under
identical conditions and measures the full observable vector.  The point is to
show that the *graph structure* of ``A`` -- not any hidden rule -- selects the
morphology, and that each morphology is separated by a specific order parameter:

  * ``dim2``      ~1 -> filaments/snakes,  ~2 -> blobs/sheets
  * ``motility``  0 -> frozen,  >0 -> alive  (predicted by ``nu``)
  * ``spin``      0 -> translating,  >0 -> rotating (cyclic non-reciprocity)
  * ``psi6``      crystal vs fluid packing
  * ``length``    2 pi / k_peak, compared with the linear-theory ``l*``
"""

import multiprocessing as mp

import numpy as np
from common import ACC, dump, save, scatter_snapshot
import matplotlib.pyplot as plt

from plife import matrices as M, model, observables as OB, stability as ST

L, N, S = 900.0, 4000, 6
T_EQ, T_MEAS = 35.0, 8.0

NAMES = ["condensate", "symmetric", "random", "gas", "chains1", "chains2", "chains3",
         "snake", "rps", "shells", "bipartite", "hub", "dimers", "triads",
         "checker", "spiral_conveyor", "patchwork", "wavefield",
         "chiral_bandpass", "rotating_conveyor", "parity_vortex"]


def build(name, seed=0):
    A = M.generate(name, S, np.random.default_rng(21))
    alpha, beta = M.random_radii(S, rng=np.random.default_rng(22))
    p = model.Params(A=A, alpha=alpha, beta=beta, L=(L, L))
    return p, model.ParticleLife(p, N=N, seed=seed)


def job(name):
    p, sim = build(name)
    sim.run(T_EQ)
    o1 = OB.observables(sim)
    sim.run(T_MEAS)
    o2 = OB.observables(sim)
    o = {k: (0.5 * (o1[k] + o2[k]) if isinstance(o1[k], float) else o1[k])
         for k in o1}
    rho = np.full(S, N / S / L**2)
    lin = ST.scan(p, rho, nk=400)
    return dict(name=name, display=M.DISPLAY_NAMES[name],
                nu=M.nonreciprocity(p.A),
                motility=o["motility"], speed=o["speed"], spin=o["spin"],
                abs_spin=o["abs_spin"], dim2=o["dim2"], psi6=o["psi6"],
                n_clusters=o["n_clusters"], mean_cluster=o["mean_cluster"],
                length_meas=o["length"], length_pred=lin.length,
                growth=lin.growth, osc_growth=lin.osc_growth,
                regime=lin.regime(), label=OB.classify(o))


if __name__ == "__main__":
    with mp.Pool(min(4, mp.cpu_count())) as pool:
        rows = pool.map(job, NAMES)

    hdr = (f"{'matrix':22s} {'nu':>5s} {'motil':>7s} {'spin':>7s} {'dim2':>5s} "
           f"{'psi6':>5s} {'l_meas':>7s} {'l_pred':>7s} {'linear':>14s}  morphology")
    print(hdr); print("-" * len(hdr))
    for r in sorted(rows, key=lambda x: -x["motility"]):
        print(f"{r['display']:22s} {r['nu']:5.2f} {r['motility']:7.4f} {r['spin']:+7.3f} "
              f"{r['dim2']:5.2f} {r['psi6']:5.2f} {r['length_meas']:7.1f} "
              f"{r['length_pred']:7.1f} {r['regime']:>14s}  {r['label']}")

    # ---------------------------------------------------------------- atlas
    ncol = 7
    nrow = int(np.ceil(len(NAMES) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.05 * ncol, 2.35 * nrow))
    for ax, name in zip(axes.ravel(), NAMES):
        p, sim = build(name)
        sim.run(T_EQ)
        scatter_snapshot(ax, sim, size=0.9, alpha=0.85)
        r = next(x for x in rows if x["name"] == name)
        ax.set_title(f"{r['display']}\n" + r"$\nu$=" + f"{r['nu']:.2f}  "
                     + f"mot={r['motility']:.3f}\n{r['label']}", fontsize=6.4)
    for ax in axes.ravel()[len(NAMES):]:
        ax.axis("off")
    fig.suptitle("E06  Structure atlas - one force law, twenty-one morphologies", fontsize=12)
    fig.tight_layout()
    save(fig, "e06_structure_atlas.png")

    # ------------------------------------------------------- order-param map
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.1))
    ax = axes[0]
    sc = ax.scatter([r["nu"] for r in rows], [r["motility"] for r in rows],
                    c=[r["dim2"] for r in rows], cmap="viridis", s=52,
                    edgecolors="#0d1117", linewidths=0.6)
    for r in rows:
        ax.annotate(r["display"][:13], (r["nu"], r["motility"]), fontsize=5.4,
                    xytext=(3, 3), textcoords="offset points", color="#8b949e")
    ax.set_xlabel(r"non-reciprocity index $\nu$")
    ax.set_ylabel("motility")
    ax.set_title("motility requires non-reciprocity", fontsize=10)
    fig.colorbar(sc, ax=ax, label=r"$D_2$", fraction=0.046)

    ax = axes[1]
    ax.scatter([r["dim2"] for r in rows], [r["psi6"] for r in rows],
               c=[r["motility"] for r in rows], cmap="magma", s=52,
               edgecolors="#0d1117", linewidths=0.6)
    for r in rows:
        ax.annotate(r["display"][:13], (r["dim2"], r["psi6"]), fontsize=5.4,
                    xytext=(3, 3), textcoords="offset points", color="#8b949e")
    ax.axvline(1.55, color="#8b949e", ls="--", lw=0.8)
    ax.set_xlabel(r"correlation dimension $D_2$")
    ax.set_ylabel(r"hexatic order $|\psi_6|$")
    ax.set_title("filaments (left) vs blobs (right)", fontsize=10)

    ax = axes[2]
    lm = [r["length_meas"] for r in rows]
    lp = [r["length_pred"] for r in rows]
    ax.scatter(lp, lm, s=52, color=ACC[0], edgecolors="#0d1117", linewidths=0.6)
    hi = max(max(lm), max(lp)) * 1.1
    ax.plot([0, hi], [0, hi], color="#8b949e", ls="--", lw=1)
    ax.set_xlabel(r"linear theory  $\ell^*=2\pi/k^*$")
    ax.set_ylabel(r"measured  $2\pi/k_{\rm peak}$")
    ax.set_title("predicted vs observed structure size", fontsize=10)

    fig.suptitle("E06b  Order parameters separate the morphologies", fontsize=12)
    fig.tight_layout()
    save(fig, "e06_order_parameters.png")
    dump(dict(rows=rows), "e06_atlas.json")
