r"""
E07 - The morphogenesis map.
============================

One two-parameter family of interaction matrices reproduces essentially the whole
Particle-Life zoo.  On the cycle of species ``0,1,...,S-1`` set

.. math::  A_{ij} = a\,\delta_{ij} \;+\; w\,\delta_{j,i+1} \;+\; b\,\delta_{j,i-1}
           \;+\; c\,(\text{everything else}),

so ``w`` is the *forward* coupling along the type cycle and ``b`` the *backward*
one.  The two coordinates have direct physical meaning:

    w + b   -> reciprocal bonding strength   (how strongly the chain holds together)
    w - b   -> non-reciprocity               (how hard each link chases the next)

and the diagonal of the map ``w = b`` is exactly the reciprocal line where the
theorem of E03 forbids motion.  Sweeping the plane and classifying each run with
the order parameters of ``observables.py`` produces a phase diagram in which
blobs, static chains, crawling snakes, rotors and gas all appear as distinct
regions -- with the reciprocal line acting as the boundary of the "alive" phase.
"""

import multiprocessing as mp

import numpy as np
from common import ACC, dump, save, scatter_snapshot
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

from plife import matrices as M, model, observables as OB, stability as ST

L, N, S = 800.0, 2400, 5
T_EQ, T_MEAS = 26.0, 6.0
A_SELF, A_OTHER = 1.0, -0.35
GRID = np.linspace(-1.0, 1.0, 9)

LABELS = ["gas", "static blobs", "crystalline clusters", "static chains / network",
          "migrating cells / chasers", "rotating cells", "crawling filaments"]
LABEL_IDX = {s: i for i, s in enumerate(LABELS)}
CMAP = ListedColormap(["#21262d", "#2d5a3d", "#39d353", "#58a6ff",
                       "#f0883e", "#a371f7", "#ff5ec4"])


def matrix(w, b):
    A = np.full((S, S), A_OTHER)
    for i in range(S):
        A[i, i] = A_SELF
        A[i, (i + 1) % S] = w
        A[i, (i - 1) % S] = b
    return A


def build(w, b, seed=0):
    A = matrix(w, b)
    alpha, beta = M.uniform_radii(S, 14.0, 42.0)
    p = model.Params(A=A, alpha=alpha, beta=beta, L=(L, L), friction=0.12)
    return p, model.ParticleLife(p, N=N, seed=seed)


def job(args):
    w, b = args
    p, sim = build(w, b)
    sim.run(T_EQ)
    o1 = OB.observables(sim)
    sim.run(T_MEAS)
    o2 = OB.observables(sim)
    o = {k: (0.5 * (o1[k] + o2[k]) if isinstance(o1[k], float) else o1[k]) for k in o1}
    lab = OB.classify(o)
    rho = np.full(S, N / S / L**2)
    lin = ST.scan(p, rho, nk=300)
    return dict(w=float(w), b=float(b), nu=M.nonreciprocity(p.A),
                motility=o["motility"], spin=o["spin"], abs_spin=o["abs_spin"],
                dim2=o["dim2"], psi6=o["psi6"], mean_cluster=o["mean_cluster"],
                length=o["length"], label=lab, idx=LABEL_IDX.get(lab, 0),
                osc_growth=lin.osc_growth, growth=lin.growth)


if __name__ == "__main__":
    jobs = [(w, b) for b in GRID for w in GRID]
    with mp.Pool(min(4, mp.cpu_count())) as pool:
        rows = pool.map(job, jobs)

    n = len(GRID)
    def field(key):
        Z = np.zeros((n, n))
        for r in rows:
            i = int(np.argmin(np.abs(GRID - r["b"])))
            j = int(np.argmin(np.abs(GRID - r["w"])))
            Z[i, j] = r[key]
        return Z

    idx, mot, spin, d2 = field("idx"), field("motility"), field("abs_spin"), field("dim2")

    print(f"{'w':>6s} {'b':>6s} {'nu':>5s} {'motil':>8s} {'|spin|':>7s} "
          f"{'dim2':>5s}  morphology")
    for r in rows:
        print(f"{r['w']:6.2f} {r['b']:6.2f} {r['nu']:5.2f} {r['motility']:8.4f} "
              f"{r['abs_spin']:7.3f} {r['dim2']:5.2f}  {r['label']}")

    # ------------------------------------------------------------------ figure
    fig = plt.figure(figsize=(13.6, 8.2))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.25, 1])

    ax = fig.add_subplot(gs[0, 0])
    im = ax.pcolormesh(GRID, GRID, idx, cmap=CMAP,
                       norm=BoundaryNorm(np.arange(-0.5, len(LABELS)), CMAP.N),
                       shading="auto")
    ax.plot([-1, 1], [-1, 1], color="w", lw=1.6, ls="--")
    ax.text(-0.92, -0.72, "reciprocal line\n(theorem: frozen)", color="w", fontsize=7)
    ax.set_xlabel("forward coupling w"); ax.set_ylabel("backward coupling b")
    ax.set_title("morphology phase diagram", fontsize=10)
    cb = fig.colorbar(im, ax=ax, ticks=range(len(LABELS)), fraction=0.046)
    cb.ax.set_yticklabels(LABELS, fontsize=6.5)

    for k, (Z, name, cmap) in enumerate([
            (mot, "motility", "magma"), (spin, r"|spin| [1/s]", "twilight"),
            (d2, r"$D_2$", "viridis")]):
        ax = fig.add_subplot(gs[0, k + 1])
        im = ax.pcolormesh(GRID, GRID, Z, cmap=cmap, shading="auto")
        ax.plot([-1, 1], [-1, 1], color="w", lw=1.4, ls="--")
        ax.set_xlabel("w"); ax.set_ylabel("b")
        ax.set_title(name, fontsize=10)
        fig.colorbar(im, ax=ax, fraction=0.046)

    picks = [(0.0, 0.0), (0.8, 0.8), (0.9, 0.0), (0.9, -0.9)]
    for col, (w, b) in enumerate(picks):
        ax = fig.add_subplot(gs[1, col])
        p, sim = build(w, b)
        sim.run(T_EQ + T_MEAS)
        scatter_snapshot(ax, sim, size=1.1)
        r = min(rows, key=lambda x: (x["w"] - w) ** 2 + (x["b"] - b) ** 2)
        ax.set_title(f"w={w:+.1f}, b={b:+.1f}   " + r"$\nu$=" + f"{r['nu']:.2f}\n"
                     f"{r['label']}  (mot={r['motility']:.3f})", fontsize=8)

    fig.suptitle("E07  Morphogenesis map: the forward/backward coupling plane", fontsize=12)
    fig.tight_layout()
    save(fig, "e07_morphogenesis_map.png")
    dump(dict(grid=GRID.tolist(), rows=rows, labels=LABELS), "e07_map.json")
