"""Figure generation for the conformational-flux pipeline."""
import os

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .common import (APO_STARTS, DT_NS, FIG, VARIANTS, amb_to_idx,
                     apo_indices, crystal_ca, idx_to_amb, kabsch_rmsd, load)
from .m2_cryptic_pocket import (OPEN_THR, POCKET_RES, collective_variables,
                                committed_crossings, opening_axis)

C = {"TEM-1": "#4a6fa5", "TEM-2": "#6aa84f", "TEM-10": "#c0392b", "TEM-52": "#e8a33d"}


def fig_overview():
    X, names = load()
    proj, cv1, _ = collective_variables()
    _, _, L = opening_axis()
    fig, ax = plt.subplots(2, 3, figsize=(16.5, 9))

    # (a) free-energy-like profile along the crystal axis, per variant
    a = ax[0, 0]
    bins = np.linspace(-3, 16, 90)
    for v in VARIANTS:
        p = proj[apo_indices(names, v)].ravel()
        h, _ = np.histogram(p, bins=bins, density=True)
        m = h > 0
        a.plot(bins[:-1][m], -np.log(h[m]), color=C[v], lw=1.8, label=v)
    a.axvline(L, color="k", ls="--", lw=1)
    a.text(L, a.get_ylim()[1] * 0.95, " 1PZO", fontsize=8, va="top")
    a.set_xlabel("position along crystal opening axis (Å)")
    a.set_ylabel("−ln P  (kT)")
    a.set_title("(a) apo landscape along the 1BTL→1PZO axis")
    a.legend(fontsize=8)

    # (b) occupancy per variant, paired by crystal start
    a = ax[0, 1]
    idx = {n: k for k, n in enumerate(names)}
    w = 0.2
    for i, v in enumerate(VARIANTS):
        y = [100 * (cv1[idx[f"{s}_{v}"]] > OPEN_THR).mean() for s in APO_STARTS]
        a.bar(np.arange(8) + (i - 1.5) * w, y, w, color=C[v], label=v)
    a.set_yscale("symlog", linthresh=1e-3)
    a.set_xticks(range(8))
    a.set_xticklabels(APO_STARTS, rotation=45, fontsize=8)
    a.set_ylabel(f"frames with mouth > {OPEN_THR} Å  (%)")
    a.set_title("(b) open-pocket occupancy, paired by crystal start")
    a.legend(fontsize=8)

    # (c) occupancy vs crossing count
    a = ax[0, 2]
    for v in VARIANTS:
        ks = apo_indices(names, v)
        x = [committed_crossings(cv1[k])[0] for k in ks]
        y = [100 * (cv1[k] > OPEN_THR).mean() for k in ks]
        a.scatter(x, y, s=46, color=C[v], label=v, edgecolor="k", lw=.5, zorder=3)
    a.set_xlabel("committed closed→open crossings per 500 ns run")
    a.set_ylabel("open occupancy (%)")
    a.set_yscale("symlog", linthresh=1e-3)
    a.set_title("(c) two independent observables agree")
    a.grid(alpha=.3)
    a.legend(fontsize=8)

    # (d) example traces
    a = ax[1, 0]
    t = np.arange(proj.shape[1]) * DT_NS
    for nm, col in [("3JYI_TEM-10", C["TEM-10"]), ("3JYI_TEM-1", C["TEM-1"])]:
        a.plot(t, cv1[idx[nm]], color=col, lw=.5, label=nm)
    a.axhline(OPEN_THR, color="k", ls="--", lw=1)
    a.set_xlabel("time (ns)")
    a.set_ylabel("L220–N276 mouth distance (Å)")
    a.set_title("(d) same crystal start, two variants")
    a.legend(fontsize=8)

    # (e) mechanism: per-residue displacement open vs closed
    a = ax[1, 1]
    sub = np.asarray(X[apo_indices(names, "TEM-10")])
    pv = proj[apo_indices(names, "TEM-10")]
    op, cl = sub[pv > 6.5], sub[pv < 1.5]
    d = np.linalg.norm(op.mean(0) - cl.mean(0), axis=-1)
    res = [idx_to_amb(i) for i in range(len(d))]
    a.fill_between(res, 0, d, color="#888", lw=0)
    for lo, hi, lab, col in [(214, 232, "H10–H11 / H11", "#c0392b"),
                             (269, 287, "H12", "#2e86c1"),
                             (161, 180, "Ω-loop", "#6aa84f")]:
        m = [(r >= lo) and (r <= hi) for r in res]
        a.fill_between(np.array(res)[m], 0, d[m], color=col, label=lab, lw=0)
    a.set_xlabel("residue (Ambler)")
    a.set_ylabel("Cα displacement, open − closed (Å)")
    a.set_title("(e) what actually moves")
    a.legend(fontsize=8)

    # (f) does the open state resemble the bound crystal?
    a = ax[1, 2]
    pr = [amb_to_idx(x_) for x_ in POCKET_RES]
    ro_ref = np.array([crystal_ca("1PZO")(x_) for x_ in POCKET_RES])
    rc_ref = np.array([crystal_ca("1BTL")(x_) for x_ in POCKET_RES])
    xs, ys, cs = [], [], []
    for k in apo_indices(names):
        for f in range(0, proj.shape[1], 300):
            fr = np.asarray(X[k][f])[pr]
            xs.append(kabsch_rmsd(fr, rc_ref))
            ys.append(kabsch_rmsd(fr, ro_ref))
            cs.append(proj[k][f])
    sc = a.scatter(xs, ys, c=cs, s=5, cmap="viridis", vmin=0, vmax=12)
    lim = [0.5, max(max(xs), max(ys)) * 1.02]
    a.plot(lim, lim, "k--", lw=1)
    a.set_xlim(lim)
    a.set_ylim(lim)
    a.set_xlabel("pocket RMSD to 1BTL (closed apo), Å")
    a.set_ylabel("pocket RMSD to 1PZO (inhibitor-bound), Å")
    a.set_title("(f) below the line = bound-like")
    plt.colorbar(sc, ax=a, label="axis position (Å)")

    fig.tight_layout()
    p = os.path.join(FIG, "msm_tem_overview.png")
    fig.savefig(p, dpi=135)
    plt.close(fig)
    return p


def fig_msm():
    from deeptime.markov import TransitionCountEstimator, reactive_flux
    from deeptime.markov.msm import MaximumLikelihoodMSM
    from .m3_msm_diagnostics import LAGS_FRAMES, _featurize, _its
    from .m4_tpt_flux import A_MAX, B_MIN, build

    X, names = load()
    apo = apo_indices(names)
    fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.6))

    # (a) ITS non-convergence
    from deeptime.clustering import KMeans
    from deeptime.decomposition import TICA
    a = ax[0]
    for tag, residues, cutoff, dim, col in [
        ("whole protein", [r for r in range(26, 291) if r not in (239, 253)],
         10.0, 4, "#c0392b"),
        ("pocket-local", list(range(214, 232)) + list(range(269, 287)),
         18.0, 2, "#2e86c1")]:
        I, J = _featurize(X, residues, cutoff)
        F = [np.linalg.norm(np.asarray(X[k])[:, I] - np.asarray(X[k])[:, J],
                            axis=-1).astype(np.float32) for k in apo]
        tica = TICA(lagtime=100, dim=max(dim, 3), scaling="kinetic_map")
        for f in F:
            tica.partial_fit(f)
        Y = np.stack([tica.fetch_model().transform(f) for f in F])
        km = (KMeans(n_clusters=100, max_iter=300, init_strategy="kmeans++",
                     n_jobs=4, fixed_seed=1)
              .fit(Y[:, ::3, :dim].reshape(-1, dim)).fetch_model())
        dt = [km.transform(Y[k, :, :dim]) for k in range(len(apo))]
        out = _its(dt, LAGS_FRAMES, n=3)
        lags = np.array([l * DT_NS for l, _, _ in out])
        for j in range(3):
            a.plot(lags, [ts[j] for _, _, ts in out], "o-", color=col, lw=1.4,
                   ms=4, alpha=1 - .25 * j, label=tag if j == 0 else None)
    a.plot([0.5, 25], [0.5, 25], "k--", lw=1, label="t = lag (unresolvable)")
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_xlabel("lag time τ (ns)")
    a.set_ylabel("implied timescale (ns)")
    a.set_title("(a) ITS test: nothing flattens → not Markovian")
    a.legend(fontsize=8)

    msm, cnt, cp = build()
    fl = reactive_flux(msm.transition_matrix,
                       np.where(cp < A_MAX)[0], np.where(cp > B_MIN)[0])

    # (b) committor along the physical coordinate
    a = ax[1]
    a.scatter(cp, fl.forward_committor, c=msm.stationary_distribution,
              cmap="magma_r", s=46, edgecolor="k", lw=.4,
              norm=matplotlib.colors.LogNorm())
    a.axhline(.5, color="k", ls="--", lw=1)
    a.set_xlabel("position along opening axis (Å)")
    a.set_ylabel("forward committor q⁺")
    a.set_title("(b) commitment happens late")
    plt.colorbar(a.collections[0], ax=a, label="stationary population")

    # (c) flux is diffuse
    a = ax[2]
    paths, caps = fl.pathways(fraction=0.9)
    frac = np.array(caps) / fl.total_flux
    a.plot(np.arange(1, len(frac) + 1), 100 * np.cumsum(frac), lw=2,
           color="#2e86c1")
    a.set_xscale("log")
    a.set_xlabel("number of pathways (ranked)")
    a.set_ylabel("cumulative share of flux (%)")
    a.set_title(f"(c) {len(paths)} routes carry 90 %; top route {100*frac[0]:.2f} %")
    a.grid(alpha=.3)

    fig.tight_layout()
    p = os.path.join(FIG, "msm_tem_kinetics.png")
    fig.savefig(p, dpi=135)
    plt.close(fig)
    return p


def main():
    os.makedirs(FIG, exist_ok=True)
    for f in (fig_overview, fig_msm):
        print("  wrote", os.path.relpath(f(), os.path.dirname(FIG)))


if __name__ == "__main__":
    main()
