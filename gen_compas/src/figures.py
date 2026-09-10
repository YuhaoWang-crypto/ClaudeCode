"""Figures for the Gen-COMPAS reproduction."""

import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402

plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 150, "font.size": 8,
    "axes.linewidth": 0.7, "axes.titlesize": 9, "figure.facecolor": "white",
})
CMAP = "viridis"


def _fel_panel(ax, f, edges, title, vmax=10.0):
    ext = [edges[0], edges[-1], edges[0], edges[-1]]
    m = np.ma.masked_invalid(f.T)
    im = ax.imshow(m, origin="lower", extent=ext, cmap=CMAP, vmin=0, vmax=vmax,
                   aspect="equal", interpolation="bilinear")
    ax.contour(0.5 * (edges[1:] + edges[:-1]), 0.5 * (edges[1:] + edges[:-1]),
               m, levels=np.arange(1, vmax, 2), colors="w", linewidths=0.4,
               alpha=0.5)
    ax.set_xlabel(r"$\phi$ (deg)")
    ax.set_ylabel(r"$\psi$ (deg)")
    ax.set_title(title)
    ax.set_xticks([-180, -90, 0, 90, 180])
    ax.set_yticks([-180, -90, 0, 90, 180])
    return im


def fig_fel(tag):
    d = np.load(os.path.join(common.RESULTS, f"analysis_{tag}.npz"))
    with open(os.path.join(common.RESULTS, f"analysis_{tag}.json")) as fh:
        j = json.load(fh)
    fel, edges, ref = d["fel"], d["edges"], d["ref"]
    has_ref = ref.size > 0

    n = 2 if has_ref else 1
    fig, axes = plt.subplots(1, n + 1, figsize=(3.9 * (n + 1), 3.2))
    axes = np.atleast_1d(axes)
    im = _fel_panel(axes[0], fel, edges,
                    "Gen-COMPAS free energy\n(from unbiased sampling only)")
    fig.colorbar(im, ax=axes[0], fraction=0.046, pad=0.03, label="kcal/mol")
    if has_ref:
        im2 = _fel_panel(axes[1], ref, edges,
                         "Reference: well-tempered\nmetadynamics on "
                         r"$(\phi,\psi)$")
        fig.colorbar(im2, ax=axes[1], fraction=0.046, pad=0.03,
                     label="kcal/mol")

    ax = axes[-1]
    if has_ref:
        both = np.isfinite(fel) & np.isfinite(ref) & (ref < 8.0)
        ax.scatter(ref[both], fel[both], s=4, alpha=0.5, color="#2b6cb0")
        lim = [0, 8]
        ax.plot(lim, lim, "k--", lw=0.7)
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_xlabel("metadynamics (kcal/mol)")
        ax.set_ylabel("Gen-COMPAS (kcal/mol)")
        ax.set_title(f"bin-by-bin, RMSE = {j.get('fel_rmse', float('nan')):.2f}")
        ax.set_aspect("equal")
    fig.tight_layout()
    p = os.path.join(common.FIGURES, f"fel_{tag}.png")
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def fig_committor(tag):
    d = np.load(os.path.join(common.RESULTS, f"analysis_{tag}.npz"))
    with open(os.path.join(common.RESULTS, f"analysis_{tag}.json")) as fh:
        j = json.load(fh)
    phi, psi, q, w = d["phi"], d["psi"], d["q"], d["weights"]
    edges = np.linspace(-180, 180, 37)

    num, _, _ = np.histogram2d(phi, psi, bins=[edges, edges], weights=q)
    cnt, _, _ = np.histogram2d(phi, psi, bins=[edges, edges])
    with np.errstate(invalid="ignore"):
        qmap = np.where(cnt >= 3, num / np.maximum(cnt, 1), np.nan)

    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.1))
    ax = axes[0]
    ext = [-180, 180, -180, 180]
    im = ax.imshow(np.ma.masked_invalid(qmap.T), origin="lower", extent=ext,
                   cmap="coolwarm", vmin=0, vmax=1, aspect="equal",
                   interpolation="nearest")
    c = 0.5 * (edges[1:] + edges[:-1])
    ax.contour(c, c, np.ma.masked_invalid(qmap.T), levels=[0.5],
               colors="k", linewidths=1.2)
    for p in j["paths"]:
        pt = p["points"]
        ax.plot([x["phi"] for x in pt], [x["psi"] for x in pt], "o-",
                ms=2.5, lw=1.0, label=f"channel {p['channel']}")
    ax.legend(fontsize=6, loc="lower left", framealpha=0.85)
    ax.set_xlabel(r"$\phi$ (deg)")
    ax.set_ylabel(r"$\psi$ (deg)")
    ax.set_title("Committor learned in conformational space\n"
                 r"(projected on $\phi,\psi$; black = separatrix)")
    fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02, label="q")

    ax = axes[1]
    tse = d["tse"]
    ax.scatter(phi[::3], psi[::3], s=1.5, c="#cbd5e0", label="all sampling")
    ax.scatter(phi[tse], psi[tse], s=5, c="#c53030",
               label=r"TSE, $|q-1/2|<0.05$")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$\phi$ (deg)")
    ax.set_ylabel(r"$\psi$ (deg)")
    ax.set_title("Transition-state ensemble")
    ax.legend(fontsize=6, loc="lower left", framealpha=0.85)
    fig.tight_layout()
    p = os.path.join(common.FIGURES, f"committor_{tag}.png")
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def fig_validation(tag):
    path = os.path.join(common.RESULTS, f"validation_{tag}.json")
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        v = json.load(fh)
    qp = np.array([r["q_pred"] for r in v["rows"]])
    qe = np.array([r["q_emp"] for r in v["rows"]])
    n = np.array([r["n"] for r in v["rows"]])

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0))
    ax = axes[0]
    err = np.sqrt(qe * (1 - qe) / np.maximum(n, 1))
    ax.errorbar(qp, qe, yerr=err, fmt="o", ms=3.5, lw=0.7, capsize=1.5,
                color="#2b6cb0")
    ax.plot([0, 1], [0, 1], "k--", lw=0.8)
    ax.set_xlabel("predicted committor")
    ax.set_ylabel("measured committor (fresh shooting)")
    ax.set_title(f"MAE = {v['mae']:.3f} (noise floor {v['expected_mae']:.3f})")
    ax.set_aspect("equal")
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)

    ax = axes[1]
    sel = np.abs(qp - 0.5) < 0.1
    ax.hist(qe[sel], bins=np.linspace(0, 1, 11), color="#c53030", alpha=0.8,
            edgecolor="w")
    ax.axvline(0.5, color="k", ls="--", lw=0.8)
    ax.set_xlabel("measured committor")
    ax.set_ylabel("count")
    ax.set_title(r"Predicted TSE ($|q-1/2|<0.1$): "
                 f"{sel.sum()} structures")
    fig.tight_layout()
    p = os.path.join(common.FIGURES, f"validation_{tag}.png")
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def fig_convergence(tag):
    with open(os.path.join(common.RESULTS, f"gencompas_{tag}_log.json")) as fh:
        L = json.load(fh)
    log = L["log"]
    d = np.load(os.path.join(common.RESULTS, f"analysis_{tag}.npz"))
    it = [r["iteration"] for r in log]

    fig, axes = plt.subplots(1, 3, figsize=(9.6, 2.8))
    ax = axes[0]
    ax.plot(it, [r["ns_cumulative"] for r in log], "o-", ms=3, color="#2b6cb0")
    ax.set_xlabel("iteration")
    ax.set_ylabel("cumulative MD (ns)")
    ax.set_title("Simulation cost")

    ax = axes[1]
    ax.plot(it, [100 * r["frac_separatrix"] for r in log], "o-", ms=3,
            color="#2f855a", label="on separatrix")
    if "committor_mae" in log[0]:
        ax2 = ax.twinx()
        ax2.plot(it, [r["committor_mae"] for r in log], "s--", ms=3,
                 color="#c05621", label="committor error")
        ax2.set_ylabel(r"$|q_{\rm pred}-q_{\rm measured}|$", color="#c05621")
        ax2.set_ylim(0, 0.55)
        ax2.tick_params(axis="y", colors="#c05621")
    ax.set_xlabel("iteration")
    ax.set_ylabel("% of shooting points with 0.2 < q < 0.8", color="#2f855a")
    ax.tick_params(axis="y", colors="#2f855a")
    ax.set_ylim(0, 100)
    ax.set_title("Hits on the separatrix")

    ax = axes[2]
    lags = d["lags"] * 0.2
    its = d["its"] * 0.2
    for k in range(its.shape[1]):
        ax.plot(lags, its[:, k], "o-", ms=3, label=f"$t_{k + 1}$")
    ax.plot(lags, lags, "k--", lw=0.8)
    ax.set_yscale("log")
    ax.set_xlabel("lag time (ps)")
    ax.set_ylabel("implied timescale (ps)")
    ax.set_title("Markov model validation")
    ax.legend(fontsize=6)
    fig.tight_layout()
    p = os.path.join(common.FIGURES, f"convergence_{tag}.png")
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def main():
    tag = sys.argv[1] if len(sys.argv) > 1 else "main"
    for fn in (fig_fel, fig_committor, fig_validation, fig_convergence):
        try:
            p = fn(tag)
            print("wrote", p)
        except Exception as e:  # keep going; a missing input is not fatal
            print(f"{fn.__name__}: {type(e).__name__}: {e}")


if __name__ == "__main__":
    main()
