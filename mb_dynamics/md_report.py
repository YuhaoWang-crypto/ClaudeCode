"""MD / R2 diagnostics and the pre-declared linker-geometry sensitivity.

    python md_report.py 0.5
Writes results/md_summary_dt{dt}.json and figures/md_diagnostics_dt{dt}.png.
"""
import glob
import json
import os
import sys

import numpy as np

import md_analysis as ma
import swv_backend as sb
from run_analysis import D_SE_NM2_S

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
FIG = os.path.join(HERE, "figures")


def variant_summary(names, beta=5.0):
    z = np.concatenate([ma.load(n)["z"].ravel() for n in names])
    m, se, _ = ma.pooled(names, beta)
    kin = sb.smoluchowski_kinetics(z, 40, D_SE_NM2_S, beta, ma.Z_MIN)
    lam = np.sort(np.linalg.eigvalsh(kin.S))[::-1]
    return dict(runs=names, obs=m, se=se, r1_slowest_relaxation_s=float(-1 / lam[1])), z


def main(dt):
    out = {}
    zs = {}
    out["default"], zs["default (0.9 nm linker, l_p 1 nm)"] = variant_summary([f"dt{dt}_s{s}" for s in (101, 102, 103, 104)])
    for tag, label in (("sens_link2", "linker 0.6 nm"), ("sens_link6", "linker 1.8 nm"),
                       ("sens_lp0.7", "linker l_p 0.7 nm"), ("sens_lp2.0", "linker l_p 2 nm")):
        names = sorted(os.path.basename(p)[:-4] for p in glob.glob(os.path.join(ma.MD_DIR, f"{tag}_dt{dt}.npz")))
        if names:
            out[tag], zs[label] = variant_summary(names)
    gates = json.load(open(os.path.join(RES, f"r2_gates_dt{dt}.json")))
    json.dump(out, open(os.path.join(RES, f"md_summary_dt{dt}.json"), "w"), indent=1, default=float)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    edges = np.arange(0.2, 11.8, 0.1)
    cen = 0.5 * (edges[1:] + edges[:-1])
    for label, z in zs.items():
        P, _ = np.histogram(z, edges, density=True)
        with np.errstate(divide="ignore"):
            ax[0].plot(cen, -np.log(P), label=label)
    ax[0].axvline(ma.Z_MIN, color="k", ls=":", lw=0.8)
    ax[0].set_xlabel("MB height z (nm)")
    ax[0].set_ylabel("PMF (k_BT, arbitrary zero)")
    ax[0].set_title("MB height free energy (CG model)")
    ax[0].legend(fontsize=7)
    for cand in gates["candidates"]:
        lags = [r["lag_ps"] / 1000 for r in cand["lags"]]
        ax[1].plot(lags, [r["timescales_ps"]["1"][0] / 1000 for r in cand["lags"]], "o-", label=f"k={cand['k']} slowest")
        ax[2].plot(lags, [max(c["max_dev"] for c in r["ck"]) for r in cand["lags"]], "o-", label=f"k={cand['k']}")
    ax[1].plot([0.01, 10], [0.01, 10], "k:", lw=0.8)
    ax[1].set_xscale("log")
    ax[1].set_xlabel("lag τ (ns)")
    ax[1].set_ylabel("slowest implied timescale (ns)")
    ax[1].set_title("R2 implied timescales (build seeds)")
    ax[1].legend(fontsize=7)
    ax[2].axhline(0.10, color="k", ls="--", lw=0.8)
    ax[2].set_xscale("log")
    ax[2].set_xlabel("lag τ (ns)")
    ax[2].set_ylabel("max |P_pred − P_emp| (2τ, 4τ, 8τ)")
    ax[2].set_title("R2 CK test on held-out seed 103")
    ax[2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, f"md_diagnostics_dt{dt}.png"), dpi=120)
    for k, v in out.items():
        print(k, {kk: round(vv, 4) for kk, vv in v["obs"].items()}, "R1 relax %.3g s" % v["r1_slowest_relaxation_s"])


if __name__ == "__main__":
    main(float(sys.argv[1]) if len(sys.argv) > 1 else 0.5)
