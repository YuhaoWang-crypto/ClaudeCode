"""Mechanistic delivery-kinetics analysis (PLAN.md Phase 4).

Uses lipidlib/kinetics.py to (1) show the LNP->protein cascade in time, (2) locate
the rate-limiting step (endosomal escape), and (3) connect Track A: map a lipid's
LiON predicted-delivery to an effective escape rate and read off expression
DYNAMICS (peak height, time-to-peak, duration) that the ML potency score alone
does not give.

Figure: results/figures/delivery_kinetics.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from lipidlib.kinetics import DeliveryParams, simulate, metrics
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE, ORANGE, GREY, GREEN = "#0072B2", "#E69F00", "#8A8A8A", "#009E73"
plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 130})


def main() -> None:
    base = DeliveryParams()
    print("baseline:", {k: round(v, 3) for k, v in metrics(base).items()})

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))

    # --- A: the cascade in time (normalised) ---
    df = simulate(base, t_end=96)
    for col, c, lab in [("L_endo", GREY, "endosomal LNP"), ("M", ORANGE, "cytosolic mRNA"),
                        ("P", BLUE, "protein")]:
        ax[0].plot(df["t"], df[col] / df[col].max(), color=c, lw=2, label=lab)
    m = metrics(base)
    ax[0].axvline(m["t_peak_h"], color=BLUE, ls="--", lw=1)
    ax[0].set_xlabel("time (h)"); ax[0].set_ylabel("normalised amount")
    ax[0].set_title(f"A  Delivery cascade\nprotein peaks ~{m['t_peak_h']:.0f} h, "
                    f"escape eff {m['escape_efficiency']*100:.1f}%")
    ax[0].legend(frameon=False, fontsize=8.5)

    # --- B: escape is rate-limiting — sweep each step, see AUC response ---
    factors = np.logspace(-1, 1, 25)  # 0.1x .. 10x
    sweeps = {
        "endosomal escape": ("k_escape", GREEN),
        "mRNA stability (1/k_mdeg)": ("k_mdeg", ORANGE),
        "uptake": ("k_uptake", GREY),
        "translation": ("k_transl", BLUE),
    }
    base_auc = metrics(base)["auc_protein_analytic"]
    for lab, (attr, c) in sweeps.items():
        aucs = []
        for f in factors:
            p = DeliveryParams()
            val = getattr(p, attr) * (1 / f if attr == "k_mdeg" else f)
            setattr(p, attr, val)
            aucs.append(metrics(p)["auc_protein_analytic"] / base_auc)
        ax[1].plot(factors, aucs, color=c, lw=2, label=lab)
    ax[1].set_xscale("log"); ax[1].set_yscale("log")
    ax[1].set_xlabel("parameter fold-change"); ax[1].set_ylabel("relative protein AUC")
    ax[1].set_title("B  Expression scales with escape,\nmRNA-stability, translation — not uptake")
    ax[1].legend(frameon=False, fontsize=8)

    # --- C: map Track-A LiON potency -> escape rate -> expression dynamics ---
    ens = ROOT / "data" / "libraries" / "combo_v1" / "combo_v1_ranked_ensemble.csv"
    curves = []
    if ens.exists():
        top = pd.read_csv(ens).sort_values("ens_mean", ascending=False)
        # pick a spread of potencies: best, median, low
        picks = [top.iloc[0], top.iloc[len(top) // 2], top.iloc[-1]]
        labs = ["top lead", "median", "weak"]
        cols = [GREEN, ORANGE, GREY]
        for row, lab, c in zip(picks, labs, cols):
            p = DeliveryParams()
            # potency -> escape rate (log-linear), capped at a physical max (~8% escape)
            p.k_escape = min(0.02 * (10 ** float(row["ens_mean"])), 0.10)
            d = simulate(p, t_end=96)
            mm = metrics(p)
            ax[2].plot(d["t"], d["P"], color=c, lw=2,
                       label=f"{lab} (z={row['ens_mean']:.2f}): peak {mm['t_peak_h']:.0f} h")
            curves.append({"lipid": f"{row['head']}+{row['linker']}+{row['tail']}",
                           "ens_mean": row["ens_mean"], **mm})
        ax[2].set_xlabel("time (h)"); ax[2].set_ylabel("protein (a.u.)")
        ax[2].set_title("C  Track-A potency → expression dynamics")
        ax[2].legend(frameon=False, fontsize=8)
        pd.DataFrame(curves).to_csv(
            ROOT / "data" / "libraries" / "combo_v1" / "combo_v1_kinetics.csv", index=False)

    fig.tight_layout()
    out = ROOT / "results" / "figures" / "delivery_kinetics.png"
    fig.savefig(out, bbox_inches="tight")
    print(f"wrote {out}")
    if curves:
        print("\nTrack-A leads -> mechanistic expression:")
        for c in curves:
            print(f"  {c['lipid']:28s} z={c['ens_mean']:+.2f}  escape={c['escape_efficiency']*100:4.1f}%  "
                  f"peak@{c['t_peak_h']:.0f}h  AUC={c['auc_protein_analytic']:.0f}")


if __name__ == "__main__":
    main()
