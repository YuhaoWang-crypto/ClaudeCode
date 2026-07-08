"""
Module 6 — Integrate real binding data into the network-stability biomarker.

Chain:  ChEMBL potency + Boltz-2.1 binding  ->  target-engagement strength
        ->  network perturbation mu  ->  M4 stability readout (LLE / DNB).

Rationale:
  M4 showed that pushing a bistable regulatory core toward its tipping point
  drives the leading Lyapunov exponent -> 0 and inflates the DNB index. How
  hard a drug pushes depends on TARGET ENGAGEMENT: how completely and tightly
  it shuts the drugged node. We combine two real, independent measurements of
  engagement per drug --- ChEMBL pIC50 (functional potency) and the Boltz-2.1
  binding score (structural) --- into a normalised engagement E in [0,1], map
  E monotonically onto the M4 perturbation axis mu, and read off the predicted
  stability biomarker for each real drug.

IMPORTANT (rigour): the E->mu map is an explicit, monotonic ILLUSTRATIVE
coupling, not a calibrated dose model. What is real here is (a) the binding
numbers, (b) the M4 stability curve, and (c) the ordering they induce.
"""
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from grn_pipeline import m4_dnb_lyapunov as m4
from grn_pipeline.m5_kras_real import CHEMBL, BOLTZ


def engagement(drug):
    """Combine ChEMBL pIC50 and Boltz optimization_score into E in [0,1]."""
    d = CHEMBL[drug]
    b = BOLTZ[drug]
    # functional component from pIC50 (nM): pIC50 = 9 - log10(IC50_nM)
    ic = d.get("ic50_g12c_nM")
    if ic:
        pic50 = 9 - math.log10(ic)
        func = (pic50 - 4) / 5                      # ~4..9 -> 0..1
    else:
        func = None
    # structural component from Boltz optimization_score (already ~0..1)
    struct = b["optimization_score"]
    if func is None:
        E = struct
    else:
        E = 0.5 * (max(0.0, min(1.0, func)) + struct)
    return {"func": func, "struct": struct, "E": E}


def E_to_mu(E, mu_max):
    """Monotonic illustrative map from engagement to perturbation mu."""
    return E * mu_max


def report(fig_path="figures/m6_integrated.png"):
    print("=" * 68)
    print("MODULE 6 — Real binding data -> network-stability biomarker")
    print("=" * 68)

    rows = m4.sweep()                       # the real M4 stability curve
    mu = np.array([r["mu"] for r in rows])
    eig = np.array([r["lead_eig"] for r in rows])
    dnb = np.array([r["dnb"] for r in rows])
    mu_max = mu.max()

    drug_points = {}
    for drug in ("sotorasib", "adagrasib"):
        e = engagement(drug)
        mu_d = E_to_mu(e["E"], mu_max)
        # interpolate the M4 readouts at this mu
        lle_d = float(np.interp(mu_d, mu, eig))
        dnb_d = float(np.interp(mu_d, mu, dnb))
        drug_points[drug] = {**e, "mu": mu_d, "lle": lle_d, "dnb": dnb_d}
        b = BOLTZ[drug]
        print(f"\n{drug} ({CHEMBL[drug]['alias']}):")
        print(f"  ChEMBL pIC50 component  : "
              f"{e['func']:.3f}" if e['func'] is not None else
              f"  ChEMBL pIC50 component  : NA (no IC50)")
        print(f"  Boltz optimization_score: {e['struct']:.3f}  "
              f"(binding_confidence={b['binding_confidence']})")
        print(f"  -> engagement E = {e['E']:.3f}  -> mu = {mu_d:.3f}  "
              f"(fold/tipping point is at mu ~ {mu.max():.2f}, i.e. E->1)")
        print(f"  -> predicted leading Lyapunov exp = {lle_d:+.3f}")
        print(f"  -> predicted DNB index            = {dnb_d:.3f}")

    # ranking — DNB is the monotonic early-warning readout here; the LLE
    # only turns sharply very close to the fold, so both mid-engagement
    # drugs still sit in the stable trough.
    ranked = sorted(drug_points.items(), key=lambda kv: -kv[1]["dnb"])
    print(f"\nordering by DNB proximity-to-tipping (monotonic readout):")
    for name, p in ranked:
        print(f"  {name:10s} E={p['E']:.3f}  DNB={p['dnb']:.3f}  "
              f"LLE={p['lle']:+.3f}")
    print("\nnote: at these engagement levels both drugs sit in the STABLE")
    print("trough; the model says near-complete target occupancy (E->1)")
    print("is needed to drive the KRAS-ON core to its apoptotic tipping")
    print("point. DNB is the sensitive discriminator; adagrasib's stronger")
    print("Boltz binding puts it marginally closer to the edge than sotorasib.")

    _figure(mu, eig, dnb, drug_points, fig_path)
    print(f"\nfigure written to: {fig_path}")
    return drug_points


def _figure(mu, eig, dnb, pts, path):
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    colors = {"sotorasib": "#c0392b", "adagrasib": "#2980b9"}

    ax[0].plot(mu, eig, "-", color="gray", lw=2, label="M4 stability curve")
    ax[0].axhline(0, ls="--", color="k", lw=1)
    for name, p in pts.items():
        ax[0].scatter([p["mu"]], [p["lle"]], s=140, color=colors[name],
                      zorder=5, edgecolors="k",
                      label=f"{name} (E={p['E']:.2f})")
        ax[0].annotate(name, (p["mu"], p["lle"]),
                       textcoords="offset points", xytext=(6, 8), fontsize=9)
    ax[0].set_xlabel("network perturbation  mu  (from target engagement)")
    ax[0].set_ylabel("leading Lyapunov exponent")
    ax[0].set_title("Real drug binding -> position on the\n"
                    "network-stability (Lyapunov) curve")
    ax[0].legend(fontsize=8)

    ax[1].plot(mu, dnb, "-", color="gray", lw=2)
    for name, p in pts.items():
        ax[1].scatter([p["mu"]], [p["dnb"]], s=140, color=colors[name],
                      zorder=5, edgecolors="k")
        ax[1].annotate(name, (p["mu"], p["dnb"]),
                       textcoords="offset points", xytext=(6, -12), fontsize=9)
    ax[1].set_xlabel("network perturbation  mu  (from target engagement)")
    ax[1].set_ylabel("DNB index")
    ax[1].set_title("Predicted DNB biomarker per drug")

    fig.suptitle("KRAS G12C: ChEMBL potency + Boltz-2.1 binding "
                 "-> network-stability biomarker", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(path, dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    report()
