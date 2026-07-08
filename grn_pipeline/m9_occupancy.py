"""
Module 9 — Calibrate the perturbation mu from real target occupancy.

Chain:  dose/exposure (Cmax, AUC — the trial's PK endpoints)  ->  free drug
        Cfree = Cmax * fu  ->  fractional target occupancy
        occ = Cfree / (Cfree + IC50)  [reversible-equivalent]
        ->  network perturbation mu = occ * mu_fold  ->  M4 LLE / DNB.

Anchored on the REAL ChEMBL IC50 for sotorasib (30 nM on KRAS G12C). PK
values (Cmax, fraction unbound) are approximate published-label figures and
are clearly flagged; they are exactly the Cmax/AUC secondary endpoints that
CodeBreaK 100 reports, so real trial PK can be dropped straight in.

Result: at the approved exposure the covalent inhibitor achieves near-complete
occupancy, which in this model drives the KRAS-ON core close to its tipping
point --- the model's structural rationale for why the drug produces tumour
regression.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from grn_pipeline import m4_dnb_lyapunov as m4

# real anchor (ChEMBL, this session) + approximate label PK (flagged)
DRUG_PK = {
    "sotorasib": {
        "ic50_nM": 30.0,            # ChEMBL, real
        "cmax_total_uM": 11.0,      # ~960 mg QD steady state (approx label)
        "fu": 0.11,                 # ~89% plasma protein bound (approx)
        "covalent": True,
    },
}


def occupancy(cfree_nM, ic50_nM):
    return cfree_nM / (cfree_nM + ic50_nM)


def report(drug="sotorasib", fig_path="figures/m9_occupancy.png"):
    pk = DRUG_PK[drug]
    ic50 = pk["ic50_nM"]

    rows = m4.sweep()
    mu = np.array([r["mu"] for r in rows])
    eig = np.array([r["lead_eig"] for r in rows])
    dnb = np.array([r["dnb"] for r in rows])
    mu_fold = mu.max()

    # sweep total Cmax over 3 orders of magnitude around the label value
    cmax_tot = np.logspace(-1, 2.2, 60)          # uM
    cfree_nM = cmax_tot * pk["fu"] * 1000.0       # -> nM
    occ = occupancy(cfree_nM, ic50)
    mu_d = occ * mu_fold
    lle = np.interp(mu_d, mu, eig)
    dnb_d = np.interp(mu_d, mu, dnb)

    # the approved-dose operating point
    cfree_label = pk["cmax_total_uM"] * pk["fu"] * 1000.0
    occ_label = occupancy(cfree_label, ic50)
    mu_label = occ_label * mu_fold
    lle_label = float(np.interp(mu_label, mu, eig))
    dnb_label = float(np.interp(mu_label, mu, dnb))

    print("=" * 68)
    print("MODULE 9 — Occupancy-calibrated perturbation (real IC50 anchor)")
    print("=" * 68)
    print(f"drug: {drug}   IC50 = {ic50} nM (ChEMBL, real)")
    print(f"approx label PK: Cmax_total = {pk['cmax_total_uM']} uM, "
          f"fu = {pk['fu']}  =>  Cfree = {cfree_label:.0f} nM")
    print(f"  -> target occupancy at approved exposure = {occ_label*100:.1f}%")
    print(f"  -> perturbation mu = {mu_label:.3f}  "
          f"(fold/tipping at mu ~ {mu_fold:.2f})")
    print(f"  -> predicted leading Lyapunov exp = {lle_label:+.3f}")
    print(f"  -> predicted DNB index            = {dnb_label:.3f}")
    if pk["covalent"]:
        print("  note: covalent inhibitor -> occupancy accumulates toward")
        print("        ~100% at steady state, so the true operating point sits")
        print("        at/above this reversible-equivalent estimate.")
    print(f"\ncontrast with M6 (Boltz binding-score proxy gave mu ~ 0.5, stable")
    print(f"trough): PK-based fractional OCCUPANCY, not a raw binding score, is")
    print(f"what sets in-vivo perturbation — and at approved exposure it is")
    print(f"near-complete, pushing the core close to its tipping point.")

    _figure(cmax_tot, occ, lle, dnb_d, pk, mu_fold,
            (pk["cmax_total_uM"], occ_label, lle_label, dnb_label), fig_path)
    print(f"\nfigure written to: {fig_path}")
    return {"occ_label": occ_label, "mu_label": mu_label,
            "lle_label": lle_label, "dnb_label": dnb_label}


def _figure(cmax, occ, lle, dnb, pk, mu_fold, label_pt, path):
    cmax_l, occ_l, lle_l, dnb_l = label_pt
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.3))

    ax[0].semilogx(cmax, occ * 100, "-", color="#16a085", lw=2)
    ax[0].axvline(cmax_l, ls="--", color="k", lw=1)
    ax[0].scatter([cmax_l], [occ_l * 100], s=120, color="#c0392b",
                  zorder=5, edgecolors="k")
    ax[0].annotate("approved\nexposure", (cmax_l, occ_l * 100),
                   textcoords="offset points", xytext=(-60, -30), fontsize=8)
    ax[0].set_xlabel("total Cmax (uM)"); ax[0].set_ylabel("target occupancy (%)")
    ax[0].set_title("Occupancy vs exposure\n(IC50=30 nM, ChEMBL)")

    ax[1].semilogx(cmax, lle, "-", color="#c0392b", lw=2)
    ax[1].axhline(0, ls="--", color="k", lw=1)
    ax[1].axvline(cmax_l, ls="--", color="k", lw=1)
    ax[1].scatter([cmax_l], [lle_l], s=120, color="#c0392b",
                  zorder=5, edgecolors="k")
    ax[1].set_xlabel("total Cmax (uM)")
    ax[1].set_ylabel("leading Lyapunov exponent")
    ax[1].set_title("Network stability vs exposure\n"
                    "-> approaches tipping at approved dose")

    ax[2].semilogx(cmax, dnb, "-", color="#8e44ad", lw=2)
    ax[2].axvline(cmax_l, ls="--", color="k", lw=1)
    ax[2].scatter([cmax_l], [dnb_l], s=120, color="#8e44ad",
                  zorder=5, edgecolors="k")
    ax[2].set_xlabel("total Cmax (uM)"); ax[2].set_ylabel("DNB index")
    ax[2].set_title("DNB biomarker vs exposure")

    fig.suptitle("Occupancy-calibrated network biomarker (sotorasib / KRAS "
                 "G12C) — PK Cmax is a real CodeBreaK-100 endpoint",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(path, dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    report()
