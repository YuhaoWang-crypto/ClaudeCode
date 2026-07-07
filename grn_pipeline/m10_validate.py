"""
Module 10 — Validate the Boltz-2.1 screen ranking against ChEMBL experiment.

For the M7 hits we pulled the REAL experimental KRAS G12C potency from ChEMBL
(this session) and test how well the Boltz screen metrics predict it.

KEY FINDING: the two Boltz metrics behave very differently.
  * binding_confidence (a POSE-confidence score) does NOT track potency
    (Spearman ~ -0.2): it over-ranked weak binders. The earlier M7 claim
    "4 analogues beat sotorasib" was based on this metric and does NOT hold up.
  * optimization_score (Boltz's affinity-oriented metric) DOES track potency
    (Spearman ~ +0.6). This is the metric to rank by.

Caveats (honest): experimental IC50s are assay-heterogeneous (cell p-ERK vs
biochemical nucleotide-exchange), covalent IC50 is time-dependent, and n is
small. The clean same-assay subset (MIAPaCa2 p-ERK, CHEMBL4357259) is flagged.
"""
import math
import numpy as np
from scipy.stats import spearmanr, pearsonr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from grn_pipeline.m7_screen import SCREEN

# Real ChEMBL experimental potency vs KRAS G12C (CHEMBL2189121), this session.
# assay 'A' = MIAPaCa2 p-ERK CHEMBL4357259 (directly comparable subset).
EXP = {
    "sotorasib_CHEMBL4535757": {"ic50_nM": 30.0,  "pIC50": 7.52, "assay": "A: MIAPaCa2 pERK"},
    "adagrasib_CHEMBL4594350": {"ic50_nM": 14.0,  "pIC50": 7.85, "assay": "NCI-H358 pERK"},
    "CHEMBL4539214":           {"ic50_nM": 11.0,  "pIC50": 7.96, "assay": "A: MIAPaCa2 pERK"},
    "CHEMBL4441771":           {"ic50_nM": 128.0, "pIC50": 6.89, "assay": "A: MIAPaCa2 pERK"},
    "CHEMBL5918612":           {"ic50_nM": 419.0, "pIC50": 6.38, "assay": "biochem nucleotide-exchange"},
    # CHEMBL4851723: only qualitative binding in ChEMBL (no numeric IC50) -> excluded
}
SAME_ASSAY = {"sotorasib_CHEMBL4535757", "CHEMBL4539214", "CHEMBL4441771"}


def _corr(names, metric_key):
    x = [SCREEN[n][metric_key] for n in names]
    y = [EXP[n]["pIC50"] for n in names]
    rho, p = spearmanr(x, y)
    return rho, p, x, y


def report(fig_path="figures/m10_validation.png"):
    names = [n for n in EXP]                      # 5 with numeric IC50
    print("=" * 68)
    print("MODULE 10 — Boltz screen ranking  vs  ChEMBL experimental potency")
    print("=" * 68)
    print(f"{'molecule':26s} {'IC50(nM)':>9} {'pIC50':>6} "
          f"{'bind_conf':>9} {'opt_score':>9}  assay")
    for n in sorted(names, key=lambda k: -EXP[k]["pIC50"]):
        s = SCREEN[n]; e = EXP[n]
        print(f"{n:26s} {e['ic50_nM']:9.0f} {e['pIC50']:6.2f} "
              f"{s['bind_conf']:9.3f} {s['opt_score']:9.3f}  {e['assay']}")

    rho_bc, p_bc, _, _ = _corr(names, "bind_conf")
    rho_os, p_os, _, _ = _corr(names, "opt_score")
    print(f"\nSpearman rank correlation vs experimental pIC50 (n={len(names)}):")
    print(f"  binding_confidence  vs pIC50 : rho = {rho_bc:+.2f}  (p={p_bc:.2f})")
    print(f"  optimization_score  vs pIC50 : rho = {rho_os:+.2f}  (p={p_os:.2f})")

    sub = [n for n in names if n in SAME_ASSAY]
    rho_bc_s, _, _, _ = _corr(sub, "bind_conf")
    rho_os_s, _, _, _ = _corr(sub, "opt_score")
    print(f"\nclean same-assay subset (MIAPaCa2 p-ERK, n={len(sub)}):")
    print(f"  binding_confidence : rho = {rho_bc_s:+.2f}   "
          f"optimization_score : rho = {rho_os_s:+.2f}")

    print("\nVERDICT:")
    print("  * binding_confidence does NOT predict potency -> the earlier M7")
    print("    '4 analogues beat sotorasib' ranking is NOT experimentally")
    print("    supported (it promoted 419 nM CHEMBL5918612 above 30 nM sotorasib).")
    print("  * optimization_score DOES rank-correlate (rho~+0.6): the correct")
    print("    Boltz metric to screen on. Only CHEMBL4539214 (11 nM) genuinely")
    print("    out-potents sotorasib (30 nM) in the same assay.")
    print("  caveats: mixed assays, covalent time-dependence, small n.")

    _figure(names, sub, rho_bc, rho_os, fig_path)
    print(f"\nfigure written to: {fig_path}")
    return {"rho_bind_conf": rho_bc, "rho_opt_score": rho_os,
            "genuinely_better": ["CHEMBL4539214"]}


def _figure(names, sub, rho_bc, rho_os, path):
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    for a, key, rho, title in [
        (ax[0], "bind_conf", rho_bc, "binding_confidence (pose)"),
        (ax[1], "opt_score", rho_os, "optimization_score (affinity)")]:
        for n in names:
            approved = SCREEN[n]["approved"]
            col = "#e67e22" if approved else "#3498db"
            mk = "*" if n in sub else "o"
            a.scatter(SCREEN[n][key], EXP[n]["pIC50"], s=170 if approved else 110,
                      c=col, marker=mk, edgecolors="k", zorder=5)
            lbl = n.split("_")[0] if approved else n.replace("CHEMBL", "…")
            a.annotate(lbl, (SCREEN[n][key], EXP[n]["pIC50"]),
                       textcoords="offset points", xytext=(6, 5), fontsize=8)
        a.set_xlabel(f"Boltz {key}")
        a.set_ylabel("experimental pIC50 (ChEMBL)")
        a.set_title(f"{title}\nSpearman rho = {rho:+.2f}")
    fig.suptitle("Boltz-2.1 screen metric vs ChEMBL potency — orange=approved, "
                 "star=same-assay subset", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(path, dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    report()
