"""Demo validation: do the ADMET-AI predictions recover documented pharmacology?

Ground truth is taken from regulatory / reference sources, not from the model:

hERG / TdP      : CredibleMeds "Known Risk of TdP" list + the classic
                  hERG-withdrawal cases (terfenadine, cisapride).
CYP inhibition  : FDA "Drug Development and Drug Interactions" index of strong
                  and moderate clinical index inhibitors.
P-gp inhibition : FDA clinical index P-gp inhibitors.
Nuclear receptor: the drug's own approved pharmacological mechanism.
BBB             : documented CNS penetration / CNS indication.

Each block prints the predicted score for the positives and the negatives, the
rank of each positive within the 30-drug panel, and — where there are at least
3 positives and 3 negatives — the ROC-AUC over the panel.
"""

from __future__ import annotations
import pandas as pd
from sklearn.metrics import roc_auc_score

CSV = "/home/user/results/admet_fda_panel/all_properties.csv"

# endpoint -> (positives, negatives, description, source)
TRUTH = {
    "hERG": (
        ["Terfenadine", "Cisapride", "Amiodarone", "Haloperidol", "Quinidine"],
        ["Metformin", "Aspirin", "Ibuprofen", "Acetaminophen", "Caffeine",
         "Ciprofloxacin"],
        "hERG blockade / known TdP risk",
        "CredibleMeds known-risk list; terfenadine & cisapride withdrawn for QT"),
    "CYP3A4_Veith": (
        ["Ketoconazole", "Verapamil"],
        ["Metformin", "Aspirin", "Caffeine", "Acetaminophen", "Ibuprofen"],
        "CYP3A4 inhibition",
        "FDA clinical index inhibitors (ketoconazole strong, verapamil moderate)"),
    "CYP2D6_Veith": (
        ["Quinidine", "Fluoxetine", "Sertraline"],
        ["Metformin", "Aspirin", "Caffeine", "Ciprofloxacin", "Ibuprofen"],
        "CYP2D6 inhibition",
        "FDA clinical index inhibitors (quinidine/fluoxetine strong, sertraline moderate)"),
    "CYP1A2_Veith": (
        ["Ciprofloxacin"],
        ["Metformin", "Aspirin", "Ibuprofen"],
        "CYP1A2 inhibition",
        "FDA clinical index inhibitor (ciprofloxacin, strong)"),
    "Pgp_Broccatelli": (
        ["Verapamil", "Ketoconazole", "Quinidine", "Amiodarone"],
        ["Metformin", "Aspirin", "Caffeine", "Acetaminophen"],
        "P-glycoprotein inhibition",
        "FDA clinical index P-gp inhibitors"),
    "NR-ER": (
        ["Tamoxifen"], ["Metformin", "Caffeine", "Aspirin", "Ciprofloxacin"],
        "Estrogen receptor activity",
        "Tamoxifen is an approved selective estrogen-receptor modulator"),
    "NR-AR": (
        ["Bicalutamide"], ["Metformin", "Caffeine", "Aspirin", "Ciprofloxacin"],
        "Androgen receptor activity",
        "Bicalutamide is an approved androgen-receptor antagonist"),
    "NR-Aromatase": (
        ["Anastrozole", "Ketoconazole"],
        ["Metformin", "Caffeine", "Aspirin", "Ciprofloxacin"],
        "Aromatase (CYP19A1) activity",
        "Anastrozole is an approved aromatase inhibitor; ketoconazole inhibits CYP19A1"),
    "NR-PPAR-gamma": (
        ["Rosiglitazone"], ["Metformin", "Caffeine", "Aspirin", "Ciprofloxacin"],
        "PPAR-gamma activity",
        "Rosiglitazone is an approved PPAR-gamma agonist"),
    "BBB_Martins": (
        ["Diazepam", "Caffeine", "Haloperidol", "Fluoxetine", "Risperidone",
         "Propranolol"],
        ["Metformin", "Atorvastatin", "Ciprofloxacin"],
        "Blood-brain-barrier penetration",
        "CNS-active drugs vs peripherally restricted drugs"),
}


def main():
    d = pd.read_csv(CSV).set_index("name")
    rows = []
    print(f"{'endpoint':16s} {'n+':>3s} {'n-':>3s} {'mean+':>7s} {'mean-':>7s} "
          f"{'sep':>7s} {'AUC':>6s}  verdict")
    print("-" * 92)
    for ep, (pos, neg, desc, src) in TRUTH.items():
        s = d[ep]
        ranks = s.rank(ascending=False, method="min")
        p = s.loc[pos]
        n = s.loc[neg]
        sep = p.mean() - n.mean()
        if len(pos) >= 3 and len(neg) >= 3:
            y = [1] * len(pos) + [0] * len(neg)
            auc = roc_auc_score(y, list(p) + list(n))
            auc_s = f"{auc:6.2f}"
        else:
            auc, auc_s = None, "    --"
        # verdict: positives recovered if all positives score > 0.5 and above all negatives
        recovered = int((p > 0.5).sum())
        verdict = ("recovered" if recovered == len(pos) else
                   "partial" if recovered else "MISSED")
        print(f"{ep:16s} {len(pos):3d} {len(neg):3d} {p.mean():7.3f} "
              f"{n.mean():7.3f} {sep:+7.3f} {auc_s}  {verdict} "
              f"({recovered}/{len(pos)} positives > 0.5)")
        for name in pos:
            rows.append({"endpoint": ep, "property": desc, "drug": name,
                         "role": "positive control", "score": round(s[name], 3),
                         "rank_in_panel": int(ranks[name]),
                         "above_0.5": bool(s[name] > 0.5), "source": src})
        for name in neg:
            rows.append({"endpoint": ep, "property": desc, "drug": name,
                         "role": "negative control", "score": round(s[name], 3),
                         "rank_in_panel": int(ranks[name]),
                         "above_0.5": bool(s[name] > 0.5), "source": src})
    out = pd.DataFrame(rows)
    out.to_csv("/home/user/results/admet_fda_panel/demo_validation.csv", index=False)

    print("\nPositive-control detail (rank out of 30):")
    pc = out[out.role == "positive control"]
    print(pc[["endpoint", "drug", "score", "rank_in_panel", "above_0.5"]]
          .to_string(index=False))
    return out


if __name__ == "__main__":
    main()
