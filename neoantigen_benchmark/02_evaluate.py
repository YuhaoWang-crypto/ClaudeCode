"""Three questions the matrix row makes claims about:
 A. Presentation: can MHCflurry separate real HLA-I epitopes from random proteome decoys? (should be yes)
 B. Immunogenicity gap: among peptides tested in T-cell assays, does the binding/presentation score
    separate immunogenic (Positive) from non-immunogenic (Negative)? And does it collapse once you
    condition on the peptide being a predicted strong binder? (the matrix's stated boundary)
 C. Leakage control: redo B on references from 2021+ only (after MHCflurry 2.0's training data)."""
import os, json, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, roc_curve, average_precision_score
from scipy.stats import mannwhitneyu
WORK = os.environ["IB_WORK"]; OUT = "results"
pairs = pd.read_parquet(f"{WORK}/pairs_scored.parquet").dropna(subset=["presentation_score"])
dec = pd.read_parquet(f"{WORK}/decoys_scored.parquet").dropna(subset=["presentation_score"])
S = {}

def auc(y, s):
    y = np.asarray(y); s = np.asarray(s)
    return round(float(roc_auc_score(y, s)), 3) if len(set(y)) == 2 and len(y) > 10 else None

# ---------- A. presentation vs decoys
pos = pairs[pairs.pos == 1]
ya = np.r_[np.ones(len(pos)), np.zeros(len(dec))]
A = {"n_epitopes": int(len(pos)), "n_decoys": int(len(dec))}
for sc in ["presentation_score", "processing_score"]:
    A[f"AUROC_{sc}"] = auc(ya, np.r_[pos[sc], dec[sc]])
A["AUROC_neg_affinity"] = auc(ya, np.r_[-pos.affinity, -dec.affinity])
# precision/recall-style: among top-scoring 2% by presentation, enrichment of true epitopes
allp = np.r_[pos.presentation_score, dec.presentation_score]
thr = np.quantile(dec.presentation_score, 0.98)
A["decoy_FPR_at_98th_decoy_pctile"] = 0.02
A["epitope_recall_at_2pct_decoy_FPR"] = round(float((pos.presentation_score > thr).mean()), 3)
A["fraction_epitopes_strong_binder_<500nM"] = round(float((pos.affinity < 500).mean()), 3)
A["fraction_decoys_strong_binder_<500nM"] = round(float((dec.affinity < 500).mean()), 3)
S["A_presentation_vs_decoys"] = A

# ---------- B. immunogenicity gap
def gap(d, label):
    r = {"n": int(len(d)), "n_pos": int(d.pos.sum()), "n_neg": int((1 - d.pos).sum())}
    if r["n_pos"] < 15 or r["n_neg"] < 15: return None
    r["AUROC_presentation_for_immunogenicity"] = auc(d.pos, d.presentation_score)
    r["AUROC_affinity_for_immunogenicity"] = auc(d.pos, -d.affinity)
    r["AUROC_processing_for_immunogenicity"] = auc(d.pos, d.processing_score)
    sb = d[d.affinity < 500]                                    # condition on predicted strong binding
    r["among_predicted_strong_binders"] = {"n": int(len(sb)), "n_pos": int(sb.pos.sum()), "n_neg": int((1 - sb.pos).sum()),
        "AUROC_presentation": auc(sb.pos, sb.presentation_score), "AUROC_affinity": auc(sb.pos, -sb.affinity),
        "pos_rate_overall": round(float(d.pos.mean()), 3), "pos_rate_among_strong_binders": round(float(sb.pos.mean()), 3)}
    return r
S["B_immunogenicity_gap_all"] = gap(pairs, "all")
S["B_by_allele"] = {al: gap(pairs[pairs.allele == al], al) for al in pairs.allele.value_counts().head(6).index}

# ---------- C. leakage control: references 2021+
recent = pairs[pairs.minyear >= 2021]
S["C_presentation_vs_decoys_note"] = "decoys are random proteome; presentation AUROC is not affected by IEDB recency"
S["C_immunogenicity_gap_refs_2021plus"] = gap(recent, "2021+")
# For presentation, also restrict positives to 2021+ (truly unseen epitopes) vs decoys
posr = recent[recent.pos == 1]
S["C_presentation_2021plus"] = {"n_epitopes": int(len(posr)),
    "AUROC_presentation": auc(np.r_[np.ones(len(posr)), np.zeros(len(dec))], np.r_[posr.presentation_score, dec.presentation_score]),
    "epitope_recall_at_2pct_decoy_FPR": round(float((posr.presentation_score > thr).mean()), 3)}

json.dump(S, open(f"{OUT}/summary.json", "w"), indent=1, default=str)

# ---------- figures
fig, ax = plt.subplots(1, 3, figsize=(17, 4.6))
# fig A: presentation separates epitopes from decoys
b = np.linspace(0, 1, 40)
ax[0].hist(dec.presentation_score, b, density=True, alpha=.6, color="#999", label=f"random decoys (n={len(dec)})")
ax[0].hist(pos.presentation_score, b, density=True, alpha=.6, color="#4575b4", label=f"real epitopes (n={len(pos)})")
ax[0].set(title=f"A. Presentation score separates epitopes from decoys\nAUROC={A['AUROC_presentation_score']}", xlabel="MHCflurry presentation score", ylabel="density"); ax[0].legend(fontsize=8)
# fig B: immunogenicity gap - score distribution pos vs neg, and ROC
g = S["B_immunogenicity_gap_all"]
ax[1].hist(pairs[pairs.pos == 0].presentation_score, b, density=True, alpha=.55, color="#d73027", label=f"non-immunogenic (T-cell Neg, n={g['n_neg']})")
ax[1].hist(pairs[pairs.pos == 1].presentation_score, b, density=True, alpha=.55, color="#1a9850", label=f"immunogenic (T-cell Pos, n={g['n_pos']})")
ax[1].set(title=f"B. Same score for immunogenicity\nAUROC={g['AUROC_presentation_for_immunogenicity']}, and {g['among_predicted_strong_binders']['AUROC_presentation']} among strong binders", xlabel="MHCflurry presentation score", ylabel="density"); ax[1].legend(fontsize=8)
# fig C: AUROC bars across tasks
tasks = ["epitope vs decoy\n(presentation)", "immunogenicity\n(all tested)", "immunogenicity\n(strong binders only)", "immunogenicity\n(refs 2021+)"]
vals = [A["AUROC_presentation_score"], g["AUROC_presentation_for_immunogenicity"], g["among_predicted_strong_binders"]["AUROC_presentation"], S["C_immunogenicity_gap_refs_2021plus"]["AUROC_presentation_for_immunogenicity"]]
cols = ["#4575b4", "#fdae61", "#d73027", "#fdae61"]
ax[2].bar(range(4), vals, color=cols); ax[2].axhline(.5, ls=":", c="k", lw=.8)
for i, v in enumerate(vals): ax[2].text(i, v + .01, f"{v:.2f}", ha="center", fontsize=9)
ax[2].set_xticks(range(4), tasks, fontsize=8); ax[2].set(title="C. MHCflurry presentation score: AUROC by task", ylabel="AUROC", ylim=(0, 1))
fig.tight_layout(); fig.savefig(f"{OUT}/fig1_presentation_vs_immunogenicity.png", dpi=130)
print(json.dumps(S, indent=1, default=str))
