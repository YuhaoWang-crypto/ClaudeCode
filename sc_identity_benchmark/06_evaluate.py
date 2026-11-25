"""Score every claim against held-out ground truth; write results/summary.json and figures."""
import json, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score, average_precision_score, f1_score, roc_curve
from common import *

Q = {"q_triana_wta": "Query A: Triana 2021 (BD Rhapsody)", "q_ts_bm": "Query B: Tabula Sapiens (10x)"}
ref_info = json.load(open(f"{OUT}/design.json"))
S = {"design": ref_info}
NOVEL = HELD_OUT_LINEAGES + ["Stroma"]

def cls(y, p):
    return {"n": int(len(y)), "accuracy": round(float((y == p).mean()), 4),
            "macro_F1": round(float(f1_score(y, p, average="macro", labels=sorted(set(y)))), 4)}

def recall_at_fpr(y, s, fpr=0.05):
    f, t, _ = roc_curve(y, s); return float(np.interp(fpr, f, t))

fig_cm, ax_cm = plt.subplots(1, 2, figsize=(15, 6.5))
fig_nv, ax_nv = plt.subplots(1, 2, figsize=(11, 4.2))
for k, (name, title) in enumerate(Q.items()):
    d = pd.read_csv(f"{OUT}/pred_scanvi_{name}.csv", index_col=0)
    c = pd.read_csv(f"{OUT}/pred_celltypist_{name}.csv", index_col=0)
    d = d.join(c.drop(columns=["fine", "coarse"]))
    d["knn_coarse"] = d.knn_fine.map(FINE2COARSE)
    r = {}
    known = d[~d.coarse.isin(NOVEL)]
    ref_fine = set(pd.unique(d.scanvi_fine)) | set(pd.unique(d.knn_fine)) | set(pd.unique(d.ct_custom_fine))
    # ---- 1. annotation of lineages the reference contains
    r["coarse_known_lineages"] = {m: cls(known.coarse, known[col]) for m, col in
        [("scANVI+scArches", "scanvi_coarse"), ("kNN on scANVI latent", "knn_coarse"),
         ("CellTypist (trained on our reference)", "ct_custom_coarse")]}
    r["coarse_all_lineages_pretrained_CellTypist"] = cls(d.coarse, d.ct_pre_coarse)
    r["per_lineage_recall_scANVI"] = {l: round(float((g.scanvi_coarse == l).mean()), 3)
                                      for l, g in known.groupby("coarse")}
    r["per_lineage_recall_CellTypist_custom"] = {l: round(float((g.ct_custom_coarse == l).mean()), 3)
                                                 for l, g in known.groupby("coarse")}
    r["per_lineage_recall_CellTypist_pretrained"] = {l: round(float((g.ct_pre_coarse == l).mean()), 3)
                                                     for l, g in d.groupby("coarse")}
    fine = known[known.fine.isin(ref_fine)]
    r["fine_subtypes_shared_with_reference"] = {
        "n_shared_fine_labels": int(fine.fine.nunique()),
        "scANVI+scArches": cls(fine.fine, fine.scanvi_fine),
        "CellTypist (trained on our reference)": cls(fine.fine, fine.ct_custom_fine),
        "per_label_recall_scANVI": {l: [int(len(g)), round(float((g.scanvi_fine == l).mean()), 3)]
                                    for l, g in fine.groupby("fine")}}
    # mature neutrophils: only 114 reference cells survive the study hold-out -> report separately
    nomn = known[~known.fine.isin(["neutrophil", "granulocyte"])]
    r["coarse_known_excluding_mature_neutrophils"] = {m: cls(nomn.coarse, nomn[col]) for m, col in
        [("scANVI+scArches", "scanvi_coarse"), ("CellTypist (trained on our reference)", "ct_custom_coarse")]}
    # ---- 2. unknown-population alert
    y = d.coarse.isin(NOVEL).values.astype(int)
    scores = {"scANVI 1-max prob": 1 - d.scanvi_maxp, "kNN distance to reference (latent)": d.knn_dist,
              "CellTypist(custom) 1-max prob": 1 - d.ct_custom_maxp}
    r["novelty"] = {"n_novel": int(y.sum()), "n_known": int((1 - y).sum()), "scores": {}}
    for sname, s in scores.items():
        e = {"AUROC": round(roc_auc_score(y, s), 3), "AUPRC": round(average_precision_score(y, s), 3),
             "prevalence": round(float(y.mean()), 3), "recall_at_5pct_FPR": round(recall_at_fpr(y, s), 3)}
        for l in NOVEL:
            m = (d.coarse == l).values | (y == 0)
            if (d.coarse == l).sum() >= 20:
                e[f"AUROC_{l}_vs_known"] = round(roc_auc_score(y[m], s[m]), 3)
        r["novelty"]["scores"][sname] = e
        f, t, _ = roc_curve(y, s); ax_nv[k].plot(f, t, label=f"{sname} (AUROC {e['AUROC']:.2f})")
    r["novelty"]["held_out_cells_were_called"] = {
        l: d.loc[d.coarse == l, "scanvi_coarse"].value_counts().head(3).to_dict() for l in NOVEL if (d.coarse == l).any()}
    ax_nv[k].plot([0, 1], [0, 1], "k:", lw=.8); ax_nv[k].set(title=f"{title}\nheld-out lineages vs known ({int(y.sum())} vs {int((1-y).sum())})",
                                                         xlabel="false-positive rate", ylabel="recall of unseen lineages")
    ax_nv[k].legend(fontsize=7.5, loc="lower right")
    cm = pd.crosstab(d.coarse, d.scanvi_coarse, normalize="index")
    cm = cm.reindex(columns=[x for x in COARSE if x in cm.columns], index=[x for x in COARSE if x in cm.index])
    im = ax_cm[k].imshow(cm.values, cmap="Blues", vmin=0, vmax=1)
    ax_cm[k].set_xticks(range(cm.shape[1]), cm.columns, rotation=60, ha="right"); ax_cm[k].set_yticks(
        range(cm.shape[0]), [f"{i}{' *' if i in NOVEL else ''} (n={int((d.coarse==i).sum())})" for i in cm.index])
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            v = cm.values[i, j]
            if v >= 0.05: ax_cm[k].text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7, color="w" if v > .5 else "k")
    ax_cm[k].set(title=f"{title}\nscANVI+scArches, row-normalised (* = lineage absent from reference)",
                 xlabel="predicted lineage", ylabel="author annotation")
    S[name] = r
fig_cm.tight_layout(); fig_cm.savefig(f"{OUT}/fig1_annotation_confusion.png", dpi=140)
fig_nv.tight_layout(); fig_nv.savefig(f"{OUT}/fig2_novelty_roc.png", dpi=140)

# ---- 3. potency (CytoTRACE 2)
fig, ax = plt.subplots(1, 2, figsize=(15, 6))
for k, (name, title) in enumerate(Q.items()):
    d = pd.read_csv(f"{OUT}/cytotrace2_{name}.csv", index_col=0).dropna(subset=["CytoTRACE2_Score"])
    hi, lo = d.potency_truth == 2, d.potency_truth == 0
    r = {"n_cells": int(len(d))}
    for sname, col in [("CytoTRACE2_Score", "CytoTRACE2_Score"), ("baseline_n_genes_detected", "n_genes_detected")]:
        r[sname] = {"spearman_vs_ordinal_potency": round(float(spearmanr(d[col], d.potency_truth)[0]), 3),
                    "AUROC_multipotent_vs_mature": round(roc_auc_score(hi[hi | lo], d.loc[hi | lo, col]), 3),
                    "AUROC_committed_progenitor_vs_mature": round(roc_auc_score((d.potency_truth == 1)[~hi], d.loc[~hi, col]), 3),
                    "AUROC_multipotent_vs_committed_progenitor": round(roc_auc_score(hi[~lo], d.loc[~lo, col]), 3)}
    r["mean_score_by_level"] = d.groupby("potency_truth").CytoTRACE2_Score.mean().round(3).to_dict()
    r["predicted_category_by_level"] = pd.crosstab(d.potency_truth, d.CytoTRACE2_Potency).to_dict()
    S[name]["potency"] = r
    order = d.groupby("fine").CytoTRACE2_Score.median().sort_values(ascending=False).index
    col = {2: "#b2182b", 1: "#ef8a62", 0: "#67a9cf"}
    bp = ax[k].boxplot([d.loc[d.fine == f, "CytoTRACE2_Score"] for f in order], vert=False, patch_artist=True,
                       showfliers=False, widths=.7)
    for patch, f in zip(bp["boxes"], order): patch.set_facecolor(col[FINE2POTENCY.get(f, 0)])
    ax[k].set_yticks(range(1, len(order) + 1), order, fontsize=7); ax[k].invert_yaxis()
    ax[k].set(title=f"{title}\nred = multipotent, orange = committed progenitor, blue = mature (author labels)",
              xlabel="CytoTRACE 2 potency score")
fig.tight_layout(); fig.savefig(f"{OUT}/fig3_cytotrace2_potency.png", dpi=140)

# ---- 4. senescence (SenCID)
d = pd.read_csv(f"{OUT}/sencid_GSE115301.csv", index_col=0)
g, ris, gfp = d.Condition2 == "Growing", d.Condition2 == "RIS", d.Condition2 == "GFP"
cols = ["SID1", "SID2", "SID3", "SID4", "SID_mean", "CDKN1A", "CDKN2A", "neg_MKI67_TOP2A"]
S["senescence_GSE115301"] = {
    "n": {"growing": int(g.sum()), "RIS_primary": int(ris.sum()), "GFP_bystander_secondary": int(gfp.sum())},
    "AUROC_RIS_vs_growing": {c: round(roc_auc_score(ris[g | ris], d.loc[g | ris, c]), 3) for c in cols},
    "AUROC_bystander_vs_growing": {c: round(roc_auc_score(gfp[g | gfp], d.loc[g | gfp, c]), 3) for c in cols},
    "fraction_called_senescent_SID>0.5": {c: {n: round(float((d.loc[m, c] > .5).mean()), 3) for n, m in
                                          [("growing", g), ("RIS", ris), ("bystander", gfp)]} for c in cols[:5]}}
fig, ax = plt.subplots(1, 5, figsize=(14, 3.2), sharey=True)
for a, c in zip(ax, cols[:5]):
    a.violinplot([d.loc[m, c] for m in (g, gfp, ris)], showmedians=True)
    a.set_xticks([1, 2, 3], ["growing", "bystander", "RIS"]); a.axhline(.5, ls=":", c="k", lw=.8)
    a.set_title(f"{c}  (AUROC RIS vs growing {S['senescence_GSE115301']['AUROC_RIS_vs_growing'][c]:.2f})", fontsize=9)
ax[0].set_ylabel("SenCID score")
fig.tight_layout(); fig.savefig(f"{OUT}/fig4_sencid.png", dpi=140)

json.dump(S, open(f"{OUT}/summary.json", "w"), indent=1, default=str)
print(json.dumps({k: v for k, v in S.items() if k != "design"}, indent=1, default=str))
