"""Subtype specificity: do inducer-specific gene programs (CellAge replicative/oncogene/stress) peak
in the matching senescence type? Plus the summary figure."""
import os, json, csv, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
SEN = os.environ["SEN_DIR"]; Q = os.environ["QB_DATA"]; OUT = "results"
R = json.load(open(f"{OUT}/scores_summary.json"))
B = pd.read_csv(f"{OUT}/bulk_scores.csv", index_col=0)
ens2sym = {}
for x in csv.DictReader(open(f"{Q}/hgnc.tsv"), delimiter="\t"):
    if x.get("ensembl_gene_id"): ens2sym[x["ensembl_gene_id"]] = x["symbol"]
# CellAge inducer-specific gene sets
ca = list(csv.DictReader(open(f"{SEN}/cellage3.tsv"), delimiter="\t"))
subsets = {}
for t in ["Replicative", "Oncogene-induced", "Stress-induced"]:
    subsets[t] = {x["Gene symbol"] for x in ca if x.get("Senescence Effect") == "Induces" and x.get("Type of senescence") == t}
# recompute bulk z from scores file is not possible (B has signature scores); reload expression
import glob
files = glob.glob(f"{SEN}/GSM*_*_geneCOUNT.txt.gz")
sm = {}
for f in files:
    base = os.path.basename(f)
    if not any(k in base for k in ["Control", "Senescence"]): continue
    sm[base.split("geneCOUNT")[0].strip("_")] = pd.read_csv(f, sep="\t", index_col=0).iloc[:, 0]
M = pd.DataFrame(sm); M.index = [ens2sym.get(i.split(".")[0]) for i in M.index]; M = M[M.index.notna()].groupby(level=0).sum()
log = np.log2(M / M.sum(0) * 1e6 + 1); z = log.sub(log.mean(1), axis=0).div(log.std(1).replace(0, 1), axis=0)
def cond(c):
    for k in ["Replicative", "Doxorubicin", "IR_Senescence", "Young_Control", "DMSO_Control", "IR_Control"]:
        if k in c: return k
    return c
cc = {c: cond(c) for c in z.columns}
pairs = {"Replicative": "Young_Control", "Doxorubicin": "DMSO_Control", "IR_Senescence": "IR_Control"}
submat = {}
for t, genes in subsets.items():
    present = [g for g in genes if g in z.index]; sig = z.loc[present].mean(0)
    submat[t] = {ind: round(float(np.mean([sig[c] for c in z.columns if cc[c] == ind]) -
                                 np.mean([sig[c] for c in z.columns if cc[c] == ctrl])), 2) for ind, ctrl in pairs.items()}
json.dump({"subtype_specific_score_diff": submat, "subset_sizes": {t: len([g for g in s if g in z.index]) for t, s in subsets.items()}},
          open(f"{OUT}/subtype.json", "w"), indent=1)

# ---------- figures
fig, ax = plt.subplots(1, 3, figsize=(18, 5))
# A: single-cell AUROC bake-off (exclude SenCID gene-set-mean misuse; annotate)
A = R["singlecell_AUROC_RAS_senescence"]
order = ["SenMayo", "Core_markers_up", "Proliferation_down", "CellAge_induces", "Fridman_up", "Reactome_senescence"]
vals = [A[k] for k in order]
cols = ["#1a9850" if v > 0.8 else "#fdae61" if v > 0.65 else "#d73027" for v in vals]
ax[0].barh(range(len(order)), vals, color=cols); ax[0].axvline(0.5, ls=":", c="k", lw=.8)
ax[0].set_yticks(range(len(order)), order, fontsize=9); ax[0].invert_yaxis()
for i, v in enumerate(vals): ax[0].text(v + .01, i, f"{v:.2f}", va="center", fontsize=9)
ax[0].set(xlabel="AUROC (RAS-senescent vs growing, single cell)", xlim=(0, 1),
          title="A. Senescence-signature bake-off\nsame data, AUROC 0.50–0.91 (SenCID model=0.89, run separately)")
# B: cross-inducer transfer heatmap
cross = pd.DataFrame(R["bulk_cross_inducer_score_difference_sen_minus_control"]).T.loc[order]
im = ax[1].imshow(cross.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
ax[1].set_xticks(range(3), ["Replicative", "Doxorubicin", "IR"], fontsize=9); ax[1].set_yticks(range(len(order)), order, fontsize=8)
for i in range(len(order)):
    for j in range(3): ax[1].text(j, i, f"{cross.values[i,j]:.2f}", ha="center", va="center", fontsize=8, color="w" if abs(cross.values[i,j])>.6 else "k")
ax[1].set_title("B. Cross-inducer transfer (bulk)\nsen − control z-score; red=up, blue=down")
# C: inter-signature agreement
corr = pd.DataFrame(R["inter_signature_spearman_singlecell"]).loc[order, order]
im2 = ax[2].imshow(corr.values, cmap="RdBu_r", vmin=-1, vmax=1)
ax[2].set_xticks(range(len(order)), order, rotation=55, ha="right", fontsize=8); ax[2].set_yticks(range(len(order)), order, fontsize=8)
for i in range(len(order)):
    for j in range(len(order)): ax[2].text(j, i, f"{corr.values[i,j]:.1f}", ha="center", va="center", fontsize=7, color="w" if abs(corr.values[i,j])>.6 else "k")
ax[2].set_title("C. Signatures disagree (per-cell Spearman)\nsame cells, correlation 0.1–0.7")
fig.tight_layout(); fig.savefig(f"{OUT}/fig1_senescence.png", dpi=130)
print("subtype specificity (diagonal should be highest if subtypes separable):")
print(pd.DataFrame(submat).T.to_string())
print("\nsubset sizes:", json.load(open(f"{OUT}/subtype.json"))["subset_sizes"])
