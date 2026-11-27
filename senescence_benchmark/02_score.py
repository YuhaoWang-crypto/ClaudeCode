"""Score senescence signatures on (A) single-cell RAS-induced senescence (GSE115301) and
(B) bulk multi-inducer senescence (Casella GSE130727: replicative/doxorubicin/IR). Compare
discrimination, inter-signature agreement, and cross-inducer transfer / subtype-specificity."""
import os, json, glob, numpy as np, pandas as pd, scanpy as sc, anndata as ad
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr
SEN = os.environ["SEN_DIR"]; Q = os.environ["QB_DATA"]; OUT = "results"
sig = {k: set(v) for k, v in json.load(open(f"{OUT}/signatures.json")).items()}
ens2sym = {}
import csv
for x in csv.DictReader(open(f"{Q}/hgnc.tsv"), delimiter="\t"):
    e = x.get("ensembl_gene_id", "")
    if e: ens2sym[e] = x["symbol"]

# ---------- A. single-cell GSE115301 (growing vs RAS-induced senescent)
cnt = pd.read_csv(f"{SEN}/../data/sen/GSE115301_Growing_Sen_10x_count.txt.gz", sep="\t", index_col=0)
meta = pd.read_csv(f"{SEN}/../data/sen/GSE115301_Growing_Sen_10x_metadata.txt.gz", sep="\t", index_col=0).loc[cnt.columns]
a = ad.AnnData(cnt.T.values.astype(np.float32), obs=meta.copy(), var=pd.DataFrame(index=cnt.index))
sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
a.obs["sen"] = (a.obs["Condition2"] == "RIS").astype(int)
g = a[a.obs.Condition2.isin(["Growing", "RIS"])].copy()
sc_scores = {}
for name, genes in sig.items():
    present = [x for x in genes if x in g.var_names]
    if len(present) < 3: continue
    sc.tl.score_genes(g, present, score_name=f"s_{name}", ctrl_size=50)
    s = g.obs[f"s_{name}"].values
    if name == "Proliferation_down": s = -s           # loss of proliferation = senescence
    sc_scores[name] = s
sc_df = pd.DataFrame(sc_scores, index=g.obs_names); sc_df["sen"] = g.obs["sen"].values
sc_df.to_csv(f"{OUT}/singlecell_scores.csv")
A = {name: round(float(roc_auc_score(sc_df.sen, sc_df[name])), 3) for name in sc_scores}

# ---------- B. Casella bulk multi-inducer
files = glob.glob(f"{SEN}/GSM*_*_geneCOUNT.txt.gz")
samples = {}
for f in files:
    base = os.path.basename(f)
    if not any(k in base for k in ["Control", "Senescence"]): continue
    s = pd.read_csv(f, sep="\t", index_col=0).iloc[:, 0]
    samples[base.split("geneCOUNT")[0].strip("_")] = s
M = pd.DataFrame(samples)
M.index = [ens2sym.get(i.split(".")[0]) for i in M.index]
M = M[M.index.notna()]; M = M.groupby(level=0).sum()
cpm = M / M.sum(0) * 1e6; logcpm = np.log2(cpm + 1)
z = logcpm.sub(logcpm.mean(1), axis=0).div(logcpm.std(1).replace(0, 1), axis=0)   # z across samples per gene
def cond(name):
    for k in ["Replicative", "Doxorubicin", "IR_Senescence", "Young_Control", "DMSO_Control", "IR_Control"]:
        if k in name: return k
    return name
colcond = {c: cond(c) for c in z.columns}
bulk = {}
for name, genes in sig.items():
    present = [x for x in genes if x in z.index]
    sc_sig = z.loc[present].mean(0)
    if name == "Proliferation_down": sc_sig = -sc_sig
    bulk[name] = sc_sig
B = pd.DataFrame(bulk); B["condition"] = [colcond[c] for c in B.index]
B.to_csv(f"{OUT}/bulk_scores.csv")
# per-inducer: mean signature score sen vs its control
pairs = {"Replicative": "Young_Control", "Doxorubicin": "DMSO_Control", "IR_Senescence": "IR_Control"}
cross = {}
for name in sig:
    cross[name] = {ind: round(float(B[B.condition == ind][name].mean() - B[B.condition == ctrl][name].mean()), 2)
                   for ind, ctrl in pairs.items()}
# inter-signature agreement (single-cell, per-cell scores)
corr = sc_df[list(sc_scores)].corr(method="spearman").round(2)
R = {"singlecell_AUROC_RAS_senescence": A,
     "bulk_cross_inducer_score_difference_sen_minus_control": cross,
     "inter_signature_spearman_singlecell": corr.to_dict(),
     "n_singlecell": {"growing": int((sc_df.sen == 0).sum()), "RIS": int((sc_df.sen == 1).sum())}}
json.dump(R, open(f"{OUT}/scores_summary.json", "w"), indent=1)
print("single-cell AUROC (RAS sen vs growing):", json.dumps(A, indent=1))
print("\ncross-inducer sen-minus-control (z):"); print(pd.DataFrame(cross).T.to_string())
print("\ninter-signature Spearman (single-cell):"); print(corr.to_string())
