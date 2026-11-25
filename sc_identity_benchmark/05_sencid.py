"""SenCID pretrained SVMs on an independent senescence dataset with experimental ground truth:
GSE115301 (Teo et al. 2019), IMR90 fibroblasts, 10x: growing vs RAS-induced senescence (RIS).

Deviations from the official pipeline (stated in the report): no DCA denoising (needs TF<2.5) and
no SID recommender (its weights are not shipped in the repository), so SID1-4 are all reported.
"""
import os, sys, numpy as np, pandas as pd, joblib, warnings
from sklearn.preprocessing import StandardScaler
from common import *
warnings.filterwarnings("ignore")
SEN = os.environ["SENCID_DIR"]; D = os.environ["SEN_DATA"]
sys.path.insert(0, SEN); from Pred import GetFeatures

counts = pd.read_csv(f"{D}/GSE115301_Growing_Sen_10x_count.txt.gz", sep="\t", index_col=0)
meta = pd.read_csv(f"{D}/GSE115301_Growing_Sen_10x_metadata.txt.gz", sep="\t", index_col=0).loc[counts.columns]
genes = pd.read_csv(f"{SEN}/resource/seneset.txt", header=None)[0]
sub = counts.reindex(genes).fillna(0.01)                 # official: absent genes get 0.01
cpm = sub / sub.sum(0) * 1e6
z = pd.DataFrame(StandardScaler().fit_transform(np.log2(cpm + 1)), index=genes, columns=counts.columns).T
out = meta.copy()
out["n_counts"] = counts.sum(0); out["sen_genes_found"] = int(genes.isin(counts.index).sum())
for i in [1, 2, 3, 4]:
    x = z.loc[:, GetFeatures(i)].fillna(0)
    out[f"SID{i}"] = joblib.load(f"{SEN}/model/SID{i}_L.pkl").predict_proba(x.values)[:, 1]
out["SID_mean"] = out[[f"SID{i}" for i in range(1, 5)]].mean(1)
# transparent baselines: single marker and proliferation-loss
lc = np.log1p(counts / counts.sum(0) * 1e4)
out["CDKN1A"] = lc.loc["CDKN1A"]; out["CDKN2A"] = lc.loc["CDKN2A"]
out["neg_MKI67_TOP2A"] = -(lc.loc["MKI67"] + lc.loc["TOP2A"])
out.to_csv(f"{OUT}/sencid_GSE115301.csv")
print(out.groupby(["Condition1", "Condition2"]).mean(numeric_only=True).round(3).T)
