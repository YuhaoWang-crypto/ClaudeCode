"""Build aligned matrices: drug sensitivity (PRISM AUC) and gene dependency (CRISPR) over shared cell lines."""
import os, numpy as np, pandas as pd
P = os.environ["PB_DATA"]
prism = pd.read_csv(f"{P}/prism_drc.csv", usecols=["depmap_id", "auc", "name", "moa", "target"]).dropna(subset=["auc", "target"])
prism = prism[~prism.target.str.contains(",")]                       # single-gene target
cr = pd.read_csv(f"{P}/CRISPRGeneEffect.csv", index_col=0)
cr.columns = [c.split(" (")[0] for c in cr.columns]
lines = sorted(set(prism.depmap_id) & set(cr.index))
cr = cr.loc[lines]
# drug x cell line AUC matrix (median across screens/replicates)
auc = prism.pivot_table(index="name", columns="depmap_id", values="auc", aggfunc="median").reindex(columns=lines)
auc = auc[auc.notna().sum(1) >= 200]                                  # drugs measured in >=200 shared lines
meta = prism.drop_duplicates("name").set_index("name")[["target", "moa"]].reindex(auc.index)
meta = meta[meta.target.isin(cr.columns)]; auc = auc.loc[meta.index]
cr.to_parquet(f"{os.environ['PB_WORK']}/cr.parquet"); auc.to_parquet(f"{os.environ['PB_WORK']}/auc.parquet")
meta.to_parquet(f"{os.environ['PB_WORK']}/meta.parquet")
print("shared cell lines", len(lines), "| drugs (single target in CRISPR, >=200 lines)", len(auc), "| genes", cr.shape[1],
      "| unique targets", meta.target.nunique())
