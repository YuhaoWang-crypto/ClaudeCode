"""Proper analysis: within cell type, and per dataset, to expose batch confounding.

Two distinct questions the recorder framework asks, which need different readouts:
  (1) Is the marker INDUCED in its source cell?  -> per-cell expression within
      a fixed cell type (KIM-1, NGAL).
  (2) Is the source cell POPULATION lost?        -> cell-type proportion
      (uromodulin's "parenchymal reserve" claim is a cell-MASS claim, so
      per-cell expression is the wrong readout for it).
"""
import warnings, json
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, cellxgene_census

GENES = ["UMOD", "LCN2", "HAVCR1", "EGF", "COL3A1", "SLC34A1"]

with cellxgene_census.open_soma(census_version="2025-11-08") as census:
    adata = cellxgene_census.get_anndata(
        census=census, organism="Homo sapiens", measurement_name="RNA",
        X_name="normalized",
        obs_value_filter=("tissue_general=='kidney' and is_primary_data==True and "
                          "disease in ['normal','chronic kidney disease',"
                          "'acute kidney failure']"),
        var_value_filter=f"feature_name in {GENES}",
        obs_column_names=["disease", "cell_type", "dataset_id", "assay"],
    )

X = adata.X
X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
obs = adata.obs.reset_index(drop=True)
gname = list(adata.var["feature_name"])
df = pd.DataFrame(X, columns=gname)
df["disease"] = obs["disease"].values
df["cell_type"] = obs["cell_type"].values
df["dataset_id"] = obs["dataset_id"].values

TAL = "kidney loop of Henle thick ascending limb epithelial cell"
PT = "epithelial cell of proximal tubule"

print("=" * 76)
print("Q1  Marker INDUCTION inside a fixed cell type (% of cells positive)")
print("=" * 76)
for ct in [PT, TAL]:
    sub = df[df.cell_type == ct]
    if not len(sub):
        continue
    print(f"\ncell type: {ct}")
    print(f"  {'gene':<9}" + "".join(f"{d[:14]:>16}" for d in
                                     ["normal", "chronic kidney", "acute kidney"]))
    for g in GENES:
        line = f"  {g:<9}"
        for d in ["normal", "chronic kidney disease", "acute kidney failure"]:
            s = sub[sub.disease == d]
            line += f"{(100*(s[g]>0).mean() if len(s) else float('nan')):>15.1f}%"
        n = [len(sub[sub.disease == d]) for d in
             ["normal", "chronic kidney disease", "acute kidney failure"]]
        print(line + f"   n={n}")

print("\n" + "=" * 76)
print("Q1b PER-DATASET check for KIM-1/NGAL induction in proximal tubule")
print("    (does the effect hold inside every dataset, or is it batch?)")
print("=" * 76)
sub = df[df.cell_type == PT]
for g in ["HAVCR1", "LCN2"]:
    print(f"\n  {g}: % positive by dataset")
    t = sub.groupby(["disease", "dataset_id"], observed=True).apply(
        lambda s: pd.Series({"pct_pos": 100*(s[g] > 0).mean(), "n": len(s)}))
    print(t.round(2).to_string())

print("\n" + "=" * 76)
print("Q2  Source-cell POPULATION loss — the uromodulin 'reserve' claim")
print("    TAL cells as % of all kidney cells, PER DATASET")
print("=" * 76)
comp = (df.groupby(["disease", "dataset_id"], observed=True)
          .apply(lambda s: pd.Series({
              "TAL_pct": 100*(s.cell_type == TAL).mean(),
              "PT_pct": 100*(s.cell_type == PT).mean(),
              "n_cells": len(s)})))
print(comp.round(2).to_string())

df.to_parquet("census_kidney_cells.parquet")
print("\nsaved per-cell table")
