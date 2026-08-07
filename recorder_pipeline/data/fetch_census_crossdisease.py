"""SLE / neurodegeneration / atherosclerosis single-cell, healthy vs disease."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, cellxgene_census

OUT = "/home/user/ClaudeCode/recorder_pipeline/data"


def pull(census, obs_filter, genes, cols=("disease", "cell_type", "dataset_id")):
    a = cellxgene_census.get_anndata(
        census=census, organism="Homo sapiens", measurement_name="RNA",
        X_name="normalized", obs_value_filter=obs_filter,
        var_value_filter=f"feature_name in {genes}",
        obs_column_names=list(cols))
    X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
    df = pd.DataFrame(X, columns=list(a.var["feature_name"]))
    for c in cols:
        df[c] = a.obs[c].values
    return df


def summarise(df, genes, group_cols):
    rows = []
    for keys, s in df.groupby(group_cols, observed=True):
        if len(s) < 50:
            continue
        rec = dict(zip(group_cols, keys if isinstance(keys, tuple) else (keys,)))
        rec["n_cells"] = len(s)
        for g in genes:
            if g not in s:
                continue
            rec[f"{g}_pct"] = round(100 * (s[g] > 0).mean(), 3)
            rec[f"{g}_mean"] = float(s[g].mean())
        rows.append(rec)
    return pd.DataFrame(rows)


with cellxgene_census.open_soma(census_version="2025-11-08") as census:

    # ---------------- 1. SLE blood, one dataset with both arms -------------
    G = ["C4A", "C4B", "C3", "CR1", "C1QA", "PADI4", "IFI44L", "ISG15", "MS4A1"]
    sle = pull(census,
               "tissue_general=='blood' and is_primary_data==True and "
               "dataset_id=='218acb0f-9f2f-4f76-b90b-15a4b7c7f629' and "
               "disease in ['normal','systemic lupus erythematosus']", G)
    print("SLE cells:", len(sle), sle.disease.value_counts().to_dict())
    s = summarise(sle, G, ["cell_type", "disease"])
    s.to_csv(f"{OUT}/census_sle.tsv", sep="\t", index=False)
    print("SLE rows:", len(s))

    # ---------------- 2. brain, AD vs normal, 8 paired datasets ------------
    G2 = ["GFAP", "AQP4", "VIM", "S100B", "CAPN1", "CAPN2"]
    ds = ['cff99df2', 'ac0c6561', '85c60876', '6c600df6', '2727d83a',
          '203025fe', '0a2d7e87', '9813a1d4']
    obs = census["census_data"]["homo_sapiens"].obs
    full = obs.read(column_names=["dataset_id"],
                    value_filter=("tissue_general=='brain' and is_primary_data==True "
                                  "and disease in ['normal','Alzheimer disease']")
                    ).concat().to_pandas()
    ids = [d for d in full.dataset_id.unique() if d[:8] in ds]
    brain = pull(census,
                 "tissue_general=='brain' and is_primary_data==True and "
                 f"disease in ['normal','Alzheimer disease'] and dataset_id in {ids}",
                 G2)
    brain["ds"] = brain.dataset_id.str[:8]
    print("brain cells:", len(brain), brain.disease.value_counts().to_dict())
    astro = brain[brain.cell_type.str.contains("astrocyte", case=False, na=False)]
    print("astrocytes:", len(astro), astro.cell_type.unique()[:5])
    b = summarise(astro, G2, ["ds", "disease"])
    b.to_csv(f"{OUT}/census_brain_astro.tsv", sep="\t", index=False)

    # ---------------- 3. atherosclerosis: PRESENCE only, no control --------
    G3 = ["MPO", "APOA1", "CD68", "LYZ", "ELANE"]
    ath = pull(census,
               "tissue_general=='vasculature' and is_primary_data==True and "
               "disease=='atherosclerosis'", G3)
    print("athero cells:", len(ath))
    a = summarise(ath, G3, ["cell_type"])
    a.to_csv(f"{OUT}/census_athero.tsv", sep="\t", index=False)

print("done")
