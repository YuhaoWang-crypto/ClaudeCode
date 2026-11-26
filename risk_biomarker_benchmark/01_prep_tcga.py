"""Build per-cancer expression + clinical tables from the TCGA pan-cancer Xena matrix."""
import os, numpy as np, pandas as pd
T = os.environ["TB_DATA"]
surv = pd.read_csv(f"{T}/survival.tsv", sep="\t", low_memory=False)
surv["stage"] = surv["ajcc_pathologic_tumor_stage"].astype(str).str.extract(r"Stage\s*(IV|III|II|I)", expand=False)
stage_num = {"I": 1, "II": 2, "III": 3, "IV": 4}
surv["stage_num"] = surv["stage"].map(stage_num)
surv["age"] = pd.to_numeric(surv["age_at_initial_pathologic_diagnosis"], errors="coerce")
surv["OS"] = pd.to_numeric(surv["OS"], errors="coerce"); surv["OS.time"] = pd.to_numeric(surv["OS.time"], errors="coerce")
surv = surv[surv["sample"].str.endswith("-01")]                       # primary tumours only
clin = surv[["sample", "cancer type abbreviation", "age", "gender", "stage_num", "OS", "OS.time"]].rename(
    columns={"cancer type abbreviation": "cancer"})
clin = clin.dropna(subset=["OS", "OS.time"]); clin = clin[clin["OS.time"] > 0]

cohorts = ["LUAD", "LUSC", "KIRC", "BRCA"]
hdr = pd.read_csv(f"{T}/exprPANCAN.gz", sep="\t", nrows=0).columns.tolist()
hset = set(hdr)
want = clin[clin.cancer.isin(cohorts)]
cols = ["sample"] + [s for s in want["sample"] if s in hset]
expr = pd.read_csv(f"{T}/exprPANCAN.gz", sep="\t", usecols=cols).set_index("sample")
expr = expr.dropna(how="any")                                         # genes measured in all kept samples
expr = expr.loc[expr.var(axis=1).sort_values(ascending=False).index[:8000]]   # 8000 most variable genes
expr.to_parquet(f"{T}/expr_sub.parquet")
clin[clin["sample"].isin(expr.columns)].to_parquet(f"{T}/clin_sub.parquet")
print("genes", expr.shape[0], "samples", expr.shape[1])
for c in cohorts:
    s = clin[(clin.cancer == c) & clin["sample"].isin(expr.columns)]
    print(f"{c}: n={len(s)}, deaths={int(s.OS.sum())}, median OS days={s['OS.time'].median():.0f}, stage known={s.stage_num.notna().mean():.2f}, age known={s.age.notna().mean():.2f}")
