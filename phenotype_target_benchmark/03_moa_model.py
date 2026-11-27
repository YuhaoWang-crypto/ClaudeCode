"""Test 2: MoA connectivity - does nearest drug by sensitivity profile share mechanism of action?
   Test 3: model matching + cross-lineage consistency - does a target's most-dependent line = most
   drug-sensitive, and does a drug's sensitivity transfer across lineages?"""
import os, json, numpy as np, pandas as pd
from scipy.stats import rankdata, spearmanr
P = os.environ["PB_DATA"]; W = os.environ["PB_WORK"]; OUT = "results"
auc = pd.read_csv(f"{W}/auc.csv", index_col=0); meta = pd.read_csv(f"{W}/meta.csv", index_col=0)
cr = pd.read_csv(f"{W}/cr_sub.csv", index_col=0); model = pd.read_csv(f"{P}/Model.csv")
lin = model.set_index("ModelID")["OncotreeLineage"]
sens = -auc.values; drugs = auc.index.values; lines = auc.columns.values
# ---- Test 2: MoA connectivity (drug-drug similarity by sensitivity profile)
Rs = np.vstack([rankdata(np.where(np.isnan(r), np.nanmedian(r), r)) for r in sens]).astype(float)
Rs = (Rs - Rs.mean(1, keepdims=True)) / (Rs.std(1, keepdims=True) + 1e-9)
sim = (Rs @ Rs.T) / Rs.shape[1]; np.fill_diagonal(sim, -9)
moa = meta.moa.reindex(drugs).values
has = pd.Series(moa).notna().values & (pd.Series(moa).map(pd.Series(moa).value_counts()) >= 3).values
nn = sim.argmax(1); same = []
for i in np.where(has)[0]:
    if has[nn[i]] and not pd.isna(moa[nn[i]]): same.append(moa[i] == moa[nn[i]])
# null: random pairing
rng = np.random.default_rng(0); idx = np.where(has)[0]
null = np.mean([moa[i] == moa[rng.choice(idx)] for i in idx for _ in range(3)])
S = {"MoA_connectivity": {"n_drugs_with_moa": int(has.sum()),
     "nearest_neighbour_shares_MoA": round(float(np.mean(same)), 3),
     "random_baseline": round(float(null), 3), "n_moa_classes": int(pd.Series(moa[has]).nunique())}}
# top-5 neighbour MoA recovery
top5 = []
for i in np.where(has)[0]:
    nbrs = np.argsort(-sim[i])[:5]; top5.append(any(moa[j] == moa[i] for j in nbrs if has[j]))
S["MoA_connectivity"]["any_of_top5_shares_MoA"] = round(float(np.mean(top5)), 3)

# ---- Test 3: model matching + cross-lineage
ri = pd.read_csv(f"{OUT}/reverse_infer.csv")
ontgt = ri[ri.target_rank <= 300].drug.values            # the subset where phenotype tracks the target
mm = []
for drug in ontgt:
    tgt = meta.target[drug]; s = pd.Series(sens[list(drugs).index(drug)], index=lines)
    dep = cr[tgt].reindex(lines)                           # gene effect (negative = dependent)
    m = s.notna() & dep.notna()
    if m.sum() < 100: continue
    # most dependent 15% lines: are they more sensitive than the rest?
    d = dep[m]; thr = d.quantile(0.15); topdep = d <= thr
    mm.append({"drug": drug, "target": tgt, "sens_dep_spearman": spearmanr(s[m], -dep[m])[0],
               "sens_in_dependent_models": s[m][topdep].mean(), "sens_in_rest": s[m][~topdep].mean()})
mm = pd.DataFrame(mm)
S["model_matching_on_target_subset"] = {"n_drugs": len(mm),
    "mean_sens_dependency_spearman": round(float(mm.sens_dep_spearman.mean()), 3),
    "dependent_models_more_sensitive_frac": round(float((mm.sens_in_dependent_models > mm.sens_in_rest).mean()), 3)}
# cross-lineage inconsistency: for each drug, variance of median sensitivity across lineages, and
# whether the most-sensitive lineage's sensitivity predicts another lineage (rank corr of lineage medians across drugs)
auc_f = auc.copy(); ln = lin.reindex(lines).values
big = pd.Series(ln).value_counts(); keep_lin = big[big >= 12].index.tolist()
linmed = {}
for L in keep_lin:
    cols = lines[pd.Series(ln).values == L]; linmed[L] = auc[cols].median(1)
LM = pd.DataFrame(linmed)                                  # drugs x lineages median AUC
cross = LM.corr(method="spearman").where(~np.eye(len(keep_lin), dtype=bool)).stack()
S["cross_lineage_consistency"] = {"n_lineages": len(keep_lin),
    "median_cross_lineage_sensitivity_corr": round(float(cross.median()), 3),
    "note": "Spearman of drug-sensitivity rankings between lineages; 1.0 = a drug's potency transfers across tissue models"}
json.dump(S, open(f"{OUT}/moa_model_summary.json", "w"), indent=1)
mm.to_csv(f"{OUT}/model_matching.csv", index=False)
print(json.dumps(S, indent=1))
