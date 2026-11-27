"""Phenotype -> direct-target reverse inference, DepMap known-answer.
For each PRISM drug, correlate its sensitivity profile across cell lines with every gene's CRISPR
dependency profile, and ask where the drug's TRUE annotated target ranks among ~18k genes.
Boundary: phenotypic-profile similarity recovers pathway/MoA but often NOT the exact physical target."""
import os, json, numpy as np, pandas as pd
from scipy.stats import rankdata
P = os.environ["PB_DATA"]; OUT = "results"
prism = pd.read_csv(f"{P}/prism_drc.csv", usecols=["depmap_id", "auc", "name", "moa", "target"]).dropna(subset=["auc", "target"])
prism = prism[~prism.target.str.contains(",")]
cr = pd.read_csv(f"{P}/CRISPRGeneEffect.csv", index_col=0); cr.columns = [c.split(" (")[0] for c in cr.columns]
lines = sorted(set(prism.depmap_id) & set(cr.index)); cr = cr.loc[lines]
auc = prism.pivot_table(index="name", columns="depmap_id", values="auc", aggfunc="median").reindex(columns=lines)
auc = auc[auc.notna().sum(1) >= 250]
meta = prism.drop_duplicates("name").set_index("name")[["target", "moa"]].reindex(auc.index)
meta = meta[meta.target.isin(cr.columns)]; auc = auc.loc[meta.index]
genes = cr.columns.values; gidx = {g: i for i, g in enumerate(genes)}

# rank-transform gene dependency per gene (Spearman via Pearson on ranks); impute gene NaNs by column mean-rank
Rg = np.vstack([rankdata(cr[g].values) for g in genes]).T.astype(np.float64)   # lines x genes (ranked within gene over lines)
Rg = (Rg - Rg.mean(0)) / (Rg.std(0) + 1e-9)
rows = []
sens = -auc.values                                   # higher = more sensitive
for i, drug in enumerate(auc.index):
    s = sens[i]; m = ~np.isnan(s)
    rs = rankdata(s[m]); rs = (rs - rs.mean()) / (rs.std() + 1e-9)
    corr = (Rg[m].T @ rs) / m.sum()                  # corr(sensitivity, dependency) per gene
    order = np.argsort(-corr)                         # descending: on-target => dependent lines sensitive => positive
    tgt = meta.target[drug]; ti = gidx[tgt]
    rank = int(np.where(order == ti)[0][0]) + 1
    rows.append({"drug": drug, "target": tgt, "moa": meta.moa[drug], "n_lines": int(m.sum()),
                 "target_corr": round(float(corr[ti]), 3), "target_rank": rank,
                 "pct_genes_beaten": round(100 * (1 - rank / len(genes)), 1), "best_gene": genes[order[0]],
                 "best_corr": round(float(corr[order[0]]), 3)})
df = pd.DataFrame(rows); df.to_csv(f"{OUT}/reverse_infer.csv", index=False)
N = len(genes)
S = {"n_drugs": len(df), "n_genes": int(N), "n_cell_lines": len(lines),
     "target_rank_median": int(df.target_rank.median()), "target_rank_mean": int(df.target_rank.mean()),
     "frac_target_is_top1": round(float((df.target_rank == 1).mean()), 3),
     "frac_target_in_top10": round(float((df.target_rank <= 10).mean()), 3),
     "frac_target_in_top100": round(float((df.target_rank <= 100).mean()), 3),
     "frac_target_in_top1pct": round(float((df.target_rank <= N * 0.01).mean()), 3),
     "frac_positive_target_corr": round(float((df.target_corr > 0).mean()), 3),
     "random_expectation_median_rank": int(N / 2),
     "frac_target_corr_gt_0.2": round(float((df.target_corr > 0.2).mean()), 3)}
# by MoA class: are some mechanisms recoverable and others not?
df["ontarget_strong"] = df.target_rank <= 100
S["best_recovered_moas"] = df.groupby("moa").agg(n=("drug", "size"), median_rank=("target_rank", "median")).query("n>=5").sort_values("median_rank").head(8).reset_index().to_dict("records")
json.dump(S, open(f"{OUT}/reverse_infer_summary.json", "w"), indent=1)
np.save(f"{os.environ['PB_WORK']}/sens.npy", sens); auc.to_csv(f"{os.environ['PB_WORK']}/auc.csv")
meta.to_csv(f"{os.environ['PB_WORK']}/meta.csv"); cr.to_csv(f"{os.environ['PB_WORK']}/cr_sub.csv")
print(json.dumps(S, indent=1))
