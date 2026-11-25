"""cis-eQTL discovery with tensorQTL in EUR_train, then replication of each lead variant in people
the scan never saw: EUR_test (same ancestry) and YRI (different ancestry)."""
import json, numpy as np, pandas as pd, torch
from scipy import stats
from statsmodels.stats.multitest import multipletests
from tensorqtl import cis
from common import *

torch.set_num_threads(4)
G, bim, E, genes, sets = load()
genes = genes[genes.index.isin(E.columns)]
def prep(name):
    m = (sets.set == name).values; C = pcs(G[m]); Y = np.apply_along_axis(int_norm, 0, E.values[m])
    return m, C, Y
m, C, Y = prep("EUR_train"); ids = E.index[m]
keep = (G[m].mean(0) / 2).clip(0, 1); keep = np.minimum(keep, 1 - keep) >= 0.05          # MAF filter in training set only
gdf = pd.DataFrame(G[m][:, keep].T, index=bim.snp[keep], columns=ids)
vdf = bim[keep].set_index("snp")[["chr", "pos"]].rename(columns={"chr": "chrom"})
pdf = pd.DataFrame(Y.T, index=E.columns, columns=ids).loc[genes.index]
ppos = pd.DataFrame({"chr": genes.chr.astype(str), "pos": genes.tss.astype(int)}, index=genes.index)
cov = pd.DataFrame(C, index=ids, columns=["PC1", "PC2", "PC3"])
res = cis.map_cis(gdf, vdf, pdf, ppos, cov, window=500_000, nperm=1000, seed=0)
res["qval"] = multipletests(res.pval_beta, method="fdr_bh")[1]
res.to_csv(f"{WORK}/cis_train.tsv", sep="\t")
eg = res[res.qval < 0.05].copy(); snp_ix = pd.Series(np.arange(len(bim)), index=bim.snp)

def assoc(name, gene, snp):
    mm, CC, YY = cache[name]; y = residualize(YY[:, E.columns.get_loc(gene)], CC); g = residualize(G[mm][:, snp_ix[snp]], CC)
    if g.std() < 1e-6: return np.nan, np.nan
    r = np.corrcoef(g, y)[0, 1]; n = len(y); t = r * np.sqrt((n - 5) / (1 - r * r))
    return r * y.std() / g.std(), 2 * stats.t.sf(abs(t), n - 5)
cache = {k: prep(k) for k in ["EUR_test", "YRI"]}
rows = []
for gene, r in eg.iterrows():
    row = {"gene": gene, "symbol": genes.symbol[gene], "lead": r.variant_id, "slope_train": r.slope, "p_train": r.pval_nominal}
    for k in cache: row[f"slope_{k}"], row[f"p_{k}"] = assoc(k, gene, r.variant_id)
    rows.append(row)
rep = pd.DataFrame(rows); rep.to_csv(f"{OUT}/cis_eqtl_replication.csv", index=False)
S = {"genes_tested": int(len(res)), "n_train": int(m.sum()), "eGenes_FDR5pct": int(len(eg))}
for k in cache:
    d = rep.dropna(subset=[f"p_{k}"]); same = np.sign(d.slope_train) == np.sign(d[f"slope_{k}"])
    S[k] = {"n_people": int(cache[k][0].sum()), "lead_variants_testable": int(len(d)),
            "lead_variants_monomorphic_or_missing": int(len(rep) - len(d)),
            "same_direction": round(float(same.mean()), 3),
            "replicated_p<0.05_same_direction": round(float(((d[f"p_{k}"] < .05) & same).mean()), 3),
            "replicated_as_fraction_of_all_eGenes": round(float(((d[f"p_{k}"] < .05) & same).sum() / len(rep)), 3),
            "median_abs_effect_ratio_vs_train": round(float((d[f"slope_{k}"].abs() / d.slope_train.abs()).median()), 3),
            "effect_size_correlation": round(float(np.corrcoef(d.slope_train, d[f"slope_{k}"])[0, 1]), 3)}
json.dump(S, open(f"{OUT}/cis_eqtl.json", "w"), indent=1); print(json.dumps(S, indent=1))
