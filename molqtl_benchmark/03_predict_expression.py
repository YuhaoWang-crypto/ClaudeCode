"""Genotype -> expression (PrediXcan-style), scored against measured expression.
 (a) published GTEx v8 elastic-net weights for LCLs (PredictDB), applied unchanged to GEUVADIS people
     (none of whom are GTEx donors);
 (b) elastic net trained here on EUR_train, applied to EUR_test and YRI."""
import json, sqlite3, warnings, numpy as np, pandas as pd
from pyliftover import LiftOver
from sklearn.linear_model import ElasticNetCV
from joblib import Parallel, delayed
from common import *
warnings.filterwarnings("ignore")
G, bim, E, genes, sets = load()
key = pd.Series(np.arange(len(bim)), index=bim.chr + ":" + bim.pos.astype(str))
groups = {"EUR_all": (sets.set != "YRI").values, "EUR_test": (sets.set == "EUR_test").values, "YRI": (sets.set == "YRI").values}
gid = pd.Series(E.columns, index=[c.split(".")[0] for c in E.columns])

# ---------- (a) pretrained weights
con = sqlite3.connect(f"{DATA}/en_Cells_EBV-transformed_lymphocytes.db")
w = pd.read_sql("select * from weights", con); ex = pd.read_sql("select * from extra", con).set_index("gene")
w[["c", "p38", "ref", "alt", "_"]] = w.varID.str.split("_", expand=True)
w = w[w.c.isin(["chr21", "chr22"])].copy()
lo = LiftOver(f"{DATA}/hg38ToHg19.over.chain.gz")
def lift(c, p):
    r = lo.convert_coordinate(c, int(p) - 1); return r[0][1] + 1 if r and r[0][0] == c else -1
w["p19"] = [lift(c, p) for c, p in zip(w.c, w.p38)]
w["ix"] = (w.c.str[3:] + ":" + w.p19.astype(str)).map(key)
rows = []
for gene, d in w.groupby("gene"):
    g0 = gene.split(".")[0]
    if g0 not in gid.index: continue
    n_model = len(d); d = d.dropna(subset=["ix"]); ix = d.ix.astype(int).values
    cnt, oth = bim.a_counted.values[ix], bim.a_other.values[ix]
    same = (d.eff_allele.values == cnt) & (d.ref_allele.values == oth); flip = (d.eff_allele.values == oth) & (d.ref_allele.values == cnt)
    ok = same | flip
    if ok.sum() == 0: continue
    dose = G[:, ix[ok]]; dose = np.where(flip[ok], 2 - dose, dose); pred = dose @ d.weight.values[ok]
    row = {"gene": g0, "symbol": ex.genename[gene], "snps_in_model": n_model, "snps_found": int(ok.sum()),
           "weight_fraction_found": float(np.abs(d.weight.values[ok]).sum() / np.abs(w.weight[w.gene == gene]).sum()),
           "gtex_cv_R2": ex.cv_R2_avg[gene], "gtex_cv_rho": ex.rho_avg[gene]}
    y = E[gid[g0]].values
    for k, m in groups.items(): row[f"r_{k}"] = np.corrcoef(pred[m], y[m])[0, 1] if pred[m].std() > 0 else np.nan
    rows.append(row)
pre = pd.DataFrame(rows); pre.to_csv(f"{OUT}/predict_pretrained.csv", index=False)

# ---------- (b) trained here on EUR_train
tr = (sets.set == "EUR_train").values; maf = G[tr].mean(0) / 2; maf_ok = np.minimum(maf, 1 - maf) >= 0.05
egenes = set(pd.read_csv(f"{OUT}/cis_eqtl_replication.csv").gene)
def fit(gene):
    c, tss = str(genes.chr[gene]), genes.tss[gene]
    ix = np.where((bim.chr.values == c) & (np.abs(bim.pos.values - tss) <= 500_000) & maf_ok)[0]
    if len(ix) < 10: return None
    X, y = G[:, ix], E[gene].values
    X = (X - X[tr].mean(0)) / X[tr].std(0)
    m = ElasticNetCV(l1_ratio=0.5, alphas=20, cv=5, max_iter=1000, tol=1e-3, eps=1e-2).fit(X[tr], y[tr]); p = m.predict(X)
    out = {"gene": gene.split(".")[0], "eGene_in_training": gene in egenes, "snps_selected": int((m.coef_ != 0).sum())}
    for k in ["EUR_test", "YRI"]: out[f"r_{k}"] = np.corrcoef(p[groups[k]], y[groups[k]])[0, 1] if p[groups[k]].std() > 0 else 0.0
    return out
own = pd.DataFrame([r for r in Parallel(n_jobs=3)(delayed(fit)(g) for g in genes.index if g in E.columns) if r])
own.to_csv(f"{OUT}/predict_selftrained.csv", index=False)

def summ(v, n): v = np.asarray(v, float); v = np.nan_to_num(v); return {"n_genes": int(len(v)), "median_r": round(float(np.median(v)), 3),
    "mean_r2_signed": round(float(np.mean(np.sign(v) * v * v)), 3), "frac_r>0.1": round(float((v > .1).mean()), 3),
    "frac_r>0.3": round(float((v > .3).mean()), 3), "frac_r<=0": round(float((v <= 0).mean()), 3), "n_people": int(n)}
S = {"pretrained_GTEx_LCL_weights": {"genes_with_model_and_measured": int(len(pre)),
        "median_fraction_of_model_snps_found": round(float((pre.snps_found / pre.snps_in_model).median()), 3),
        "gtex_reported_median_cv_rho": round(float(pre.gtex_cv_rho.median()), 3),
        **{k: summ(pre[f"r_{k}"], m.sum()) for k, m in groups.items()},
        "spearman_reported_rho_vs_observed_r_EUR": round(float(pre[["gtex_cv_rho", "r_EUR_all"]].corr(method="spearman").iloc[0, 1]), 3)},
     "trained_here_on_EUR_train": {
        "all_genes": {k: summ(own[f"r_{k}"], groups[k].sum()) for k in ["EUR_test", "YRI"]},
        "training_eGenes_only": {k: summ(own.loc[own.eGene_in_training, f"r_{k}"], groups[k].sum()) for k in ["EUR_test", "YRI"]}}}
both = pre.merge(own, on="gene", suffixes=("_pre", "_own"))
S["same_genes_pretrained_vs_trained_here_EUR_test"] = {"n_genes": int(len(both)),
    "median_r_pretrained": round(float(both.r_EUR_test_pre.median()), 3), "median_r_trained_here": round(float(both.r_EUR_test_own.median()), 3)}
json.dump(S, open(f"{OUT}/predict_expression.json", "w"), indent=1); print(json.dumps(S, indent=1))
