"""Does colocalization uniquely identify the validated causal gene at a locus?"""
import json, numpy as np, pandas as pd
OUT = "results"
R = pd.read_csv(f"{OUT}/coloc_all.csv")
loci = json.load(open(f"{OUT}/causal_loci.json"))
THR = 0.8
S = {"n_loci": len(loci), "n_coloc_tests": int(len(R)), "per_locus": {}, "threshold_PP_H4": THR}
rows = []
for g in loci:
    d = R[R.locus == g]
    if len(d) == 0: continue
    # take best tissue per gene
    best = d.sort_values("PP_H4", ascending=False).drop_duplicates("gene_id")
    coloc_genes = best[best.PP_H4 > THR]
    causal = best[best.is_causal == 1]
    causal_pp = float(causal.PP_H4.max()) if len(causal) else 0.0
    causal_colocs = bool(causal_pp > THR)
    # rank of causal gene by PP_H4
    ranked = best.sort_values("PP_H4", ascending=False).reset_index(drop=True)
    rank = int(ranked.index[ranked.is_causal == 1][0]) + 1 if (ranked.is_causal == 1).any() else None
    rows.append({"locus": g, "n_genes_tested": int(best.gene_id.nunique()),
                 "n_genes_colocalizing": int(len(coloc_genes)), "causal_PP_H4": round(causal_pp, 3),
                 "causal_colocalizes": causal_colocs, "causal_is_top": rank == 1 if rank else False,
                 "causal_rank": rank, "top_gene": ranked.gene_symbol.iloc[0], "top_PP_H4": round(float(ranked.PP_H4.iloc[0]), 3)})
pl = pd.DataFrame(rows); pl.to_csv(f"{OUT}/per_locus.csv", index=False)
coloc_loci = pl[pl.n_genes_colocalizing > 0]
S["summary"] = {
    "loci_with_any_coloc": int((pl.n_genes_colocalizing > 0).sum()),
    "mean_genes_colocalizing_per_locus": round(float(coloc_loci.n_genes_colocalizing.mean()), 2),
    "loci_where_causal_colocalizes": int(pl.causal_colocalizes.sum()),
    "loci_where_causal_is_TOP_gene": int(pl.causal_is_top.sum()),
    "loci_where_coloc_UNIQUELY_identifies_causal(1 gene & its the causal)": int(((pl.n_genes_colocalizing == 1) & pl.causal_colocalizes).sum()),
    "frac_coloc_loci_with_multiple_genes": round(float((coloc_loci.n_genes_colocalizing > 1).mean()), 2),
    "causal_gene_median_rank": float(pl.causal_rank.dropna().median())}
json.dump(S, open(f"{OUT}/summary.json", "w"), indent=1)
print(pl.to_string(index=False)); print("\n", json.dumps(S["summary"], indent=1))
