"""Score the GWAS: genomic inflation, hit counts, and whether the lead signal sits on the known cloned gene."""
import os, json, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import chi2
W = os.environ["QTL_WORK"] + "/rice"
traits = open(f"{W}/traits.txt").read().split("\n")
# cloned genes reported for these traits in this panel (Zhao et al. 2011); approximate Mb positions
KNOWN = {"Amylose content": ("Waxy", 6, 1.77), "Seed length": ("GS3", 3, 16.73),
         "Seed width": ("qSW5/GW5", 5, 5.36), "Pericarp color": ("Rc", 7, 6.07)}
WIN = 0.5  # Mb
out = {}; fig, ax = plt.subplots(2, 4, figsize=(17, 6.5))
for i, t in enumerate(traits, 1):
    r = {}
    for m, col in [("naive_linear", "p_wald"), ("LMM_kinship", "p_wald")]:
        d = pd.read_csv(f"{W}/{'lm' if m == 'naive_linear' else 'lmm'}_{i}.assoc.txt", sep="\t")
        p = d[col].clip(lower=1e-300); bonf = 0.05 / len(d)
        e = {"n_snps": int(len(d)), "lambda_GC": round(float(np.median(chi2.isf(p, 1)) / chi2.ppf(.5, 1)), 2),
             "n_bonferroni_hits": int((p < bonf).sum()),
             "n_chromosomes_with_hits": int(d.loc[p < bonf, "chr"].nunique())}
        if t in KNOWN:
            g, c, mb = KNOWN[t]
            near = (d.chr == c) & ((d.ps / 1e6 - mb).abs() <= WIN)
            # rank of the known locus among independent 1-Mb bins ordered by their best p-value
            bins = d.assign(bin=d.chr.astype(str) + "_" + (d.ps // 1e6).astype(int).astype(str), p=p).groupby("bin").p.min().sort_values()
            lead = d.loc[p.idxmin()]
            e.update({"known_gene": g, "best_p_within_0.5Mb_of_gene": float(p[near].min()),
                      "gene_region_passes_bonferroni": bool(p[near].min() < bonf),
                      "rank_of_gene_bin_among_1Mb_bins": int(np.where(bins.index == f"{c}_{int(mb)}")[0][0]) + 1,
                      "lead_snp": f"chr{int(lead.chr)}:{lead.ps/1e6:.2f}Mb",
                      "lead_snp_is_at_gene": bool(lead.chr == c and abs(lead.ps / 1e6 - mb) <= WIN)})
            a = ax[0 if m == "naive_linear" else 1, list(KNOWN).index(t)]
            x = d.ps + d.chr.map(d.groupby("chr").ps.max().cumsum().shift(fill_value=0))
            a.scatter(x, -np.log10(p), s=2, c=np.where(d.chr % 2, "#4575b4", "#91bfdb"))
            a.scatter(x[near], -np.log10(p[near]), s=6, c="#d73027", label=f"±0.5 Mb of {g}")
            a.axhline(-np.log10(bonf), ls=":", c="k", lw=.8); a.set_xticks([])
            a.set_title(f"{t} — {m.replace('_', ' ')}\nλ={e['lambda_GC']}, Bonferroni hits={e['n_bonferroni_hits']}", fontsize=9)
            a.legend(fontsize=7, loc="upper right"); a.set_ylabel("-log10 p")
        r[m] = e
    out[t] = r
fig.tight_layout(); fig.savefig("results/fig2_rice_gwas.png", dpi=130)
json.dump(out, open("results/rice_gwas.json", "w"), indent=1); print(json.dumps(out, indent=1))
