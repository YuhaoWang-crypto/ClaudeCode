"""Response-eQTL (reQTL) discovery: compare eGenes across stimulus conditions (same donors, same cell).
A baseline-only (resting) eQTL scan misses variants whose regulatory effect appears only under challenge."""
import os, json, numpy as np, pandas as pd
from statsmodels.stats.multitest import multipletests
R = os.environ["RQ_DATA"]; OUT = "results"
conds = {"naive": "mono_naive", "LPS": "mono_LPS", "R848": "mono_R848", "IAV": "mono_IAV", "Pam3": "mono_Pam3",
         "mac_naive": "mac_naive", "mac_IFNg": "mac_IFNg"}
def load(k):
    d = pd.read_csv(f"{R}/{k}.permuted.tsv.gz", sep="\t").dropna(subset=["p_beta"])
    d = d.drop_duplicates("molecular_trait_id")
    d["q"] = multipletests(d.p_beta, method="fdr_bh")[1]
    return d.set_index("molecular_trait_id")
D = {k: load(k) for k in conds}
def egenes(k, q=0.05): return set(D[k].index[D[k].q < q])
pairs = [("naive", "LPS"), ("naive", "R848"), ("naive", "IAV"), ("naive", "Pam3"), ("mac_naive", "mac_IFNg")]
S = {"n_genes_tested": {k: int(len(D[k])) for k in conds}, "n_eGenes_q0.05": {k: len(egenes(k)) for k in conds}, "comparisons": {}}
for base, stim in pairs:
    eb, es = egenes(base), egenes(stim)
    gained = es - eb; shared = es & eb; lost = eb - es
    # among shared, do lead effect sizes differ? (context-dependent magnitude)
    sh = list(shared); bb = D[base].loc[sh, "beta"].abs(); sb = D[stim].loc[sh, "beta"].abs()
    S["comparisons"][f"{base}->{stim}"] = {
        "n_eGenes_baseline": len(eb), "n_eGenes_stimulated": len(es),
        "response_gained_eGenes": len(gained), "pct_of_stim_eGenes_baseline_invisible": round(100 * len(gained) / len(es), 1),
        "shared_eGenes": len(shared), "baseline_only_lost": len(lost),
        "median_abs_beta_shared_baseline": round(float(bb.median()), 3), "median_abs_beta_shared_stim": round(float(sb.median()), 3)}
    # save gained eGene list with stim lead variant for the tabix confirmation stage
    g = D[stim].loc[sorted(gained), ["variant", "chromosome", "position", "pvalue", "beta", "q"]].reset_index()
    g.to_csv(f"{OUT}/gained_{base}_to_{stim}.csv", index=False)
json.dump(S, open(f"{OUT}/egene_overlap.json", "w"), indent=1)
print(json.dumps(S, indent=1))
