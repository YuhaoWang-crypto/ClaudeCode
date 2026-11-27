"""Confirm response-gained eGenes are true reQTLs, not threshold/power artifacts:
look up each stimulated lead variant in the NAIVE full summary stats (same gene) and check it is
genuinely null at baseline. A real reQTL: strong in stimulus, absent at baseline."""
import os, sys, subprocess, numpy as np, pandas as pd
from concurrent.futures import ThreadPoolExecutor
OUT = "results"
base_url = {"naive": "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000024/QTD000409/QTD000409.all.tsv.gz",
            "mac_naive": "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000001/QTD000001/QTD000001.all.tsv.gz"}
comp = sys.argv[1]; base = comp.split("_to_")[0] if "_to_" in comp else "naive"
g = pd.read_csv(f"{OUT}/gained_{comp}.csv")
g = g.sample(min(250, len(g)), random_state=0)
URL = base_url[base]
def lookup(row):
    c, p, gene, var = str(row.chromosome), int(row.position), row.molecular_trait_id, row.variant
    try:
        out = subprocess.run(["tabix", URL, f"{c}:{p}-{p}"], capture_output=True, text=True, timeout=90).stdout
    except Exception: return None
    for line in out.strip().split("\n"):
        if not line: continue
        f = line.split("\t")
        # columns: molecular_trait_id ... variant ... pvalue beta ... (match gene + variant)
        if row.molecular_trait_id in line and (var in line):
            # find pvalue/beta by header positions: use known eQTL Catalogue all.tsv columns
            return {"gene": gene, "variant": var, "beta_stim": row.beta, "p_stim": row.pvalue,
                    "beta_naive": float(f[9]), "p_naive": float(f[8])}
    return {"gene": gene, "variant": var, "beta_stim": row.beta, "p_stim": row.pvalue, "beta_naive": np.nan, "p_naive": np.nan}
rows = [r for r in ThreadPoolExecutor(max_workers=16).map(lookup, [r for _, r in g.iterrows()]) if r]
d = pd.DataFrame(rows).dropna(subset=["p_naive"])
d["abs_beta_stim"] = d.beta_stim.abs(); d["abs_beta_naive"] = d.beta_naive.abs()
d.to_csv(f"{OUT}/confirm_{comp}.csv", index=False)
res = {"comparison": comp, "n_checked": int(len(d)),
       "median_p_naive": float(d.p_naive.median()), "median_abs_beta_stim": round(float(d.abs_beta_stim.median()), 3),
       "median_abs_beta_naive": round(float(d.abs_beta_naive.median()), 3),
       "frac_naive_nonsignificant_p>1e-4": round(float((d.p_naive > 1e-4).mean()), 3),
       "frac_confirmed_reQTL_naiveP>0.01_and_betadrop2x": round(float(((d.p_naive > 0.01) & (d.abs_beta_naive < d.abs_beta_stim / 2)).mean()), 3),
       "median_effect_fold_change_stim_over_naive": round(float((d.abs_beta_stim / d.abs_beta_naive.clip(lower=0.01)).median()), 2)}
import json; json.dump(res, open(f"{OUT}/confirm_{comp}.json", "w"), indent=1)
print(json.dumps(res, indent=1))
