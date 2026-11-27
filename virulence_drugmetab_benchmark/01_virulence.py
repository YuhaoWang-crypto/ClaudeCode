"""Virulence known-answer: do invasive Klebsiella isolates carry hypervirulence (hvKp) loci more than
carriage/environmental isolates? Markers: aerobactin (iuc), yersiniabactin (ybt), salmochelin (iro),
colibactin (clb), rmpA regulator. (Russo 2018: aerobactin is the strongest hvKp marker.)"""
import os, json, re, numpy as np, pandas as pd
from scipy.stats import fisher_exact
V = os.environ["VIRU_DIR"]; OUT = "results"
src = pd.read_csv(f"{V}/source_sel.tsv", sep="\t"); src["genome_id"] = src.genome_id.astype(str)
vir = pd.read_csv(f"{V}/virulence.tsv", sep="\t"); vir["genome_id"] = vir.genome_id.astype(str)
vir["g"] = vir["gene"].astype(str).str.lower(); vir["p"] = vir["product"].astype(str).str.lower()
loci = {"aerobactin_iuc": r"\biuc[abcd]\b|\biuta\b|aerobactin",
        "yersiniabactin_ybt": r"\bybt|\birp[12]\b|\bfyua\b|yersiniabactin",
        "salmochelin_iro": r"\biro[bcdn]\b|salmochelin",
        "colibactin_clb": r"\bclb|colibactin|pks island",
        "rmpA_regulator": r"\brmpa", }
flags = pd.DataFrame(index=sorted(set(vir.genome_id)))
for name, pat in loci.items():
    hit = vir[vir.g.str.contains(pat, regex=True) | vir.p.str.contains(pat, regex=True)].genome_id.unique()
    flags[name] = flags.index.isin(hit).astype(int)
flags["n_virulence_genes"] = vir.groupby("genome_id").size().reindex(flags.index).fillna(0).astype(int)
d = src.set_index("genome_id").join(flags, how="inner")
d["invasive"] = (d.cat == "invasive").astype(int)
d.to_csv(f"{OUT}/virulence_flags.csv")
S = {"n_genomes": int(len(d)), "by_category_n": d.cat.value_counts().to_dict(), "marker_prevalence_by_category": {}, "invasive_vs_rest_enrichment": {}}
for m in loci:
    S["marker_prevalence_by_category"][m] = {c: round(float(d[d.cat == c][m].mean()), 3) for c in ["invasive", "mucosal", "carriage", "environmental"]}
    a = int(((d.invasive == 1) & (d[m] == 1)).sum()); b = int(((d.invasive == 1) & (d[m] == 0)).sum())
    c_ = int(((d.invasive == 0) & (d[m] == 1)).sum()); dd = int(((d.invasive == 0) & (d[m] == 0)).sum())
    orr, p = fisher_exact([[a, b], [c_, dd]], alternative="greater")
    S["invasive_vs_rest_enrichment"][m] = {"odds_ratio": round(float(orr), 2), "p": float(p),
        "prev_invasive": round(a / (a + b), 3), "prev_noninvasive": round(c_ / (c_ + dd), 3)}
# convergence (hvKp score = number of loci) vs category
d["hvkp_score"] = d[list(loci)].sum(1)
S["mean_hvkp_loci_by_category"] = {c: round(float(d[d.cat == c].hvkp_score.mean()), 2) for c in ["invasive", "mucosal", "carriage", "environmental"]}
S["mean_total_virulence_genes_by_category"] = {c: int(d[d.cat == c].n_virulence_genes.mean()) for c in ["invasive", "mucosal", "carriage", "environmental"]}
json.dump(S, open(f"{OUT}/virulence_summary.json", "w"), indent=1)
print(json.dumps(S, indent=1))
