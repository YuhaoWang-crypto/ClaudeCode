"""Frame the real finding: AMR gene detection has high PPV (gene present -> resistant, reliable rule-in)
but low sensitivity / NPV (gene absent != susceptible -> missed resistance, the clinical warning)."""
import os, json, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
A = os.environ["AMR_DATA"]; OUT = "results"
res = pd.read_csv(f"{OUT}/amr_concordance.csv")
pheno = pd.read_csv(f"{A}/pheno_sel.tsv", sep="\t"); geno = pd.read_csv(f"{A}/genotype.tsv", sep="\t")
g = geno[geno.source == "NDARO"].copy(); g["gene"] = g.gene.astype(str)
# NPV per drug: among genotype-predicted-susceptible, what fraction are actually resistant (missed)
import re
# reuse prediction: load from a quick recompute consistent with 01
S = json.load(open(f"{OUT}/summary.json"))
pd_rows = S["per_drug"]
df = pd.DataFrame(pd_rows)
df["missed_resistance_rate_VME"] = df["very_major_error_missedR"]
# known-answer: carbapenemase -> meropenem
carb = set(g[g.gene.str.contains("KPC|NDM|VIM|IMP|OXA-48|OXA-181|OXA-232", case=False, regex=True)].genome_id.astype(str))
mero = pheno[pheno.antibiotic == "meropenem"].copy(); mero["gid"] = mero.genome_id.astype(str)
known = {"blaKPC+_meropenem_R_rate": round(float(mero[mero.gid.isin(set(g[g.gene.str.contains("KPC",case=False)].genome_id.astype(str)))].R.mean()), 3),
         "carbapenemase+_R_rate": round(float(mero[mero.gid.isin(carb)].R.mean()), 3),
         "carbapenemase-_R_rate": round(float(mero[~mero.gid.isin(carb)].R.mean()), 3),
         "frac_of_R_without_carbapenemase": round(float((~mero[mero.R==1].gid.isin(carb)).mean()), 3)}
S["carbapenem_known_answer"] = known; json.dump(S, open(f"{OUT}/summary.json", "w"), indent=1)

fig, ax = plt.subplots(1, 2, figsize=(14, 5))
d = df.sort_values("PPV", ascending=True)
y = np.arange(len(d)); h = 0.4
ax[0].barh(y + h/2, d.PPV, h, color="#1a9850", label="PPV (gene present → resistant)")
ax[0].barh(y - h/2, d.sensitivity, h, color="#d73027", label="sensitivity (resistant genomes flagged)")
ax[0].set_yticks(y, d.drug, fontsize=8); ax[0].axvline(0.9, ls=":", c="#999", lw=.8)
ax[0].set(xlabel="fraction", title="A. Resistance-gene detection:\nreliable rule-IN (green PPV), unreliable rule-OUT (red sensitivity)", xlim=(0,1)); ax[0].legend(fontsize=8, loc="lower right")
# B: carbapenem known-answer
labels = ["blaKPC+", "any\ncarbapenemase+", "carbapenemase−"]
vals = [known["blaKPC+_meropenem_R_rate"], known["carbapenemase+_R_rate"], known["carbapenemase-_R_rate"]]
ax[1].bar(range(3), vals, color=["#1a9850","#1a9850","#d73027"])
for i,v in enumerate(vals): ax[1].text(i, v+0.02, f"{v:.0%}", ha="center", fontsize=11)
ax[1].set_xticks(range(3), labels, fontsize=9); ax[1].set(ylabel="meropenem-resistant rate", ylim=(0,1.05),
    title=f"B. Known-answer: carbapenemase → carbapenem resistance\nbut {known['frac_of_R_without_carbapenemase']:.0%} of resistant isolates carry NO carbapenemase\n(resistance via porin loss / efflux — missed by gene panels)")
fig.tight_layout(); fig.savefig(f"{OUT}/fig1_amr.png", dpi=130)
print(json.dumps({"per_drug_PPV_vs_sensitivity": df[["drug","PPV","sensitivity","specificity","very_major_error_missedR","AUROC"]].round(3).to_dict("records"), "carbapenem": known}, indent=1))
