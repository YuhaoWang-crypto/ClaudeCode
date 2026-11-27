"""Does the PAV flag mark a distinct class? Effect-size inflation + consequence breakdown, among
genome-wide-significant cis-pQTLs. (Epitope artifacts inflate apparent effect because the reagent
simply stops binding.)"""
import os, json, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT = "results"
d = pd.read_csv(f"{OUT}/cis_leads_annotated.tsv", sep="\t")
d["p"] = d.p.astype(float); d["abs_beta"] = d.beta.abs()
sig = d[d.p < 5e-8].copy()
S = {"n_cis_leads_annotated": int(len(d)), "n_genome_wide_sig": int(len(sig)),
     "pct_sig_lead_is_PAV_in_cognate": round(float(sig.pav_in_cognate.mean()) * 100, 1),
     "consequence_breakdown_sig": sig.most_severe.value_counts().head(12).to_dict()}
pav = sig[sig.pav_in_cognate == 1]; nonpav = sig[sig.pav_in_cognate == 0]
from scipy.stats import mannwhitneyu
S["effect_size"] = {"median_abs_beta_PAV": round(float(pav.abs_beta.median()), 3),
    "median_abs_beta_nonPAV": round(float(nonpav.abs_beta.median()), 3),
    "mannwhitney_p": float(mannwhitneyu(pav.abs_beta, nonpav.abs_beta, alternative="greater")[1]) if len(pav) > 3 else None,
    "frac_PAV_abs_beta_gt_1": round(float((pav.abs_beta > 1).mean()), 3),
    "frac_nonPAV_abs_beta_gt_1": round(float((nonpav.abs_beta > 1).mean()), 3)}
S["n_PAV_flagged"] = int(len(pav)); S["n_nonPAV"] = int(len(nonpav))
S["flagged_examples"] = pav.sort_values("abs_beta", ascending=False)[["gene", "lead_rsid", "most_severe", "cognate_consequences", "abs_beta", "p"]].head(15).to_dict("records")
json.dump(S, open(f"{OUT}/analysis.json", "w"), indent=1, default=str)

fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
ax[0].hist([nonpav.abs_beta, pav.abs_beta], bins=np.linspace(0, max(sig.abs_beta.max(), 1.5), 25),
           label=[f"non-PAV cis-pQTL (n={len(nonpav)})", f"PAV in cognate gene (n={len(pav)})"], color=["#4575b4", "#d73027"])
ax[0].set(xlabel="|effect size| of cis-pQTL lead (SD/allele)", ylabel="proteins",
          title=f"Effect-size inflation of PAV-flagged cis-pQTLs\nmedian |β|: PAV {pav.abs_beta.median():.2f} vs non-PAV {nonpav.abs_beta.median():.2f}")
ax[0].legend(fontsize=8)
cons = sig.most_severe.value_counts()
pav_mask = sig.set_index("gene").pav_in_cognate
vc = sig.groupby("most_severe").size().sort_values(ascending=False).head(8)
colors = ["#d73027" if c in {"missense_variant","stop_gained","splice_acceptor_variant","splice_donor_variant","frameshift_variant","inframe_deletion"} else "#4575b4" for c in vc.index]
ax[1].barh(range(len(vc)), vc.values, color=colors); ax[1].set_yticks(range(len(vc)), vc.index, fontsize=8); ax[1].invert_yaxis()
ax[1].set(xlabel="number of proteins (genome-wide-sig cis-pQTL)", title="Lead-variant consequence\n(red = protein-altering → candidate epitope artifact)")
fig.tight_layout(); fig.savefig(f"{OUT}/fig1_pqtl_epitope.png", dpi=130)
print(json.dumps(S, indent=1, default=str))
