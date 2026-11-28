import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
RES="/home/user/ClaudeCode/protein_design_benchmark/results"
sizes=[("esm2_t12_35M_UR50D",35),("esm2_t30_150M_UR50D",150),("esm2_t33_650M_UR50D",650)]
data={m:json.load(open(f"{RES}/esm_{m}.json")) for m,_ in sizes}
bl=json.load(open(f"{RES}/blosum62.json"))
assays=[("BLAT_ECOLX_Stiffler_2015","β-lactamase\n(enzyme activity)"),
        ("SPG1_STRSG_Olson_2014","GB1\n(IgG-Fc binding)"),
        ("GFP_AEQVI_Sarkisyan_2016","GFP\n(fluorescence)")]
col={"35":"#cfe0c3","150":"#8bb174","650":"#3f6634","BLOSUM":"#b0592c"}

fig,axes=plt.subplots(1,2,figsize=(14,6),gridspec_kw={"width_ratios":[1.25,1]})
# Panel A: grouped bars
ax=axes[0]; x=np.arange(len(assays)); w=0.2
vals={m:[data[m][a]["spearman_masked_marginal"] for a,_ in assays] for m,_ in sizes}
blv=[bl[a] for a,_ in assays]
ax.bar(x-1.5*w,vals["esm2_t12_35M_UR50D"],w,label="ESM-2 35M",color=col["35"])
ax.bar(x-0.5*w,vals["esm2_t30_150M_UR50D"],w,label="ESM-2 150M",color=col["150"])
ax.bar(x+0.5*w,vals["esm2_t33_650M_UR50D"],w,label="ESM-2 650M",color=col["650"])
ax.bar(x+1.5*w,blv,w,label="BLOSUM62 (trivial)",color=col["BLOSUM"])
ax.axhline(0,c="#333",lw=0.8)
ax.set_xticks(x); ax.set_xticklabels([l for _,l in assays],fontsize=10)
ax.set_ylabel("Zero-shot Spearman ρ  (ESM score vs measured DMS)",fontsize=11)
ax.set_title("A. Zero-shot variant-effect prediction, by assay & model size",fontsize=12,loc="left")
ax.legend(fontsize=9,ncol=2,loc="upper right"); ax.set_ylim(-0.02,0.85)
# annotate the GFP failure
gi=2
ax.annotate("PLM worse than\ntrivial BLOSUM →\nreal boundary",(gi+1.5*w,blv[gi]),
            (gi-0.1,0.42),fontsize=8.5,ha="center",color="#842",
            arrowprops=dict(arrowstyle="->",color="#842"))
# Panel B: scaling curves
ax=axes[1]
pr=[35,150,650]
for a,lab in assays:
    ys=[data[m][a]["spearman_masked_marginal"] for m,_ in sizes]
    ax.plot(pr,ys,"o-",label=lab.replace("\n"," "),lw=2,ms=7)
    ax.axhline(bl[a],ls=":",lw=1,color="#999")
ax.set_xscale("log")
from matplotlib.ticker import NullLocator, FixedLocator, FixedFormatter
ax.xaxis.set_minor_locator(NullLocator())
ax.xaxis.set_major_locator(FixedLocator(pr)); ax.xaxis.set_major_formatter(FixedFormatter(["35M","150M","650M"]))
ax.set_xlabel("ESM-2 parameters",fontsize=11)
ax.set_ylabel("Spearman ρ",fontsize=11)
ax.set_title("B. Does scaling help? (dotted = BLOSUM62 floor)",fontsize=12,loc="left")
ax.legend(fontsize=9,loc="center right"); ax.grid(alpha=0.25)
plt.tight_layout(); plt.savefig(f"{RES}/fig1_esm_dms.png",dpi=135,bbox_inches="tight")
print("saved",f"{RES}/fig1_esm_dms.png")
print("\nSummary (masked-marginal ρ):")
print(f"{'assay':28s}{'35M':>7}{'150M':>7}{'650M':>7}{'BLOSUM':>8}")
for a,_ in assays:
    print(f"{a:28s}"+"".join(f"{data[m][a]['spearman_masked_marginal']:7.3f}" for m,_ in sizes)+f"{bl[a]:8.3f}")
