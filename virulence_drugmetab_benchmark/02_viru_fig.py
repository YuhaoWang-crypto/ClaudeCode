import json, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT="results"; S=json.load(open(f"{OUT}/virulence_summary.json"))
loci=["aerobactin_iuc","yersiniabactin_ybt","salmochelin_iro"]; labels=["aerobactin\n(iuc)","yersiniabactin\n(ybt)","salmochelin\n(iro)"]
cats=["invasive","mucosal","carriage","environmental"]; colors=["#d73027","#fc8d59","#91bfdb","#4575b4"]
fig,ax=plt.subplots(1,2,figsize=(13,4.8))
x=np.arange(len(loci)); w=0.2
for i,c in enumerate(cats):
    vals=[S["marker_prevalence_by_category"][m][c] for m in loci]
    ax[0].bar(x+(i-1.5)*w, vals, w, label=c, color=colors[i])
ax[0].set_xticks(x,labels,fontsize=9); ax[0].set(ylabel="fraction of genomes carrying locus",
    title="A. Hypervirulence siderophore loci by isolation source\n(invasive > carriage > environmental; OR 1.5–1.7, p≤2e-8)")
ax[0].legend(fontsize=8,title="source")
hv=S["mean_hvkp_loci_by_category"]
ax[1].bar(range(4),[hv[c] for c in cats],color=colors)
for i,c in enumerate(cats): ax[1].text(i,hv[c]+0.02,f"{hv[c]:.2f}",ha="center",fontsize=10)
ax[1].set_xticks(range(4),cats,fontsize=9,rotation=20); ax[1].set(ylabel="mean # hvKp loci per genome",
    title="B. Hypervirulence gene burden gradient\n(monotonic with invasiveness)")
fig.tight_layout(); fig.savefig(f"{OUT}/fig1_virulence.png",dpi=130); print("fig done")
