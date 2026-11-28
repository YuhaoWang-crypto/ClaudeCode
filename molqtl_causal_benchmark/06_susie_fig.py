import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT="results"; d=pd.read_csv(f"{OUT}/susie_coloc_raw.csv")
scen=["shared_single","two_distinct","shared_plus_extra"]
lab=["shared\nsingle variant\n(truth: SHARED)","two distinct\nvariants\n(truth: NOT shared)","shared signal +\nextra 2nd signal in trait1\n(truth: SHARED)"]
fig,ax=plt.subplots(1,2,figsize=(13,4.8))
abf=[d[d.scenario==s].abf_H4.mean() for s in scen]; sus=[d[d.scenario==s].susie_H4.mean() for s in scen]
x=np.arange(3); w=0.35
ax[0].bar(x-w/2,abf,w,label="coloc.abf (single-variant)",color="#91bfdb")
ax[0].bar(x+w/2,sus,w,label="coloc.susie (multi-signal)",color="#1a9850")
for i,(a,s) in enumerate(zip(abf,sus)): ax[0].text(i-w/2,a+.02,f"{a:.2f}",ha="center",fontsize=8); ax[0].text(i+w/2,s+.02,f"{s:.2f}",ha="center",fontsize=8)
ax[0].set_xticks(x,lab,fontsize=8); ax[0].axhline(0.8,ls=":",c="k",lw=.8); ax[0].set(ylabel="mean PP.H4 (colocalization)",ylim=(0,1.1),
    title="A. coloc.abf vs coloc.susie (real GEUVADIS LD, known truth)")
ax[0].legend(fontsize=8,loc="upper center")
# B: shared-call rate in the multi-signal scenario
sp=d[d.scenario=="shared_plus_extra"]
rates=[(sp.abf_H4>0.8).mean(),(sp.susie_H4>0.8).mean()]
ax[1].bar([0,1],rates,color=["#91bfdb","#1a9850"])
for i,v in enumerate(rates): ax[1].text(i,v+.01,f"{v:.0%}",ha="center",fontsize=12)
ax[1].set_xticks([0,1],["coloc.abf","coloc.susie"],fontsize=10)
ax[1].set(ylabel="% correctly calling SHARED (H4>0.8)",ylim=(0,1),
    title="B. When a locus has multiple signals,\nSuSiE-coloc recovers the shared signal more often")
fig.tight_layout(); fig.savefig(f"{OUT}/fig2_susie_coloc.png",dpi=130); print("done",abf,sus)
