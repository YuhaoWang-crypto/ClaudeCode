import json, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT="results"
dg=pd.read_csv(f"{OUT}/diagnosis.csv"); R=pd.read_csv(f"{OUT}/coloc_all.csv"); loci=json.load(open(f"{OUT}/causal_loci.json"))
# status per locus
def status(r):
    if not r.causal_has_cis_eQTL: return 0   # no detectable eQTL
    if r.causal_PPH4<=0.5: return 1           # eQTL but no coloc
    return 2                                   # colocalizes
dg["status"]=dg.apply(status,axis=1)
order=dg.sort_values(["status","locus"])
fig,ax=plt.subplots(1,3,figsize=(18,5))
cmap={0:"#d73027",1:"#fdae61",2:"#1a9850"}; labs={0:"no detectable cis-eQTL",1:"eQTL present but no coloc\n(causal variant acts via protein/function)",2:"colocalizes"}
y=np.arange(len(order))
ax[0].barh(y,[1]*len(order),color=[cmap[s] for s in order.status])
for i,(_,r) in enumerate(order.iterrows()): ax[0].text(0.02,i,f"{r.locus}  (eQTL minP {r.causal_eQTL_minP}, {r.n_coloc_0_5 if 'n_coloc_0_5' in r else r['n_coloc_0.5']} genes coloc)",va="center",fontsize=8)
ax[0].set_yticks([]); ax[0].set_xticks([]); ax[0].set(title="A. Can eQTL-coloc recover the validated causal gene?\n15 LDL drug-target / Mendelian loci")
from matplotlib.patches import Patch
ax[0].legend(handles=[Patch(color=cmap[k],label=labs[k]) for k in [2,1,0]],fontsize=8,loc="lower right")
# B: SORT1 multi-gene ambiguity
s=R[(R.locus=="SORT1")].sort_values("PP_H4",ascending=False).drop_duplicates("gene_id").head(8)
ax[1].barh(range(len(s)),s.PP_H4,color=["#1a9850" if i else "#4575b4" for i in s.is_causal])
ax[1].set_yticks(range(len(s)),s.gene_symbol,fontsize=8); ax[1].invert_yaxis(); ax[1].axvline(0.8,ls=":",c="k",lw=.8)
ax[1].set(xlabel="coloc PP.H4 (eQTL vs LDL GWAS)",title="B. Even where coloc fires, it lights up MANY genes\nSORT1 locus: 3 genes colocalize (blue=causal SORT1)",xlim=(0,1))
# C: breakdown
cnt=dg.status.value_counts().reindex([2,1,0]).fillna(0)
ax[2].bar([0,1,2],cnt.values,color=[cmap[2],cmap[1],cmap[0]])
for i,v in enumerate(cnt.values): ax[2].text(i,v+0.1,int(v),ha="center",fontsize=12)
ax[2].set_xticks([0,1,2],["colocalizes","eQTL,\nno coloc","no\neQTL"],fontsize=9)
ax[2].set(ylabel="# of 15 validated causal genes",title="C. eQTL-coloc recovers 2/15 causal genes\n(half lack eQTL; most of rest act via protein, not expression)",ylim=(0,9))
fig.tight_layout(); fig.savefig(f"{OUT}/fig1_causal_gene.png",dpi=130); print("done")
