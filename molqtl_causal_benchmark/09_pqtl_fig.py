import json, csv, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
# eQTL causal-gene PP.H4 from diagnosis.csv
eqtl={}
for r in csv.DictReader(open("results/diagnosis.csv")):
    eqtl[r["locus"]]=float(r["causal_PPH4"])
# pQTL results
pq={}
for r in csv.DictReader(open("results/pqtl_coloc.csv")):
    pq[r["gene"]]=dict(h4=float(r["PP_H4"]),h3=float(r["PP_H3"]),h2=float(r["PP_H2"]),
                       minp=float(r["pqtl_minP"]))
genes=[g for g in ["PCSK9","APOE","APOB","LDLR","HMGCR","CETP","SORT1"] if g in pq]
eq=[eqtl.get(g,np.nan) for g in genes]
pp=[pq[g]["h4"] for g in genes]

fig,axes=plt.subplots(1,2,figsize=(14,6),gridspec_kw={"width_ratios":[1.35,1]})
# --- Panel A: eQTL vs pQTL PP.H4 grouped bars ---
ax=axes[0]; x=np.arange(len(genes)); w=0.38
b1=ax.bar(x-w/2,eq,w,label="cis-eQTL coloc (expression)",color="#9db4c0")
b2=ax.bar(x+w/2,pp,w,label="cis-pQTL coloc (protein)",color="#c1666b")
ax.axhline(0.8,ls="--",c="#444",lw=1); ax.text(len(genes)-0.5,0.82,"PP.H4=0.8",fontsize=9,ha="right")
ax.set_xticks(x); ax.set_xticklabels(genes,fontsize=11)
ax.set_ylabel("Colocalization PP.H4 with LDL-C GWAS",fontsize=11)
ax.set_ylim(0,1.08); ax.legend(fontsize=10,loc="upper right")
ax.set_title("A. eQTL vs pQTL colocalization per causal gene",fontsize=12,loc="left")
for xi,g in zip(x,genes):
    for bx,val in [(xi-w/2,eqtl.get(g,np.nan)),(xi+w/2,pq[g]["h4"])]:
        if val>=0.8: ax.text(bx,val+0.01,"✓",ha="center",fontsize=13,color="#2a7",weight="bold")
# annotate APOB H3
i=genes.index("APOB")
ax.annotate("strong cis-pQTL\nbut PP.H3≈1\n(distinct variants/\nepitope artifact)",
            (i+w/2,0.05),(i+w/2,0.42),fontsize=8,ha="center",
            arrowprops=dict(arrowstyle="->",color="#852"))
# --- Panel B: pQTL cis-signal strength (minP) ---
ax=axes[1]
mp=[-np.log10(pq[g]["minp"]) for g in genes]
colors=["#c1666b" if pq[g]["h4"]>=0.8 else ("#d9a441" if mp_i>8 else "#9db4c0") for g,mp_i in zip(genes,mp)]
ax.barh(np.arange(len(genes)),mp,color=colors)
ax.axvline(8,ls="--",c="#444",lw=1); ax.text(8.2,len(genes)-0.5,"cis-pQTL\n-log10P=8",fontsize=8)
ax.set_yticks(np.arange(len(genes))); ax.set_yticklabels(genes,fontsize=11); ax.invert_yaxis()
ax.set_xlabel("cis-pQTL signal strength  -log10(min P)",fontsize=11)
ax.set_title("B. Does the protein even have a cis-pQTL?",fontsize=12,loc="left")
for i,g in enumerate(genes):
    tag="coloc✓" if pq[g]["h4"]>=0.8 else ("H3" if pq[g]["h3"]>0.5 else "no cis-pQTL" if mp[i]<8 else "")
    ax.text(mp[i]+1,i,tag,va="center",fontsize=8,color="#333")
plt.tight_layout(); plt.savefig("results/fig3_pqtl_coloc.png",dpi=135,bbox_inches="tight")
print("saved results/fig3_pqtl_coloc.png")

# summary table
print("\ngene   eQTL_H4  pQTL_H4  pQTL_minP   verdict")
for g in genes:
    v=("pQTL rescues" if pq[g]["h4"]>=0.8 and eqtl.get(g,0)<0.5 else
       "eQTL only" if eqtl.get(g,0)>=0.8 and pq[g]["h4"]<0.5 else
       "pQTL H3 (artifact?)" if pq[g]["h3"]>0.5 and pq[g]["minp"]<1e-8 else
       "neither (no cis-pQTL)")
    print(f"{g:6s} {eqtl.get(g,float('nan')):6.2f}  {pq[g]['h4']:6.3f}  {pq[g]['minp']:.1e}  {v}")
