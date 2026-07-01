#!/usr/bin/env python3
"""Compile all evidence into a single master table (CSV + Markdown) and a
summary figure. Combines: existing-binder re-validation, new-library screen,
round-2 maturation (forced + un-forced), and developability."""
import json, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- gather final wet-lab panel ----
# New peptide leads (round-2, with un-forced binding_confidence = honest metric)
unforced={r['id']:r for r in json.load(open("results/_unforced_round2.json"))}
r2dev={d['lib_id']:d for d in json.load(open("design/round2_developability.json"))}

PANEL=[
 # (wetlab_id, source_id, seq, forced_bind, unforced_bind, iptm_unforced, tier, class)
 ("TZP-P1","R2_008","SLSTLENELSTLENEISTIENEFSTG",0.375,0.369,0.923,"A","de novo helix (Phe-anchor)"),
 ("TZP-P2","R2_002","SLSTLENELSTLENEISTIENEWSTG",0.367,0.368,0.914,"A","de novo helix (Trp-anchor)"),
 ("TZP-P3","R2_035","SLSTLENELSTLENEISTIENEWSTGLENEISTG",0.409,0.409,0.952,"A","de novo helix (extended)"),
 ("TZP-P4","R2_010","SLSTLENELSTLENEISTIENEYSTG",0.356,0.353,0.948,"A","de novo helix (Tyr-anchor)"),
 ("TZP-P5","R2_007","SLSTWENELSTLENEISTIENEWSTG",0.345,0.337,0.923,"B","de novo helix (di-Trp)"),
 ("TZP-P6","R2_034","SLSTLENELSTLENEISTIEG",0.322,0.330,0.940,"B","de novo helix (21-mer, minimal)"),
 ("TZP-B1","spec_13","MVSKTFTVTLPTGATLTVTVTYDYETNVLTVKVTLPQTGKSDEATFKVAAGVTVETILPGVGLVVKAHYDAEKNTITITLECPALGATFTFTLNLS",0.200,0.199,0.714,"A","existing mini-binder (best; orthogonal scaffold)"),
]

import importlib.util
_spec=importlib.util.spec_from_file_location("dev","design/developability.py")
_dev=importlib.util.module_from_spec(_spec); _spec.loader.exec_module(_dev)

rows=[]
for wid,sid,seq,fb,ub,ip,tier,cls in PANEL:
    dev=r2dev.get(sid) or _dev.analyze(seq)
    rows.append({
      "wetlab_id":wid,"source_id":sid,"tier":tier,"class":cls,"length":len(seq),
      "boltz_bind_conf_epitope_forced":fb,"boltz_bind_conf_UNFORCED":ub,"iptm_unforced":ip,
      "net_charge":dev.get("net_charge"),"gravy":dev.get("gravy"),
      "aromatic_frac":dev.get("aromatic_frac"),"cys":dev.get("cys"),
      "n_liabilities":dev.get("n_liabilities",0),"liabilities":";".join(dev.get("liabilities",[])) or "none",
      "sequence":seq})

rows.sort(key=lambda r:(r['tier'],-r['boltz_bind_conf_UNFORCED']))
with open("results/MASTER_panel.csv","w",newline="") as fh:
    w=csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

with open("results/MASTER_panel.md","w") as fh:
    fh.write("| Wet-lab ID | Src | Tier | Class | Len | bind(forced) | **bind(un-forced)** | ipTM | charge | GRAVY | liab | Sequence |\n")
    fh.write("|---|---|---|---|---|---|---|---|---|---|---|---|\n")
    for r in rows:
        fh.write(f"| {r['wetlab_id']} | {r['source_id']} | {r['tier']} | {r['class']} | {r['length']} | "
                 f"{r['boltz_bind_conf_epitope_forced']:.2f} | **{r['boltz_bind_conf_UNFORCED']:.2f}** | {r['iptm_unforced']:.2f} | "
                 f"{r['net_charge']:+d} | {r['gravy']:.2f} | {r['n_liabilities']} | `{r['sequence']}` |\n")
print("wrote results/MASTER_panel.csv and .md with", len(rows), "candidates")

# ---- figure: maturation progression (binding_confidence, un-forced) ----
fig,ax=plt.subplots(1,2,figsize=(11,4.2))
# left: maturation
labels=["decoy\n(field median)","spec_13\n(best existing)","NB094\n(new lead, R1)","R2_002/008\n(matured)","R2_035\n(matured,ext)"]
vals=[0.02,0.199,0.177,0.369,0.409]
colors=["#bbb","#4C78A8","#54A24B","#E4572E","#B00020"]
ax[0].bar(labels,vals,color=colors)
ax[0].set_ylabel("Boltz-2 binding_confidence (un-forced)")
ax[0].set_title("Affinity-maturation progression")
ax[0].axhline(0.2,ls='--',c='gray',lw=0.8); ax[0].text(0,0.21,"existing-best level",fontsize=8,color='gray')
for i,v in enumerate(vals): ax[0].text(i,v+0.008,f"{v:.2f}",ha='center',fontsize=9)
# right: forced vs unforced for round2 (shows robustness)
u=json.load(open("results/_unforced_round2.json"))
forced={"R2_035":0.409,"R2_008":0.375,"R2_002":0.367,"R2_010":0.356,"R2_007":0.345,"R2_034":0.322,"NB094":0.186,"spec_13":0.200}
ids=[x['id'] for x in u]
fx=[forced[i] for i in ids]; ux=[x['bind'] for x in u]
ax[1].scatter(fx,ux,c=['#B00020' if i.startswith('R2') else '#4C78A8' for i in ids],s=60)
for i,idn in enumerate(ids): ax[1].annotate(idn,(fx[i],ux[i]),fontsize=7,xytext=(3,3),textcoords='offset points')
lim=[0,0.45]; ax[1].plot(lim,lim,ls='--',c='gray',lw=0.8)
ax[1].set_xlabel("binding_confidence (epitope-forced)"); ax[1].set_ylabel("binding_confidence (un-forced)")
ax[1].set_title("Forced vs un-forced (on diagonal = robust)")
plt.tight_layout(); plt.savefig("results/figure_maturation.png",dpi=140)
print("wrote results/figure_maturation.png")
EOF_GUARD=1
