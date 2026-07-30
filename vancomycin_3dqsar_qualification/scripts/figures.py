import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
OUT="/tmp/claude-0/-home-user-ClaudeCode/0e1c2284-3690-5d68-a96a-480ebe908c1c/scratchpad"
R=json.load(open(f"{OUT}/analysis_results.json"))
labels=["Parent\nVancomycin B","RRT 0.87\n(7S,26R)-Vanco B","RRT 0.75\n(26R)-DMe-DeLI"]
keys=["Parent","RRT087","RRT075"]
COL=["#2b6cb0","#38a169","#dd6b20"]

# ---- Figure 3: toxicophore / functional-group inventory comparison ----
fgs=list(R["Parent"]["functional_groups"].keys())
# keep the toxicologically-relevant ones
show=["phenol","aromatic_ring","aromatic_chloride(ArCl)","aryl_ether","benzylic_alcohol",
      "primary_aliphatic_amine","secondary_aliphatic_amine(NHMe/NHR)","primary_amide",
      "secondary_amide","carboxylic_acid"]
mat=np.array([[R[k]["functional_groups"][f] for f in show] for k in keys])
fig,ax=plt.subplots(figsize=(11,6))
im=ax.imshow(mat,aspect="auto",cmap="Blues",vmin=0,vmax=6)
ax.set_xticks(range(len(show))); ax.set_xticklabels(show,rotation=40,ha="right",fontsize=9)
ax.set_yticks(range(3)); ax.set_yticklabels([l.replace("\n"," ") for l in labels],fontsize=10)
for i in range(3):
    for j in range(len(show)):
        v=mat[i,j]
        # highlight the changed cell (RRT075 amines)
        chg = (i==2 and show[j] in ("primary_aliphatic_amine","secondary_aliphatic_amine(NHMe/NHR)"))
        ax.text(j,i,str(v),ha="center",va="center",fontsize=11,
                color=("red" if chg else ("white" if v>=4 else "black")),
                fontweight=("bold" if chg else "normal"))
ax.set_title("Toxicophore / functional-group inventory (counts) — red = only difference vs parent",fontsize=11)
plt.colorbar(im,label="count",shrink=0.7)
plt.tight_layout(); plt.savefig(f"{OUT}/fig3_toxicophore_inventory.png",dpi=140); plt.close()
print("saved fig3")

# ---- Figure 4: 3D shape descriptors + key physchem ----
d3=["RadiusOfGyration","Asphericity","SpherocityIndex","Eccentricity","NPR1","NPR2"]
x=np.arange(len(d3)); w=0.26
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(13,5))
for i,k in enumerate(keys):
    vals=[R[k]["descr3d"][p] for p in d3]
    ax1.bar(x+(i-1)*w,vals,w,label=labels[i].replace("\n"," "),color=COL[i])
ax1.set_xticks(x); ax1.set_xticklabels(d3,rotation=30,ha="right",fontsize=9)
ax1.set_title("3D whole-molecule shape descriptors (single low-energy conformer)")
ax1.legend(fontsize=8)
# physchem panel
pc=["MW","TPSA","HBD","HBA","RotBonds","cLogP_Crippen"]
xp=np.arange(len(pc))
for i,k in enumerate(keys):
    vals=[R[k]["physchem"][p] for p in pc]
    ax2.bar(xp+(i-1)*w,vals,w,label=labels[i].replace("\n"," "),color=COL[i])
ax2.set_xticks(xp); ax2.set_xticklabels(pc,rotation=30,ha="right",fontsize=9)
ax2.set_yscale("symlog")
ax2.set_title("Physicochemical descriptors (symlog scale)")
ax2.legend(fontsize=8)
plt.tight_layout(); plt.savefig(f"{OUT}/fig4_descriptors.png",dpi=140); plt.close()
print("saved fig4")

# ---- Figure 5: read-across similarity metrics ----
mcs=json.load(open(f"{OUT}/mcs_core_alignment.json"))
metrics=["ECFP4\nTanimoto","MACCS\nTanimoto","Common heavy\natoms (frac)"]
def frac(k):
    if k=="Parent": return 1.0
    name={"RRT087":"RRT 0.87 (7S,26R)-Vanco B","RRT075":"RRT 0.75 (26R)-DMe-DeLI"}[k]
    return mcs[name]["core_atoms_matched"]/R["Parent"]["physchem"]["HeavyAtoms"]
fig,ax=plt.subplots(figsize=(8,5))
x=np.arange(3); w=0.26
for i,k in enumerate(keys):
    vals=[R[k]["ECFP4_Tanimoto_vs_parent"],R[k]["MACCS_Tanimoto_vs_parent"],frac(k)]
    b=ax.bar(x+(i-1)*w,vals,w,label=labels[i].replace("\n"," "),color=COL[i])
    ax.bar_label(b,fmt="%.2f",fontsize=8)
ax.set_xticks(x); ax.set_xticklabels(metrics)
ax.set_ylim(0,1.15); ax.set_ylabel("similarity to parent (1.0 = identical)")
ax.set_title("2D structural similarity & shared-skeleton read-across metrics")
ax.legend(fontsize=8,loc="lower center")
plt.tight_layout(); plt.savefig(f"{OUT}/fig5_similarity.png",dpi=140); plt.close()
print("saved fig5")
