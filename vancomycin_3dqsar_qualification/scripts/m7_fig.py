import json, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
OUT="/tmp/claude-0/-home-user-ClaudeCode/0e1c2284-3690-5d68-a96a-480ebe908c1c/scratchpad"
M=json.load(open(f"{OUT}/m7_results.json"))
comp=list(M["per_compound"].keys())
short=["Parent\nVancomycin B","RRT 0.87\n(7S,26R)","RRT 0.75\n(DMe-DeLI)"]

fig,(ax1,ax2)=plt.subplots(1,2,figsize=(13,4.6),gridspec_kw={"width_ratios":[1.05,1]})

# Panel A: concordance matrix (2 methodologies x 3 compounds), green = negative
axm=ax1
grid=np.zeros((2,3))
labels=[["","",""],["","",""]]
for j,c in enumerate(comp):
    r=M["per_compound"][c]
    grid[0,j]= 0 if not r["expert_rulebase_alerts"] else 1
    grid[1,j]= 0 if r["stat_Ames_call"]=="negative" else 1
    labels[0][j]="NEGATIVE\n(no alert)" if grid[0,j]==0 else "ALERT"
    labels[1][j]=f"NEGATIVE\np={r['stat_Ames_prob']:.2f}" if grid[1,j]==0 else f"POSITIVE\np={r['stat_Ames_prob']:.2f}"
axm.imshow(grid,cmap=matplotlib.colors.ListedColormap(["#2e7d32","#c0392b"]),vmin=0,vmax=1,aspect="auto")
axm.set_xticks(range(3)); axm.set_xticklabels(short,fontsize=9)
axm.set_yticks([0,1]); axm.set_yticklabels(["Methodology 1\nExpert rule-based\n(Benigni-Bossa)","Methodology 2\nStatistical QSAR\n(Ames, Hansen N6512)"],fontsize=9)
for i in range(2):
    for j in range(3):
        axm.text(j,i,labels[i][j],ha="center",va="center",color="white",fontsize=9,fontweight="bold")
axm.set_title("A.  ICH M7 two-methodology concordance\n(both negative -> Class 5)",fontsize=11,loc="left")

# Panel B: Ames QSAR probability vs threshold + AD
axp=ax2
probs=[M["per_compound"][c]["stat_Ames_prob"] for c in comp]
ad=[M["per_compound"][c]["stat_AD_maxTanimoto"] for c in comp]
thr=M["methodology_2"]["threshold_youden"]
x=np.arange(3)
b=axp.bar(x,probs,color=["#2b6cb0","#38a169","#dd6b20"],width=0.55)
axp.axhline(thr,ls="--",color="k",lw=1.2)
axp.text(2.5,thr+0.02,f"decision threshold {thr:.2f}",ha="right",fontsize=8)
axp.axhline(0.5,ls=":",color="grey",lw=0.8)
axp.set_xticks(x); axp.set_xticklabels(short,fontsize=9)
axp.set_ylim(0,1); axp.set_ylabel("P(Ames positive)")
axp.set_title(f"B.  Statistical Ames QSAR (5-fold AUC={M['methodology_2']['cv_auc']:.3f})",fontsize=11,loc="left")
for i,(p,a) in enumerate(zip(probs,ad)):
    axp.text(i,p+0.03,f"{p:.2f}",ha="center",fontsize=9,fontweight="bold")
    axp.text(i,0.04,f"AD sim\n{a:.2f}",ha="center",fontsize=7,color="white")
axp.text(0.5,0.90,"All below threshold -> negative\n(note: AD similarity low -> low-confidence;\nexpert rulebase is the primary M7 call)",
         transform=axp.transAxes,fontsize=7.5,ha="center",va="top",
         bbox=dict(boxstyle="round",fc="#fff3cd",ec="#e0a800"))
plt.tight_layout(); plt.savefig(f"{OUT}/fig7_m7.png",dpi=140,bbox_inches="tight"); print("saved fig7_m7.png")
