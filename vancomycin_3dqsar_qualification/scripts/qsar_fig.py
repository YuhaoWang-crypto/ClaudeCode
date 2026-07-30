import json, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
OUT="/tmp/claude-0/-home-user-ClaudeCode/0e1c2284-3690-5d68-a96a-480ebe908c1c/scratchpad"
R=json.load(open(f"{OUT}/qsar_results.json"))
A=json.load(open(f"{OUT}/qsar_augment.json"))
EP=[r["endpoint"] for r in R["cv_metrics"]]
auc=[r["CV_AUC"] for r in R["cv_metrics"]]
cols=["Parent (Vancomycin B)","RRT 0.87 (7S,26R)-Vanco B","RRT 0.75 (26R)-DMe-DeLI"]
COL=["#2b6cb0","#38a169","#dd6b20"]

fig=plt.figure(figsize=(14,9))
gs=fig.add_gridspec(2,2,height_ratios=[1.15,1],hspace=0.55,wspace=0.25)

# Panel A: QSAR predicted probabilities heatmap (3 molecules x 12 endpoints)
axA=fig.add_subplot(gs[0,:])
M=np.array([[R["predictions"][c][e] for e in EP] for c in cols])
im=axA.imshow(M,aspect="auto",cmap="RdYlGn_r",vmin=0,vmax=1)
axA.set_xticks(range(len(EP))); axA.set_xticklabels(EP,rotation=40,ha="right",fontsize=9)
axA.set_yticks(range(3)); axA.set_yticklabels(["Parent\n(measured+in-set)","RRT 0.87\n(=parent FP)","RRT 0.75\n(QSAR pred.)"],fontsize=9)
for i in range(3):
    for j in range(len(EP)):
        axA.text(j,i,f"{M[i,j]:.2f}",ha="center",va="center",fontsize=8,
                 color="white" if (M[i,j]>0.65 or M[i,j]<0.2) else "black")
axA.set_title("A.  12-endpoint Tox21 QSAR — predicted activity probability (RandomForest, Morgan r2)",fontsize=11,loc="left")
plt.colorbar(im,ax=axA,label="P(active)",shrink=0.8)

# Panel B: CV-AUC per endpoint
axB=fig.add_subplot(gs[1,0])
order=np.argsort(auc)
axB.barh(np.array(EP)[order],np.array(auc)[order],color="#4a6fa5")
axB.axvline(np.mean(auc),ls="--",color="k",lw=1)
axB.set_xlim(0.5,1.0); axB.set_xlabel("5-fold OOF ROC-AUC")
axB.set_title(f"B.  Model performance (mean AUC={np.mean(auc):.3f})",fontsize=11,loc="left")
axB.tick_params(labelsize=8)

# Panel C: delta RRT075 vs held-out parent (honest, both external to identical rows)
axC=fig.add_subplot(gs[1,1])
hp=A["heldout_parent_pred"]; h7=A["heldout_pred_rrt075"]
delta=[h7[e]-hp[e] for e in EP]
c=["#c0392b" if d>0 else "#2e7d32" for d in delta]
axC.barh(EP,delta,color=c)
axC.axvline(0,color="k",lw=0.8)
axC.set_xlabel("Δ P(active)  =  RRT 0.75  −  parent")
axC.set_title("C.  Impurity − parent (held-out; RRT 0.87 Δ≡0)",fontsize=11,loc="left")
axC.tick_params(labelsize=8)
mx=max(0.15,max(abs(d) for d in delta)*1.3); axC.set_xlim(-mx,mx)
for i,d in enumerate(delta):
    axC.text(d+(0.005 if d>=0 else -0.005),i,f"{d:+.2f}",va="center",
             ha="left" if d>=0 else "right",fontsize=7)

plt.savefig(f"{OUT}/fig6_qsar.png",dpi=140,bbox_inches="tight"); print("saved fig6_qsar.png")
print("max |delta| RRT075 vs parent (held-out):",round(max(abs(d) for d in delta),3))
