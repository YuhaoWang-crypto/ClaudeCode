import json, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT="results"; S=json.load(open(f"{OUT}/drugmetab_summary.json"))
P=pd.DataFrame(S["prevalence_fraction"]).T
enz_lab={"GUS":"β-glucuronidase\n(SN-38/NSAID)","azoR":"azoreductase\n(sulfasalazine)","nfs":"nitroreductase\n(metronidazole)",
         "BSH":"bile salt hydrolase\n(deconjugation)","tyrDC":"tyrosine decarb.\n(L-DOPA)","cgr":"glycoside reductase\n(digoxin)"}
P=P[["GUS","azoR","nfs","BSH","tyrDC","cgr"]]
fig,ax=plt.subplots(figsize=(10,6))
im=ax.imshow(np.log1p(P.values),cmap="YlOrRd",aspect="auto")
ax.set_xticks(range(len(P.columns)),[enz_lab[c] for c in P.columns],fontsize=8)
ax.set_yticks(range(len(P.index)),P.index,fontsize=8)
for i in range(P.shape[0]):
    for j in range(P.shape[1]):
        v=P.values[i,j]; ax.text(j,i,f"{v:.1f}" if v>0 else "·",ha="center",va="center",fontsize=8,color="k" if v<4 else "w")
ax.set_title("Drug-metabolizing enzyme genes across gut taxa (copies per genome)\ntyrDC L-DOPA → Enterococcus-specific; nitroreductase broad; cgr/digoxin not captured by annotation",fontsize=10)
fig.colorbar(im,label="log(1+copies/genome)",shrink=.7)
fig.tight_layout(); fig.savefig(f"{OUT}/fig2_drugmetab.png",dpi=130); print("done")
