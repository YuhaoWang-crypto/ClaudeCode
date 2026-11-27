import json, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT = "results"
ri = pd.read_csv(f"{OUT}/reverse_infer.csv"); S1 = json.load(open(f"{OUT}/reverse_infer_summary.json")); S2 = json.load(open(f"{OUT}/moa_model_summary.json"))
N = S1["n_genes"]
fig, ax = plt.subplots(1, 3, figsize=(18, 4.8))
# A: true-target rank distribution
ax[0].hist(ri.target_rank, bins=np.linspace(0, N, 40), color="#4575b4", alpha=.85)
ax[0].axvline(N / 2, ls="--", c="k", lw=1, label=f"random expectation (≈{N//2})")
ax[0].axvline(ri.target_rank.median(), ls="-", c="#d73027", lw=2, label=f"observed median ({int(ri.target_rank.median())})")
ax[0].set(xlabel="genome-wide rank of the drug's TRUE physical target\n(by sensitivity–dependency correlation, 1 = best of 18,443)",
          ylabel="number of drugs", title=f"A. Physical-target recovery fails\n({S1['n_drugs']} drugs, only {S1['frac_target_in_top100']*100:.1f}% put true target in top 100)")
ax[0].legend(fontsize=8)
# B: MoA connectivity
mc = S2["MoA_connectivity"]
ax[1].bar([0, 1], [mc["nearest_neighbour_shares_MoA"], mc["random_baseline"]], color=["#1a9850", "#999"], width=.6)
ax[1].set_xticks([0, 1], ["nearest drug by\nsensitivity profile", "random\npairing"]); ax[1].set_ylim(0, 1)
for i, v in enumerate([mc["nearest_neighbour_shares_MoA"], mc["random_baseline"]]): ax[1].text(i, v + .02, f"{v:.2f}", ha="center", fontsize=11)
ax[1].set(ylabel="fraction sharing mechanism of action (MoA)", title=f"B. Mechanism class IS recoverable\n({mc['n_drugs_with_moa']} drugs, {mc['n_moa_classes']} MoA classes)")
# C: the on-target subset
sub = ri.sort_values("target_rank").head(20)
ax[2].barh(range(len(sub)), sub.target_corr, color="#4575b4")
ax[2].set_yticks(range(len(sub)), [f"{d[:16]}→{t}" for d, t in zip(sub.drug, sub.target)], fontsize=7)
ax[2].invert_yaxis(); ax[2].set(xlabel="sensitivity–dependency correlation", title="C. The minority that works: genuinely on-target drugs\n(top 20 of 531 by target rank)")
fig.tight_layout(); fig.savefig(f"{OUT}/fig1_phenotype_target.png", dpi=130)
json.dump({**S1, **S2}, open(f"{OUT}/summary.json", "w"), indent=1)
print("fig written")
