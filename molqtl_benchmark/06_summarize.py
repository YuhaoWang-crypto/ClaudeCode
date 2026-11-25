"""Merge all experiments into results/summary.json and draw the figures."""
import json, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
S = {"cis_eqtl": json.load(open("results/cis_eqtl.json")), "predict_expression": json.load(open("results/predict_expression.json")),
     "finemap_simulation": pd.read_csv("results/finemap_sim.csv").to_dict("records"),
     "coloc_simulation": pd.read_csv("results/coloc_sim.csv").to_dict("records"),
     "sort1_real_locus": pd.read_csv("results/sort1_coloc.csv").to_dict("records")}
json.dump(S, open("results/summary.json", "w"), indent=1)

# fig 1: replication of lead eQTL effects + expression prediction
rep = pd.read_csv("results/cis_eqtl_replication.csv"); pre = pd.read_csv("results/predict_pretrained.csv"); own = pd.read_csv("results/predict_selftrained.csv")
fig, ax = plt.subplots(1, 4, figsize=(18, 4.3))
for a, k, t in [(ax[0], "EUR_test", "unseen Europeans (n=112)"), (ax[1], "YRI", "Yoruba (n=89)")]:
    d = rep.dropna(subset=[f"slope_{k}"]); ok = (d[f"p_{k}"] < .05) & (np.sign(d.slope_train) == np.sign(d[f"slope_{k}"]))
    a.scatter(d.slope_train, d[f"slope_{k}"], s=14, c=np.where(ok, "#1a9850", "#d73027")); a.plot([-2, 2], [-2, 2], "k:", lw=.8)
    a.axhline(0, c="#bbb", lw=.5); a.axvline(0, c="#bbb", lw=.5)
    a.set(title=f"Lead cis-eQTL effect in {t}\ngreen = replicated ({ok.sum()}/{len(rep)} eGenes)", xlabel="effect in discovery (261 Europeans)", ylabel="effect in held-out group")
b = np.linspace(-.4, 1, 29)
ax[2].hist(pre.r_EUR_all.dropna(), b, alpha=.7, label=f"Europeans (median r={pre.r_EUR_all.median():.2f})", color="#4575b4")
ax[2].hist(pre.r_YRI.dropna(), b, alpha=.6, label=f"Yoruba (median r={pre.r_YRI.median():.2f})", color="#d73027")
ax[2].set(title=f"Published GTEx LCL weights, {len(pre)} genes\npredicted vs measured expression", xlabel="r per gene", ylabel="genes"); ax[2].legend(fontsize=8)
ax[3].hist(own.r_EUR_test, b, alpha=.7, label=f"unseen Europeans (median r={own.r_EUR_test.median():.2f})", color="#4575b4")
ax[3].hist(own.r_YRI, b, alpha=.6, label=f"Yoruba (median r={own.r_YRI.median():.2f})", color="#d73027")
ax[3].set(title=f"Elastic net trained here, all {len(own)} genes", xlabel="r per gene", ylabel="genes"); ax[3].legend(fontsize=8)
fig.tight_layout(); fig.savefig("results/fig1_eqtl_replication_prediction.png", dpi=130)

# fig 2: fine-mapping + coloc simulations
A = pd.read_csv("results/finemap_sim.csv"); B = pd.read_csv("results/coloc_sim.csv")
fig, ax = plt.subplots(1, 3, figsize=(18, 4.8), gridspec_kw={"width_ratios": [1, 1, 1]})
ms = list(dict.fromkeys(A.method)); cols = ["#4575b4", "#91bfdb", "#fdae61", "#d73027"]; x = np.arange(3)
for j, m in enumerate(ms):
    d = A[A.method == m].sort_values("n_causal")
    ax[0].bar(x + (j - 1.5) * .2, d.coverage, .19, color=cols[j], label=m); ax[1].bar(x + (j - 1.5) * .2, d.power, .19, color=cols[j])
ax[0].axhline(.95, ls=":", c="k", lw=.8); ax[0].set_xticks(x, ["1 causal", "2 causal", "3 causal"]); ax[1].set_xticks(x, ["1 causal", "2 causal", "3 causal"])
ax[0].set(ylim=(0, 1.45)); ax[0].legend(fontsize=7.5, loc="upper left"); ax[0].set_title("SuSiE: 95% credible sets containing a causal variant", fontsize=10)
ax[1].set(ylim=(0, 1)); ax[1].set_title("SuSiE: causal variants captured in any credible set", fontsize=10)
order = ["same causal variant, strong eQTL", "same causal variant, weak eQTL (4% of variance)", "different variants, r2 < 0.1",
         "different variants, r2 0.3-0.7", "different variants, r2 0.8-0.95", "GWAS signal only, no eQTL"]
B = B.set_index("scenario").reindex(order); y = np.arange(len(B))
ax[2].barh(y, B["called_shared_H4_gt_0.8"], color=["#1a9850"] * 2 + ["#d73027"] * 4); ax[2].set_yticks(y, [f"{s} (n={int(n)})" for s, n in zip(B.index, B.runs)], fontsize=8)
ax[2].invert_yaxis(); ax[2].set(xlim=(0, 1)); ax[2].set_title("coloc: fraction called 'shared' (PP.H4 > 0.8)\ngreen = truly shared, red = truly not shared", fontsize=10)
fig.tight_layout(); fig.savefig("results/fig2_finemap_coloc_simulation.png", dpi=130)

# fig 3: SORT1 locus
d = pd.read_csv("results/sort1_coloc.csv"); p = d.pivot(index="gene", columns="tissue", values="PP_H4_shared")[["liver", "muscle", "adipose", "blood", "LCL"]].loc[["SORT1", "PSRC1", "CELSR2"]]
fig, ax = plt.subplots(figsize=(6.5, 3)); im = ax.imshow(p.values, cmap="Greens", vmin=0, vmax=1)
ax.set_xticks(range(5), p.columns); ax.set_yticks(range(3), p.index)
for i in range(3):
    for j in range(5): ax.text(j, i, f"{p.values[i, j]:.2f}", ha="center", va="center", color="w" if p.values[i, j] > .5 else "k")
ax.set_title("1p13 locus x LDL cholesterol: coloc PP.H4 by gene and tissue", fontsize=10)
fig.tight_layout(); fig.savefig("results/fig3_sort1_coloc.png", dpi=130)
print(json.dumps(S["predict_expression"])); print(pd.DataFrame(S["finemap_simulation"]).to_string()); print(pd.DataFrame(S["coloc_simulation"]).to_string())
