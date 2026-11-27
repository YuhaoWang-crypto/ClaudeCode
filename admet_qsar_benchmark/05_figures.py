import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT = "results"
q = json.load(open(f"{OUT}/qsar_splits.json")); ct = json.load(open(f"{OUT}/curve_transfer.json")); fc = json.load(open(f"{OUT}/fair_compare.json"))
fig, ax = plt.subplots(1, 3, figsize=(18, 4.8))
# A. random vs scaffold (applicability domain), RF vs RF
eps = list(q.keys()); met = ["spearman" if q[e]["random_split"].get("spearman") is not None else "AUROC" for e in eps]
rand = [q[e]["random_split"][m] for e, m in zip(eps, met)]; scaf = [q[e]["scaffold_split"][m] for e, m in zip(eps, met)]
x = np.arange(len(eps))
ax[0].bar(x - .2, rand, .38, label="random split", color="#4575b4"); ax[0].bar(x + .2, scaf, .38, label="scaffold split (new chemotypes)", color="#d73027")
ax[0].set_xticks(x, [e.split("_")[0] for e in eps], rotation=35, ha="right", fontsize=8); ax[0].axhline(0.5, ls=":", c="#999", lw=.7)
ax[0].set(ylabel="Spearman (reg) / AUROC (clf)", title="A. Applicability domain: same model, new chemotypes\n(RF fingerprint, leakage-free)", ylim=(0, 1)); ax[0].legend(fontsize=8)
# B. learning curve
for name, c in [("Lipophilicity_AstraZeneca", "#4575b4"), ("Solubility_AqSolDB", "#1a9850")]:
    cu = ct["learning_curve"][name]; xs = sorted(int(k) for k in cu); ys = [cu[str(k)] for k in xs]
    ax[1].plot(xs, ys, "o-", color=c, label=name.split("_")[0])
ax[1].set_xscale("log"); ax[1].set(xlabel="number of client training molecules", ylabel="Spearman on held-out test",
    title="B. Can few client points learn SAR?\n(fixed held-out test, growing training set)", ylim=(-0.05, 1)); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
# C. ADMET-AI vs RF on identical canonical test
eps2 = list(fc.keys()); rf = [fc[e]["RF_fingerprint"] for e in eps2]; ai = [fc[e]["ADMET_AI"] for e in eps2]
x = np.arange(len(eps2))
ax[2].bar(x - .2, rf, .38, label="RF fingerprint (baseline)", color="#91bfdb"); ax[2].bar(x + .2, ai, .38, label="ADMET-AI (pretrained D-MPNN)", color="#762a83")
ax[2].set_xticks(x, [f"{e.split('_')[0]}\n({fc[e]['metric']})" for e in eps2], rotation=35, ha="right", fontsize=7.5)
ax[2].set(ylabel="metric on TDC test split", title="C. Shipped tool vs simple baseline\n(same TDC scaffold test)", ylim=(0, 1)); ax[2].legend(fontsize=8)
fig.tight_layout(); fig.savefig(f"{OUT}/fig1_admet_qsar.png", dpi=130)

S = {"applicability_domain_random_vs_scaffold": q, "learning_curve": ct["learning_curve"],
     "cross_endpoint_transfer": ct["cross_endpoint_transfer_spearman"], "endpoint_label_correlation": ct["endpoint_label_correlation_on_shared_molecules"],
     "shipped_tool_vs_baseline": fc}
json.dump(S, open(f"{OUT}/summary.json", "w"), indent=1)
print("figure + summary written")
