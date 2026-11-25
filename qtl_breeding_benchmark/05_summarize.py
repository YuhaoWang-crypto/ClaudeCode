"""Collect the three experiments into results/summary.json and draw figures 1 and 3."""
import json, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
S = {}
# ---- cross
ka = pd.read_csv("results/cross_known_answer.csv"); dv = pd.read_csv("results/cross_discovery_validation.csv")
dv["replicated"] = (dv.p_valid < 0.05) & dv.same_direction
S["cross_known_answer"] = [{"trait": r.trait, "gene": r.known_gene, "gene_chr": int(r.known_chr), "gene_Mb": r.known_Mb,
    "peak_chr": int(r.peak_chr), "peak_marker_Mb": round(r.peak_marker_Mb, 2), "LOD": round(r.peak_LOD, 1),
    "interval_Mb": [round(r.ci_lo_Mb, 1), round(r.ci_hi_Mb, 1)], "gene_inside_interval": bool(r.ci_lo_Mb <= r.known_Mb <= r.ci_hi_Mb),
    "best_LOD_on_any_other_chr": round(r.second_best_chr_LOD, 1)} for r in ka.itertuples()]
S["cross_insulin_discovery_validation"] = {"splits": 20, "n_discovery": 272, "lod_threshold_5pct": round(float(dv.threshold[0]), 2),
    "splits_with_any_significant_qtl": int(dv.split.nunique()), "significant_qtl_calls": int(len(dv)),
    "calls_replicated_in_unseen_half": int(dv.replicated.sum()),
    "mean_variance_explained_discovery": round(float(dv.r2_disc.mean()), 3),
    "mean_variance_explained_validation": round(float(dv.r2_valid.mean()), 3),
    "calls_by_chr": dv.chr.astype(str).value_counts().to_dict()}
full = pd.read_csv("results/cross_insulin_fullscan.csv")
fig, ax = plt.subplots(1, 2, figsize=(13, 4), gridspec_kw={"width_ratios": [2.2, 1]})
x = np.arange(len(full)); ax[0].plot(x, full.lod, lw=.9, c="#4575b4")
b = full.groupby("chr", sort=False).size().cumsum(); ax[0].set_xticks(b - full.groupby("chr", sort=False).size() / 2, b.index, fontsize=7)
for v in b[:-1]: ax[0].axvline(v, c="#ddd", lw=.5)
ax[0].axhline(dv.threshold[0], ls=":", c="k", lw=.8); ax[0].set(xlabel="chromosome", ylabel="LOD",
    title=f"Insulin, all 544 F2 mice (dotted: 5% threshold for a 272-mouse half, LOD {dv.threshold[0]:.2f})")
ax[1].scatter(dv.r2_disc * 100, dv.r2_valid * 100, c=np.where(dv.replicated, "#1a9850", "#d73027"), s=30)
ax[1].plot([0, 10], [0, 10], "k:", lw=.8); ax[1].set(xlabel="% variance explained, discovery half", ylabel="% variance explained, unseen half",
    title=f"{len(dv)} significant calls in 20 splits\ngreen = replicated ({int(dv.replicated.sum())})", xlim=(0, 10), ylim=(0, 10))
fig.tight_layout(); fig.savefig("results/fig1_cross_qtl.png", dpi=130)
# ---- rice GWAS
S["rice_gwas"] = json.load(open("results/rice_gwas.json"))
# ---- genomic selection
g = pd.read_csv("results/genomic_selection.csv")
t = g.groupby(["data", "scenario", "model", "trait"]).r.agg(["mean", "min", "max", "size"]).round(3).reset_index()
S["genomic_selection"] = t.to_dict("records"); t.to_csv("results/genomic_selection_summary.csv", index=False)
fig, ax = plt.subplots(1, 2, figsize=(14, 4.6), gridspec_kw={"width_ratios": [1, 1.5]})
w = g[g.data == "wheat"]; envs = sorted(w.trait.unique())
same = [w[(w.trait == e) & (w.model == "GBLUP") & (w.scenario.str.contains("same"))].r.mean() for e in envs]
bb = [w[(w.trait == e) & (w.model == "BayesB")].r.mean() for e in envs]
diff = [w[(w.trait == e) & (w.scenario.str.contains("different"))].r.mean() for e in envs]
k = np.arange(4)
for off, v, lab, c in [(-.27, same, "GBLUP, same environment", "#4575b4"), (0, bb, "BayesB, same environment", "#91bfdb"), (.27, diff, "GBLUP, trained in this env, scored in the other three", "#d73027")]:
    ax[0].bar(k + off, v, .26, label=lab, color=c)
ax[0].set_xticks(k, envs); ax[0].axhline(0, c="k", lw=.6); ax[0].legend(fontsize=7.5); ax[0].set(ylabel="prediction accuracy r (unseen lines)", title="Wheat grain yield, 599 lines")
r = g[g.data == "rice"]; tr = ["Flowering time at Arkansas", "Plant height", "Seed length", "Amylose content"]
sc = [("unseen lines, random split (whole panel r)", "GBLUP, random split, r over whole panel", "#4575b4"),
      ("unseen lines, random split (r within subpopulation)", "GBLUP, random split, r within subpopulation", "#91bfdb"),
      ("unseen subpopulation", "GBLUP, whole subpopulation unseen", "#d73027"),
      ("10 top markers chosen in training fold", "10 top SNPs, r over whole panel", "#fdae61"),
      ("10 top markers chosen in training fold (r within subpopulation)", "10 top SNPs, r within subpopulation", "#fee090")]
for j, (s, lab, c) in enumerate(sc):
    ax[1].bar(k + (j - 2) * .16, [r[(r.trait == x) & (r.scenario == s)].r.mean() for x in tr], .155, label=lab, color=c)
ax[1].set_xticks(k, [x.replace(" at Arkansas", "") for x in tr]); ax[1].legend(fontsize=7.5); ax[1].axhline(0, c="k", lw=.6)
ax[1].set(ylabel="prediction accuracy r", title="Rice Diversity Panel, 413 accessions")
fig.tight_layout(); fig.savefig("results/fig3_genomic_selection.png", dpi=130)
json.dump(S, open("results/summary.json", "w"), indent=1)
print(json.dumps(S["cross_known_answer"]), json.dumps(S["cross_insulin_discovery_validation"]))
print(t.to_string())
