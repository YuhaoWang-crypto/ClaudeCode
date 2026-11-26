"""Figures + combined summary for the risk-stratification row."""
import json, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT = "results"
prs = pd.read_csv(f"{OUT}/prs_scores.csv", index_col=0); psum = json.load(open(f"{OUT}/prs_summary.json"))
R = pd.read_csv(f"{OUT}/prognostic_splits.csv"); TR = pd.read_csv(f"{OUT}/prognostic_transfer.csv")
nul = pd.read_csv(f"{OUT}/prognostic_null.csv")
order = ["EUR", "SAS", "AMR", "EAS", "AFR"]; col = {"EUR": "#4575b4", "SAS": "#91bfdb", "AMR": "#fee090", "EAS": "#fc8d59", "AFR": "#d73027"}
traits = ["T2D", "Height", "CAD", "LDL"]
name2pgs = {"T2D": "PGS000020", "LDL": "PGS002605", "Height": "PGS000297", "CAD": "PGS000337"}

fig, ax = plt.subplots(1, 4, figsize=(18, 4.3))
for a, t in zip(ax, traits):
    eur = prs[prs.super_pop == "EUR"][name2pgs[t]] if name2pgs[t] in prs else None
    colname = name2pgs[t]
    mu = prs[prs.super_pop == "EUR"][colname].mean(); sd = prs[prs.super_pop == "EUR"][colname].std()
    for sp in order:
        z = (prs[prs.super_pop == sp][colname] - mu) / sd
        a.hist(z, np.linspace(-5, 5, 45), density=True, histtype="step", lw=2, color=col[sp], label=sp)
    a.axvline(0, c="#bbb", lw=.6)
    a.set(title=f"{t} PRS (public weights)\nz vs European distribution", xlabel="PRS (SD from European mean)", ylabel="density")
    a.legend(fontsize=7, title="1000G ancestry")
fig.suptitle(f"Same public PRS weights, different ancestries — partial genome (chr {psum['chromosomes_used']}); a raw score is uninterpretable without an ancestry-matched reference", fontsize=10)
fig.tight_layout(); fig.savefig(f"{OUT}/fig1_prs_ancestry.png", dpi=130)

# prognostic panel figure
fig, ax = plt.subplots(1, 2, figsize=(14, 4.6), gridspec_kw={"width_ratios": [1.3, 1]})
ag = R.groupby("cancer").mean(numeric_only=True).loc[["LUAD", "LUSC", "KIRC", "BRCA"]]
x = np.arange(len(ag)); w = .2
for i, (c, lab, cc) in enumerate([("c_train", "training (apparent)", "#999"), ("c_test_panel", "panel, held-out", "#4575b4"),
                                   ("c_test_clinical", "age+stage, held-out", "#1a9850"), ("c_test_both", "panel+age+stage", "#d73027")]):
    ax[0].bar(x + (i - 1.5) * w, ag[c], w, label=lab, color=cc)
ax[0].axhline(.5, ls=":", c="k", lw=.8); ax[0].set_xticks(x, ag.index); ax[0].set(ylim=(0.45, 0.88), ylabel="C-index (OS)",
    title="Expression prognostic panel: apparent vs held-out, and vs clinical variables"); ax[0].legend(fontsize=8)
piv = TR.groupby(["train", "test"]).c_index.mean().unstack().loc[["LUAD", "LUSC", "KIRC", "BRCA"], ["LUAD", "LUSC", "KIRC", "BRCA"]]
im = ax[1].imshow(piv.values, cmap="RdBu_r", vmin=0.35, vmax=0.65)
ax[1].set_xticks(range(4), piv.columns); ax[1].set_yticks(range(4), piv.index); ax[1].set(xlabel="applied to (test cohort)", ylabel="trained on")
for i in range(4):
    for j in range(4):
        v = piv.values[i, j]; ax[1].text(j, i, "—" if np.isnan(v) else f"{v:.2f}", ha="center", va="center", fontsize=9, color="w" if abs(v-.5) > .1 else "k")
ax[1].set_title("Locked panel transferred across cancers\n(diagonal blank = same-cohort, see left)")
fig.tight_layout(); fig.savefig(f"{OUT}/fig2_prognostic.png", dpi=130)

S = {"prs_ancestry": psum,
     "prognostic_panel": {"per_cancer_mean": R.groupby("cancer").mean(numeric_only=True).round(3).to_dict("index"),
        "cross_cancer_transfer_cindex": TR.groupby(["train","test"]).c_index.mean().round(3).reset_index().to_dict("records"),
        "null_permuted_cindex": nul.to_dict("records")}}
json.dump(S, open(f"{OUT}/summary.json", "w"), indent=1, default=str)
print("figures written")
