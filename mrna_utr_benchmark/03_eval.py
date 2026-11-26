import os, json, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
OUT = "results"
rt = pd.read_csv(f"{OUT}/pred_random_test.csv"); de = pd.read_csv(f"{OUT}/pred_designed.csv")
def metrics(y, p):
    r = pearsonr(y, p)[0]; return {"pearson_r": round(float(r), 3), "r2": round(float(r**2), 3),
        "spearman": round(float(spearmanr(y, p)[0]), 3), "n": int(len(y))}
S = {}
S["random_test_held_out"] = {"CNN_Optimus": metrics(rt.rl, rt.pred_cnn), "biophysical_baseline": metrics(rt.rl, rt.pred_baseline)}
hum = de[de.library == "human_utrs"]
S["human_5UTRs_endogenous"] = {"CNN_Optimus": metrics(hum.rl, hum.pred_cnn), "biophysical_baseline": metrics(hum.rl, hum.pred_baseline)}
# read-depth filtered human (higher-quality measurements)
humq = hum[hum.r_total >= hum.r_total.median()]
S["human_5UTRs_high_read_depth"] = {"CNN_Optimus": metrics(humq.rl, humq.pred_cnn), "n_note": "top 50% by reads"}
# designed "target" UTRs (model-designed high-MRL): does predicted-high actually measure high?
for lib in ["sergii_utrs", "step_worst_to_best_no_uatg", "snv"]:
    g = de[de.library == lib]
    if len(g) > 50: S[f"designed_{lib}"] = metrics(g.rl, g.pred_cnn)
# uAUG mechanism (known answer): uAUG lowers MRL; does each predictor capture it?
for nm, d in [("random_test", rt), ("human", hum)]:
    has = d.uAUG_count >= 1; S.setdefault("uAUG_effect", {})[nm] = {
        "measured_MRL_no_uAUG": round(float(d[~has].rl.mean()), 2), "measured_MRL_with_uAUG": round(float(d[has].rl.mean()), 2),
        "CNN_MRL_no_uAUG": round(float(d[~has].pred_cnn.mean()), 2), "CNN_MRL_with_uAUG": round(float(d[has].pred_cnn.mean()), 2)}
json.dump(S, open(f"{OUT}/summary.json", "w"), indent=1)

fig, ax = plt.subplots(1, 3, figsize=(17, 5))
def hexp(a, y, p, title, r2):
    a.hexbin(p, y, gridsize=45, cmap="viridis", mincnt=1, bins="log")
    a.set(xlabel="predicted MRL (CNN)", ylabel="measured MRL", title=f"{title}\nPearson r²={r2}")
hexp(ax[0], rt.rl, rt.pred_cnn, f"Random library, held-out test (n={len(rt)})", S["random_test_held_out"]["CNN_Optimus"]["r2"])
hexp(ax[1], hum.rl, hum.pred_cnn, f"Real human 5'UTRs (n={len(hum)}) — never trained on", S["human_5UTRs_endogenous"]["CNN_Optimus"]["r2"])
tasks = ["random test\n(CNN)", "random test\n(baseline)", "human UTRs\n(CNN)", "human UTRs\n(baseline)"]
vals = [S["random_test_held_out"]["CNN_Optimus"]["r2"], S["random_test_held_out"]["biophysical_baseline"]["r2"],
        S["human_5UTRs_endogenous"]["CNN_Optimus"]["r2"], S["human_5UTRs_endogenous"]["biophysical_baseline"]["r2"]]
ax[2].bar(range(4), vals, color=["#4575b4", "#91bfdb", "#d73027", "#fdae61"])
for i, v in enumerate(vals): ax[2].text(i, v + .01, f"{v:.2f}", ha="center", fontsize=9)
ax[2].set_xticks(range(4), tasks, fontsize=8); ax[2].set(ylabel="variance explained (r²)", title="MRL prediction: random library vs endogenous human UTRs", ylim=(0, 1))
fig.tight_layout(); fig.savefig(f"{OUT}/fig1_mrl.png", dpi=130)
print(json.dumps(S, indent=1))
