"""Evaluate SpliceAI (pretrained, zero-shot) on ClinVar splice variants.
Known-answer: pathogenic canonical splice variants should score high, benign low.
Boundary: deep-intronic pathogenic variants (non-canonical) are the hard case the matrix flags."""
import os, re, json, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, roc_curve
V = os.environ["VB_DATA"]; OUT = "results"
sel = pd.read_parquet(f"{V}/variants_sel.parquet")
sel["key"] = sel.chrom + ":" + sel.pos.astype(str) + ":" + sel.ref + ":" + sel.alt

rows = {}
for l in open(f"{V}/spliceai_out.vcf"):
    if l.startswith("#"): continue
    f = l.rstrip("\n").split("\t"); key = f"{f[0]}:{f[1]}:{f[3]}:{f[4]}"
    m = re.search(r"SpliceAI=([^;]+)", f[7])
    if not m: continue
    for ann in m.group(1).split(","):
        p = ann.split("|")
        if len(p) < 10 or p[2] == ".": continue
        ds = [float(x) for x in p[2:6]]; dp = [int(x) for x in p[6:10]]
        dmax = max(ds); dist = min(abs(d) for d, s in zip(dp, ds))   # distance of the strongest event
        dist_strong = abs(dp[int(np.argmax(ds))])
        rows[key] = max(rows.get(key, (-1, 0))[0:1] + (0,), (dmax, dist_strong)) if key in rows else (dmax, dist_strong)
sa = pd.DataFrame([(k, v[0], v[1]) for k, v in rows.items()], columns=["key", "delta_max", "dist_to_event"])
d = sel.merge(sa, on="key", how="left")
d["scored"] = d.delta_max.notna(); d["delta_max"] = d.delta_max.fillna(0.0)
d.to_csv(f"{OUT}/spliceai_scored.csv", index=False)

def auc(pos, neg):
    y = np.r_[np.ones(len(pos)), np.zeros(len(neg))]; s = np.r_[pos, neg]
    return round(float(roc_auc_score(y, s)), 3) if len(pos) > 5 and len(neg) > 5 else None
ben = d[d.cls == "benign"]
S = {"n_variants": int(len(d)), "n_scored_by_spliceai": int(d.scored.sum()),
     "groups": {g: int((d.grp == g).sum()) for g in d.grp.unique()}}
pc = d[d.grp == "path_canonical_splice"]; pi = d[d.grp == "path_intronic_noncanonical"]
for thr in [0.2, 0.5, 0.8]:
    S[f"recall_at_delta_{thr}"] = {"path_canonical_splice": round(float((pc.delta_max >= thr).mean()), 3),
        "path_intronic_noncanonical": round(float((pi.delta_max >= thr).mean()), 3),
        "benign_false_positive_rate": round(float((ben.delta_max >= thr).mean()), 3)}
S["AUROC"] = {"canonical_splice_vs_benign": auc(pc.delta_max, ben.delta_max),
              "intronic_noncanonical_vs_benign": auc(pi.delta_max, ben.delta_max),
              "all_pathogenic_vs_benign": auc(d[d.cls == "path"].delta_max, ben.delta_max)}
S["median_delta"] = {g: round(float(d[d.grp == g].delta_max.median()), 3) for g in d.grp.unique()}
# distance dependence for pathogenic variants
pth = d[d.cls == "path"].copy(); pth["distbin"] = pd.cut(pth.dist_to_event.fillna(999), [-1, 2, 10, 50, 500, 1e9],
    labels=["0-2 (canonical)", "3-10", "11-50", "51-500", ">500"])
S["pathogenic_detection_by_distance_to_splice_event"] = {str(k): {"n": int(len(g)), "median_delta": round(float(g.delta_max.median()), 3),
    "pct_detected_delta>=0.5": round(float((g.delta_max >= .5).mean()) * 100, 1)} for k, g in pth.groupby("distbin", observed=True)}
json.dump(S, open(f"{OUT}/summary.json", "w"), indent=1)

fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
b = np.linspace(0, 1, 30); labels = {"path_canonical_splice": ("#b2182b", "pathogenic canonical splice"),
    "path_intronic_noncanonical": ("#ef8a62", "pathogenic intronic (non-canonical)"),
    "benign_synonymous": ("#67a9cf", "benign synonymous"), "benign_other": ("#999", "benign other")}
for g, (c, lab) in labels.items():
    ax[0].hist(d[d.grp == g].delta_max, b, density=True, histtype="step", lw=2, color=c, label=lab)
ax[0].axvline(0.5, ls=":", c="k", lw=.8, label="SpliceAI 0.5 threshold")
ax[0].set(title="SpliceAI delta score by ClinVar class (chr 1/2/17)", xlabel="SpliceAI delta_max", ylabel="density"); ax[0].legend(fontsize=7.5)
dd = S["pathogenic_detection_by_distance_to_splice_event"]; xs = list(dd.keys()); ys = [dd[k]["pct_detected_delta>=0.5"] for k in xs]; ns = [dd[k]["n"] for k in xs]
ax[1].bar(range(len(xs)), ys, color="#b2182b"); ax[1].set_xticks(range(len(xs)), [f"{x}\n(n={n})" for x, n in zip(xs, ns)], fontsize=8)
ax[1].set(title="Pathogenic splice variants detected (delta>=0.5)\nby distance to the predicted splice event (bp)", ylabel="% detected", ylim=(0, 100))
fig.tight_layout(); fig.savefig(f"{OUT}/fig1_spliceai.png", dpi=130)
print(json.dumps(S, indent=1))
