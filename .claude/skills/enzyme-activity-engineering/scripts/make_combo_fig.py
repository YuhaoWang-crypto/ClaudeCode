import json, itertools
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

recs = json.load(open("data/combo_ranked.json"))
byset = {frozenset(r["pos"]): r for r in recs}
SITES = [168, 193, 231, 248, 257]
LAB = {168: "L168K", 193: "H193R", 231: "F231Y", 248: "A248K", 257: "V257A"}

# pairwise min atom distances (recomputed)
byres = {}
for line in open("data/8WT9.pdb"):
    if line.startswith("ATOM") and line[21] == "A":
        r = int(line[22:26]); byres.setdefault(r, []).append(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])])
def mind(a, b):
    A = np.array(byres[a]); B = np.array(byres[b])
    return float(np.linalg.norm(A[:, None, :] - B[None, :, :], axis=-1).min())

INK = "#12141c"; FG = "#e8eaf0"
plt.rcParams.update({"figure.facecolor": INK, "axes.facecolor": INK, "savefig.facecolor": INK,
                     "text.color": FG, "axes.labelcolor": FG, "xtick.color": FG, "ytick.color": FG})
fig, (axE, axB) = plt.subplots(1, 2, figsize=(15, 6.6), gridspec_kw={"width_ratios": [1, 1.25]})

# --- Panel A: pairwise epistasis matrix (lower tri) + distance (upper tri) ---
n = len(SITES)
M = np.full((n, n), np.nan)
for i, a in enumerate(SITES):
    for j, b in enumerate(SITES):
        if i > j:  # lower: epistasis of the double
            M[i, j] = byset[frozenset([a, b])]["epistasis"]
        elif i < j:  # upper: distance (put NaN, annotate separately)
            M[i, j] = np.nan
im = axE.imshow(M, cmap="PuOr", vmin=-0.6, vmax=0.6)
for i, a in enumerate(SITES):
    for j, b in enumerate(SITES):
        if i > j:
            axE.text(j, i, f"{M[i,j]:+.2f}", ha="center", va="center", fontsize=11,
                     color="#e8eaf0")
        elif i < j:
            d = mind(a, b)
            col = "#a3e635" if d < 12 else "#5b6474"
            axE.text(j, i, f"{d:.0f}Å", ha="center", va="center", fontsize=10, color=col)
        else:
            axE.text(j, i, "—", ha="center", va="center", color="#5b6474")
axE.set_xticks(range(n)); axE.set_xticklabels([LAB[s] for s in SITES], rotation=30, ha="right", fontsize=10)
axE.set_yticks(range(n)); axE.set_yticklabels([LAB[s] for s in SITES], fontsize=10)
axE.set_title("Pairwise interaction\nlower ▽ = ESM epistasis (mut-context − additive)\nupper △ = min atom distance",
              fontsize=10.5, color=FG)
cb = fig.colorbar(im, ax=axE, fraction=0.046, pad=0.03); cb.ax.tick_params(colors=FG)
cb.set_label("epistasis (log-lik units)", color=FG, fontsize=9)

# --- Panel B: recommended combination panel ranked by ESM joint ---
panel = ["L168K+F231Y", "L168K+A248K", "L168K+H193R",
         "L168K+H193R+F231Y", "L168K+F231Y+A248K", "L168K+H193R+F231Y+V257A",
         "L168K+F231Y+A248K+V257A", "L168K+H193R+F231Y+A248K+V257A"]
pr = [next(r for r in recs if r["key"] == k) for k in panel]
pr.sort(key=lambda r: r["esm_joint"])
y = np.arange(len(pr))
divcol = {1: "#6b7280", 2: "#38bdf8", 3: "#a3e635", 4: "#fbbf24"}
cols = [divcol[r["ndiv"]] for r in pr]
axB.barh(y, [r["esm_joint"] for r in pr], color=cols, edgecolor="#0a0c12")
for yi, r in zip(y, pr):
    axB.text(r["esm_joint"] + 0.08, yi, f"{r['esm_joint']:+.1f}  (epi {r['epistasis']:+.2f})",
             va="center", fontsize=8.5, color=FG)
axB.set_yticks(y); axB.set_yticklabels([r["key"].replace("+", "+\n") if r["n"] >= 4 else r["key"] for r in pr],
                                       fontsize=8.5, fontfamily="monospace")
axB.set_xlabel("epistasis-aware ESM joint score (higher = more favored)", fontsize=10)
axB.set_title("Recommended combination panel", fontsize=11, color=FG)
axB.set_xlim(0, 10)
from matplotlib.patches import Patch
axB.legend(handles=[Patch(color=divcol[k], label=f"{k} mechanism{'s' if k>1 else ''}") for k in (1,2,3,4)],
           fontsize=8.5, facecolor="#1c1f2b", edgecolor="#333", labelcolor=FG, loc="lower right",
           title="activity routes covered", title_fontsize=8.5)
for ax in (axE, axB):
    for sp in ax.spines.values(): sp.set_color("#333")
fig.suptitle("IS621 — combination-mutation analysis (five spatially-independent sites → near-additive, mild positive synergy)",
             fontsize=13, y=1.0)
fig.tight_layout(rect=[0, 0, 1, 0.96]); fig.savefig("figures/fig4_combinations.png", dpi=150, bbox_inches="tight")
print("wrote figures/fig4_combinations.png")
