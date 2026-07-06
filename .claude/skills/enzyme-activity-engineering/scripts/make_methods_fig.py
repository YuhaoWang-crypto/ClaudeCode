import json, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

rows = json.load(open("data/ranked_mutations.json"))
esm = json.load(open("data/esm_scores.json"))

# focus set: top-5 + a few notable, keep order by consensus
FOCUS = ["L168K", "H193R", "F231Y", "A248K", "V257A", "M268L", "S209A", "K219E", "S45A", "F231W"]
by = {r["mutation"]: r for r in rows}
sel = [by[m] for m in FOCUS if m in by]

methods = ["ESM-1v", "ESM2", "ProteinMPNN", "ESM-IF"]
# recompute per-method percentile the same way as rank.py for the heatmap
def build(name, key=None):
    if name in ("ESM-1v", "ESM2"):
        src = esm["esm1v_t33_650M_UR90S_1" if name == "ESM-1v" else "esm2_t33_650M_UR50D"]["scores"]
        d = {(int(p), mut): v for p, rec in src.items() for mut, v in rec["muts"].items()}
    return d

# use raw_llr stored in rows for coloring (sign & magnitude) + consensus percentile
INK = "#12141c"; FG = "#e8eaf0"
plt.rcParams.update({"figure.facecolor": INK, "axes.facecolor": INK, "savefig.facecolor": INK,
                     "text.color": FG, "axes.labelcolor": FG, "xtick.color": FG, "ytick.color": FG})

M = np.array([[ (by[m]["raw_llr"][meth] if by[m]["raw_llr"][meth] is not None else np.nan)
                for meth in methods] for m in FOCUS if m in by])
labels = [m for m in FOCUS if m in by]

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13, 6), gridspec_kw={"width_ratios": [3, 1]})
vmax = np.nanmax(np.abs(M))
im = ax.imshow(M, cmap="RdYlGn", vmin=-vmax, vmax=vmax, aspect="auto")
ax.set_xticks(range(len(methods))); ax.set_xticklabels(methods, rotation=20, ha="right", fontsize=11)
ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=11)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        v = M[i, j]
        if not np.isnan(v):
            ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=8.5,
                    color="black" if abs(v) < vmax*0.6 else "white")
ax.set_title("Per-method mutation score (log-likelihood ratio, mut vs WT)\ngreen = favored over wild type", fontsize=11)
cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02); cb.ax.tick_params(colors=FG)

cons = [by[m]["consensus_pct"] for m in labels]
y = np.arange(len(labels))[::-1]
ax2.barh(y, cons, color="#a3e635")
ax2.set_yticks(y); ax2.set_yticklabels(labels, fontsize=10)
ax2.set_xlim(0.9, 1.0); ax2.set_xlabel("consensus percentile", fontsize=10)
ax2.set_title("Cross-method\nconsensus rank", fontsize=11)
for yi, c in zip(y, cons):
    ax2.text(c-0.001, yi, f"{c:.3f}", ha="right", va="center", fontsize=9, color="black")
fig.suptitle("IS621 — cross-method comparison (ESM-1v · ESM2 · ProteinMPNN · ESM-IF)", fontsize=13, y=1.02)
fig.tight_layout(); fig.savefig("figures/fig3_methods.png", dpi=150, bbox_inches="tight"); plt.close(fig)
print("wrote figures/fig3_methods.png")
