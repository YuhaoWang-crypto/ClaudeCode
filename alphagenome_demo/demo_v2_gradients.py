"""AlphaGenome v2 demo: input gradients (predict_interval_gradients).

One backward pass gives a per-base attribution score for a gene's predicted
expression over the whole input window (here 131 kb), instead of the 3 forward
passes per base that saturation mutagenesis needs.

Target: APOL4 RNA-seq in transverse colon. The script reports where the
attribution concentrates and cross-checks it against the ISM scores from
demo_v2_ism.py (v2_ism_scores.csv) on the same 80 bp.

Requires ALPHAGENOME_API_KEY.

    python alphagenome_demo/demo_v2_gradients.py
"""

import os
import pathlib
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client
from alphagenome.visualization import plot_components

OUT_DIR = pathlib.Path(__file__).parent
pd.set_option("display.width", 200)

GENE_ID = "ENSG00000100336"  # APOL4
CONDITION = "UBERON:0001157"  # transverse colon

model = dna_client.create(os.environ["ALPHAGENOME_API_KEY"], timeout=60)

t0 = time.time()
out = model.predict_interval_gradients(
    GENE_ID,
    requested_output=dna_client.OutputType.RNA_SEQ,
    sequence_length=dna_client.SEQUENCE_LENGTH_100KB,
    condition=CONDITION,
)
print(f"Gradients returned in {time.time() - t0:.1f} s")
gene = out.target_gene.interval
print("Target gene:", out.gene_name, gene, "| biosample:", out.biosample_name)
print("Input interval:", out.interval, "| gradients shape:", out.gradients.shape)
print("scalar_score (log of predicted expression):", round(out.scalar_score, 3))
print("Tracks aggregated in the objective:")
print(out.track_metadata[["name", "strand"]].to_string())

grad = out.gradients  # (L, 4): gradient x one-hot input, non-zero at the reference base
score = grad.sum(axis=1)  # per-position attribution
start = out.interval.start
pos1 = np.arange(start + 1, start + 1 + len(score))  # 1-based coordinates
np.save(OUT_DIR / "v2_gradients_APOL4.npy", score.astype(np.float32))

# ------------------------------------------------- where is the attribution?
abs_score = np.abs(score)
total = abs_score.sum()
tss = gene.end if gene.strand == "-" else gene.start  # gene start in transcription direction
print(f"\nGene TSS (annotated gene boundary): {gene.chromosome}:{tss}")
in_gene = (pos1 > gene.start) & (pos1 <= gene.end)
print(f"Share of |attribution| inside the gene body ({gene.width} bp = "
      f"{gene.width / len(score):.0%} of the window): {abs_score[in_gene].sum() / total:.0%}")
for w in [100, 500, 2000, 10000]:
    near = np.abs(pos1 - tss) <= w
    print(f"  within +/-{w:>5} bp of the TSS: {abs_score[near].sum() / total:.1%}")
order = np.argsort(-abs_score)
for frac in [0.1, 0.5, 0.9]:
    n = int(np.searchsorted(np.cumsum(abs_score[order]) / total, frac)) + 1
    print(f"  {frac:.0%} of |attribution| is carried by the top {n} bases ({n / len(score):.2%})")

# Top attribution hotspots: 20 bp bins.
BIN = 20
nb = len(score) // BIN
binned = abs_score[: nb * BIN].reshape(nb, BIN).sum(axis=1)
hot = pd.DataFrame({"bin_start": start + 1 + np.arange(nb) * BIN, "abs_attr": binned})
hot["share"] = hot.abs_attr / total
hot["dist_to_tss"] = hot.bin_start + BIN // 2 - tss
hot["in_gene"] = (hot.bin_start > gene.start) & (hot.bin_start <= gene.end)
print(f"\nTop 10 hotspots ({BIN} bp bins):")
print(hot.sort_values("abs_attr", ascending=False).head(10).round(4).to_string(index=False))

# ------------------------------------------------- cross-check against ISM
ism_csv = OUT_DIR / "v2_ism_scores.csv"
ism_interval = None
if ism_csv.exists():
    ism = pd.read_csv(ism_csv).groupby("position")["RNA_SEQ"].mean()
    # ISM gives the effect of *mutating* a base (negative = base is needed);
    # the gradient gives the contribution of *having* it, so flip the sign.
    ism_importance = -ism
    g = pd.Series(score[ism.index.to_numpy() - (start + 1)], index=ism.index)
    r = np.corrcoef(g, ism_importance)[0, 1]
    rho = g.rank().corr(ism_importance.rank())
    print(f"\nGradient vs ISM over {len(ism)} bp (chr22:{ism.index.min()}-{ism.index.max()}):")
    print(f"  Pearson r = {r:.3f}, Spearman rho = {rho:.3f}")
    cmp = pd.DataFrame({"gradient": g, "ism_importance": ism_importance})
    print("  10 positions with the largest ISM effect:")
    print(cmp.reindex(cmp.ism_importance.abs().sort_values(ascending=False).index).head(10).round(3).to_string())
    top_g = set(g.abs().sort_values(ascending=False).head(6).index)
    top_i = set(ism_importance.abs().sort_values(ascending=False).head(6).index)
    print(f"  overlap of top-6 positions by |score|: {len(top_g & top_i)} / 6")
    ism_interval = genome.Interval("chr22", int(ism.index.min()) - 1, int(ism.index.max()))

# ------------------------------------------------------------------ figures
fig, ax = plt.subplots(figsize=(14, 3.2))
ax.plot(pos1, score, lw=0.4, color="tab:blue")
ax.axvspan(gene.start, gene.end, color="orange", alpha=0.15, label=f"{out.gene_name} gene body")
ax.axvline(tss, color="red", lw=0.8, ls="--", label="TSS")
ax.set_xlabel(f"{gene.chromosome} position")
ax.set_ylabel("gradient x input")
ax.set_title(f"Input gradients of {out.gene_name} RNA-seq, {out.biosample_name} (AlphaGenome v2, 131 kb window)")
ax.legend(loc="upper left")
fig.tight_layout()
fig.savefig(OUT_DIR / "v2_gradients_overview.png", dpi=150, bbox_inches="tight")

if ism_interval is not None:
    sl = slice(ism_interval.start - start, ism_interval.end - start)
    ism_full = pd.read_csv(ism_csv)
    ism_mat = np.zeros((ism_interval.width, 4), dtype=np.float32)
    ref_idx = grad[sl].__abs__().argmax(axis=1)
    ism_mat[np.arange(ism_interval.width), ref_idx] = (-ism_full.groupby("position")["RNA_SEQ"].mean()).to_numpy()
    plot_components.plot(
        [
            plot_components.SeqLogo(scores=grad[sl].astype(np.float32), scores_interval=ism_interval, ylabel="gradient x input"),
            plot_components.SeqLogo(scores=ism_mat, scores_interval=ism_interval, ylabel="ISM (-mean LFC)"),
        ],
        interval=ism_interval,
        fig_width=16,
    )
    plt.savefig(OUT_DIR / "v2_gradients_vs_ism.png", dpi=150, bbox_inches="tight")
print("\nSaved v2_gradients_overview.png and v2_gradients_vs_ism.png")
