"""AlphaGenome v2 demo: score_variant with the recommended variant scorers.

`predict_variant` returns raw REF/ALT tracks; `score_variant` instead reduces
them server-side to one number per (track) or per (gene, track), using scorers
such as "log2 fold change of summed signal in a 501 bp window" or "log fold
change of mean RNA-seq signal over each gene's exons".

Scores chr22:36201698 A>C in transverse colon with all recommended scorers.
Requires ALPHAGENOME_API_KEY.

    python alphagenome_demo/demo_v2_score.py
"""

import os
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client, variant_scorers

OUT_DIR = pathlib.Path(__file__).parent
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 30)
pd.set_option("display.max_colwidth", 45)

model = dna_client.create(os.environ["ALPHAGENOME_API_KEY"], timeout=60)

interval = genome.Interval(chromosome="chr22", start=35677410, end=36725986)
variant = genome.Variant(
    chromosome="chr22", position=36201698, reference_bases="A", alternate_bases="C"
)
CONDITION = "UBERON:0001157"  # transverse colon

recommended = variant_scorers.get_recommended_scorers(
    dna_client.Organism.HOMO_SAPIENS.to_proto()
)
names = {v: k for k, v in variant_scorers.RECOMMENDED_VARIANT_SCORERS.items()}
print(f"{len(recommended)} recommended scorers:", [names[s] for s in recommended])

# The server caps scorers per request, so score in chunks.
chunk = dna_client.MAX_VARIANT_SCORERS_PER_REQUEST
results = []
for i in range(0, len(recommended), chunk):
    batch = recommended[i : i + chunk]
    scores = model.score_variant(interval, variant, CONDITION, variant_scorers=batch)
    for scorer, ad in zip(batch, scores):
        print(f"  {names[scorer]:28s} AnnData {ad.shape} (rows x tracks)")
        tidy = variant_scorers.tidy_scores([ad])
        if tidy is not None and len(tidy):
            tidy.insert(0, "scorer", names[scorer])
            results.append(tidy)

df = pd.concat(results, ignore_index=True)
df.to_csv(OUT_DIR / "v2_variant_scores.csv", index=False)
print("\nTidy table:", df.shape)
print("Columns:", list(df.columns))

# v2 returns raw scores only (no quantile_score as in v1), so scores are
# comparable within a scorer but not across scorers.
df["abs_score"] = df["raw_score"].abs()
summary = df.groupby("scorer", sort=False).agg(
    n_scores=("raw_score", "size"),
    min_raw=("raw_score", "min"),
    max_raw=("raw_score", "max"),
)
print("\nPer-scorer summary (raw scores):")
print(summary.round(3).to_string())

LFC_SCORERS = ["RNA_SEQ", "SC_GENE_EXPRESSION", "CAGE", "ATAC", "DNASE", "SC_ATAC", "CHIP_HISTONE", "CHIP_TF", "POLYADENYLATION", "SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS", "METHYLATION", "CONTACT_MAPS"]
cols = ["gene_name", "histone_mark", "transcription_factor", "cpg_position", "track_name", "track_strand", "imputation_mode", "raw_score"]
for name in LFC_SCORERS:
    sub = df[df.scorer == name].sort_values("abs_score", ascending=False).head(4)
    keep = [c for c in cols if c in sub.columns and sub[c].notna().any()]
    print(f"\n{name}: top scores")
    print(sub[keep].round(3).to_string(index=False))

# Figure: gene-level RNA-seq effects and the most affected TF ChIP tracks.
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
rna = df[df.scorer == "RNA_SEQ"]
gene = rna.groupby("gene_name")["raw_score"].mean()
gene = gene.reindex(gene.abs().sort_values(ascending=False).index).head(10)[::-1]
axes[0].barh(gene.index, gene.values, color=["tab:red" if v < 0 else "tab:blue" for v in gene.values])
axes[0].axvline(0, color="black", lw=0.8)
axes[0].set_xlabel("RNA_SEQ gene log fold change (mean over tracks)")
axes[0].set_title("Genes most affected (GeneMaskLFCScorer)")

tf = df[df.scorer == "CHIP_TF"].sort_values("abs_score", ascending=False).head(12)[::-1]
tf_labels = tf["transcription_factor"].astype(str) + " (" + tf["imputation_mode"].astype(str).str.replace("imputed_", "imp. ") + ")"
axes[1].barh(range(len(tf)), tf["raw_score"], color=["tab:red" if v < 0 else "tab:blue" for v in tf["raw_score"]])
axes[1].set_yticks(range(len(tf)), tf_labels, fontsize=8)
axes[1].axvline(0, color="black", lw=0.8)
axes[1].set_xlabel("CHIP_TF log2 fold change, 501 bp window")
axes[1].set_title("TF binding tracks most affected (CenterMaskScorer)")
fig.suptitle(f"score_variant: {variant.chromosome}:{variant.position} A>C, transverse colon (AlphaGenome v2)")
fig.tight_layout()
out_png = OUT_DIR / "v2_variant_scores.png"
fig.savefig(out_png, dpi=150, bbox_inches="tight")
print(f"\nSaved {out_png}")
