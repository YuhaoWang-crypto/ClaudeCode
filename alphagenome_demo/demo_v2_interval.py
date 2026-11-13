"""AlphaGenome v2 demo: score_interval (gene-level scores without a variant).

`score_interval` reduces the predicted tracks for an interval to one number
per (gene, track), e.g. mean predicted RNA-seq signal over each gene. No
variant is involved, so this is "what does the model think is expressed here,
in this tissue?".

Scores the same 1 Mb chr22 window in several tissues and compares the
predicted tissue pattern of each gene. Requires ALPHAGENOME_API_KEY.

    python alphagenome_demo/demo_v2_interval.py
"""

import os
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client, interval_scorers, variant_scorers

OUT_DIR = pathlib.Path(__file__).parent
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 30)
OT = dna_client.OutputType

model = dna_client.create(os.environ["ALPHAGENOME_API_KEY"], timeout=60)
interval = genome.Interval(chromosome="chr22", start=35677410, end=36725986)

# Pick tissues by name among conditions with *observed* RNA-seq.
TISSUES = [
    "transverse colon",
    "liver",
    "heart left ventricle",
    "gastrocnemius medialis",
    "cerebellum",
    "frontal cortex",
    "spleen",
    "upper lobe of left lung",
    "thyroid gland",
    "testis",
]
supported = model.supported_conditions()
observed = supported[(supported.rna_seq == "observed") & (supported.biosample_type == "tissue")]
conditions = {}
for tissue in TISSUES:
    hit = observed[observed.biosample_name == tissue]
    if len(hit):
        conditions[tissue] = hit.condition.iloc[0]
    else:
        print(f"(no observed RNA-seq condition named {tissue!r}, skipped)")
print("Conditions:", conditions)

Agg = interval_scorers.IntervalAggregationType
scorers = {
    "RNA_SEQ": interval_scorers.RECOMMENDED_INTERVAL_SCORERS["RNA_SEQ"],
    "CAGE": interval_scorers.GeneMaskScorer(OT.CAGE, 200_001, Agg.MEAN),
    "DNASE": interval_scorers.GeneMaskScorer(OT.DNASE, 200_001, Agg.MEAN),
}
print("Scorers:", {k: str(v) for k, v in scorers.items()})

tables = []
for tissue, condition in conditions.items():
    scores = model.score_interval(interval, condition, interval_scorers=list(scorers.values()))
    for name, ad in zip(scorers, scores):
        if not tables or name not in {t["scorer"].iloc[0] for t in tables}:
            print(f"  {name}: AnnData {ad.shape} (genes x tracks)")
        tidy = variant_scorers.tidy_scores([ad])
        tidy.insert(0, "scorer", name)
        tidy.insert(1, "tissue", tissue)
        tables.append(tidy)
    print(f"done: {tissue}")

df = pd.concat(tables, ignore_index=True)
df.to_csv(OUT_DIR / "v2_interval_scores.csv", index=False)
print("\nTidy table:", df.shape, "| columns:", list(df.columns))

# Gene x tissue matrix from the unstranded GTEx polyA+ track, which exists for
# every tissue; where GTEx has no data for a tissue the track is imputed, and
# those columns are marked with "*".
rna = df[df.scorer == "RNA_SEQ"]
gtex = rna[rna.track_name.str.contains("gtex")].copy()
mode = gtex.groupby("tissue").imputation_mode.first()
print("GTEx polyA+ track per tissue:", mode.to_dict())
gtex["tissue"] = [t if mode[t] == "observed" else t + " *" for t in gtex.tissue]
conditions = {(t if mode[t] == "observed" else t + " *"): c for t, c in conditions.items()}
mat = gtex.pivot_table(index="gene_name", columns="tissue", values="raw_score", aggfunc="mean")
mat = mat[[t for t in conditions if t in mat.columns]]
coding = set(df.loc[df.gene_type == "protein_coding", "gene_name"])
mat = mat[mat.index.isin(coding)]
mat = mat.loc[mat.max(axis=1).sort_values(ascending=False).index]
print(f"\nPredicted mean RNA-seq signal per gene ({len(mat)} protein-coding genes):")
print(mat.round(3).to_string())

print("\nTissue with the highest predicted expression, per gene:")
best = pd.DataFrame({"top_tissue": mat.idxmax(axis=1), "top": mat.max(axis=1), "median": mat.median(axis=1)})
best["fold_over_median"] = best.top / best["median"]
print(best.round(3).sort_values("fold_over_median", ascending=False).to_string())

# Figure: row-scaled heatmap (log2 relative to each gene's median tissue).
rel = np.log2((mat + 1e-3).div((mat + 1e-3).median(axis=1), axis=0))
fig, ax = plt.subplots(figsize=(8, 0.33 * len(rel) + 2))
vmax = float(np.abs(rel.values).max())
im = ax.imshow(rel.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
ax.set_xticks(range(rel.shape[1]), rel.columns, rotation=45, ha="right")
ax.set_yticks(range(rel.shape[0]), rel.index, fontsize=8)
fig.colorbar(im, ax=ax, shrink=0.6, label="log2(signal / gene median across tissues)")
ax.set_title(f"score_interval RNA_SEQ, {interval.chromosome}:{interval.start}-{interval.end}\n(AlphaGenome v2, GTEx polyA+ track; * = imputed)")
out_png = OUT_DIR / "v2_interval_scores.png"
fig.savefig(out_png, dpi=150, bbox_inches="tight")
print(f"\nSaved {out_png}")
