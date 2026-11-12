"""AlphaGenome v2 demo: in silico saturation mutagenesis (score_ism_variants).

Mutates every base in an 80 bp window around chr22:36201698 to each of the
three alternative bases (240 variants) and scores each one in transverse
colon, to see which positions the model considers important and whether the
quick-start variant (A>C) is special among its neighbours.

Requires ALPHAGENOME_API_KEY.

    python alphagenome_demo/demo_v2_ism.py
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
from alphagenome.interpretation import ism
from alphagenome.models.v2 import dna_client, variant_scorers
from alphagenome.visualization import plot_components

OUT_DIR = pathlib.Path(__file__).parent
pd.set_option("display.width", 200)

model = dna_client.create(os.environ["ALPHAGENOME_API_KEY"], timeout=60)

VARIANT_POS = 36201698  # 1-based
CONDITION = "UBERON:0001157"  # transverse colon
GENE = "APOL4"
ISM_WIDTH = 80

center = genome.Interval("chr22", VARIANT_POS - 1, VARIANT_POS)
# 131 kb of sequence context is enough to contain all of APOL4 (16 kb) and is
# much cheaper than 1 Mb for 240 forward passes.
interval = center.resize(dna_client.SEQUENCE_LENGTH_100KB)
ism_interval = center.resize(ISM_WIDTH)

REC = variant_scorers.RECOMMENDED_VARIANT_SCORERS
scorers = {"RNA_SEQ": REC["RNA_SEQ"], "SPLICE_SITES": REC["SPLICE_SITES"], "DNASE": REC["DNASE"]}

t0 = time.time()
results = model.score_ism_variants(
    interval,
    ism_interval,
    CONDITION,
    variant_scorers=list(scorers.values()),
    progress_bar=False,
)
print(f"{len(results)} variants scored in {time.time() - t0:.0f} s")

# One scalar per variant per scorer.
rows = []
for per_variant in results:
    v = per_variant[0].uns["variant"]
    row = {"variant": str(v), "position": v.position, "ref": v.reference_bases, "alt": v.alternate_bases}
    for name, ad in zip(scorers, per_variant):
        if "gene_name" in ad.obs:  # gene-centric scorer: take the APOL4 row
            ad = ad[ad.obs["gene_name"] == GENE]
        x = np.asarray(ad.X, dtype=float)
        if name == "RNA_SEQ":  # gene LFC on the unstranded GTEx track
            track = ad.var["name"].str.contains("gtex").to_numpy()
            row[name] = float(x[:, track].mean())
        elif name == "SPLICE_SITES":  # largest change over donor/acceptor
            row[name] = float(x.max()) if x.size else 0.0
        else:  # mean log2FC over DNase tracks
            row[name] = float(x.mean())
    rows.append((v, row))
variants = [v for v, _ in rows]
df = pd.DataFrame([r for _, r in rows]).sort_values(["position", "alt"]).reset_index(drop=True)
df.to_csv(OUT_DIR / "v2_ism_scores.csv", index=False)

ref_seq = "".join(df.groupby("position")["ref"].first())
print(f"Reference {ism_interval} (+ strand): {ref_seq}")

for name in scorers:
    print(f"\n{name}: 8 strongest variants by |score|")
    top = df.reindex(df[name].abs().sort_values(ascending=False).index).head(8)
    print(top[["variant", name]].round(3).to_string(index=False))

focal = df[(df.position == VARIANT_POS)]
print(f"\nAll three substitutions at {VARIANT_POS}:")
print(focal[["variant"] + list(scorers)].round(3).to_string(index=False))
for name in scorers:
    rank = int(df[name].abs().rank(ascending=False)[focal.index[focal.alt == "C"][0]])
    print(f"  A>C rank by |{name}|: {rank} / {len(df)}")

# Per-position summary: mean effect of the three substitutions.
pos = df.groupby("position")[list(scorers)].mean()
print("\nPositions with the largest mean |RNA_SEQ| effect:")
print(pos.reindex(pos["RNA_SEQ"].abs().sort_values(ascending=False).index).head(10).round(3).to_string())

# Figure: sequence logos of ISM contribution scores (mean-centred per position).
components = []
for name, label in [
    ("RNA_SEQ", f"{GENE} RNA-seq LFC"),
    ("SPLICE_SITES", f"{GENE} splice sites"),
    ("DNASE", "DNase log2FC"),
]:
    by_variant = dict(zip(df["variant"], df[name]))
    matrix = ism.ism_matrix(
        [by_variant[str(v)] for v in variants], variants, multiply_by_sequence=True
    )
    components.append(
        plot_components.SeqLogo(scores=matrix, scores_interval=ism_interval, ylabel=label)
    )
plot_components.plot(
    components,
    interval=ism_interval,
    fig_width=16,
    annotations=[
        plot_components.VariantAnnotation(
            [genome.Variant("chr22", VARIANT_POS, "A", "C")], alpha=0.6
        )
    ],
)
out_png = OUT_DIR / "v2_ism.png"
plt.savefig(out_png, dpi=150, bbox_inches="tight")
print(f"\nSaved {out_png}")
