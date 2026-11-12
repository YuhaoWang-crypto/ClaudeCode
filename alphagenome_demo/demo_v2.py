"""AlphaGenome v2 API demo: one variant, bulk tissue vs. individual cell types.

v2 replaces `ontology_terms=[...]` with a single `condition=` and adds
composite "tissue_celltype" conditions (e.g. UBERON:0001157_CL:0011108 =
transverse colon colonocyte), most of which are imputed rather than observed.

Same variant as demo.py (chr22:36201698 A>C). Requires ALPHAGENOME_API_KEY.

    python alphagenome_demo/demo_v2.py
"""

import os
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client

OUT_DIR = pathlib.Path(__file__).parent

model = dna_client.create(os.environ["ALPHAGENOME_API_KEY"], timeout=60)

interval = genome.Interval(chromosome="chr22", start=35677410, end=36725986)
variant = genome.Variant(
    chromosome="chr22",
    position=36201698,
    reference_bases="A",
    alternate_bases="C",
)
zoom = interval.resize(2**15)  # 32 kb around the variant
lo, hi = zoom.start - interval.start, zoom.end - interval.start

# What v2 can be conditioned on, and how each output is backed per condition.
supported = model.supported_conditions()
print(f"{len(supported)} supported human conditions")
print("RNA-seq backing:", supported["rna_seq"].value_counts().to_dict())

TISSUE = "UBERON:0001157"  # transverse colon
CELL_TYPES = [
    "CL:0011108",  # colonocyte
    "CL:0019031",  # intestine goblet cell
    "CL:0009010",  # transit amplifying cell
    "CL:0000186",  # myofibroblast cell
    "CL:0000235",  # macrophage
    "CL:0000236",  # B cell
    "CL:0000624",  # CD4-positive, alpha-beta T cell
    "CL:0002144",  # capillary endothelial cell
]
conditions = [TISSUE] + [f"{TISSUE}_{cl}" for cl in CELL_TYPES]

rows = []
for condition in conditions:
    out = model.predict_variant(
        interval,
        variant,
        requested_outputs=[dna_client.OutputType.RNA_SEQ],
        condition=condition,
    )
    ref, alt = out.reference.rna_seq, out.alternate.rna_seq
    ref_sum = ref.values[lo:hi].sum(axis=0)
    alt_sum = alt.values[lo:hi].sum(axis=0)
    meta = ref.metadata
    for i in range(len(meta)):
        rows.append(
            {
                "condition": condition,
                "biosample": meta["biosample_name"].iloc[i],
                "assay": meta["Assay title"].iloc[i],
                "endedness": meta["endedness"].iloc[i],
                "strand": meta["strand"].iloc[i],
                "mode": meta["imputation_mode"].iloc[i],
                "ref": float(ref_sum[i]),
                "alt": float(alt_sum[i]),
                "log2fc": float(np.log2(alt_sum[i] / ref_sum[i])),
            }
        )
    print(f"done: {condition} ({len(meta)} tracks)")

df = pd.DataFrame(rows)
df.to_csv(OUT_DIR / "v2_celltype_effects.csv", index=False)

# The affected gene is on the minus strand; compare one track type across
# conditions (polyA+ RNA-seq, paired-ended, minus strand).
sel = df[
    (df.assay == "polyA plus RNA-seq")
    & (df.endedness == "paired")
    & (df.strand == "-")
].copy()
sel["label"] = sel["biosample"].str.replace("transverse colon ", "")
print(f"\npolyA+ RNA-seq (paired, - strand), summed within {zoom}:")
print(sel[["label", "mode", "ref", "alt", "log2fc"]].round(3).to_string(index=False))

fig, ax = plt.subplots(figsize=(8, 4))
colors = ["dimgrey" if c == TISSUE else "tab:red" for c in sel["condition"]]
ax.barh(sel["label"], sel["log2fc"], color=colors)
ax.invert_yaxis()
ax.axvline(0, color="black", lw=0.8)
ax.set_xlabel("log2(ALT / REF), polyA+ RNA-seq, minus strand, 32 kb window")
ax.set_title(f"{variant.chromosome}:{variant.position} A>C in transverse colon (AlphaGenome v2)")
out_png = OUT_DIR / "v2_celltype_effects.png"
fig.savefig(out_png, dpi=150, bbox_inches="tight")
print(f"\nSaved {out_png}")
