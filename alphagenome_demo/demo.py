"""AlphaGenome demo: predict the effect of a single variant on RNA-seq.

Uses the official quick-start example (chr22:36201698 A>C, colon tissue).
Requires the ALPHAGENOME_API_KEY environment variable.

    pip install alphagenome
    python alphagenome_demo/demo.py

Behind a TLS-intercepting proxy, also set GRPC_DEFAULT_SSL_ROOTS_FILE_PATH
to the proxy CA bundle (the client talks gRPC to gdmscience.googleapis.com).
"""

import os
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from alphagenome.data import genome
from alphagenome.models.v1 import dna_client
from alphagenome.visualization import plot_components

OUT_DIR = pathlib.Path(__file__).parent

model = dna_client.create(os.environ["ALPHAGENOME_API_KEY"], timeout=60)

# 1 Mb window on chr22 (hg38) and a single-nucleotide variant inside it.
interval = genome.Interval(chromosome="chr22", start=35677410, end=36725986)
variant = genome.Variant(
    chromosome="chr22",
    position=36201698,
    reference_bases="A",
    alternate_bases="C",
)
ONTOLOGY = "UBERON:0001157"  # transverse colon

outputs = model.predict_variant(
    interval=interval,
    variant=variant,
    ontology_terms=[ONTOLOGY],
    requested_outputs=[dna_client.OutputType.RNA_SEQ],
)

ref, alt = outputs.reference.rna_seq, outputs.alternate.rna_seq
print("Tracks returned:")
print(ref.metadata.to_string())
print("Prediction shape (bp x tracks):", ref.values.shape)

# Numeric summary in a 32 kb window centred on the variant.
zoom = ref.interval.resize(2**15)
lo = zoom.start - ref.interval.start
hi = zoom.end - ref.interval.start
ref_sum = ref.values[lo:hi].sum(axis=0)
alt_sum = alt.values[lo:hi].sum(axis=0)
print(f"\nSummed RNA-seq signal within {zoom}:")
for name, strand, r, a in zip(ref.metadata["name"], ref.metadata["strand"], ref_sum, alt_sum):
    print(f"  {name} ({strand}): REF={r:.1f}  ALT={a:.1f}  log2(ALT/REF)={np.log2(a / r):+.3f}")

# REF vs ALT overlay plot.
plot_components.plot(
    [
        plot_components.OverlaidTracks(
            tdata={"REF": ref, "ALT": alt},
            colors={"REF": "dimgrey", "ALT": "red"},
        ),
    ],
    interval=zoom,
    annotations=[plot_components.VariantAnnotation([variant], alpha=0.8)],
)
out_png = OUT_DIR / "variant_effect_rna_seq.png"
plt.savefig(out_png, dpi=150, bbox_inches="tight")
print(f"\nSaved plot to {out_png}")
