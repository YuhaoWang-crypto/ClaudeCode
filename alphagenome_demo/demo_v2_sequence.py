"""AlphaGenome v2 demo: predict_sequence (predict from a raw DNA string).

predict_interval / predict_variant look the sequence up in the reference
genome server-side. predict_sequence takes any ACGTN string of a supported
length, so it can score edits that are not a single variant: scrambled or
masked elements, designed sequences, reverse complements, random DNA.

Experiments on the 131 kb window around APOL4 (transverse colon, RNA-seq):
  1. reference string vs predict_interval        (should be identical)
  2. scramble the 200 bp gradient hotspot found by demo_v2_gradients.py
  3. scramble a 200 bp control region 5 kb away  (should do little)
  4. mask the hotspot with N
  5. reverse-complement the whole window         (strands should swap)
  6. random DNA with the same base composition   (should be silent)

The reference sequence is fetched from the UCSC REST API.
Requires ALPHAGENOME_API_KEY.

    python alphagenome_demo/demo_v2_sequence.py
"""

import json
import os
import pathlib
import urllib.request

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client

OUT_DIR = pathlib.Path(__file__).parent
pd.set_option("display.width", 200)
OT = dna_client.OutputType
CONDITION = "UBERON:0001157"  # transverse colon

# Same window predict_interval_gradients used for APOL4 (0-based, half-open).
interval = genome.Interval("chr22", 36139270, 36270342)
GENE = genome.Interval("chr22", 36189123, 36204840, "-")  # APOL4
HOTSPOT = (36201810, 36202010)  # strongest gradient hotspot
CONTROL = (36206810, 36207010)  # same size, 5 kb away, outside the gene

url = (
    "https://api.genome.ucsc.edu/getData/sequence?genome=hg38;"
    f"chrom={interval.chromosome};start={interval.start};end={interval.end}"
)
with urllib.request.urlopen(url) as resp:
    reference = json.load(resp)["dna"].upper()
assert len(reference) == interval.width == dna_client.SEQUENCE_LENGTH_100KB
print(f"Fetched {len(reference)} bp of hg38 {interval} from UCSC")

rng = np.random.default_rng(0)
COMP = str.maketrans("ACGTN", "TGCAN")


def scramble(seq, region):
    a, b = region[0] - interval.start, region[1] - interval.start
    return seq[:a] + "".join(rng.permutation(list(seq[a:b]))) + seq[b:]


def mask(seq, region):
    a, b = region[0] - interval.start, region[1] - interval.start
    return seq[:a] + "N" * (b - a) + seq[b:]


sequences = {
    "reference": reference,
    "hotspot scrambled": scramble(reference, HOTSPOT),
    "control scrambled": scramble(reference, CONTROL),
    "hotspot masked (N)": mask(reference, HOTSPOT),
    "reverse complement": reference.translate(COMP)[::-1],
    "random DNA": "".join(rng.permutation(list(reference))),
}

model = dna_client.create(os.environ["ALPHAGENOME_API_KEY"], timeout=60)
preds = {}
for name, seq in sequences.items():
    out = model.predict_sequence(
        seq, requested_outputs=[OT.RNA_SEQ], condition=CONDITION, interval=interval
    )
    preds[name] = out.rna_seq
    print(f"done: {name}")

meta = preds["reference"].metadata
minus = int(np.flatnonzero(meta.name.str.contains("total RNA-seq-encode-paired") & (meta.strand == "-"))[0])
plus = int(np.flatnonzero(meta.name.str.contains("total RNA-seq-encode-paired") & (meta.strand == "+"))[0])
gtex = int(np.flatnonzero(meta.name.str.contains("gtex"))[0])
print("Tracks used:", meta.name[minus], "(-/+),", meta.name[gtex], "| modes:",
      meta.imputation_mode[minus], meta.imputation_mode[gtex])

# 1. predict_sequence(reference) vs predict_interval.
by_interval = model.predict_interval(interval, requested_outputs=[OT.RNA_SEQ], condition=CONDITION).rna_seq
a, b = preds["reference"].values, by_interval.values
print("\n1. reference string vs predict_interval:")
print(f"   max |difference| = {np.abs(a - b).max():.2e} (max signal {b.max():.2f}), "
      f"Pearson r = {np.corrcoef(a.ravel(), b.ravel())[0, 1]:.6f}")

# 2-6. APOL4 signal: sum over the gene body.
g0, g1 = GENE.start - interval.start, GENE.end - interval.start
rows = []
for name, td in preds.items():
    v = td.values
    if name == "reverse complement":
        # Map back to reference coordinates: reverse positions, swap strands.
        v = v[::-1]
        m, u = v[g0:g1, plus].sum(), v[g0:g1, gtex].sum()
    else:
        m, u = v[g0:g1, minus].sum(), v[g0:g1, gtex].sum()
    rows.append({"sequence": name, "total RNA-seq (gene strand)": m, "GTEx polyA+ (unstranded)": u})
res = pd.DataFrame(rows).set_index("sequence")
ref_row = res.loc["reference"]
for c in list(res.columns):
    res[f"log2FC {c.split(' (')[0]}"] = np.log2(res[c] / ref_row[c])
print("\nPredicted APOL4 signal summed over the gene body:")
print(res.round(3).to_string())
res.to_csv(OUT_DIR / "v2_sequence_edits.csv")

# 5. Strand equivariance of the reverse complement, base by base.
rc = preds["reverse complement"].values[::-1]
r_same = np.corrcoef(a[:, minus], rc[:, plus])[0, 1]
r_wrong = np.corrcoef(a[:, minus], rc[:, minus])[0, 1]
print("\n5. reverse complement, per-base correlation with the reference minus-strand track:")
print(f"   RC plus-strand track (expected match): r = {r_same:.4f}")
print(f"   RC minus-strand track (expected mismatch): r = {r_wrong:.4f}")

# Figure: gene-strand RNA-seq around the gene for each sequence.
zoom = (GENE.start - 3000, GENE.end + 6000)
z0, z1 = zoom[0] - interval.start, zoom[1] - interval.start
x = np.arange(zoom[0], zoom[1])
fig, axes = plt.subplots(len(preds), 1, figsize=(13, 1.5 * len(preds)), sharex=True, sharey=True)
for ax, (name, td) in zip(axes, preds.items()):
    v = td.values
    y = v[::-1][z0:z1, plus] if name == "reverse complement" else v[z0:z1, minus]
    ax.fill_between(x, y, color="dimgrey" if name == "reference" else "tab:red", lw=0)
    ax.axvspan(*HOTSPOT, color="orange", alpha=0.5)
    ax.axvspan(*CONTROL, color="tab:blue", alpha=0.3)
    ax.set_ylabel(name, rotation=0, ha="right", va="center", fontsize=9)
axes[0].set_title("predict_sequence: APOL4-strand total RNA-seq, transverse colon "
                  "(orange = hotspot, blue = control region)")
axes[-1].set_xlabel("chr22 position (reference coordinates)")
fig.tight_layout()
out_png = OUT_DIR / "v2_sequence_edits.png"
fig.savefig(out_png, dpi=150, bbox_inches="tight")
print(f"\nSaved {out_png}")
