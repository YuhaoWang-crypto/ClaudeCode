"""AlphaGenome v2 demo: single-cell and methylation outputs.

For one variant (chr22:36201698 A>C) and several transverse-colon cell types,
request three v2-only outputs in a single call per condition:

  SC_GENE_EXPRESSION  gene-level CPM per single-cell protocol (genes x tracks)
  SC_ATAC             base-resolution single-cell ATAC tracks
  METHYLATION         base-resolution CpG methylation probability tracks

Requires ALPHAGENOME_API_KEY.

    python alphagenome_demo/demo_v2_sc.py
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
OT = dna_client.OutputType

model = dna_client.create(os.environ["ALPHAGENOME_API_KEY"], timeout=60)

interval = genome.Interval(chromosome="chr22", start=35677410, end=36725986)
variant = genome.Variant(
    chromosome="chr22",
    position=36201698,
    reference_bases="A",
    alternate_bases="C",
)
var_idx = variant.start - interval.start  # 0-based index of the variant base
LOCAL = 250  # +/- bp around the variant for ATAC / methylation summaries

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

gene_rows, track_rows, profiles = [], [], {}
for condition in conditions:
    out = model.predict_variant(
        interval,
        variant,
        requested_outputs=[OT.SC_GENE_EXPRESSION, OT.SC_ATAC, OT.METHYLATION],
        condition=condition,
    )
    ref, alt = out.reference, out.alternate
    label = (
        ref.sc_atac.metadata["biosample_name"].iloc[0].replace("transverse colon ", "")
    )

    # --- single-cell gene expression: genes x protocol tracks, in CPM ---
    ge_ref, ge_alt = ref.sc_gene_expression, alt.sc_gene_expression
    if condition == TISSUE:
        print("SC_GENE_EXPRESSION:", ge_ref.values.shape, "(genes x tracks)")
        print(ge_ref.track_metadata.to_string())
        print(ge_ref.gene_metadata.head(5).to_string())
    for g in range(len(ge_ref.gene_metadata)):
        gm = ge_ref.gene_metadata.iloc[g]
        for t in range(len(ge_ref.track_metadata)):
            tm = ge_ref.track_metadata.iloc[t]
            gene_rows.append(
                {
                    "condition": condition,
                    "label": label,
                    "gene_id": gm.get("gene_id"),
                    "gene_name": gm.get("gene_name", gm.get("gene_id")),
                    "track": tm["name"],
                    "mode": tm.get("imputation_mode"),
                    "ref_cpm": float(ge_ref.values[g, t]),
                    "alt_cpm": float(ge_alt.values[g, t]),
                }
            )

    # --- base-resolution tracks: summarise +/- LOCAL bp around the variant ---
    sl = slice(var_idx - LOCAL, var_idx + LOCAL + 1)
    for kind, r, a in [
        ("sc_atac", ref.sc_atac, alt.sc_atac),
        ("methylation", ref.methylation, alt.methylation),
    ]:
        for t in range(len(r.metadata)):
            rv, av = r.values[sl, t], a.values[sl, t]
            track_rows.append(
                {
                    "condition": condition,
                    "label": label,
                    "output": kind,
                    "assay": r.metadata["Assay title"].iloc[t],
                    "track": r.metadata["name"].iloc[t],
                    "mode": r.metadata["imputation_mode"].iloc[t],
                    "ref_mean": float(rv.mean()),
                    "alt_mean": float(av.mean()),
                    "ref_at_variant": float(r.values[var_idx, t]),
                    "alt_at_variant": float(a.values[var_idx, t]),
                    "max_abs_delta": float(np.abs(av - rv).max()),
                }
            )
    profiles[label] = {
        "atac_ref": ref.sc_atac.values[sl].mean(axis=1),
        "atac_alt": alt.sc_atac.values[sl].mean(axis=1),
        "meth_ref": ref.methylation.values[sl].mean(axis=1),
        "meth_alt": alt.methylation.values[sl].mean(axis=1),
    }
    print(f"done: {condition}")

genes = pd.DataFrame(gene_rows)
tracks = pd.DataFrame(track_rows)
genes.to_csv(OUT_DIR / "v2_sc_gene_expression.csv", index=False)
tracks.to_csv(OUT_DIR / "v2_sc_atac_methylation.csv", index=False)

# Gene-level effect: average CPM over protocol tracks, then log2 fold change.
PSEUDO = 0.1
g = genes.groupby(["label", "gene_name"], sort=False)[["ref_cpm", "alt_cpm"]].mean()
g["log2fc"] = np.log2((g.alt_cpm + PSEUDO) / (g.ref_cpm + PSEUDO))
lfc = g["log2fc"].unstack("label")[list(profiles)]
top = lfc.abs().max(axis=1).sort_values(ascending=False).head(8).index
print("\nSC gene expression: genes with the largest |log2FC| (mean over protocols)")
print(lfc.loc[top].round(2).to_string())
print("\nREF CPM of those genes (mean over protocols)")
print(g["ref_cpm"].unstack("label")[list(profiles)].loc[top].round(1).to_string())

t = tracks.groupby(["output", "label"], sort=False)[
    ["ref_mean", "alt_mean", "ref_at_variant", "alt_at_variant", "max_abs_delta"]
].mean()
t["log2fc_mean"] = np.log2(t.alt_mean / t.ref_mean)
print(f"\nSC ATAC and methylation within +/-{LOCAL} bp (mean over tracks)")
print(t.round(4).to_string())

# Figure: gene heatmap + local ATAC / methylation profiles for two cell types.
fig = plt.figure(figsize=(14, 8))
gs = fig.add_gridspec(2, 2, width_ratios=[1.1, 1])
ax = fig.add_subplot(gs[:, 0])
vmax = float(np.abs(lfc.loc[top].values).max())
im = ax.imshow(lfc.loc[top].values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
ax.set_xticks(range(len(lfc.columns)), lfc.columns, rotation=45, ha="right")
ax.set_yticks(range(len(top)), top)
ax.set_title("SC gene expression: log2(ALT/REF)")
fig.colorbar(im, ax=ax, shrink=0.6)

x = np.arange(-LOCAL, LOCAL + 1)
show = ["transverse colon", "macrophage", "colonocyte"]
palette = dict(zip(show, ["dimgrey", "tab:red", "tab:blue"]))
for row, (key, title) in enumerate(
    [("atac", "SC ATAC (mean over tracks)"), ("meth", "CpG methylation probability")]
):
    ax = fig.add_subplot(gs[row, 1])
    for name in show:
        ax.plot(x, profiles[name][f"{key}_ref"], color=palette[name], label=f"{name} REF")
        ax.plot(x, profiles[name][f"{key}_alt"], color=palette[name], ls="--", label=f"{name} ALT")
    ax.axvline(0, color="orange", lw=0.8)
    ax.set_title(title)
    ax.set_xlabel("bp from variant")
    if row == 0:
        ax.legend(fontsize=7, ncol=3)
fig.suptitle(f"{variant.chromosome}:{variant.position} A>C, AlphaGenome v2 single-cell and methylation outputs")
fig.tight_layout()
out_png = OUT_DIR / "v2_sc_methylation.png"
fig.savefig(out_png, dpi=150, bbox_inches="tight")
print(f"\nSaved {out_png}")
