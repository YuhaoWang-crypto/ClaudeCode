# How to use this

Two capabilities came out of this work. This is what each takes in, what it
gives back, and a runnable demo for both.

---

# 1. The virtual cell model (`ContextTransfer`)

## What it does

Predicts what silencing a gene does to a cell line **you have no perturbation
data for**, using only that line's non-targeting control cells.

## Input / output contract

| | |
|---|---|
| **In** | control cells of your cell line — an `.h5ad`, an `AnnData`, or a mean-expression vector — plus a list of gene symbols to knock down |
| **Out** | a `Prediction`: post-perturbation profile, effect vs control, and a ranked differential-expression list per knockdown |
| **Not needed** | any perturbation measured in your cell line. That is the point. |
| **Gene matching** | by symbol (`gene_name` / `gene_symbol` / `symbol` in `var`) or by ENSG in `var_names`; raw counts are fine, normalisation is handled |

```python
from virtualcell.predict import VirtualCell

vc  = VirtualCell.from_atlas()                       # trains on 4 CRISPRi lines
out = vc.predict("my_controls.h5ad", ["TP53", "MYC", "RPL5"])

out.expression        # (3, 6642) predicted profile, log1p CP10K
out.effect            # (3, 6642) change from control
out.seen              # which knockdowns existed in training
out.top_genes("TP53", 20)
out.to_frame()        # long table: knockdown, gene, control, predicted, effect, log2fc
```

Pass `VirtualCell.from_atlas(exclude="K562")` to drop a line from training —
that is how to get an honest zero-shot estimate for a line the atlas contains.

## Demo 1 — predict

```bash
python demos/demo1_predict.py
```

Predicts five knockdowns in K562 with K562 excluded from training, so nothing
about it is memorised. Output:

```
knockdown seen?     |effect|   genes moved   on-target
RPL5      yes           7.04          1078       -1.65
SF3B1     yes           3.37           182       -0.30
PSMA1     yes           6.97          1079       -0.45
MYC       yes           5.23           624       -0.64
ACOT12    yes           0.71             1         nan
```

It separates a ribosome subunit that guts the cell (1,078 genes moved) from a
metabolic gene that does nothing (1 gene) — and the ranking for SF3B1, a
splicing factor, is led by SNHG1, GAS5, ZFAS1 and SNHG6, the snoRNA host genes
that splicing inhibition disrupts, plus the DDIT3/GADD45A stress response. That
is correct biology the model was never told.

`nan` on-target for ACOT12 is honest: that gene is not in the measured 6,642, so
there is nothing to report.

Writes `demos/out/demo1_predictions.csv.gz` — 33,210 rows, one per
knockdown × gene.

## Demo 2 — validate

```bash
python demos/demo2_validate.py --target K562 --n-knockdowns 300
python demos/demo2_validate.py --target RPE1        # check it is not one lucky line
```

Predicts a held-out line from its controls, then scores against what was
actually measured, with the Challenge's metrics and the baselines that matter:

```
model                               discrim   DE@100      MAE  direction   score  balanced
control (delta=0)                     0.503    0.004   0.0431      0.000   0.023    +0.006
global mean [challenge baseline]      0.503    0.053   0.0463      0.660   0.000    +0.000
naive transfer                        0.728    0.262   0.0536      0.839   0.225    +0.172
ContextTransfer (ours)                0.681    0.282   0.0438      0.851   0.218    +0.218
```

How to read it, honestly:

- **0.500 discrimination is chance.** Both transfer models are far above it; the
  two trivial baselines are at it, as they must be.
- **The comparison that counts is against `naive transfer`**, not against
  predicting nothing. It already uses cross-context information.
- Ours wins DE overlap, direction and MAE; naive wins discrimination.
- **Naive transfer is worse than baseline on MAE**; ours is not. The Challenge
  enforces minimum thresholds on every metric, so failing one matters.
- The two score columns disagree on purpose: `score` is the leaderboard's, which
  clips a sub-baseline metric to zero and so forgives naive's MAE failure;
  `balanced` does not clip. Report both.

## Limits worth knowing before trusting it

- Recovers the **transferable** part of a response — ~43% of effect variance
  across four lines. It does not predict the context-specific remainder.
- For a gene perturbed in **none** of the training lines it can say which genes
  respond (2.7× baseline DE overlap) but **cannot tell that knockdown from
  another** (discrimination stays at chance).
- CRISPRi knockdown only. Cross-modality transfer to CRISPRa fails outright.
- Trained on essential genes in four cancer/immortalised lines.

---

# 2. Reproducing Stack (Arc's foundation model)

## What is and is not reproducible

| | |
|---|---|
| ❌ Pre-training | 149M cells, H100 80GB for 2–3 days, 320 GB RAM |
| ❌ Alignment | H100, 400 GB RAM |
| ❌ Regenerating *Perturb Sapiens* | 1.34 TB |
| ✅ Inference from released weights | runs on **CPU**; 2.61 GB checkpoint |
| ✅ Verifying *Perturb Sapiens* | ~6 GB per perturbation file |

## Input / output contract

`stack-generation` predicts by analogy — give it a context population where the
conditions were measured and a query population of control cells:

```bash
python -m virtualcell.patch_stack -- stack-generation \
  --checkpoint bc_large_aligned.ckpt \
  --base-adata context.h5ad   `# perturbed cells, labelled by condition` \
  --test-adata  test.h5ad     `# control cells of the target population` \
  --genelist basecount_1000per_15000max.pkl \
  --gene-name-col gene_name --split-column gene \
  --num-steps 5 --output-dir out
```

Out: one `.h5ad` of predicted **cells** per condition, raw counts, on Stack's
15,012-gene universe. Build the inputs with
`python -m virtualcell.prep_stack_inputs --target k562`.

**The `patch_stack` wrapper is required.** Released `arc-stack` 0.1.3 crashes
partway through generation with `ValueError: Quantiles must be in the range
[0, 1]` — reproducible on the project's own tutorial. The unmasking schedule
computes `unmask_rate = 1 − mask_rate/f`, but the next two lines re-mask on the
model's own logits, so `f` is data-dependent and can fall below `mask_rate`,
making the rate negative. It is data-dependent, so it does not fire
deterministically. See `virtualcell/patch_stack.py`.

## Measured cost on 4 CPU cores

| workload | time |
|---|---|
| official demo, 15 conditions × 2,848 query cells | **96m51s** |
| this benchmark, 4 conditions × 500 query cells | **4m36s** (~69 s per condition) |

Cost tracks query-cell count strongly. Keep the query set small.

## Verified against the paper

| claim | field | result |
|---|---|---|
| 28 tissues | `tissue_in_publication` | **28** ✅ |
| 40 cell classes | `broad_cell_class` | **40** ✅ |
| 201 perturbations | file listing | 203 files (91 cytokine + 112 drug), one is `PBS` vehicle |

The tissue count needs the right column: the CellxGene ontology `tissue` field
has 75 values and the methods' ≥1000-cell filter gives 48. Checking against
either would wrongly suggest the claim fails.

**Caveat when applying Stack here.** `Stack-Large-Aligned` was post-trained on
chemical and cytokine perturbations. Asked for CRISPRi knockdowns it produces
perturbation-specific responses (pairwise effect correlation 0.26–0.44) but
**no on-target knockdown** — the silenced gene sits at ≈0 where CRISPRi gives
≈ −2. Genetic perturbation is within its stated in-context scope but is not the
modality it was aligned on, so a weak score here should not be read as a verdict
on the model.

---

## Getting the data

Both capabilities need the CRISPRi atlas:

```bash
cd virtualcell/datapackage && ./fetch_data.sh     # ~16 GB, keeps ~1.2 GB
python -m virtualcell.data                        # expect 4 lines, 6642 genes, 2053 knockdowns
```

For Stack, additionally `virtualcell/prep_singlecell.py` re-fetches the raw
single-cell files and subsamples them (34 GB → 1.5 GB), because Stack consumes
cells where this package works on pseudobulk.
