# Data package — provenance and contents

Everything here is derived. The primary sources are 16 GB and are not vendored;
`fetch_data.sh` rebuilds them from the original repositories in one command.

```bash
./fetch_data.sh                      # ~16 GB downloaded, ~1.2 GB retained
cd ../.. && python -m virtualcell.data
# expect: 4 cell lines, 6642 shared genes, 2053 shared perturbations
```

## Primary sources

| Dataset | Cell lines | Repository | Accession | Retrieved |
|---|---|---|---|---|
| Replogle et al., *Cell* 185:2559 (2022), genome-scale Perturb-seq | K562, RPE1 | Figshare+ | [article 20029387](https://plus.figshare.com/articles/dataset/_Mapping_information-rich_genotype-phenotype_landscapes_with_genome-scale_Perturb-seq_Replogle_et_al_2022_processed_Perturb-seq_datasets/20029387) | 2026-08-12 |
| Nadig et al., *Nat Genet* 57:1228 (2025), matched essential-gene CRISPRi | Jurkat, HepG2 | NCBI GEO | [GSE264667](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264667) | 2026-08-12 |
| `cell-eval` (metric definitions ported, not redistributed) | — | PyPI | `cell-eval==0.8.2` | 2026-08-12 |

Licences follow the originals: Replogle via Figshare+ (CC BY 4.0), GEO data per
NCBI terms, `cell-eval` per its own licence. Nothing here relicenses them.

Files actually downloaded, with sizes as served:

| File | Size | Source id |
|---|---|---|
| `K562_essential_raw_bulk.h5ad` | 80 MB | figshare file 35773070 |
| `rpe1_raw_bulk.h5ad` | 95 MB | figshare file 35775581 |
| `K562_gwps_raw_bulk.h5ad` | 375 MB | figshare file 35774443 |
| `GSE264667_hepg2_raw_singlecell_01.h5ad` | 5.6 GB | GEO supplementary |
| `GSE264667_jurkat_raw_singlecell_01.h5ad` | 9.4 GB | GEO supplementary |

## Derived contents

### `tables/perturbations.csv`

One row per knockdown × cell line, over the 2,053 knockdowns shared by all four
lines.

| Column | Meaning |
|---|---|
| `cell_line`, `knockdown` | context and targeted gene |
| `n_cells` | cells behind that pseudobulk |
| `effect_l2` | L2 norm of the effect, log1p CP10K space |
| `relative_effect` | `effect_l2` over that line's median, so lines are comparable |
| `n_de_genes` | FDR < 0.05 against the control-replicate null |
| `on_target_residual` | fraction of the target transcript remaining after CRISPRi; blank where the gene is not measured or is near-silent |
| `mean_cross_line_r` | mean pairwise correlation of that knockdown's effect across the four lines — how transferable it is |

### `tables/genes.csv`

One row per gene (6,642), with baseline expression and response spread in each
line, plus `shared_response_fraction`: the fraction of that gene's response
variance surviving averaging across lines. Low values mark genes whose response
is context-specific.

### `tables/hyperparameters.json`

The switch settings cross-validation chose per fold, on source lines only.

### `model_artifacts.npz`

Per held-out line (`K562/`, `RPE1/`, `HepG2/`, `Jurkat/`):

| Key | Shape | Meaning |
|---|---|---|
| `source_names` | (3,) | which lines were trained on |
| `context_weights` | (3,) | weights derived from control profiles alone |
| `perturbations` | (n,) | knockdowns in the source bank |
| `reliability` | (n,) | signal fraction, from cross-line reproducibility |
| `knockdown_residual` | (n,) | per-gene on-target efficiency; NaN where unmeasured |
| `global_effect` | (6642,) | generic across-perturbation response |
| `program_basis` | (50, 6642) | leading response programs |

Plus `genes`, `symbols`, `shared_perturbations`.

The full source consensus (≈2,400 × 6,642 per fold) is **not** included — it is
63 MB per fold and is exactly reproducible from the inputs. Rebuild it with
`SourceBank.build(sources, weights)`.

### `results/`

Raw benchmark output: per-fold metrics for every model, both regimes, plus the
context ablation and the stratified analysis. `*.txt` are the rendered tables.

## Regenerating

```bash
python -m virtualcell.run --regime context     # leave-one-cell-line-out
python -m virtualcell.run --regime double      # perturbation unseen too
python -m virtualcell.run --ablation           # score vs number of contexts
python -m virtualcell.analysis                 # stratified by effect strength
python -m virtualcell.export                   # regenerate this package
python -m virtualcell.figures
```

## Caveats that travel with the data

- **Depth.** 45–120 cells and ~11–14k UMI per perturbation, against the
  challenge's ~1,000 cells and >50k UMI on 10x Flex. Absolute scores here are
  not comparable to leaderboard numbers.
- **Pseudobulk.** `n_de_genes` comes from a control-replicate empirical null,
  not from a single-cell test, because there are no single cells at this level.
- **Essential genes.** All four screens target essential genes, so the
  perturbation panel is not representative of the transcriptome at large. This
  is not an abstract caveat: the panel shares **0 of the 300** Virtual Cell
  Challenge 2026 targets, and the genome-wide K562 arm — the same
  `K562_gwps_raw_bulk.h5ad` listed above, 9,866 knockdowns — had to be
  substituted to reach 272 of them. See `../VCC2026.md`.
- **One lab, one design.** The shared library and pipeline are what make these
  four lines comparable; they also mean shared systematic error is invisible
  here in a way it would not be across independent studies.
