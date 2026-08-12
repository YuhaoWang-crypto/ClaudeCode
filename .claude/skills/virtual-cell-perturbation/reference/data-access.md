# Getting four matched CRISPRi contexts

Context generalisation cannot be tested on one cell line. These four share a
CRISPRi library design, an essential-gene target set, and the Weissman-lab
processing conventions, so a knockdown means the same thing in all of them —
which removes the batch confound that makes cross-study perturbation
comparisons unreliable.

After harmonisation: **6,642 genes and 2,053 knockdowns shared by all four.**

| Line | Origin | Perturbations | Control replicates | Median cells/perturbation |
|---|---|---|---|---|
| K562 | chronic myeloid leukaemia | 2,057 | 109 | 110 |
| RPE1 | retinal pigment epithelium | 2,393 | 130 | 69 |
| HepG2 | hepatocellular carcinoma | 2,393 | 56 | 44 |
| Jurkat | T-cell leukaemia | 2,393 | 55 | 81 |

## K562 and RPE1 — Replogle et al., Cell 2022

Figshare+ article **20029387**. Enumerate files through the API rather than
scraping the page:

```bash
curl -s "https://api.figshare.com/v2/articles/20029387/files?page_size=100"
```

Take the **pseudobulk** files, not the single-cell ones. The single-cell K562
genome-wide matrix alone is 66 GB; the pseudobulk equivalent is 375 MB and is
all this task needs.

| File | Size | ndownloader id |
|---|---|---|
| `K562_essential_raw_bulk.h5ad` | 80 MB | 35773070 |
| `rpe1_raw_bulk.h5ad` | 95 MB | 35775581 |
| `K562_gwps_raw_bulk.h5ad` | 375 MB | 35774443 |

```bash
curl -sSL --retry 4 -o K562_essential_raw_bulk.h5ad \
  https://ndownloader.figshare.com/files/35773070
```

**Schema.** `X` is mean UMI per cell (row sums ~14k, *not* normalised to a
constant). `obs` is indexed by `gene_transcript` — `10023_ZC3H18_P1P2_ENSG…` —
so the target symbol is `name.split("_")[1]` and controls are named
`non-targeting`. `var` is indexed by ENSG with a `gene_name` column.

Three traps:

- `core_control` flags only the designated control set. Other rows are also
  non-targeting and are usable controls; catch them by target name.
- `num_cells_filtered` is **NaN** for exactly those extra control rows. Fall
  back to `num_cells_unfiltered`, or the variance model silently produces NaN
  and every DE test returns nothing.
- A gene targeted by two guide sets appears as two rows (`P1P2`, `P1`). Pool
  them by cell count before cross-line comparison.

## Jurkat and HepG2 — Nadig et al., Nat Genet 2025

GEO **GSE264667**. These are the matched essential-gene CRISPRi screens that
turn a two-context problem into a four-context one.

| File | Size |
|---|---|
| `GSE264667_hepg2_raw_singlecell_01.h5ad` | 5.6 GB |
| `GSE264667_jurkat_raw_singlecell_01.h5ad` | 9.4 GB |

```bash
curl -sSL --retry 4 -o GSE264667_hepg2_raw_singlecell_01.h5ad \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl/GSE264667_hepg2_raw_singlecell_01.h5ad
```

Only single-cell versions are published, so they must be pseudobulked. `X` is a
**dense** float32 cells × genes matrix of raw counts — 145k × 9,624 for HepG2 —
which is why the file is that large and why it must be streamed in row chunks
rather than loaded. `obs` categoricals are stored as codes plus a
`__categories` group and need rebuilding by hand; the useful columns are `gene`
(target symbol, `non-targeting` for controls) and `gem_group`.

`virtualcell/prep_nadig.py` does this in ~50 s per file:

```bash
python -m virtualcell.prep_nadig GSE264667_hepg2_raw_singlecell_01.h5ad hepg2.npz
```

It accumulates one row per targeted gene and one row per gem group over
non-targeting cells only — the control replicates are what later give an
empirical null for differential expression.

**Disk strategy.** 15 GB of downloads against a typical ~30 GB allowance, with
the repository and outputs also competing. Download one file, pseudobulk it to
~37 MB, delete the raw file, then move to the next. Do not fetch both first.

## Normalisation

Every pseudobulk row is scaled to counts-per-10k and `log1p`-ed — the space the
challenge metrics operate in. Do this *after* pooling guide-set rows, not
before, or the pooling is weighted wrongly.

## Datasets deliberately not used

- **scBaseCount** (Arc Virtual Cell Atlas) — observational only, no
  perturbations, so it cannot supply transfer signal.
- **Tahoe-100M** — chemical rather than genetic perturbations, but it is the
  only public resource with enough distinct cell lines to fit a real context
  encoder. The obvious next step for anyone with GPU budget.
- **scPerturb** — harmonises across labs and technologies, which reintroduces
  the batch confound the four lines above avoid.
- **VCC 2025 H1 hESC** — behind challenge registration. Challenge 2 explicitly
  permits its reuse, so a real entry should add it as a fifth context.
