---
name: kinome-crispr-library
description: Build and QC a pooled CRISPR-Cas9 knockout library for a set of genes (e.g. the kinome). Subsets a published library (Brunello/TKOv3) to your gene list, reports coverage with HGNC alias recovery and on-target stats, identifies gaps, and designs de-novo SpCas9 guides for gap genes from Ensembl canonical CDS. Use when someone asks to design, assemble, or check coverage of a targeted CRISPR knockout library or guide set.
---

# Kinome / targeted CRISPR-Cas9 knockout library

Two steps that were validated end-to-end on the human kinome (512 protein kinases
→ 98.8% covered by Brunello, 6 gaps filled de-novo). Reusable for any gene set.

## Step 1 — Subset a published library to your genes

Prefer starting from a validated pre-made library, not designing from scratch.

```bash
# gene list (one symbol per line), or --pkinfam to auto-fetch the human kinome
python3 scripts/subset_library.py \
  --genes my_genes.txt \
  --library broadgpp-brunello-library-contents.txt \
  --guides-per-gene 4 \
  --resolve-aliases \
  --outdir out/
```

- **Library file** (referenced by name, not bundled — it is large): the Broad GPP
  **Brunello** contents TSV (`broadgpp-brunello-library-contents.txt`, ~9 MB,
  from Addgene/Broad GPP). Any library with gene/spacer/score columns works;
  column names auto-detect, or pass `--gene-col/--spacer-col/--score-col`.
- `--pkinfam` fetches the authoritative human kinome from UniProt pkinfam instead
  of a `--genes` file.
- `--resolve-aliases` catches gene-symbol drift (e.g. `GRK2`→`ADRBK1`,
  `MAP3K20`→`ZAK`) via HGNC before declaring a gene uncovered.

**Outputs:** `out/selected_guides.tsv` (gene, library_symbol, spacer, score, GC,
source) and `out/coverage_stats.json` (covered / direct / via_alias / gaps /
median on-target). The printed "gaps" list is the input to Step 2.

## Step 2 — Fill gaps with de-novo guides

For genes with no library guide, design fresh SpCas9 (NGG) spacers from the
Ensembl canonical CDS, early-exon window, with the standard hard filters
(no poly-T, GC 30–70%, no homopolymer≥5, no BsmBI/BbsI site; cut sites ≥8 nt apart).

```bash
# put the gap symbols (from Step 1) one per line into gaps.txt
python3 scripts/design_guides_from_cds.py \
  --genes gaps.txt --out out/gap_guides.tsv \
  --guides-per-gene 6 --window 0.05 0.65
```

**Output:** `out/gap_guides.tsv` with real spacers. `rs2_score` and `offtarget`
are left `NA`/`pending` — score them with the **crispr-guide-scoring** skill
before ordering.

## Notes

- Exclude pseudogenes from the target set (they have no functional product).
- Concatenate Step 1 + Step 2 guides, add ~250 non-targeting + controls, then run
  scoring. See the **crispr-guide-scoring** and **crispr-screen-scale** skills for
  the next stages.
- Both scripts need only Python 3 (stdlib) + network (UniProt/Ensembl/HGNC REST).
