---
name: gismo-essentiality
description: >-
  Gene essentiality from deep mammalian evolution: GISMO (gene-loss frequency
  across 462 placental mammals, weighted by divergence from human) and
  GISMO-mis (fixed missense/synonymous ratio vs a mammalian consensus), from
  Liao et al., Nature Genetics 2026 (doi:10.1038/s41588-026-02783-9). Use when
  asked how constrained / essential / dispensable a gene is; to prioritise
  Mendelian, recessive or neurodevelopmental candidate genes; to rescue short
  genes that LOEUF cannot score; to flag X-linked LOEUF artefacts; to compare
  paralogues; or to annotate pathway / network nodes with evolutionary
  constraint. Provides an offline score table for 19,686 genes, a faithful
  re-implementation, and a reproduction of the paper's numbers.
---

# GISMO / GISMO-mis — mammalian gene-loss & missense constraint

## What the two scores are (one line each)

- **GISMO** = how often a gene's ORF is lost (TOGA `one2zero`) across 462
  placental mammals, each loss weighted by that species' divergence time from
  human (w = d/D_max). Reported value = upper 97.5% bound of a Bernoulli
  simulation with p̂ estimated per gene *per mammalian order*. **Lower = more
  constrained.** Captures heterozygous **and homozygous** selection.
- **GISMO-mis** = mean over species of (missense / synonymous) codon
  differences against a **mammalian consensus** sequence (not human).
  **Lower = more constrained.** A dN/dS-like, codon-level constraint.
- Deciles: **1 = most constrained**. Use scores continuously when possible.

Both are nearly independent of CDS length (ρ ≈ −0.12 vs −0.50 for LOEUF): that
is the core advantage over gnomAD LOEUF/pLI.

## Quick use

```bash
python3 -m gismo lookup KRAS CFTR MECP2 OR1A1     # offline; adds LOEUF if fetched
python3 -m gismo fetch                             # Zenodo inputs + gnomAD/gene sets -> data/gismo/
python3 -m gismo build [N_ITER]                    # recompute scores (~4 min at 2000 sims)
python3 -m gismo benchmark                         # paper-vs-recomputed table + figure
python3 -m grn_pipeline.m23_gismo_paralogs         # essential paralogue in each pathway fiber
python3 -m grn_pipeline.m24_gismo_atlas            # annotate all M14 pathway genes -> gismo/data/m14_pathway_gismo.tsv
```

```python
from gismo.core import load_scores          # 19,686 genes: symbol, gene_id, chrom,
d = load_scores()                            # gismo, gismo_decile, n_loss, gismo_mis,
                                             # gismo_mis_decile, mean_mis, mean_syn,
                                             # QC: gc_cds, gc_cds_pct, gismo_hq(_decile),
                                             #     loss_n50_shift, n_loss_hq
from gismo.interpret import lookup, classify # rule-based tags (below)
```

Committed table: `gismo/data/gismo_scores.tsv.gz` (derived from CC-BY 4.0
data — cite the paper). Raw inputs are never committed; they cache in
`$GISMO_DATA` (default `data/gismo/`, git-ignored).

## How to interpret (paper Table 1 + rules validated here)

| Situation | Read it as | Evidence in this repo |
|---|---|---|
| GISMO or GISMO-mis decile 1–2 | deep-time essential; Mendelian / lethal-KO candidate | 34% of CRISPR core-essential genes sit in GISMO-mis d1 |
| GISMO d1–4 **and** LOEUF ≥ 0.6 | **recessive-like** (kept across mammals, heterozygous LoF tolerated in humans) | 13.0% OMIM-AR vs 5.9% background; AD unchanged (3.6%) |
| LOEUF high but CDS < 1 kb and GISMO d1–2 | LOEUF is underpowered (short gene); trust GISMO | e.g. KRAS: LOEUF 1.24, GISMO d2, mouse-lethal |
| LOEUF < 0.35 but GISMO d9–10, **X-linked** | LOEUF likely inflated (hemizygosity / mapping) | X is 16.9% of these vs 3.8% of all genes (OR 5.3) |
| GISMO d9–10 and GISMO-mis d≥8 | lineage-dispensable / fast-evolving (olfactory, immune) | 61% of olfactory receptors in GISMO d10 |
| GISMO d8–10 **and** CDS GC in top 20% (`gc_cds_pct` ≥ 0.8) | "loss" may be GC-rich assembly dropout — don't call dispensable | such genes are mouse-lethal 16.0% vs 5.5% for other d8–10 genes; pHaplo>0.86 12.0% vs 3.2% |
| GISMO decile improves ≥3 using only high-quality assemblies (`gismo_hq_decile`) | losses sit in fragmented genomes; treat GISMO as uncertain | 2.7% of genes |
| Recessive gene discovery | use **GISMO**, not GISMO-mis | AR genes: GISMO enriched at d1–4, GISMO-mis *depleted* (3% in d1) |
| Missense / protein-level questions | use **GISMO-mis** (tracks AlphaMissense, |ρ| = 0.85) | — |
| Picking the essential member of a paralogue family | **GISMO-mis or LOEUF, not GISMO** | M23: LOEUF 24/29, GISMO-mis 20/27, GISMO 16/29 ≈ random |

Tags come from `gismo.interpret.classify` and are ⚠️ heuristics — the paper
asks users to calibrate thresholds against known gene sets for each use.

## Rigor ledger (what is reproduced vs not)

- ✅ **GISMO-mis recomputed exactly** (max |Δ| = 0 vs authors' table; deciles
  99.994% identical).
- ✅ **GISMO**: authors' final scores are not in the Zenodo archive (they are
  Supplementary Table 2). Validated by reproducing reported numbers: 6
  never-lost genes (exact list), median 18 losses/gene, all 10 Spearman
  correlations within ±0.015, short-gene medians within 20 bp, Fig. 4c,d
  decile distributions within 1.7 percentage points.
- ⚠️ **Enrichment odds ratios** (e.g. mouse-lethal OR 8.21): the Methods do not
  define them. No single definition (d1 vs rest, d1 vs d10, d1–2 vs d9–10,
  logistic per decile) matches all reported ORs; the decile *distributions*
  (what the figures show) do match. Report ORs with their definition.
- ⚠️ Not reproducible from public data: UK Biobank CNV rates (Fig. 2g,h). LDSC
  heritability is summarised from the authors' shipped outputs, not re-run.
- Full table: `reference/reproduction.md`. Paper content in depth:
  `reference/paper-notes.md`.

## Gotchas found while re-implementing (real traps)

- **GISMO is GC-biased (not reported in the paper).** ρ(CDS GC, GISMO) = +0.19
  genome-wide and +0.37 among mouse-lethal genes; lethal genes in the top GC
  quintile land in GISMO d9–10 19.8% of the time vs ~1% in low-GC quintiles.
  Restricting to high-N50 assemblies does not remove it (ρ +0.16). GISMO-mis is
  unaffected (ρ = 0.02). A global GC correction does not improve gene-set
  prediction, so flag per gene instead (rule above). Examples: SHH, IRS2,
  SMAD6, KEAP1. Gene losses in general concentrate in fragmented assemblies
  (losing species have ~half the average contig N50).
- Curated pathway gene lists are biased toward famous essential genes: compare
  against other mouse-lethal genes, not just the genome (M24: same GISMO
  decile, p = 0.28).

- The LDSC files `partitioned-heritability/inputs/decile/d*/gene-set-file`
  look like GISMO deciles but are **GISMO-mis** deciles.
- GISMO-mis in the paper/downstream = `mean.comb` (per-species ratio, then
  mean). The R script also writes `mis.cal` = mean(mis)/mean(syn); its deciles
  agree only 89% — do not use it as "GISMO-mis".
- Paper formula divides by Σw, code divides by N: same ranking, different scale.
  Absolute GISMO values from other sources may differ by this constant.
- The R script reads `gen.dist.withassembly.txt`, which lacks `lin.num`; the
  species table that works is `genestat_26-07-2023.txt` (same rows + order).
- "1,130 losses per species" is TOGA's `lost + uncertain_loss`; the binary
  `one2zero` count used inside GISMO has median 1,424.
- AlphaMissense correlation is reported as +0.84; with both on their native
  scales it is −0.85 (low GISMO-mis = high pathogenicity).
- Mouse gene lists are matched by upper-casing symbols (approximate orthology).
- Paralogues: TOGA can call a gene lost when reads map to its paralogue, and
  all members of a family are usually retained, so GISMO rarely separates them.
