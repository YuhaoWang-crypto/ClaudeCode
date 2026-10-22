# Reproduction — paper vs recomputed

Generated from `gismo/data/benchmark.json` (`python3 -m gismo benchmark`).
GISMO built with 2,000 simulations (paper: 10,000); reference data: gnomAD v2.1.1
(LOEUF/MOEUF/CDS), Collins 2022 pHaplo/pTriplo, AlphaMissense gene-level,
gnomAD benchmark gene lists, OMIM AR (1,183) / AD (709) lists from macarthur-lab/gene_lists.

## Descriptive and correlations

| metric | paper | recomputed | n |
|---|---|---|---|
| `rho_gismo_loeuf` | 0.38 | 0.385 | 17265 |
| `rho_mis_loeuf` | 0.56 | 0.545 | 15910 |
| `rho_gismo_moeuf` | 0.33 | 0.32 | 17265 |
| `rho_mis_moeuf` | 0.61 | 0.602 | 15910 |
| `rho_gismo_phaplo` | -0.32 | -0.331 | 16456 |
| `rho_mis_phaplo` | -0.44 | -0.445 | 15274 |
| `rho_gismo_ptriplo` | -0.32 | -0.324 | 16456 |
| `rho_mis_ptriplo` | -0.5 | -0.502 | 15274 |
| `rho_mis_alphamissense` | 0.84 | -0.849 | 11287 |
| `rho_gismo_cds` | -0.12 | -0.117 | 17265 |
| `rho_mis_cds` | -0.12 | -0.111 | 15910 |
| `rho_loeuf_cds` | -0.5 | -0.504 | 17382 |
| `median_cds_top15_gismo` | 1502 | 1509.0 | 2590 |
| `median_cds_top15_mis` | 1302 | 1284.0 | 2387 |
| `median_cds_loeuf_lt_0.35` | 2267 | 2262.0 | 2853 |
| `median_losses_per_gene` | 18 | 18.0 |  |
| `median_losses_per_species` | 1130 | 1130.5 |  |
| `median_total_missense_per_species` | 610603 | 610602.917 |  |
| `median_total_synonymous_per_species` | 1084983 | 1084982.583 |  |
| `median_mean_missense` | 27.7 | 26.8 |  |
| `median_mean_synonymous` | 53.4 | 53.4 |  |
| `n_zero_loss_genes` | 6 | 6 |  |

AlphaMissense: paper reports +0.84; on native scales the correlation is negative
(low GISMO-mis ↔ high pathogenicity), magnitude matches.

## Enrichment odds ratios (definition not given in the paper)

| set / metric | paper | d1 vs rest | d1 vs d10 | d1–2 vs d9–10 | logistic (d1 vs d10 scale) |
|---|---|---|---|---|---|
| mouse_het_lethal_gismo | 3.43 | 1.35 | 8.52 | 3.03 | 2.94 |
| mouse_het_lethal_mis | 5.78 | 2.59 | 9.38 | 6.25 | 7.15 |
| mouse_lethal_gismo | 8.21 | 2.04 | 17.06 | 7.34 | 6.36 |
| mouse_lethal_mis | 13.1 | 2.53 | 10.93 | 6.48 | 6.99 |
| cell_essential_gismo | 4.34 | 1.87 | 117.59 | 55.7 | 7.04 |
| cell_essential_mis | 10.7 | 5.16 | 51.23 | 20.84 | 16.54 |
| recessive_gismo | 3.37 | 1.38 | 14.79 | 5.01 | 3.46 |
| recessive_mis | 0.91 | 0.29 | 0.56 | 0.61 | 0.71 |
| dominant_gismo | 9.92 | 1.7 | 4.96 | 3.05 | 2.6 |
| dominant_mis | 4.38 | 1.72 | 4.4 | 3.92 | 3.61 |

Direction and ordering reproduce (GISMO-mis > GISMO for essential/lethal sets;
GISMO enriched and GISMO-mis depleted for recessive genes), magnitudes depend on the
unstated definition. The figure-level quantity (share of a gene set per decile) matches:

| Fig. 4 series | max abs. difference (percentage points) |
|---|---|
| recessive_gismo | 1.3 |
| recessive_mis | 0.56 |
| dominant_gismo | 1.65 |
| dominant_mis | 0.5 |

(paper values read off the bar chart, ±1 pp)

## Share of each gene set per decile (recomputed; decile 1 = most constrained)

| set | metric | d1 | d2 | d3 | d4 | d5 | d6 | d7 | d8 | d9 | d10 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| mouse_het_lethal | GISMO | 13.0 | 14.8 | 13.0 | 12.5 | 10.7 | 10.7 | 7.8 | 8.3 | 7.8 | 1.6 |
| mouse_het_lethal | GISMO-mis | 21.8 | 15.6 | 14.3 | 12.4 | 10.0 | 8.9 | 7.3 | 3.5 | 3.8 | 2.4 |
| mouse_lethal | GISMO | 16.0 | 15.2 | 14.2 | 12.1 | 11.5 | 10.1 | 8.1 | 7.0 | 4.5 | 1.3 |
| mouse_lethal | GISMO-mis | 17.7 | 15.3 | 13.3 | 12.0 | 11.6 | 9.2 | 7.5 | 6.0 | 4.8 | 2.5 |
| cell_essential | GISMO | 16.8 | 15.1 | 15.3 | 14.4 | 12.7 | 10.7 | 9.8 | 4.5 | 0.5 | 0.1 |
| cell_essential | GISMO-mis | 33.9 | 15.0 | 10.6 | 9.8 | 8.8 | 8.3 | 5.9 | 5.0 | 1.8 | 0.8 |
| nonessential | GISMO | 2.4 | 2.6 | 2.9 | 4.3 | 5.1 | 6.7 | 9.0 | 14.6 | 21.7 | 30.8 |
| nonessential | GISMO-mis | 4.1 | 3.6 | 3.8 | 5.2 | 5.8 | 8.7 | 13.1 | 10.2 | 15.0 | 30.6 |
| olfactory | GISMO | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.3 | 0.8 | 3.3 | 34.9 | 60.7 |
| olfactory | GISMO-mis | 0.0 | 0.0 | 0.0 | 0.9 | 5.4 | 20.5 | 39.3 | 30.4 | 3.6 | 0.0 |
| recessive | GISMO | 13.1 | 13.0 | 14.4 | 13.0 | 12.0 | 10.7 | 9.7 | 8.6 | 4.6 | 0.9 |
| recessive | GISMO-mis | 3.3 | 5.9 | 10.5 | 12.5 | 12.6 | 15.1 | 11.8 | 13.5 | 9.1 | 5.8 |
| dominant | GISMO | 15.6 | 14.3 | 9.7 | 12.1 | 9.1 | 9.1 | 9.3 | 10.7 | 6.8 | 3.3 |
| dominant | GISMO-mis | 15.7 | 15.1 | 12.9 | 12.3 | 11.1 | 8.2 | 8.2 | 8.4 | 4.5 | 3.7 |

## Short genes

- GISMO-intolerant (top 15%) & LOEUF-tolerant (>0.35): median CDS 1278 bp (n=1873)
- GISMO-intolerant & LOEUF-intolerant: median CDS 2236 bp (n=714)

## Recessive rule (GISMO d1–2; both_constrained = LOEUF<0.35, gismo_only = LOEUF≥0.6)

| class | n | % OMIM AR | % OMIM AD |
|---|---|---|---|
| both_constrained | 979 | 2.66 | 9.4 |
| gismo_only | 1914 | 11.86 | 4.02 |
| other | 14372 | 6.25 | 3.69 |

## X chromosome

LOEUF<0.35 but GISMO decile 9–10: n=124, 16.94% X-linked vs 3.82% of all genes (Fisher OR 5.27, p=9.1e-09).

## Heritability (authors' shipped S-LDSC outputs, mean enrichment per decile)

- GISMO (650 trait IDs in archive; paper used 276 independent + 5 disease traits): d1 1.228, d2 1.173, d3 1.146, d4 1.105, d5 1.101, d6 1.127, d7 1.042, d8 1.008, d9 1.019, d10 0.939
- GISMO-mis (650 trait IDs in archive; paper used 276 independent + 5 disease traits): d1 1.4, d2 1.315, d3 1.296, d4 1.268, d5 1.275, d6 1.154, d7 1.104, d8 1.101, d9 1.104, d10 0.943

## Not reproduced

- UK Biobank exome CNV deletion/duplication per kb (Fig. 2g,h): individual-level data.
- S-LDSC re-run (needs GWAS summary statistics + 1000G reference panels).
- GTEx tissue-breadth, sHet, phylogenetic GISMO, DDD trio analyses.
