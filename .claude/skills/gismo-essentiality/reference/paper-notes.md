# Paper notes — Liao et al. 2026, Nature Genetics (Technical Report)

"The landscape of gene loss and missense variation across the mammalian tree
informs on gene essentiality." Liao C, Ye R, Ivankovic F, Fu JM, Walters R,
Lowther C, Valkanas E, Churchhouse C, Samocha KE, Lindblad-Toh K, Karlsson E,
Hiller M, Talkowski ME, Neale BM. doi:10.1038/s41588-026-02783-9.
Received 18 Jul 2024, accepted 9 Sep 2026, published 20 Oct 2026. CC-BY 4.0.
Code + data: Zenodo doi:10.5281/zenodo.20547722 (GitHub cliao5/GISMO v1.0.0).

## Motivation
- gnomAD LOEUF/pLI quantify human LoF depletion (>800k individuals) but are
  length-biased: short genes have too few possible LoF sites, so their
  constraint cannot be estimated. ~1,000 genes have zero observed LoF.
- Deep evolutionary conservation (orthologs across ~100 Myr) is a
  length-independent alternative. 462 species ≈ 10% of placental mammals.

## Data
- Orthologs: **TOGA** (Kirilenko et al. Science 2023) projecting human
  GENCODE v38 onto 462 mammal assemblies (Zoonomia + others). Per species the
  highest contig-N50 assembly is used.
- Gene loss = TOGA classifies **all** transcripts as lost → ORF disrupted
  (frameshift, stop, splice-site, large deletion; multiple disruptions needed).
  Stored as orthology class `one2zero`. Other classes: one2one, one2many,
  many2one, many2many.
- Codon alignments: TOGA multiple codon alignments (MACSE v2), one-to-one
  orthologs; a species is excluded for a gene if it has >4 predicted orthologs;
  a gene needs a 1:1 ortholog in ≥75% of species.
- Divergence times: TimeTree.

## Descriptive results (Fig. 1)
- Median gene has 18 loss events across species; median species has ~1,130
  likely losses (TOGA lost + uncertain).
- Only 6 genes are never lost: NAP1L1, ZC3H14, MEIS2, STRIP1, CCT5, ATP6V0D1
  (DNA replication, poly(A) binding, development).
- 12 genes lost in every species except human (human-specific; mostly
  pseudogenes / single-exon).
- Median per-gene mean missense 27.7, synonymous 53.4 (vs consensus); per
  species median 610,603 missense and 1,084,983 synonymous differences.

## Methods (formulas)
GISMO: p_i = Σ_j I_ij w_j / Σ_j w_j, w_j = d_j / D_max.
CI: Ĩ_ijk ~ Bernoulli(p̂_i) (code: p̂ per gene per order), k = 1..10,000; the
95th/97.5th percentile is the score (code: `upper.ci`, quantile 0.975).
Phylogenetic GISMO (tree-aware) released too (Supp. Table 11) but weaker
correlations, so main analyses use plain GISMO. Forward weighting (distant
losses count more) beat reverse weighting.

GISMO-mis: consensus codon per position (ties keep several codons). Gaps in
either sequence skipped. Ambiguous N: if any resolution is synonymous, call
synonymous (conservative). With several reference codons: missense only if
missense against all. GISMO-mis_i = mean_j m_ij / s_ij.

Coverage: GISMO scores 97.7% of GENCODE v38 genes (highest of all metrics).

## Main results
- Essentiality gene sets (decile 1 vs others, ORs as reported):
  mouse het-lethal 3.43 / 5.78; mouse homozygous KO lethal 8.21 / 13.1;
  cell-essential 4.34 / 10.70 (GISMO / GISMO-mis). Least constrained deciles
  hold nonessential and olfactory genes. Frequently lost genes are expressed in
  fewer GTEx tissues. No gene-family-size bias.
- Correlations (Spearman): LOEUF 0.38 / 0.56; MOEUF 0.33 / 0.61; pHaplo
  −0.32 / −0.44; pTriplo −0.32 / −0.50; AlphaMissense (GISMO-mis) 0.84;
  sHet weaker for GISMO (−0.28 vs −0.61) — consistent with GISMO capturing
  homozygous selection. Beat conservation scores (mean GERP-like evolutionary
  rate, phylogenetic substitution frequency).
- CNVs (UK Biobank exomes, n = 197,306): most constrained deciles have fewest
  deletions/duplications; GISMO tracks deletions better than pHaplo.
- Length: ρ(CDS) = −0.12 for both vs −0.50 LOEUF. Top-15% intolerant median
  CDS: GISMO 1,502 bp (2,897 genes), GISMO-mis 1,302 (2,502) vs LOEUF<0.35
  2,267 (2,968). GISMO-intolerant & LOEUF-tolerant genes enrich for sperm
  immotility, absent acrosome, actin cytoskeleton (fertility).
- X chromosome: LOEUF-constrained but GISMO-least-constrained genes are
  enriched on X; no Mendelian/essential enrichment among them → LOEUF likely
  overestimates constraint there (hemizygosity, sex-specific selection,
  short-read artefacts).
- Heritability (S-LDSC, 276 independent UKB traits + SCZ, BIP, ASD, ADHD, CAD;
  100 kb windows; baseline v2.2): linear enrichment across deciles, lowest in
  least constrained. Rare-variant (fecundity-reducing) traits — autism,
  schizophrenia, NDD — strongest in the most constrained deciles.
- Recessive genes (1,183 OMIM AR): enriched in most constrained GISMO deciles
  (OR 3.37), not GISMO-mis (OR 0.91). Dominant (709 OMIM AD): 9.92 / 4.38.
  DDD (14,000 trios): confident/strong recessive genes follow the Mendelian
  GISMO distribution; limited/moderate ones are LOEUF-constrained (likely het
  selection, not recessive). Some definitive DDD recessive genes are in the
  least constrained GISMO deciles (corneal, thoracic, skeletal) → human-
  essential but mammal-dispensable.

## Recommended usage (Table 1)
| Application | Metric | Usage | Threshold |
|---|---|---|---|
| Mendelian disease–gene prioritisation | both | decile or continuous | deciles 1–2 |
| Recessive disease–gene discovery | GISMO | decile | deciles 1–4 |
| Pathway-level analyses | both | continuous (ranked enrichment) | — |
| Cross-mammalian conservation | both | continuous | — |
| Lineage-specific evolution | GISMO | continuous | — |
| Duplicated genes | GISMO | continuous | — (scores where LOEUF/pLI fail from mapping) |

GISMO-mis is also informative for polygenic traits in broader windows (top
20–40%). Combine: "a gene is much more likely recessive if GISMO is low and
LOEUF is high".

## Limitations (authors')
- Mammalian conservation ≠ human-specific selection or recent innovation.
- GISMO mixes homozygous + heterozygous selection; pair with LOEUF/GISMO-mis.
- Species sampling: survivorship / accessibility bias.
- Non-phylogenetic: one ancestral loss shared by a clade counts many times
  (mitigated by TOGA's strict multi-disruption calls; phylo version released).
- Recent human-specific duplications may inherit the paralogue's score.
