# Analysis Plan

## 1. From reads to a guide-count matrix

1. **Demultiplex** by sample index.
2. **Trim** to the constant scaffold anchor and extract the 20-nt spacer.
3. **Map** spacers to the library (exact or ≤1-mismatch) — `MAGeCK count`,
   `mageck-count`, or a custom exact matcher. Discard reads not matching any
   designed guide.
4. Produce a **guide × sample count matrix**. Record mapping rate (should be
   ≥ 65–75%) and per-sample total reads.

## 2. QC gates before calling hits

- **Mapping rate** ≥ 65%; **zero-count guides** low (most guides detected).
- **Gini index** of counts per sample (plasmid/T0 tight, < ~0.1; increasing
  Gini over time reflects real selection).
- **Replicate correlation** (Pearson/Spearman on log-counts or LFCs) ≥ ~0.8.
- **PCA / clustering** separates timepoints/conditions, not batches.
- **Positive controls drop out:** essential kinases (PLK1, CDK1, AURKB…) and
  Hart CEGv2 pan-essentials deplete; NTC and NEGv2 stay flat.
- **ROC-AUC / precision-recall** for known essential vs. non-essential
  benchmark sets — the headline "did the screen work" metric.

## 3. Statistical analysis (pick by screen type)

| Screen type | Method | Notes |
|-------------|--------|-------|
| **Dropout / essentiality** | **MAGeCK (RRA)**, **MAGeCK-MLE**, **BAGEL2** (Bayes Factor) | BAGEL2 gives calibrated essentiality via CEGv2/NEGv2 priors |
| **Positive/negative selection vs. control** | **MAGeCK test**, **drugZ** | drugZ is tuned for chemogenomic/synthetic-lethal (drug vs. DMSO) |
| **Multi-condition / time-course** | **MAGeCK-MLE**, **JACKS**, **CRISPhieRmix** | Model design matrix; share info across guides |

- **Normalize** to total reads or to non-targeting/median; MAGeCK's
  median-ratio normalization is a good default. Use NTCs to build the null.
- **Correct copy-number bias** for essentiality (CERES / **Chronos**), so
  amplified-region cut-toxicity isn't misread as dependency.
- Report **gene-level LFC, FDR, and rank**, plus per-guide consistency (do the
  6 guides agree? flag genes carried by 1 outlier guide).

## 4. Interpreting kinome hits

- Cross-reference hits against **DepMap/CCLE** dependency data for the same
  lineage (is this a known common vs. selective dependency?).
- Map hits onto **kinase families / pathways** (Manning groups, KEGG, Reactome)
  to find enriched signaling nodes rather than isolated genes.
- For therapeutic angle: intersect hits with **druggable kinases** and existing
  inhibitors (ChEMBL mechanism/bioactivity, clinical-trial landscape) to
  prioritize follow-up.

## 5. Validation of top hits

- **Multiple-guide concordance** (built in) is the first filter.
- Individual-guide arrayed KO / competition assays; **cDNA rescue** to prove
  on-target; where a tool compound exists, pharmacological confirmation.

## Outputs

- `results/guide_counts.tsv` (small — can link).
- `results/gene_summary.tsv` (MAGeCK/BAGEL gene-level — small, can link).
- QC report (HTML/PDF). Raw FASTQ and full count BAMs are large and
  **referenced by name** in the summary, not linked.
