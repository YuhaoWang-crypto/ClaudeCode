# Analysis Command Cheat-Sheet (MAGeCK / BAGEL2 / drugZ)

Ready-to-adapt commands stitching together `06_analysis.md` and the dual-line
analysis of `07/08`. Placeholders marked `<...>`. Assumes `mageck`, `BAGEL2`
(`BAGEL.py`), and `drugz` (`drugz.py`) are installed.

---

## 0. Input conventions

- **Library file** `library.csv`: three columns `sgRNA_id,sequence,gene`.
- **Samples** (per `templates/sample_sheet.csv`):
  `plasmid, T0_1..3, DMSO_1..3, DRUG_1..3`.

---

## 1. Counting: FASTQ → count matrix

```bash
mageck count \
  --list-seq library.csv \
  --sample-label plasmid,T0_1,T0_2,T0_3,DMSO_1,DMSO_2,DMSO_3,DRUG_1,DRUG_2,DRUG_3 \
  --fastq plasmid.fq T0_1.fq T0_2.fq T0_3.fq \
          DMSO_1.fq DMSO_2.fq DMSO_3.fq DRUG_1.fq DRUG_2.fq DRUG_3.fq \
  -n kinome_screen
# → kinome_screen.count.txt (guide×sample matrix) + .countsummary.txt (QC)
```

Check `.countsummary.txt`: mapping rate (≥65%), zero-count guides, Gini index.

---

## 2. Line A — Essentiality/dropout (vehicle vs. T0)

### 2a. BAGEL2 (primary, Bayes Factor)

```bash
# fold-change relative to T0 reference columns
BAGEL.py fc -i kinome_screen.count.txt -o kinome_fc -c T0_1,T0_2,T0_3

# Bayes Factor with CEGv2/NEGv2 priors
BAGEL.py bf -i kinome_fc.foldchange -o results/essentiality_bagel.tsv \
  -e CEGv2.txt -n NEGv2.txt -c DMSO_1,DMSO_2,DMSO_3
# higher BF = more essential; calibrate threshold with CEG/NEG precision-recall
```

### 2b. MAGeCK RRA (secondary, gene-level LFC/rank)

```bash
mageck test -k kinome_screen.count.txt \
  -t DMSO_1,DMSO_2,DMSO_3 -c T0_1,T0_2,T0_3 \
  -n results/essentiality_mageck
```

> Copy-number correction: pass essentiality results through **Chronos/CERES** so
> amplicon cut-toxicity isn't misread as dependency (`03` §5, `06` §3).

---

## 3. Line B — Chemogenomic/synthetic-lethal (drug vs. vehicle)

### drugZ (primary)

```bash
python drugz.py \
  -i kinome_screen.count.txt \
  -o results/synlethal_drugz.tsv \
  -c DMSO_1,DMSO_2,DMSO_3 \
  -x DRUG_1,DRUG_2,DRUG_3
# normZ < 0 & small FDR → sensitizer/synthetic lethal; normZ > 0 → resistance
# paired by default; add --unpaired for non-paired replicates
```

### MAGeCK MLE (secondary, multi-condition design matrix)

`design_matrix.txt`:

```
Samples     baseline   drug
DMSO_1      1          0
DMSO_2      1          0
DMSO_3      1          0
DRUG_1      1          1
DRUG_2      1          1
DRUG_3      1          1
```

```bash
mageck mle -k kinome_screen.count.txt -d design_matrix.txt \
  -n results/synlethal_mle
# read the 'drug' beta: negative = sensitizer, positive = resistance
```

---

## 4. Interpretation & validation checklist

1. **Pass QC gates first** (mapping rate, Gini, replicate corr ≥0.8, CEG/NEG ROC-AUC).
2. **Essentiality:** BF > threshold + 6-guide concordance → dependency kinase; map
   to Manning families/pathways.
3. **Synthetic lethal:** drugZ sensitizer hits + known categories appearing in the
   expected direction (`08` §4).
4. **Cross-reference** DepMap/CCLE same-lineage dependencies.
5. **Validate:** multi-guide concordance (built in) + arrayed single-guide KO +
   **cDNA rescue** to prove on-target.

---

## 5. Output files (small linkable; large referenced by name)

| File | Content |
|------|---------|
| `kinome_screen.count.txt` | guide×sample count matrix (small, linkable) |
| `results/essentiality_bagel.tsv` | essentiality BF/gene-level (small, linkable) |
| `results/synlethal_drugz.tsv` | synthetic-lethal drugZ results (small, linkable) |
| `results/qc_report.html` | QC report |
| raw FASTQ, full count BAMs | large — **referenced by name** in summaries, not linked |
