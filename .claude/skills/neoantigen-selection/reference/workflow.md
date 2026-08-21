# The workflow this mirrors: what is public, what is not, what you must supply

## The publicly disclosed steps

The individualized neoantigen therapy class of product — intismeran autogene
(formerly mRNA-4157/V940; Moderna/Merck), and the same shape of pipeline behind
autogene cevumeran (BioNTech/Genentech) — is described publicly as:

1. obtain tumor tissue and matched normal blood
2. sequence tumor and normal DNA; obtain tumor RNA expression
3. call somatic mutations present only in the tumor
4. determine the patient's HLA type
5. predict, algorithmically, which mutant peptides are presented on HLA-I/HLA-II
   and may be recognized by T cells
6. select up to 34 patient-specific neoantigens
7. encode them in one individualized synthetic mRNA, delivered as an LNP

Everything in that list is a *shape*, not a specification. The scoring function,
its training data, the thresholds, the construct rules and the manufacturing
release criteria are proprietary and are **not** reproduced by this skill.
See Moderna's public description:
<https://www.modernatx.com/media-center/all-media/blogs/advancing-fight-against-cancer>

## What is genuinely open, and used here

| need | open source used | note |
|---|---|---|
| somatic variants (demo) | cBioPortal REST, TCGA PanCanAtlas | real calls with tumor read counts |
| somatic variants (real patient) | your own tumor/normal caller | the `clairs-somatic` skill covers ONT/HiFi/Illumina tumor-normal calling |
| tumor RNA | cBioPortal RSEM (demo) / your Salmon-RSEM TPM | gene-level in the demo, transcript-level preferred |
| reference proteome | UniProt reviewed human (~20.4k canonical seqs) | mutant protein reconstruction + self k-mer novelty gate |
| MHC-I / MHC-II presentation | NetMHCpan-4.1 / NetMHCIIpan via the IEDB cloud REST API | no licence, no local install; `mhcflurry` swaps in locally |
| T-cell-recognition prior | IEDB `query-api` positive T-cell assays | also the benchmark's ground truth |

## What is *not* open, and what you must supply

**HLA typing.** There is no open endpoint that returns a TCGA sample's
four-digit class-I type — the PanImmune OptiType calls are controlled access.
For a real patient this comes from OptiType / xHLA / HLA-HD run on the
**normal** sample (never the tumor: LOH at the HLA locus corrupts tumor-derived
typing). The demo therefore declares a common class-I haplotype and labels the
whole run `[assumed]` on that point. Everything downstream inherits this
assumption; it is the single largest error source in the pipeline.

**Tumor purity and copy number.** Needed for a real CCF. From ABSOLUTE, FACETS,
Sequenza or PureCN. The demo fixes purity at 0.7 and CN at 2 and says so.

**Allele-specific expression.** The mutant-allele fraction in RNA, not just the
gene's TPM, is what catches nonsense-mediated decay and allelic silencing. It
needs the patient's own aligned RNA-seq, which cBioPortal does not expose.

**Frameshift / neo-ORF and fusion neoantigens.** These need transcript-level
(nucleotide) annotation to translate the novel reading frame. The pipeline
does not invent them: such variants land in `peptides_skipped` with the reason,
and can be re-introduced by supplying a translated `neo_orf` column. In a
melanoma exome this is a small fraction of candidates; in an MSI-high tumor it
is the majority of the good ones, so do not skip it there.

## Where each output goes in a real program

| output | used for |
|---|---|
| `gate_waterfall.csv` | the "why only N candidates" conversation; catches a bad expression matrix immediately |
| `ranked.csv` | the full ordered candidate list with all per-feature columns |
| `selected.csv` | the payload, with `why_selected` per slot |
| `coverage.csv` | HLA spread — the hedge against allele-loss escape |
| `minigenes.csv` | the 25-mers in construct order |
| `junction_scan.csv` | junction epitopes created by concatenation; anything `flagged` is a non-tumor decoy epitope |
| `construct.fasta` | the concatemer protein and its codon-optimized CDS |
| `REPORT.md` | the labelled write-up |

## Adjacent skills

- `clairs-somatic` — tumor/normal BAM → somatic VCF (the step upstream of this one)
- `immunogenicity-multimodel` — independent DTU/ImmunoGeNN and IEDB cross-check on the final construct
- `codon-optimize-qc` — production-grade codon optimization and gene-synthesis QC
- `lnp-delivery-kinetics` — LNP uptake / endosomal escape / expression ODE modelling of the delivered mRNA
- `synbio-cassette-designer` — if the payload needs cell-state-restricted expression rather than a plain CDS
