# Guide Selection — Choosing Optimal sgRNAs

Goal: for each kinase gene, pick a small set of sgRNAs that (a) cut early and
in all transcripts to guarantee loss of function, (b) cut efficiently
(high on-target activity), and (c) are specific (minimal off-target risk,
covered in `03_offtarget_prediction.md`).

## 1. Nuclease and PAM

- **SpCas9**, PAM = `NGG`, 20-nt spacer. Default choice: best-characterized
  on/off-target models and the largest set of validated pre-made libraries.
- Consider **enhanced-specificity variants** (eSpCas9(1.1), HiFi Cas9) at the
  protein/mRNA level to reduce off-target cutting; they do not change guide
  design but improve editing precision.

## 2. Where to cut (positional rules)

1. **Target early constitutive coding exons** — restrict the cut to the
   5′-most exons shared by **all** protein-coding transcripts (constitutive),
   typically within the **5–65% CDS** window. Early frameshift → NMD / truncation.
2. **Avoid the very first exon / start codon region** where re-initiation at
   downstream ATGs can rescue function; and avoid the last exon (escapes NMD).
3. **Prefer exons encoding the kinase domain** (or upstream of it) so an
   in-frame indel still disrupts catalytic function.
4. **Avoid alternative/skipped exons** — a guide in a cassette exon misses
   isoforms lacking it.
5. Ensure the 20-mer maps **uniquely** to the intended locus (see off-target doc).

## 3. On-target efficiency scoring

Rank candidate spacers by a modern activity model and keep the top scorers:

| Model | Notes |
|-------|-------|
| **Rule Set 3 (Sanson et al. 2018)** | Current Broad standard; basis of the Brunello library. **Primary ranker.** |
| **DeepSpCas9 (Kim et al. 2019)** | Deep-learning on-target; strong independent cross-check |
| **Azimuth / Rule Set 2 (Doench et al. 2016)** | Widely used, well-validated fallback |

Take the **union rank** (e.g., mean percentile across Rule Set 3 + DeepSpCas9)
so a guide isn't chosen on a single model's quirk.

## 4. Hard sequence filters (reject if any fail)

- **Poly-T:** no `TTTT` in the spacer (premature Pol III termination).
- **GC content:** keep 30–70% (reject <25% or >80%).
- **Homopolymers:** no run ≥5 of any base; avoid extreme repeats.
- **Cloning-site collisions:** no internal **BsmBI (`CGTCTC`) / Esp3I** or
  **BbsI (`GAAGAC`)** sites (or their reverse complements) — they break
  Golden-Gate assembly. Match this to your vector's cloning enzyme.
- **Common SNPs:** reject spacers or PAMs overlapping high-frequency dbSNP
  variants in your cell line's background (guide won't match the actual genome).
- **U6 start:** if using the U6 promoter, a leading `G` is preferred; if the
  spacer doesn't start with G, prepend one (21-nt) or pick a `G`-starting site.

## 5. Guides per gene and redundancy

- **6 sgRNA/gene** is the recommended default for a boutique kinome screen —
  more redundancy than genome-wide libraries (which use 4) because a focused
  library is cheap to scale. Use **8/gene** if you want maximal confidence and
  can afford the cells/reads.
- Spread the 6 guides across **≥2 constitutive exons** where possible, so a
  single mis-annotated exon doesn't sink a gene.
- De-duplicate: never include two guides with identical or near-identical
  (≤2 nt offset) cut sites.

## 6. Controls (build into the library)

| Control type | Count | Purpose |
|--------------|-------|---------|
| **Non-targeting (NTC)** | ~200–250 | Null distribution / normalization baseline |
| **Intergenic "safe" cutting** | ~50 | Controls for cut-toxicity (DNA damage independent of gene KO) |
| **Positive essential (kinase)** | e.g. PLK1, CDK1, AURKB, EEF2K... | Expected strong dropout → confirms screen works |
| **Pan-essential (non-kinase)** | Hart CEGv2 subset | Benchmark for essentiality classifier (ROC/PR) |
| **Non-essential reference** | Hart NEGv2 subset | Negative benchmark for classifier calibration |

Intergenic cutting controls are important in a **kinase** screen specifically:
kinases include DNA-damage nodes, and you want to separate real dependency
from generic double-strand-break toxicity.

## 7. Reuse a validated library vs. design fresh

**Strongly prefer starting from a pre-validated library and subsetting to the
kinome**, then top up:

- **Brunello** (human, 4 sgRNA/gene, Rule Set 2/3, minimized off-target) —
  subset to kinome symbols → instant, empirically-validated core.
- **TKOv3** (Hart lab, ~71k guides genome-wide) — subset alternative.
- Add 2 extra fresh guides/gene (to reach 6) designed with the rules above via
  **CRISPOR**, **GuideScan2**, or **FlashFry** to fill gaps and boost redundancy.

Designing entirely from scratch is fine but redundant with well-benchmarked
public libraries; the hybrid approach gives validated guides plus redundancy.

## 8. Tooling summary

| Task | Tool(s) |
|------|---------|
| End-to-end candidate design + scores + off-target | **CRISPOR** (web/CLI), **GuideScan2** |
| Batch/genome-scale design | **FlashFry**, **GuideScan2** (precomputed) |
| On-target rescoring | Rule Set 3 (`rs3`), DeepSpCas9 |
| Constitutive-exon annotation | Ensembl BioMart / GENCODE basic + APPRIS principal isoform |

## Output of this stage

A per-gene table (`data/selected_guides.tsv`, generated by your run — this is a
working file, referenced by name) with columns:
`gene, transcript, exon, spacer(20nt), PAM, strand, cut_site, rs3_score,
deepcas9_score, GC, cfd_specificity, chosen(bool)`.
