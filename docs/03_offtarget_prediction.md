# Off-Target Prediction

Every candidate spacer is scored for the risk that Cas9 cuts elsewhere in the
genome. In a pooled screen, off-target cutting creates false dependencies
(a "hit" that is really a bystander gene) and adds cut-toxicity noise, so
filtering on specificity is as important as on-target efficiency.

## 1. What drives off-target cutting

- **Sequence similarity** of a genomic site to the spacer (mismatches tolerated,
  especially PAM-distal / 5′ end).
- **Mismatch position and identity** — PAM-proximal (seed, ~10–12 nt) mismatches
  are far more disruptive than PAM-distal ones. rG:dT wobble mismatches are
  often tolerated.
- **Non-canonical PAMs** — SpCas9 also cuts weakly at `NAG`/`NGA`; include them
  when enumerating.

## 2. Enumerate candidate off-target sites

Exhaustively find all genomic sites within an edit distance of each spacer:

- **Cas-OFFinder** — exhaustive, alignment-free enumeration up to N mismatches
  (and bulges), GPU-accelerated. Ground truth for "what sites exist."
- **GuideScan2** — precomputed, genome-wide specificity with an efficient
  index; returns off-target counts by mismatch and an aggregate specificity.
- **CRISPOR** — convenient for candidate sets; reports MIT + CFD specificity
  and lists individual off-targets with annotations (exon/intron/intergenic).

Search parameters: enumerate sites with **≤4 mismatches** to `NGG`, and
additionally check `NAG`/`NGA` PAMs at ≤3 mismatches.

## 3. Score specificity

Two complementary aggregate scores (0–1, higher = more specific):

- **CFD specificity score** (Doench 2016) — position/identity-weighted; the
  **primary** off-target metric.
- **MIT specificity score** (Hsu 2013) — legacy but still reported; use as a
  secondary check.

Compute an aggregate genome-wide specificity per guide (as GuideScan2/CRISPOR
do) — it summarizes the whole off-target profile into one number.

## 4. Filtering rules (reject a guide if it fails)

Apply in order, hard first:

1. **No perfect (0-mismatch) match anywhere except the intended site.**
2. **No 1-mismatch off-target that lands in a CDS/exon** of any gene
   (especially another kinase — that would confound the screen).
3. **No 2-mismatch off-target in a seed region + exon** combination.
4. **CFD specificity ≥ 0.2** (raise to ≥ 0.3–0.5 if you can afford to be strict;
   a focused library can afford strictness because you have many candidates per
   gene).
5. Prefer guides whose only off-targets are **intergenic/intronic** and ≥3
   mismatches away.

Then, among survivors, rank by on-target score (see `02_guide_design.md`) and
pick the top 6/gene.

## 5. Cell-line-aware refinement (recommended)

Off-target risk depends on the **actual genome you screen in**, not the
reference:

- Use the cell line's variant calls (if available, e.g., from CCLE/DepMap or
  your own WGS) to mask sites where SNPs create or destroy off-target matches.
- Copy-number amplifications create **cut-toxicity** artifacts (many cuts in an
  amplified region kill cells regardless of gene function). Flag guides whose
  intended or off-target cut sites fall in amplicons (Munoz 2016 / CERES-style
  copy-number bias) — the analysis stage (`06_analysis.md`) also corrects for
  this via CERES/Chronos.

## 6. Empirical validation (optional, usually not needed for pooled screens)

If a specific high-value hit must be de-risked for off-target effects:

- **GUIDE-seq**, **CIRCLE-seq / CHANGE-seq**, or **DISCOVER-seq** to measure
  genome-wide cutting empirically for individual guides.
- For screen hits, the cleaner validation is **multiple independent guides
  agreeing** (built in via 6 guides/gene) plus **cDNA rescue**.

## Output of this stage

A per-guide off-target table (large, genome-wide — **referenced by name**,
not linked, in the summary): `data/offtarget_report.tsv.gz` with
`spacer, cfd_specificity, mit_specificity, n_off_0mm, n_off_1mm, n_off_2mm,
n_off_exonic, worst_offtarget_locus, pass(bool)`.
