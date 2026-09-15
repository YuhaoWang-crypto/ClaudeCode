# Real inputs

This directory is empty on purpose. The three things needed to turn the
demonstration into a reproduction all live behind access control or an author
request. Put them here and the pipeline switches over without code changes.

## 1. Trial proteome (blocking)

| | |
|---|---|
| Accession | `OMIX008341`, file `OMIX008341-01` |
| Repository | China National Center for Bioinformation |
| URL | https://ngdc.cncb.ac.cn/omix/release/OMIX008341 |
| Contents | "OLINK data for IPF patient plasma", Olink Explore 3072 |
| Format | XLSX, 3.9 MB |
| Access | **Controlled.** Submit a request; approval is not automatic. |

Registered against the 2025 *Nature Medicine* rentosertib trial paper rather
than the 2026 aging-clock paper. Confirm on receipt that it carries all four
visits (baseline, weeks 2, 4, 12) and the treatment assignment, not baseline
only. Without arm labels and visit numbers nothing downstream runs.

Save as `data/OMIX008341-01.xlsx`.

## 2. Reference cohort (blocking)

The paper calibrates age acceleration against 55,319 UK Biobank participants
with Olink Explore 3072. This needs an approved UK Biobank application, which
takes months.

Any large Olink cohort with age and sex can substitute. Doing so makes the
age-acceleration residuals non-comparable with the paper's, because the
reference distribution is what defines the zero point. Say so in any result
produced that way.

Save as `data/reference_npx.parquet` plus `data/reference_meta.csv` with
columns `age`, `sex`, and `died_10y` if you intend to fit mortality clocks.

## 3. Clock weights (partly obtainable)

Drop one CSV per clock into `data/weights/`, named after the clock, with
columns `protein` and `coefficient` and an optional `(Intercept)` row. Then
pass `--weight-dir data/weights`.

| Clock | Status | Where |
|---|---|---|
| `OrganAge_chrono` | obtainable | Goeminne 2025 supplementary Table S1A/S1C |
| `OrganAge_mortality` | obtainable | same, mortality coefficients, no intercept |
| `PAC` | obtainable | `github.com/kuo-lab-uchc/PAC`, coefficients inside `pac_proteomic_age.R` |
| `PAOPAC` | unverified | `github.com/41way5/Organ-PAC` ships a training script and no README |
| `ProtAge` | on request | not in `github.com/miargentieri/proteomic-age-ukb`; email the authors |
| `ipfP3GPT` | unavailable | Precious3GPT base model is on HuggingFace, the IPF fine-tune is not |

Any clock without a CSV falls back to a surrogate, and `run_all` reports how
many of each it used.

## Identifier namespaces

The likeliest failure is not access, it is naming. Olink exports use assay
names, the organAging tables use UniProt accessions, UK Biobank uses its own
identifiers. A total mismatch raises. A **partial** match does not, and it
silently scores a clock on a fraction of its proteins, which produces
plausible and wrong numbers. `m3_clocks.load_linear_weights` prints how many
coefficients matched; check it against the clock's declared `n_features`
every single time.
