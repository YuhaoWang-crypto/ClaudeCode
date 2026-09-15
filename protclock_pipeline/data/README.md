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

## 3. Clock weights (mostly obtainable)

Install the paper's own library and three of the six clocks run on their real
published weights with no further work:

```bash
git clone https://github.com/Insilico-org/proteoclock
cd proteoclock && pip install -e . && pip install --upgrade scikit-posthocs
```

The upgrade is required. `setup.py` pins `scikit-posthocs==0.11`, which
imports `multipletests` from `statsmodels.sandbox.stats.multicomp`, removed
in statsmodels 0.15. Without it the package raises on import. pip will warn
that the pin is violated; the package works anyway.

| Clock | Status | Where |
|---|---|---|
| `PAC` | **real weights** | proteoclock `kuo_2024` (Gompertz, needs age) |
| `OrganAge_chrono` | **real weights** | proteoclock `goeminne_2025_full_chrono` |
| `OrganAge_mortality` | **real weights** | proteoclock `goeminne_2025_full_mortality` |
| `PAOPAC` | blocked on platform | `github.com/JackieHanLab/PAOPAC`; `model.bin` is on Releases but the extension is a Windows-only `.pyd` for Python 3.9 |
| `ipfP3GPT` | restricted | code at `osf.io/457w8`; proteoclock ships only the feature order. Weights are usable solely inside the UK Biobank Research Analysis Platform |
| `ProtAge` | on request | not in `github.com/miargentieri/proteomic-age-ukb`; email the authors |

To force a specific coefficient set instead, drop a CSV per clock into
`data/weights/` with columns `protein` and `coefficient` plus an optional
`(Intercept)` row, and pass `--weight-dir data/weights`. That overrides
proteoclock. Any clock with neither falls back to a surrogate, and `run_all`
reports how many of each it used.

## 4. Supplementary tables (open, fetched automatically)

`m8_supplementary` downloads the workbook from static-content.springer.com on
first use and caches it here as `supplementary_tables.xlsx`. It is
git-ignored, being publisher-hosted content. No credentials are needed even
though the article itself is paywalled.

## Identifier namespaces

The likeliest failure is not access, it is naming. Olink exports use assay
names, the organAging tables use UniProt accessions, UK Biobank uses its own
identifiers. A total mismatch raises. A **partial** match does not, and it
silently scores a clock on a fraction of its proteins, which produces
plausible and wrong numbers. `m3_clocks.load_linear_weights` prints how many
coefficients matched; check it against the clock's declared `n_features`
every single time.
