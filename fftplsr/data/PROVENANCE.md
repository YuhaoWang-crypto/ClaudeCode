# Provenance of the vendored PylRS datasets

Everything in this directory is a format conversion — no values were altered,
rounded, filtered or re-derived.

## Source

Hu, Y. *et al.* **Machine learning-guided evolution of pyrrolysyl-tRNA synthetase
for improved incorporation efficiency of diverse noncanonical amino acids.**
*Nature Communications* **16** (2025). doi:[10.1038/s41467-025-61952-2](https://doi.org/10.1038/s41467-025-61952-2)
· [PMC12274524](https://pmc.ncbi.nlm.nih.gov/articles/PMC12274524/)

Code and data repository: <https://github.com/zjuhaoran/FPFORCOM> (MIT License).
Corresponding author contact listed there: yuhaoran@zju.edu.cn

The article is open access under CC BY 4.0; the repository is MIT. Attribution is
retained here and in `fftplsr/datasets.py`.

## What was converted

| This directory | Converted from | Transformation |
|---|---|---|
| `IFRS.fasta` | `input/seq_IFRS.txt` | wrapped at 60 columns, FASTA header added |
| `Com1-IFRS.fasta` | `input/seq_Com1.txt` | wrapped at 60 columns, FASTA header added |
| `trainset1_ifrs_singles.csv` | `input/data.xlsx`, sheet `Trainset1` | `.xlsx` → `.csv` |
| `trainset2_ifrs_combos.csv` | `input/data.xlsx`, sheet `Trainset2` | `.xlsx` → `.csv` |
| `trainset3_com1_singles.csv` | `input/data.xlsx`, sheet `Trainset3` | `.xlsx` → `.csv`; `Unnamed: 2` renamed `Source` |
| `trainset4_com1_combos.csv` | `input/data.xlsx`, sheet `Trainset4` | `.xlsx` → `.csv` |
| `measured_ifrs_panel.csv` | `data/IFRS-mutants-data.json` | EnzymeML → 2-column table |
| `measured_com1_panel.csv` | `data/Com1-IFRS-mutants-data.json` | EnzymeML → 2-column table |

The two panels are extracted from EnzymeML records, where each measured variant
appears as a `proteins[]` entry (giving the variant name) and a matching
`measurements[0].species_data[]` entry whose **`initial`** field holds the relative
activity. The non-synthetase entry `sfGFP-3BrF` (the reporter) carries no
measurement and is excluded.

Reproduce the conversion with `scripts/convert_paper_data.py`.

## Units and meaning

`Fitness` is **relative stop-codon-suppression efficiency**, measured as
`Fluorescence / OD600^t` on a sfGFP reporter bearing the noncanonical amino acid
3-bromophenylalanine (3BrF), normalised so the parent enzyme on that background
equals 1.0. Values are therefore fold-changes, not absolute rates, and are only
comparable **within** a background.

Two backgrounds, which must not be mixed:

- **IFRS** (454 aa) — parent of `trainset1`, `trainset2`, `measured_ifrs_panel`.
- **Com1-IFRS** (454 aa) — parent of `trainset3`, `trainset4`, `measured_com1_panel`.
  It differs from IFRS at 7 positions: `D2N/V31I/T56P/R61K/H62Y/T122S/S193R`
  (read off the two sequences; the article's narrative mentions fewer).

## Overlap between tables

`trainset2` ⊃ `trainset1` (13 of its 38 rows), and both are subsets of
`measured_ifrs_panel`. `trainset4` is a subset of `measured_com1_panel`. The
held-out splits in `fftplsr/m3_baselines.py` are built by set difference, so a
variant is never scored against a model that was trained on it.

Caveat on `measured_com1_panel` minus `trainset4`: those 21 variants are the ones
the authors' own model nominated and then assayed, so they are a
high-fitness, narrow-range sample (0.89–2.75). R² computed on them suffers
restriction of range and should not be read as a general accuracy figure — rank
metrics are reported alongside it for that reason.
