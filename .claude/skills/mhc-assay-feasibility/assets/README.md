# Installing the feasibility pipeline

Paths resolve relative to the parent of `scripts/`:

```
<project>/
  scripts/     <- assets/scripts/*
  data/        <- _iedb_popcov/ (fetched once), candidate_peptides.tsv (optional)
  results/     <- created by the modules
```

```bash
mkdir -p myproject && cd myproject
cp -r <skill>/assets/scripts .
pip install numpy
python scripts/c5_donor_feasibility.py     # seconds — run this first
```

`c5` downloads and caches the IEDB population-coverage package into
`data/_iedb_popcov/` on first run.

## Pointing it at a different study

Everything study-specific is a constant at the top of a module, deliberately —
these are properties of the study, not of the pipeline, so they are stated with
their evidence rather than buried in a config:

| What to change | Where |
|---|---|
| The alleles under study | `ALLELES` in `c1_benchmark.py`, `c4_prescreen.py`, `c5_donor_feasibility.py` |
| Peptide length | `LENGTH` in `c1_benchmark.py` |
| Candidate peptides | `data/candidate_peptides.tsv` with columns `id`, `peptide` |
| Assumed binding prior | `ASSUMED_PREVALENCE` in `c3_calibrate.py` |
| How much sensitivity to keep | `MIN_SENSITIVITY` in `c3_calibrate.py` |
| Class II system + panel | `ALLELE`, `COVA`, `PANEL` in `c6_nanopeptide_classII.py` |
| Donor cohort sizes | `WANT_DONORS` in `c5_donor_feasibility.py` |

## What to re-run when

| Change | Re-run |
|---|---|
| New candidate peptides, same alleles | `c4`, `c7` |
| New alleles | everything — the calibration is per allele |
| New assumed prior or sensitivity floor | `c3` (analysis only; the scoring cache is reused), `c4`, `c7` |
| Donor requirement changed | `c5`, `c7` |

`c3` is the long pole (~30–60 min for ~13k peptides). It appends to
`results/c3_benchmark_scored.tsv` as it goes and resumes from it, so an
interrupted run costs only the unscored remainder. That cache is regenerable and
belongs in `.gitignore`.

## The one thing not to skip

`c6` runs the system's own positive control. If the restriction element and its
reference epitope do not reproduce, stop — nothing downstream is interpretable.
