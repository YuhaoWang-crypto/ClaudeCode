---
name: crispr-screen-scale
description: Compute the scale and sequencing depth of a pooled CRISPR screen — cells to transduce, coverage to maintain at every bottleneck, gDNA per sample, and reads per sample — from library size, target coverage, MOI, replicates, and arms/timepoints. Use when planning a CRISPR screen's cell numbers, gDNA input, or NGS read depth, or when someone asks "how many cells / how much sequencing" for a pooled screen.
---

# Pooled CRISPR screen scale & sequencing-depth calculator

Turns a handful of biological parameters into the numbers that keep library
representation intact at every bottleneck (transduction, selection, passaging,
gDNA PCR template, reads). Prints all assumptions so the output is auditable.

## Run

```bash
python3 scripts/coverage_calculator.py \
  --genes 518 --guides-per-gene 6 --controls 400 \
  --coverage 1000 --moi 0.3 \
  --replicates 3 --arms 2 --timepoints 2
```

## Key parameters

| Flag | Meaning | Typical |
|------|---------|---------|
| `--genes`, `--guides-per-gene`, `--controls` | library size | — |
| `--coverage` | cells (and reads) per guide, held at every step | 500 min, **1000** for dropout |
| `--moi` | functional MOI for transduction | 0.3 (≈1 guide/cell) |
| `--reads-per-guide` | sequencing depth target | 500–1000 |
| `--ploidy-pg` | gDNA mass per cell | 6.6 pg (human diploid) |
| `--arms`, `--timepoints`, `--replicates` | sample accounting | e.g. 2 × 2 × 3 |

## What it returns

- **Integrants needed** and **cells to transduce** (Poisson-corrected for MOI).
- **Cells to maintain** at every passage/timepoint (the non-negotiable).
- **gDNA per sample** (sized to preserve coverage in the PCR) + number of PCR reactions.
- **Reads per sample** and **total reads** across all samples.

## Interpreting

The single rule: coverage must hold at **every** bottleneck, not on average.
Under-sampling anywhere (too few cells at a passage, too little gDNA, too few
reads) permanently loses guides to drift and creates false hits. Use the printed
"maintain ≥ N cells" figure as the floor for every step of the screen.
