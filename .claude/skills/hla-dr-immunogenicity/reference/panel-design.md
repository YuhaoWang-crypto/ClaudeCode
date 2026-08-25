# Designing an HLA-DR panel to a stated coverage target

## Why not use the standard reference set

The DR subset of the IEDB class-II reference set (~15 molecules) is widely
quoted as "representative". Measured against the IEDB allele-frequency tables
it reaches **84.2 %** weighted US/EU DRB1 phenotypic coverage. If the
requirement is 95–98 %, it does not meet it, and quoting it as though it does
is the single most common defect in these assessments.

A greedy build to a stated target reaches **97.3 %** with 21 DRB1 alleles, plus
4 DRB3/4/5 molecules for presentation breadth — 25 molecules total. Coverage is
strongly diminishing: state the target, build to it, and publish the coverage
curve so the reader can see where the knee is.

## The coverage model, and the renormalisation most implementations miss

Single-locus Hardy–Weinberg phenotypic frequency over the IEDB DRB1
allele-frequency tables:

```
cov = 1 - (1 - f)^2,   f = sum(panel allele frequencies) / N
```

**`N` is the population's total DRB1 allele frequency when that total exceeds
1, and 1 otherwise.** This is the part that is easy to get wrong. Several
curated IEDB populations sum above 1 — Europe sums to **1.040** — and skipping
the renormalisation overstates coverage. On a 20-allele panel it inflated
Europe by ~0.8 points here, which was enough to change which alleles the greedy
build selected.

```python
class CoverageModel:
    def __init__(self, tables, weights):
        self.norm = {p: max(sum(t.values()), 1.0) for p, t in tables.items()}

    def population(self, panel, pop):
        drb1 = [a for a in panel if a.startswith("HLA-DRB1")]
        f = min(sum(self.tables[pop].get(a, 0.0) for a in drb1) / self.norm[pop], 1.0)
        return 1.0 - (1.0 - f) ** 2
```

**Verify against the IEDB CLI before trusting any coverage number you publish.**
`calculate_population_coverage.py -c II` on the same allele list should agree to
two decimals. If it does not, the discrepancy is in the normalisation or in
which population tables were summed — not in the Hardy–Weinberg step.

## DRB3/4/5 — report separately, never double-count

DRB3, DRB4 and DRB5 carry **no frequencies in the DRB1-locus tables**. They are
a second DR molecule expressed alongside DRB1 on many haplotypes, so they add
real presentation breadth, but they contribute nothing to the DRB1 coverage
arithmetic. Including them in the coverage sum is a straightforward
double-count.

Handle them as this pipeline does: add them to the prediction panel, report the
molecule count as "21 DRB1 + 4 DRB3/4/5", and keep the coverage percentage
strictly a DRB1 number.

## Population weighting

The headline is a weighted composite, not a single population:

- **US composite**, weighted by US Census 2020 race/ethnicity shares across the
  IEDB "United States Caucasoid / Hispanic / Black / Asian / Amerindian /
  Mestizo" tables.
- **Europe**, the IEDB curated table.
- Combined 50:50 by default (`panel.us_eu_split`).

All of it lives in `config/config.yaml`. Change the shares there rather than in
code, and re-run M2 — the greedy selection depends on the weights, so a
different market means a different panel, not just a different percentage.

## Widening a panel later

`m3_binding_prediction.py` is resumable: it reads any existing
`results/m3_binding_long.tsv` and fetches only the (sequence, allele, head)
combinations missing from it. Adding an allele therefore costs one allele's
worth of prediction, not a full re-run. Note that the calibration layer
(M10–M13) *is* panel-dependent and does need re-running when the panel changes.
