# MHC-I epitope report — NY-ESO-1 (CTAG1B, P78358)

A second real target, chosen to **validate the pipeline against a known
epitope**. Produced live through the integrated stack (IEDB NetMHCpan_el) — no
token, no local license.

## Target & panel

- **Protein**: NY-ESO-1 / CTAG1B, 180 aa. A cancer-testis antigen and one of the
  most-studied targets in cancer immunotherapy (vaccines, TCR-T cell therapy).
  Here the goal is the **opposite** of the asparaginase case: *find* the strong
  T-cell epitopes rather than remove them.
- **HLA-I panel**: A\*02:01, A\*01:01, A\*03:01, A\*24:02, B\*07:02, B\*35:01.
- **Model**: IEDB NetMHCpan_el, 9-mers. 1032 peptide×allele predictions.

## Validation against the literature ✅

The canonical NY-ESO-1 epitope **SLLMWITQC (157–165)** is a well-characterized
**HLA-A\*02:01-restricted** CD8 epitope (used in NY-ESO-1 TCR-T therapies). The
pipeline recovers exactly that restriction:

| allele | %Rank of SLLMWITQC | call |
|--------|--------------------|------|
| **HLA-A\*02:01** | **0.63** | binder (top 0.6%) |
| HLA-A\*03:01 | 9.6 | none |
| HLA-A\*24:02 | 19 | none |
| HLA-A\*01:01 | 32 | none |
| HLA-B\*07:02 | 38 | none |
| HLA-B\*35:01 | 40 | none |

Correct HLA-A\*02:01 specificity, and a binder-grade rank on the right allele
only — a clean sanity check that the cloud pipeline reproduces known biology.
(0.63 sits just over the strict 0.5 strong/weak line; it is unambiguously an
A\*02:01 ligand and a non-binder everywhere else.)

## Epitope load

17 strong + 28 weak MHC-I binders (`results_nyeso/mhc1_all.csv`). Top epitopes
(`results_nyeso/mhc1_consensus.csv`):

| peptide | alleles | best %Rank |
|---------|---------|-----------|
| APRGPHGGA | B\*07:02 | 0.06 |
| MPFATPMEA | B\*07:02, B\*35:01 | 0.12 |
| GPESRLLEF | B\*07:02, B\*35:01 | 0.13 |
| SLAQDAPPL | A\*02:01 | 0.14 |
| RLLEFYLAM | A\*02:01 | 0.21 |
| FATPMEAEL | A\*02:01, B\*35:01 | 0.29 |

The C-terminal region (~150–170), which contains SLLMWITQC and several other
A\*02:01/B\*35:01 binders, is the epitope-dense stretch — consistent with why
this region dominates NY-ESO-1 immunotherapy design.

Per-residue landscape (known epitope marked): `figures/nyeso1_mhc1_landscape.png`.

## Reproduce

```python
from pipeline import integrate, aggregate
alleles = ["HLA-A*02:01","HLA-A*01:01","HLA-A*03:01",
           "HLA-A*24:02","HLA-B*07:02","HLA-B*35:01"]
res  = integrate.run_iedb_mhci("examples/nyeso1.fasta", alleles)
cons = aggregate.consensus([res], top_n=15)
```
