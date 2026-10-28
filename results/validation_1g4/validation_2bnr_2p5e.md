# Retrospective check: 2BNR -> 2P5E

Peptide: `SLLMWITQC`

Scanned **494 single mutants** of 2BNR with no knowledge of the engineered variant, then looked up where the known substitutions ranked.

| mutation | loop | rank | percentile | rank at position | predicted ddG | verdict | pep contacts | MHC contacts |
|---|---|---|---|---|---|---|---|---|
| `D:Q51T` | A-CDR2 | 253/494 | top 51% | 14/19 | +0.72 | reject-specificity | 0 | 26 |
| `D:S52P` | A-CDR2 | 73/494 | top 15% | 9/19 | +0.15 | reject-specificity | 0 | 3 |
| `D:S53W` | A-CDR2 | 1/494 | top 0% | 1/19 | -0.73 | reject-specificity | 0 | 8 |
| `D:T95L` | A-CDR3 | 101/494 | top 20% | 5/19 | +0.25 | below-threshold | 3 | 5 |
| `D:G97D` | A-CDR3 | 59/494 | top 12% | 13/19 | +0.07 | below-threshold | 5 | 3 |
| `E:G50A` | B-CDR2 | 312/494 | top 63% | 16/19 | +0.90 | reject-specificity | 0 | 10 |
| `E:A51I` | B-CDR2 | 72/494 | top 15% | 4/19 | +0.14 | reject-specificity | 0 | 3 |
| `E:G52Q` | B-CDR2 | 55/494 | top 11% | 12/19 | +0.05 | reject-specificity | 0 | 3 |
| `E:I53T` | B-CDR2 | 401/494 | top 81% | 14/19 | +1.13 | reject-specificity | 0 | 22 |
| `E:T70I` | B-HV4 | 159/494 | top 32% | 3/19 | +0.47 | reject-specificity | 0 | 8 |
| `E:V95L` | B-CDR3 | 215/494 | top 44% | 3/19 | +0.63 | below-threshold | 6 | 9 |

**Outcome:** FAILED: none of the 11 experimentally validated mutations would have been shortlisted; median rank was the 20th percentile. 3 of 11 ranked in the top 3 of the 19 substitutions available at their own position (D:S53W 1/19, E:T70I 3/19, E:V95L 3/19). So there is some signal in which substitution to choose at a given position, but not enough to clear the threshold for calling a mutation an improvement at all.

Not reachable by the scan (the position fell below the contact threshold or sits outside the paratope): `D:R72S`, `D:S96L`, `D:S99T`

