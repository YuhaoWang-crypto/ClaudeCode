# Retrospective check: 3QDG -> 4L3E

Peptide: `ELAGIGILTV`

Scanned **456 single mutants** of 3QDG with no knowledge of the engineered variant, then looked up where the known substitutions ranked.

| mutation | loop | rank | percentile | rank at position | predicted ddG | verdict | pep contacts | MHC contacts |
|---|---|---|---|---|---|---|---|---|
| `D:D26Y` | A-CDR1 | 112/456 | top 25% | 5/19 | +0.14 | reject-specificity | 0 | 14 |
| `E:L98W` | B-CDR3 | 117/456 | top 26% | 1/19 | +0.16 | below-threshold | 4 | 5 |

**Outcome:** FAILED: none of the 2 experimentally validated mutations would have been shortlisted; median rank was the 26th percentile. 1 of 2 ranked in the top 3 of the 19 substitutions available at their own position (E:L98W 1/19). So there is some signal in which substitution to choose at a given position, but not enough to clear the threshold for calling a mutation an improvement at all.

