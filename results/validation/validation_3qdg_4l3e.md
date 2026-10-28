# Retrospective check: 3QDG -> 4L3E

Peptide: `ELAGIGILTV`

Scanned **456 single mutants** of 3QDG with no knowledge of the engineered variant, then looked up where the known substitutions ranked.

| mutation | loop | rank | percentile | rank at position | verdict | pep contacts | MHC contacts |
|---|---|---|---|---|---|---|---|
| `D:D26Y` | A-CDR1 | 216/456 | top 47% | 15/19 | reject-specificity | 0 | 14 |
| `E:L98W` | B-CDR3 | 199/456 | top 44% | 5/19 | below-threshold | 4 | 5 |

**Outcome:** FAILED: none of the 2 experimentally validated mutations would have been shortlisted. Median rank was the 47th percentile, which is indistinguishable from picking at random. The score should not be trusted to enrich for real affinity gains.

