# Defining the Target Gene Set

See [`../data/kinome_targets.md`](../data/kinome_targets.md) for the full
target-set definition, the Manning-classification group counts, authoritative
source URLs, and reproducible fetch commands (HGNC → symbols → Ensembl
transcript IDs).

## Summary

- **Core target set:** the **518 human protein kinases** (Manning et al. 2002).
- **Recommended default:** 518 protein kinases **+ ~20 lipid kinases ≈ 538**.
- **Extensions:** ~635 (add metabolic small-molecule kinases) or ~750–900
  (add pseudokinases and regulatory subunits) if the biological question calls
  for it.

## Why pull from an authority instead of hardcoding

Gene symbols drift and the kinome boundary is a judgment call (pseudokinases,
atypical PIKKs, metabolic kinases). Fetch from **HGNC gene group 1258
("Protein kinases")** and map to **Ensembl/GENCODE** transcripts so the guide
designer can target **constitutive early exons** per `02_guide_design.md`.
The exact per-gene guide list is a generated working file
(`data/selected_guides.tsv`), referenced by name.
