# Unified aptamer dataset schema (Apta-Index ∪ UTexas ∪ literature)

Normalize every source into these columns (`seed_dataset.csv` uses them; `ingest.py`
maps source columns onto them). One row = one aptamer with a known target.

| column | type | notes |
|---|---|---|
| `apt_id` | str | stable id, e.g. `UTX_0123`, `APTAGEN_045`, `LIT_tba` |
| `name` | str | descriptive name |
| `chem` | enum | DNA \| RNA \| 2F-RNA \| LNA \| L-DNA \| L-RNA \| XNA \| chimeric |
| `sequence` | str | 5'→3', uppercase, ACGT/ACGU only (strip mods to a parallel column) |
| `mods` | str | modification notes (2'-F pyrimidines, 3'-invdT, PEG, ...) |
| `length` | int | nt |
| `target_name` | str | e.g. "Thrombin", "EGFR", "Nucleolin" |
| `target_type` | enum | protein \| cell \| peptide \| tissue \| small-molecule \| other |
| `target_seq` | str | target protein/peptide sequence if applicable (for ML features) |
| `target_uniprot` | str | accession if protein |
| `kd_value` | float | numeric affinity |
| `kd_unit` | enum | pM \| nM \| uM \| mM \| fM \| CFU/mL |
| `kd_nM` | float | **normalized** affinity in nM (derived; use this for ML) |
| `buffer` | str | binding buffer / conditions |
| `specificity` | str | cross-reactivity notes |
| `source` | enum | aptagen \| utexas \| literature |
| `doi` | str | citation |
| `label` | int | 1 = experimentally confirmed binder (default); 0 reserved for negatives |

## Sources & access (verified July 2026)
- **UTexas Aptamer Database** (~1,475 curated, **downloadable public dataset**, 1990–2022):
  https://sites.utexas.edu/aptamerdatabase/ — the primary training source. NAR 2024, D351.
- **Aptagen Apta-Index** (~800 entries): web search/filter only, **no API/CSV/bulk export**;
  commercial DB — respect its terms of use; use for spot-checking, not bulk scraping.
- **Literature** (bioRxiv/PubMed MCP): fill gaps, recent aptamers with Kd.

## Negatives (there is no public "non-binder" set)
Generate synthetic negatives by **target-shuffling**: pair each aptamer with a random
*different* target (assumed non-binding). Keep a held-out positive set per target for
evaluation. Mark generated rows `label=0`, `source=synthetic-negative`.

## Affinity normalization (`kd_nM`)
pM→×0.001, nM→×1, uM→×1000, mM→×1e6, fM→×1e-6. CFU/mL cannot be converted → leave `kd_nM` null.
