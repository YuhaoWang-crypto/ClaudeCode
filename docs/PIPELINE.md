# Generic target evidence pipeline (Problem B)

One command mines all three predicted/measured sources for any human target and runs
the molecular + fragment + protein-layer enrichment analyses. Engine:
`lipidlib/targetpipe.py`; CLI: `analysis/run_target.py`.

## Usage

```bash
python analysis/run_target.py --uniprot P43220 --name GLP1R
python analysis/run_target.py --uniprot Q04609 --name PSMA \
    --drugclip-decoy P00533 --ppi-decoys P04406,P00533
```

Outputs per target:
- `data/targets/<NAME>/` — `*_chembl_ligands.csv`, `*_drugclip_predicted.csv`,
  `*_humanppi_partners.csv`, `*_surface_partners.csv`
- `results/figures/<NAME>_molecular.png`, `<NAME>_humanppi.png`
- `results/reports/<NAME>.md` — a one-page summary

It degrades gracefully: a target with no DrugCLIP coverage (e.g. ASGR1) or sparse
ChEMBL data still completes, skipping the sections it can't compute.

### What the engine does (`lipidlib/targetpipe.py`)
| function | purpose |
|---|---|
| `resolve_chembl_target(uniprot)` | UniProt → ChEMBL target id |
| `mine_chembl / mine_drugclip / mine_humanppi` | the three source miners |
| `molecular_enrichment(fg, ref, decoy)` | nearest-neighbour Tanimoto vs decoy (MWU, EF) |
| `fragment_enrichment(fg, decoy, actives)` | substructure containment vs decoy (Fisher) |
| `ppi_surface_enrichment(target, decoys)` | cell-surface fraction vs decoys + genome |

## Cross-target results so far

| target | ChEMBL sm | DrugCLIP (MW) | molec. enrichment | PPI surface enrichment |
|---|---|---|---|---|
| **GLP1R** | 693 | 173 (257) | none (EF 0×) | **yes** — 36.7% vs 17–20% decoys (p≈1e-10) |
| **PSMA** | 308 | 414 (298) | none (EF 0×; MWU p=2e-7 but nothing ≥0.35) | no — 16.8%, n.s. |
| **ASGR1** | 1 | none | n/a | (162 partners; ChEMBL/DrugCLIP too sparse) |

Two observations that generalise:
- **DrugCLIP is a fragment library** — hit MW ≈ 250–300 Da for *every* target checked
  (GLP1R, PSMA, EGFR, CDK2, thrombin, AR, matriptase). It will systematically look
  un-enriched against drug-sized (MW ≈ 500–600) measured actives on whole-molecule
  or fragment metrics.
- **The PPI surface-enrichment test discriminates**: GLP1R is significantly
  surface-enriched, PSMA is not — so a positive result is meaningful, not automatic.

---

# Pipeline validation (does it enrich when it should?)

`analysis/validate_pipeline.py`. Feeds the *same* `molecular_enrichment` routine a
foreground known to share chemistry with the reference — a held-out 40% split of a
target's own measured actives — plus DrugCLIP hits and an unrelated target's actives,
all scored against a 60% reference split of PSMA actives.

| foreground vs PSMA reference | median NN-Tanimoto | EF @ 0.35 | MWU p |
|---|---|---|---|
| **held-out PSMA actives (positive control)** | **0.776** | **1674×** | 9.6e-82 |
| PSMA DrugCLIP hits | 0.157 | 0.0× | 1.0 |
| EGFR actives (decoy floor) | 0.168 | 1.0× | — |

**Conclusion:** the machinery detects genuine chemical similarity emphatically
(1674× enrichment for held-out actives), so the GLP1R/PSMA DrugCLIP negatives are
**real, not artefacts** — DrugCLIP's fragment hits simply don't share chemotype with
the drug-sized measured actives. Figure: `results/figures/pipeline_validation.png`.
