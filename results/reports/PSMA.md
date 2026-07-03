# PSMA (Q04609) — evidence report

- **ChEMBL** (CHEMBL1892): 465 actives (308 small-molecule).
- **DrugCLIP**: 414 predicted hits (MW median 298).
- **humanPPI**: 1456 predicted partners (244 cell-surface).

## Molecular cross-comparison (DrugCLIP vs measured)
- whole-molecule NN-Tanimoto: target median 0.162 vs decoy 0.150; MWU p=2e-07; EF@0.35=0.00×
- fragment-in-active containment: 1/414 vs decoy 0/874; Fisher p=0.32, OR=inf
- MW mismatch: DrugCLIP 298 vs actives 503 Da

## Protein-layer enrichment (humanPPI)
- PSMA cell-surface partners: 244/1456 = 16.8%
- vs P04406: 17.1% (p=0.6); vs P00533: 20.3% (p=1); vs genome ~25%: p=1
