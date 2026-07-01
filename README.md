# Tirzepatide Binder Design & Evaluation

Computational campaign to (1) uniformly re-evaluate 90 pre-existing designed
binders against **Tirzepatide** (39-aa GIP/GLP-1 dual agonist), and (2) design +
screen new **peptide binders**, using **Boltz-2.1** (co-folding + binding module)
plus local biophysical analysis. Output is a wet-lab-ready candidate panel.

## TL;DR

- Your 5 yellow-highlighted picks = **exactly the top-5** by an independent
  composite score. Well-calibrated selection.
- Independent Boltz-2 re-validation: **`design_spec_13` is the single best
  existing binder** (cross-validated on 3 signals), ahead of the yellow set.
- Consensus epitope mapped: tirzepatide **C-terminal amphipathic face
  (F22/V23/L26/I27) + polyproline tail**.
- New de novo peptide **NB094** matched the best existing binder; one round of
  **affinity maturation ~doubled** its predicted binding, giving **6 clean,
  liability-free peptide leads** (`TZP-P1…P6`) that hold up under the strict
  **un-forced** test.
- Deliverables: ranked panel, SPPS-ready sequences, humanized **IgG1/IgG4-Fc
  fusion** constructs, predicted complex structures.

## Start here
- **`report/EVALUATION_REPORT.md`** — full narrative, methods, honest caveats.
- **`results/MASTER_panel.md`** / `.csv` — the recommended wet-lab panel.
- **`results/wetlab_constructs.fasta`** — peptides + Fc-fusion sequences to order.
- **`results/figure_maturation.png`** — maturation progression + robustness.

## Layout
```
data/        parsed input spreadsheet, rankings, top-20 payload
design/      peptide generators, developability, Fc-fusion builder, libraries
analysis/    ranking + master-compile scripts
results/     screen outputs, master panel, FASTA, figures, lead CIF structures
report/      EVALUATION_REPORT.md
```

## Boltz-2 jobs run (all via Boltz API remote MCP)
| Job | id | n | purpose |
|---|---|---|---|
| protein_screen | prot_scr_RngI9Lruq9Zj0bgrzIB3 | 20 | re-validate top-20 existing binders |
| protein_screen (epitope) | prot_scr_mxFEVbn9tjG988gzzjCt | 100 | new library round 1 |
| protein_screen (epitope) | prot_scr_lbhTypFeQABzEqrYY7w9 | 35 | NB094 affinity maturation |
| protein_screen (un-forced) | prot_scr_7ZbIFOwLhyJbWwfhUpy8 | 12 | honest binding check, leads |
| protein_screen (un-forced) | prot_scr_Rv6ZsCz6SeOinuGZaDSN | 8 | honest binding check, matured |

Total ≈ 175 Boltz-2 predictions (~$4.4 est). Reproduce metrics via the scripts
in `analysis/` and `design/`.

## Key caveat
Tirzepatide is a small, heavily modified, partly-disordered peptide — a hard
target for structure-based affinity. Scores are **relative triage signals**, not
affinity guarantees. Confirm binding empirically (SPR/BLI); the Fc-avidity
format and further maturation are the levers to improve affinity.
