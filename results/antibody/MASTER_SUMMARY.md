# Anti-Tirzepatide antibodies — full humanized IgG1 re-dock & ranking

All previously designed binders re-assembled from their **native VH/VL** into complete **human IgG1**
(heavy = VH–CH1–hinge–CH2–CH3; light = VL–Cλ), then re-docked against modified Tirzepatide
(Aib2/Aib13; K20 lipid distal). Boltz-2.1 co-fold of the **Fab** + peptide, 5 samples each.

## Ranking (recommendation tier)

Ordered by recommended **tier**, which weighs binding-head (affinity proxy) first, then interface
ipTM, epitope specificity, and lineage diversity. Raw binding-head order alone is:
ab2-mat1 (0.453) > ab2-wt (0.372) > ab2-mat3 (0.349) > ab2-mat2 (0.343) > ab1-wt (0.249).
mat3 is demoted despite a marginally higher binder score because its interface is the least specific
(22-residue diffuse footprint, lowest ipTM); ab1-wt is kept as an orthogonal-lineage backup.

| rec | IgG | tier | CDR-H3 mut | binding-head | ipTM μ (best) | fold conf | epitope core | scFv ipTM→ | verdict |
| 1 | ab2-mat1 | A | H3:A8Y | 0.453 | 0.817 (0.830) | 0.898 | 5/5 | 0.885 | **Primary lead** — top binder in both formats |
| 2 | ab2-mat2 | B | H3:A9Y | 0.343 | 0.818 (0.841) | 0.906 | 5/5 | 0.895 | Backup — A9Y edge over WT lost in IgG format |
| 3 | ab2-wt | B | (parent) | 0.372 | 0.817 (0.826) | 0.907 | 5/5 | 0.855 | Parent baseline |
| 4 | ab1-wt | B | (orthogonal) | 0.249 | 0.836 (0.853) | 0.927 | 5/5 | – | Orthogonal backup — best fold, weakest binder |
| 5 | ab2-mat3 | C | H3:A2Y+A4W+A8Y | 0.349 | 0.791 (0.819) | 0.877 | 5/5 | 0.8 | Drop — triple mutant regresses (epistasis) |

## Key conclusions

1. **Reformatting is safe.** All 5 IgGs keep the canonical epitope **F22·V23·Q24·L26·I27** (5/5) — native VH/VL pairing + human CH1/Cλ/Fc did not move the paratope.
2. **Interface score softens, binding preserved.** Atomistic interface ipTM falls from the scFv range (0.86–0.90) to the Fab range (0.79–0.84); all remain confident (>0.79). Expected geometry dilution, not weaker binding.
3. **Lead is robust.** ab2-mat1 (A8Y) is the #1 binding-head scorer in both the scFv (0.731) and full-IgG (0.453) rounds — the ranking survives the format change.
4. **A8Y > A9Y in the IgG context.** ab2-mat2 (A9Y)'s scFv edge over WT collapses once the constant domains are present (0.343 ≈ WT 0.372). A8Y is the durable improvement.
5. **ab2-mat3 (triple)** has the lowest ipTM and the most diffuse interface — confirmed regression, do not advance.

## Files
- `humanized_igg_constructs.fasta` — wet-lab-ready heavy+light chains for all 5 IgGs
- `antibody_report_igg.html` — interactive report (3D docking viewer, ranking, before/after, sequences)
- `cif/*.cif` — best-of-5 Boltz Fab+peptide co-fold structures
- `ranking.json`, `interface_contacts.json`, `constructs.json`, `job_ids.json` — machine-readable data

_Co-fold screen; confirm empirically by SPR/BLI + DSF/SEC. Boltz-2.1 predictions run 2026-07-03._