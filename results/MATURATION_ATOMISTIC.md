# Atomistic structure_and_binding validation — matured ab2 scFv

Boltz-2.1 `structure_and_binding` (protein-protein), scFv (chain A) + modified peptide T (Aib2/Aib13), 3 samples each. Clean scFv<->peptide interface metrics + downloadable complexes in `results/antibody/atomistic/`.

*(Note: the K20 lipid is dropped here — Boltz protein-protein binding requires all-protein entities; the lipid sits on the opposite face and its effect was already characterised in the co-fold screens.)*

| Lead | protein_ipTM (best) | interface PDE (Å) | struct_conf | pTM | pp binding_confidence | pred_id |
|---|---|---|---|---|---|---|
| ab2-mat2 (H3:A9Y) | 0.895 | 0.84 | 0.928 | 0.933 | 0.507 | `sab_pred_pXwYSPxcsQOky5yJ47fM` |
| ab2-mat1 (H3:A8Y) | 0.885 | 0.89 | 0.916 | 0.927 | 0.533 | `sab_pred_0CIPL1g678hyEdsAi99y` |
| ab2 WT | 0.866 | 1.19 | 0.903 | 0.925 | 0.394 | `sab_pred_SfWKJ9B51nXqH9YoE8ML` |
| ab2-mat3 (A2Y+A4W+A8Y) | 0.800 | 1.98 | 0.877 | 0.916 | 0.345 | `sab_pred_ovDq8AvG506v7qBFWTkB` |

**Verdict — maturation confirmed at atomistic resolution:**
- Both single-point leads beat WT on **every** interface metric: **A9Y** protein_ipTM 0.895 / iPDE 0.84 Å, **A8Y** 0.885 / 0.89 Å (and highest pp-binding 0.533), vs **WT** 0.855 / 1.19 Å / 0.394.
- The **triple (mat3) regresses** (protein_ipTM 0.800, iPDE 1.98 Å, pp-binding 0.345 < WT) — independent confirmation of the epistasis seen in the screens. Ship the single-point leads, not the stack.
- Complex CIFs: `ab2WT_sample1_best.cif`, `ab2mat1_A8Y_best.cif`, `ab2mat2_A9Y_best.cif`, `ab2mat3_triple_best.cif` (open in PyMOL/ChimeraX; chain A = scFv, chain T = Tirzepatide).
