# Boltz-2 job log — IL-6 miniprotein campaign (2026-08-21)

Target for all design runs: `https://files.rcsb.org/download/1P9M.cif`, chain B (IL-6),
`crop_residues: "all"`, binder `custom_protein`, designed length `55..85`,
rules: `excluded_amino_acids: ["C"]`, `max_hydrophobic_fraction: 0.42`,
`excluded_sequence_motifs: ["NXS","NXT"]`.

## Design runs (120 designs each, $0.05/design)

| Job ID | Epitope | `epitope_residues` (0-indexed into the 186-residue entity) | Note |
|---|---|---|---|
| `prot_des_FHH27c3CTUeVhCW5g6gC` | site I | 11,14,32,39,44,47,51,52,53,56,150,153,156,157,158,160,161 | **wrong** — indexed into resolved residues only |
| `prot_des_ByQiu8waaBaPAw4DWR2a` | site II | 0,5,8,9,11,12,15,88,89,91,92,95,96,99,102,103,106 | **wrong** — same bug |
| `prot_des_TIvb388MQVUZWeqrMcuC` | site I | 32,35,56,63,68,71,75,76,77,80,174,177,180,181,182,184,185 | corrected (`= mature position + 2`) |
| `prot_des_whXVXpevA2FLJYXGGwhR` | site II | 21,26,29,30,32,33,36,112,113,115,116,119,120,123,126,127,130 | corrected |

## Independent confirmation co-folds

`boltz_start_structure_and_binding`, model `boltz-2.1`, `num_samples: 5`,
target = mature IL-6 (183 aa) with automatic MSA, binder with `msa: {type: empty}`,
`binding: protein_protein_binding`. No design-time template is provided.

| Prediction ID | Design | Target | ipTM (5 samples) | binding score |
|---|---|---|---|---|
| `sab_pred_B6SkhqhNC0YFLCaXshm2` | IL6-S2-01 | IL-6 | 0.950, 0.955, 0.948, 0.951, 0.950 | 0.765 |
| `sab_pred_y80UgAdJWjjGvqSdM8s5` | IL6-S1-01 | IL-6 | 0.895, 0.891, 0.897, 0.870, 0.885 | 0.211 |
| `sab_pred_XLs3z1kNIg8LsrR4ADYe` | IL6-S1-02 | IL-6 | 0.854, 0.769, 0.852, 0.891, 0.862 | 0.0027 |
| `sab_pred_mmP1L32OqvFwehCVbg1i` | IL6-S1-03 | IL-6 | (run, not tabulated) | — |
| `sab_pred_yDcfKK6cdv9EtMwZ2haS` | IL6-S2-02 | IL-6 | (run, not tabulated) | — |
| `sab_pred_nEgKWxS8AFJHQbplfjkk` | site II design from the buggy run | IL-6 | 0.949, 0.943, 0.929, 0.946, 0.949 | 0.00029 |
| `sab_pred_c6Xiq1FGrAVqCPUkE9pA` | IL6-S1-01 | **CLEC12A ectodomain (off-target control)** | 0.745, 0.813, 0.755, 0.765, 0.847 | 0.00067 |

CLEC12A construct: UniProt Q5QGZ9 residues 65–265 (C-type lectin ectodomain).

## Reference interfaces measured from 1P9M (same code path as the designs)

| Interface | buried surface | interface residues | atom pairs < 4 Å | H-bonds | salt bridges |
|---|---|---|---|---|---|
| IL-6 · IL-6Rα (site I) | 1454 Å² | 17 / 19 | 71 | 3 | 3 |
| IL-6 · gp130 (site II) | 1352 Å² | 17 / 19 | 42 | 5 | 0 |
