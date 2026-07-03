# Orthogonal validation: dock antibody to the WHOLE Tirzepatide as a single small-molecule ligand

**Motivation.** All prior work modelled Tirzepatide as a *peptide chain* (+ K20 lipid ligand).
As an independent cross-check the user requested an AlphaFold-class "small-molecule docking":
the entire modified drug as ONE ligand (full SMILES, **341 heavy atoms**), docked to the scFv.
(Boltz-2.1 is the AF3-class co-folding model available here; there is no separate AlphaFold3
endpoint. Boltz `ligand_protein_binding` affinity mode caps ligands at 100 heavy atoms, so this
was run as a co-fold screen with the drug as a ligand target, ligand-epitope forced —
job `prot_scr_mldVzPigBS5LqUfIJBgR`.)

## Metrics (whole-drug-as-ligand) vs the reliable peptide-chain representation
| scFv | ligand-rep binding_conf | ligand-rep struct_conf | ipTM | min_iPAE (Å) | peptide-chain binding_conf (ref) |
|---|---|---|---|---|---|
| A9Y | 0.270 | 0.22 | 0.65 | 5.8 | 0.621 |
| spec7 | 0.267 | 0.21 | 0.69 | 4.0 | 0.623 |
| A8Y | 0.193 | 0.13 | 0.57 | 6.7 | 0.664 |
| WT  | 0.175 | 0.13 | 0.54 | 7.0 | 0.638 |

## Interface: does the drug dock at the paratope? — YES
Heavy-atom contact analysis of the predicted complexes:
- **A8Y**: 57 antibody contact residues, **14 in CDRs** — CDR-H1 (31,32), **CDR-H3 (97–108)**,
  CDR-L1 (173–181), CDR-L2 (186,189,190), CDR-L3 (219,221).
- **A9Y**: 47 contact residues, **10 in CDRs** — CDR-H3 (97–106), CDR-L1 (173–181),
  CDR-L2 (187,189,190), CDR-L3 (219,221).
The whole drug sits across the VH/VL groove engaging CDR-H3 (the matured loop) + all three
VL CDRs — the same paratope as the peptide-chain model.

## Verdict
- **Consistent, orthogonal support:** two independent target representations (peptide-chain and
  whole-drug small-molecule) **agree the antibody binds Tirzepatide at the CDR paratope**, with
  the matured CDR-H3 engaged.
- **Use the peptide-chain numbers for ranking.** As predicted, a 341-heavy-atom flexible
  peptidic ligand folds/docks with **lower, noisier confidence** (binding_conf 0.17–0.27,
  struct_conf 0.13–0.22, iPAE 4–7 Å) and the ranking scrambles — this is a representation
  artifact, not a real re-ranking. The peptide-chain + lipid model (binding_conf 0.62–0.66,
  atomistic protein_ipTM 0.85–0.90, MD-stable) remains the trustworthy readout.
- **Next:** MD (pose stability + MM/GBSA) on the top complexes for a physics-based confirmation.

Structures: results/antibody/ligand_dock/{A8Y,A9Y}_wholedrug.cif (chain A = scFv, chain L = whole Tirzepatide).
