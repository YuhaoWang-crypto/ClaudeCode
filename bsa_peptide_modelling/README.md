# Modelling of two 16-mer biotinylated peptides binding to BSA

Replication of project **CCB100725-YW01** ("Modelling of binding of a 16-amino-acid
biotinylated peptide to BSA") for two new peptide sequences, using the identical
computational protocol.

| | Sequence | Change vs. original | N-cap | C-cap |
|---|---|---|---|---|
| **Original** (CCB100725-YW01) | `CFAGTPSILMLAGGGS` | — | Biotin (BTN) @ res 1 | Amide (NH2) @ res 16 |
| **Peptide A** | `CFAGTPSILKKNGGGS` | M10K, L11K, A12N (+2 charge, central) | Biotin (BTN) @ res 1 | Amide (NH2) @ res 16 |
| **Peptide B** | `KKAGTPSILMLAGGGS` | C1K, F2K (+2 charge, N-terminal) | Biotin (BTN) @ res 1 | Amide (NH2) @ res 16 |

Target: **Bovine serum albumin (BSA)**, 583 aa, 66.5 kDa (mature sequence, chain A).

## Method (identical to CCB100725-YW01)

1. **Structure** — Boltz-2.1 co-folding of the BSA + peptide complex, 5 sampled models
   per peptide. Biotin is represented by the CCD component `BTN` substituted at peptide
   residue 1 and C-terminal amidation by CCD `NH2` at residue 16 (exactly as in the
   original order's input). Automatic MSA.
2. **Binding energy & interface** — PRODIGY (contact-based ΔG / Kd predictor, 5.5 Å
   contact cutoff, 25 °C). The BTN and NH2 caps are stripped before PRODIGY (which
   requires standard residues), so the interface analysis reflects peptide residues 2–15 —
   the same behaviour noted in the original report.
3. **Interface analysis** — per-model ranking of BSA residues by number of peptide
   contacts, per-model ranking of peptide residues, cross-model consensus and hotspot bins.
4. **Figures** — ray-traced PyMOL cartoons: BSA (chain A) in gray cartoon (15%
   transparency), peptide (chain B) in orange cartoon + sticks, and BSA residues within
   5 Å of the peptide shown as marine-blue sticks (matching the CCB_BSA_v2 style).

The pipeline was validated two ways: (1) re-running PRODIGY on the original CCB100725-YW01
AF3 structures reproduces that report's ΔG, Kd, contact counts and interface tables exactly;
(2) re-running the *original* peptide through the full Boltz-2.1 + PRODIGY pipeline agrees
with the AF3-based numbers to within ~0.3 kcal/mol on the mean ΔG, with both predictors
ranking the same pose strongest on the same hydrophobic BSA patch — see
`validation_original_AF3_vs_Boltz/VALIDATION.md`.

## Key results

| Peptide | Best model ΔG (kcal/mol) | Best Kd (M) | ΔG range (5 models) | Max contacts |
|---|---|---|---|---|
| Peptide A `CFAGTPSILKKNGGGS` | **−11.1** (Model 4) | 6.8×10⁻⁹ | −7.3 … −11.1 | 92 |
| Peptide B `KKAGTPSILMLAGGGS` | **−9.8** (Model 4) | 6.2×10⁻⁸ | −8.1 … −9.8 | 81 |
| Original `CFAGTPSILMLAGGGS` | −12.4 (Model 5) | 8.6×10⁻¹⁰ | −8.1 … −12.4 | 102 |

- Both new peptides are predicted to retain **sub-micromolar** BSA binding.
- **Peptide A**: its best-scoring model re-docks into the *same* canonical BSA pocket as the
  original peptide (390–410 / 540–548: GLN393, LEU397, ARG409, GLU540, MET547…), while the
  other models sample an alternative groove around residues 208–353.
- **Peptide B**: all five models converge on a distinct, more N-terminal/charged BSA region
  (~114–189: LYS114/LYS116, GLU125/GLU140, ARG144/HIS145, ARG185), recruiting charged BSA
  residues to complement its two new N-terminal lysines. Most reproducible of the three.

> All ΔG/Kd values are Boltz-2.1 + PRODIGY **predictions** and should be confirmed
> experimentally (e.g. BLI/SPR with the biotinylated peptides on streptavidin sensors).

## Deliverables

```
Project_Report_BSA_two_peptides.pdf   Full formatted report (12 pp), CCB format
comparison.json                       Cross-peptide summary + interpretation
<peptide>/
  structures/model_{1..5}.cif         Raw Boltz-2.1 complexes (BTN/NH2 present)
  structures/model_{1..5}.pdb         Same, PDB format (chain A = BSA, B = peptide)
  prodigy/model_{1..5}/output.csv     PRODIGY ΔG/Kd + contact-type breakdown
  prodigy/model_{1..5}/ic.csv         PRODIGY interface contact list
  analysis.json                       Parsed rankings / consensus / hotspots
  boltz_confidence.json               Boltz pLDDT / pTM / ipTM per sample
  figures/model_{1..5}.png            Ray-traced PyMOL cartoon of each model
  figure.png                          Best-model complex figure (used in report)
scripts/                              cif2pdb.py, analyze.py, render_pymol.py, generate_report.py
```

### Reproducing

```bash
# structures come from Boltz-2.1 (structure_and_binding, no binding block, num_samples=5)
python scripts/cif2pdb.py <sample>.cif model_n.pdb          # CIF -> PDB
python scripts/analyze.py <label> <dir_with_model_1..5.pdb> <outdir>   # strip caps, PRODIGY, rank
python scripts/render.py  model_n.pdb figure.png "title" <iface_resnums>
python scripts/generate_report.py Project_Report_BSA_two_peptides.pdf
```

Requires `prodigy-prot` (+ `freesasa`), `biopython`, `reportlab`, `matplotlib`.

## Notes / caveats

- **Biotin representation.** As in the original protocol, biotin occupies peptide
  position 1 (the residue there is replaced by the BTN cap). For Peptide B, position 1 is a
  lysine, so the modeled construct is `Biotin–K·AGTPSILMLAGGG–NH2` — i.e. the *first* of the
  two N-terminal lysines is represented by the biotin group. If the intended construct is
  instead "both lysines present **plus** an added biotin cap", the model should be rerun with
  a 17-position peptide.
- The original CCB100725-YW01 job included a Boltz protein–protein *binding-affinity* pass;
  that module timed out for this larger BSA complex and is not needed here (PRODIGY provides
  the affinity), so it was omitted. The structure prediction is otherwise identical.
