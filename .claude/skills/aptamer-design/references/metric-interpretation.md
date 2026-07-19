# Interpreting the metrics (and the mandatory caveats)

## What the numbers mean
- **ipTM** (0–1): confidence in the *relative docking pose* of the two chains.
  > 0.8 = model is confident about the geometry. It is **not** a probability of binding
  and **not** an affinity.
- **pTM**: confidence in the overall predicted fold.
- **complex_iplddt**: per-atom confidence at the interface.
- **complex_ipde**: predicted distance error at the interface (lower is better).

## What the numbers do NOT mean
1. **No Kd / IC50.** Boltz's affinity head covers small-molecule ligands only. There is
   no predicted affinity for aptamers. Ranking uses confidence as a *proxy*.
2. **Optimistic for nucleic acids.** Predictors are trained mostly on protein complexes;
   ipTM for a flexible single-stranded ligand overstates certainty.
3. **Single conformation.** Aptamers are ensembles; Mg²⁺/ion-dependent tertiary structure
   and alternative folds are not modeled.
4. **No specificity by default.** A high on-target ipTM says nothing about cross-reactivity
   until you run paralog/off-target and decoy controls.

## Honesty labeling (apply to every claim)
- ✅ = measured or taken from an experimental structure (e.g. footprint from a PDB complex).
- ⚠️ = predicted / hypothesis (all ipTM/pTM/pLDDT values, all designed sequences,
  all "likely binds" statements).

## Red flags when ranking
- Tiny interface (< ~10 residues) with "100%-on-site" → small-denominator artifact.
- On-target ipTM not clearing the scrambled-decoy baseline → likely non-specific stickiness.
- Species mismatch between a structure-derived footprint and the modeled construct →
  residue numbering must be re-mapped by alignment.

## One-paragraph disclaimer to paste into deliverables
> These candidates are a computationally prioritized starting pool for experimental
> validation (SELEX / binding assays), not validated binders. Interface-confidence scores
> (ipTM/pLDDT) are relative ranking signals, not binding affinities. No sequence here has
> been experimentally tested; specificity, affinity, and mechanism require wet-lab confirmation.
