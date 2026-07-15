# Methodology — artificial allosteric biosensors

Detailed companion to `SKILL.md`. Source: Guo, Smutok, Lee, … Baker & Alexandrov,
*Artificial allosteric protein switches with machine-learning-designed
receptors*, Nature Biotechnology (2026), doi:10.1038/s41587-026-03081-9.

## 1. The conceptual shift

Classical single-protein biosensors need a receptor that undergoes a large
**ligand-induced conformation change**, mechanically coupled to a reporter
(fluorescent protein, enzyme). This paper shows that requirement is unnecessary:
small, rigid **ML/de-novo binders** that show *no global conformation change*
still make efficient switches. The mechanism is **entropic**, not mechanical —
ligand binding lowers the conformational entropy of the chimera, which is
thermodynamically coupled to the reporter's catalytic activity (HDX-MS + 19F-NMR
in the paper; MD was inconclusive). Circular permutation deliberately *raises*
the receptor's conformational entropy (it inserts a long unstructured loop before
the last helix), suppressing reporter activity in the OFF state and giving room
for ligand binding to switch it ON — at the cost of a ~4-fold affinity penalty.

## 2. The design recipe (implemented in `design.py`)

Given a **receptor** R (a binder) and a **reporter** enzyme:

1. **Permutation site** — pick a residue inside a loop of R (`structure.py`
   derives loops from `annotate_sse`). Remove it → new N/C termini.
2. **Circular permutation** — join the *native* N/C termini with a flexible
   Gly/Ser linker: `cpR = R[site+1:] + GSlinker + R[:site]`.
3. **Insertion** — graft cpR into a permissive surface loop of the reporter,
   one glycine each side: `chimera = reporter[:q] + G + cpR + G + reporter[q:]`.
   For TEM-1 β-lactamase the paper uses **position 253** (a surface coil);
   positions **196/197** are used for logic gates.
4. **Focused screen** — because ML binders are small, the permutation library is
   "fewer than ten variants"; each is assayed for ligand-dependent activity.

## 3. Reporters demonstrated

| Reporter | Readout | Notes |
|---|---|---|
| TEM-1 β-lactamase | colorimetric (nitrocefin, 486 nm) | primary; also confers ampicillin resistance in E. coli |
| LuxSit Pro / NanoLuc (de-novo luciferases) | luminescent | shows *fully synthetic* switches (designed receptor + designed reporter) |
| PQQ-glucose dehydrogenase | electrochemical | bioelectrode; quantifies steroid hormones, LOD < 0.5 nM |

## 4. Receptor classes demonstrated

Small molecules (17-OHP, cortisol via NTF2-like binders; colchicine via
anticalin), peptides (BCL-11, C-peptide minibinders), proteins (VirB8, MDM2,
Gal-3 via FN3con). **The recipe is analyte-agnostic** — this is the headline.

## 5. Logic gates and enhancements

- **YES gate** — duplicate the receptor (e.g. TEM-1 loops 41 and 197): >20-fold
  dynamic-range increase in the paper.
- **AND gate** — two *orthogonal* receptors at 41 and 197: activity rises only
  when *both* ligands are present.
- **Auxiliary binding domain** — a second binder fused by a flexible linker
  raises local ligand concentration and partly offsets the CP affinity penalty
  (Gal-3 example: ~7-fold affinity gain).
- **Linker tuning** — lengthening the receptor↔reporter linker with Gly lowers
  dynamic range but raises kcat; Ser rigidification reverts it.

## 6. Performance parameters (wet-lab ground truth)

- **Dynamic range** `DR = kobs(saturating ligand) / kobs(no ligand)` — the linear
  phase of the absorbance trace is fit to `kobs`. **This is the only real
  measure of switchability.**
- **Kd** — fit `kobs` vs [ligand] to the E+S⇌ES explicit binding solution
  (paper eq. 1).
- **Latency, kcat** — activation kinetics and catalytic rate.

## 7. The honesty-labeling contract

Every claim in this skill's outputs is one of:

- **✅ rigorous** — deterministic sequence construction (round-trip verified),
  structure-derived loop/motif annotations, real geometric measurements on a
  model, or numbers returned by the predictor.
- **⚠️ hypothesis** — anything standing in for a measurement: the switch proxy,
  "active site intact ⇒ catalytically ON", "apo→holo confidence shift ⇒
  allostery", or any claim a construct "works". These require the bench.

Negative and partial in-silico results are reported, not hidden. Structure
prediction ranks and de-risks a small library; it does not replace the kobs
titration that defines a biosensor.
