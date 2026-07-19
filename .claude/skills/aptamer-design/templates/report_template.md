# Aptamer design against <TARGET> (<GENE>)

**Target:** <organism> <protein> — UniProt **<ACC>**
**Use case:** <diagnostic probe | therapeutic antagonist>
**Design date:** <YYYY-MM-DD>
**Deliverable:** structure-informed, prioritized DNA/RNA aptamer shortlist for experimental validation.

> **Caveat — read first.** Computationally prioritized starting pool for SELEX / binding
> assays, **not** validated binders. ipTM/pLDDT are relative confidence signals, **not**
> affinities. Nothing here has been experimentally tested.

## 1. Target rationale
- Domain architecture, folded functional module targeted, construct crop (residues).
- Ligand-binding footprint (PDB <id>, residues within 4.5 Å) — note species/numbering mapping.
- Electrostatics / basic or heparin groove (polyanion docking site) — and its specificity risk.
- Prior aptamers in the literature? (PubMed / bioRxiv check.)

## 2. Methods
- Candidate generation: scaffolds sampled, length/GC constraints, ViennaRNA fold filter.
- Co-folding: Boltz-2.1, target crop + aptamer, auto-MSA, N samples, scrambled decoy control.
- Ranking: composite (confidence + fold + chemistry + specificity[+on-target]); see rank_candidates.py.

## 3. Ranked candidates
| # | ID | seq (5'→3') | chem | scaffold | ipTM | score | grade | notes |
|---|----|-------------|------|----------|------|-------|-------|-------|
|   |    |             |      |          |      |       |       |       |

Decoy baseline ipTM: <value>. Candidates must clear this by a clear margin.

## 4. Specificity counter-screen
On-target vs paralog ipTM for finalists.

## 5. Files
- ranked CSV, predicted complex .cif structures, secondary-structure diagrams, this report.

## 6. Caveats & limitations
(Paste the disclaimer from references/metric-interpretation.md.)

## 7. Experimental validation plan
1. Synthesize leads (DNA oligos / 2'-F,2'-OMe RNA) + scramble & random-pool controls.
2. Primary binding: BLI / SPR / MST vs recombinant target ectodomain; rank by apparent Kd.
3. (Antagonist) ligand-competition assay to confirm mechanism.
4. Doped-SELEX affinity maturation around validated leads + counter-selection on paralogs.
5. Cell-based functional readout.
6. Truncation + 3'-inverted-dT / PEGylation; re-confirm affinity.
