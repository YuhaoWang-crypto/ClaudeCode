# MD + MM/GBSA cross-validation of the Boltz static affinity

**Method.** For 5 docked complexes (binder chain A + fully-modified Tirzepatide
chain T) we ran explicit-solvent MD on Modal GPUs (OpenMM, ff14SB/TIP3P, 4 fs HMR,
0.2 ns NVT + 0.2 ns NPT + **5 ns production**, ~49–79 k atoms), then computed
single-trajectory **chain-split MM/GBSA** (ΔG ≈ ⟨E_complex⟩ − ⟨E_binder⟩ − ⟨E_target⟩,
GBn2 implicit, entropy omitted) and **pose stability** (binder Cα-RMSD after
superposing the target; inter-chain heavy-atom contact retention). Trajectories
are PBC-unwrapped and the target is min-imaged into the binder's periodic image
before scoring. `mdscreen` run `run_7c04089573`; engine code added:
`mdscreen/binding_pp.py` + `complexpp`/`rescorepp` entrypoints.

| Rank by ΔG | System | MM/GBSA ΔG (kcal/mol) | binder RMSD (nm) | contact retention | ⟨contacts⟩ | ΔG per contact |
|---|---|---|---|---|---|---|
| 1 | design_spec_7 | −80.0 ± 7.0 | 0.26 | 1.00 | 325 | −0.246 |
| 2 | design_spec_4 | −74.2 ± 6.1 | 0.33 | 0.93 | 262 | −0.283 |
| 3 | **ab2mat1 (H3:A8Y)** | −47.9 ± 5.3 | 0.35 | 0.98 | 150 | −0.319 |
| 4 | ab2 WT | −38.9 ± 4.0 | 0.42 | 1.27 | 101 | −0.384 |
| 5 | ab2mat2 (H3:A9Y) | −37.5 ± 3.4 | 0.50 | 1.24 | 101 | −0.372 |

## What MD confirms
1. **Every Boltz pose is dynamically stable over 5 ns** — binder RMSD 0.26–0.50 nm
   and contact retention ~0.9–1.3 (contacts even grow slightly as interfaces relax).
   The Boltz-predicted binding modes are physically plausible, not artifacts. Good
   sanity check on the static docking.
2. **Within the ab2 lineage (fair, same-scaffold comparison), maturation helps:**
   **A8Y improves MM/GBSA by ~+9 kcal/mol over WT** (−47.9 vs −38.9) — corroborating
   the Boltz un-forced ranking and the atomistic protein_ipTM/pp-binding where A8Y
   was the strongest single point. A9Y ≈ WT here (its gain was smaller and within noise).

## The important, honest caveat
**MM/GBSA ranks the two FAILED wet-lab designs (spec_7, spec_4) as the "strongest"
binders — above every ab2 lead.** This must not be read as "spec_7 is good." It is a
textbook limitation, and it *reinforces* the earlier evaluation:

- Single-trajectory MM/GBSA (no entropy) is dominated by **interface size**. spec_7
  buries the largest interface (325 contacts, 252-aa scFv) → most negative ΔG. On a
  size-normalised basis (**ΔG per contact**) the ranking inverts and the ab2 family is
  actually the most efficient per contact (WT −0.384, A9Y −0.372, A8Y −0.319 vs
  spec_7 −0.246). Raw ΔG across *different scaffolds* is confounded and unreliable.
- **Both independent structure-based affinity methods — Boltz (ipTM / binding_confidence)
  AND MD MM/GBSA — agree the spec_7/spec_4 interfaces are energetically favourable and
  stable.** Yet they failed experimentally. That is exactly the point of the prior-design
  analysis: their failure is **not** an affinity or pose problem, so no affinity score
  (static or dynamic) can predict it. The cause is **developability — an N-glycosylation
  sequon in CDR-H2 (spec_7 `NGS@55`) + high liability load** — which glycosylates the
  paratope in mammalian expression and abolishes binding.

## Bottom line
- MD/MM-GBSA is a **useful supplement** where it is valid: it confirms pose stability
  for all leads and, within the ab2 scaffold, **independently ranks A8Y > WT**, adding a
  third line of evidence for **ab2-mat1 (H3:A8Y)** as the lead.
- It is **not** a substitute for developability filtering: two orthogonal physics methods
  both miss the spec_7 glycosylation failure. Use MM/GBSA only **within a scaffold** and
  always pair affinity scores with liability screening (CDR glycosites, oxidation,
  deamidation) — the actual determinant of the prior wet-lab outcome.

*Caveats: 5 ns single trajectory, implicit-solvent end-state MM/GBSA, no entropy →
ranking signal within a congeneric series, not absolute affinity. Aib2/Aib13 modelled
as Ala (PDBFixer); K20 lipid absent (protein–protein MD) — both on the face opposite the
mapped epitope.*

---

## 20 ns confirmation MD (final physical check, leads only)
Job `run_633ecd9108` — A8Y & A9Y scFv + modified peptide, 20 ns (4× longer), from the validated
atomistic complexes.

| lead | ΔG 5 ns | **ΔG 20 ns (kcal/mol)** | binder RMSD (nm) | contact retention |
|---|---|---|---|---|
| ab2-mat1 (A8Y) | −47.9 | **−53.1 ± 3.4** | 0.48 | 0.93 |
| ab2-mat2 (A9Y) | −37.5 | **−41.6 ± 4.6** | 0.66 | 1.02 |

Both complexes remain **stably bound over 20 ns** (binder RMSD < 0.7 nm, contact retention ≥ 0.93),
ΔG is consistent with (slightly stronger than) the 5 ns estimate, and **A8Y > A9Y holds**. This is the
physics-based final confirmation that the matured leads form stable complexes with the real modified drug.
