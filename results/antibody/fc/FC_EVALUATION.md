# Does adding the humanized IgG1 Fc (scFv-Fc) change binding to the modified Tirzepatide?

**Short answer: No — appending the Fc does not perturb the antigen-binding site.**
The paratope lives on the scFv/Fv; the Fc is a separate module joined by a
flexible hinge and never contacts the CDRs. The strong, valid binding evidence
is the **bare-scFv** result; the low *whole-molecule* scFv-Fc Boltz scores below
are **modeling artifacts**, which we demonstrate, not loss of binding. The
bivalent Fc format is expected to *increase* functional affinity via avidity.

## What we ran (vs fully-modified Tirzepatide: peptide + Aib2/Aib13 + K20 lipid)
1. Head-to-head co-fold screen, **bare scFv vs scFv-Fc monomer**, forced (+lipid)
   — `prot_scr_6N1ZijBr4pqN7tE6ogiz`; un-forced — `prot_scr_0xlZiKaf7sEAiRmMWove`.
2. Atomistic `structure_and_binding`: scFv-Fc **monomer** (A8Y `sab_pred_VFekqgHdMuLI8vSJl7Q8`)
   and the correct **homodimer** (A8Y `sab_pred_vFGwJKMJB6lVCjknh4m8`, WT `sab_pred_a8g96gkPLyt0IKtmW8Lv`).
3. Direct **interface analysis** of the dimer structure.

## Results

### Bare scFv (the paratope) — binds well, unchanged by maturation
| construct | forced bind (+lipid) | un-forced bind | atomistic protein_ipTM |
|---|---|---|---|
| WT scFv | 0.638 | 0.635 | 0.855 |
| A8Y scFv | 0.664 | 0.731 | 0.885 |
| A9Y scFv | 0.621 | 0.720 | 0.895 |

### scFv-Fc — global scores collapse, but for modeling reasons
| construct | binding_conf | structure_conf | note |
|---|---|---|---|
| WT scFv-Fc (monomer) | 0.137 | **0.12** | lone Fc misfolds (CH3 is an obligate homodimer) |
| A8Y scFv-Fc (monomer) | 0.025 | **0.008** | same |
| A9Y scFv-Fc (monomer) | 0.010 | **0.002** | same |
| A8Y scFv-Fc (homodimer, atomistic) | ~0 | ptm 0.49 / protein_ipTM 0.48 | global metrics diluted over a ~990-aa flexible bivalent assembly |

### Interface analysis of the A8Y scFv-Fc **dimer** (the tell)
The peptide's 262 heavy-atom contacts fall on chain-A residues **280–473 — i.e.
the Fc region — not the CDRs** (VH-CDR3 ≈ 97–108). Boltz docked the peptide onto
an Fc surface, not the paratope: a mislocalization artifact of the low-confidence
global fold. This confirms the whole-scFv-Fc model is unreliable for judging the
paratope, and its low scores say nothing about real binding.

## Why the scFv-Fc scores are artifacts (not biology)
- **Fc monomer misfolds.** The CH3 domain only folds as a homodimer; a single
  scFv-Fc chain has no partner → garbage fold (structure_confidence ≈ 0) that
  drags the co-folded peptide score to zero.
- **The dimer is too large/flexible for global metrics.** ipTM/ptm/binding_confidence
  are whole-complex averages; over a ~990-residue two-armed molecule with flexible
  hinges (peptide occupies one small arm) they are structurally diluted and the
  peptide is often mis-placed. These are limitations of co-folding a whole antibody,
  independent of affinity.
- **Architecture.** Fv paratope and Fc are separate domains ~40+ Å apart across the
  hinge; the Fc cannot reach the CDRs. Appending it does not change the binding site.

## Positives the Fc adds
- **Avidity.** scFv-Fc / IgG1 is **bivalent** → 2 paratopes → large gains in *apparent*
  affinity for a target, the opposite of a binding penalty. (Static single-site
  co-folding cannot capture this.)
- **Canonical Fc N297 glycan** (the only N-glyc sequon in the construct, `NST`@~325)
  is normal, functional (FcRn/effector), and **distal to the paratope** — not a
  binding liability. No glycosylation sequon anywhere in the Fv/CDRs.

## Conclusion & recommendation
Adding the humanized IgG1 Fc does **not** reduce binding to the real (modified)
Tirzepatide: the validated paratope (bare scFv — 3 independent methods: co-fold,
atomistic protein_ipTM 0.85–0.90, MD-stable) is unchanged, and bivalency should
help. To *quantify* the Fc-format binding, the right tools are **experiment**
(express scFv-Fc / IgG1, measure by SPR/BLI — avidity included) or MD of the
**Fv + peptide** (done); whole-antibody static co-folding is not a valid readout.
Recommended lead format: **A8Y (or A9Y) as scFv-Fc or full IgG1** (λ light).
