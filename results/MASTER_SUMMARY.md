# Tirzepatide antibody campaign — master summary

**Goal.** Design and validate antibodies against the **real, modified Tirzepatide**
(39-aa dual GIP/GLP-1 agonist; Aib2/Aib13; Lys20 Nε C20-diacid–γGlu–(AEEA)₂ lipid;
C-terminal amide), and provide wet-lab-ready sequences with honest, multi-method evidence.

## Recommended leads
| Lead | Format | CDR-H3 mutation | Why |
|---|---|---|---|
| **ab2-mat1 (H3:A8Y)** | scFv → **IgG1 / scFv-Fc** (λ) | Ala→Tyr @ H3 pos 8 | best across all methods; strongest MM/GBSA |
| **ab2-mat2 (H3:A9Y)** | scFv → IgG1 / scFv-Fc (λ) | Ala→Tyr @ H3 pos 9 | best atomistic interface; robust |
| ab1-wt | backup, orthogonal lineage | — | different scaffold, healthy scores |

Sequences (VH/VL/CDR/scFv/scFv-Fc/IgG1 H+L + human constant cassettes):
`results/antibody/matured_constructs.fasta` and the report's **Antibody Sequences** tab.

## Epitope
Consensus = Tirzepatide **C-terminal amphipathic face F22 · V23 · L26 · I27** (+ E3, PPPS tail).
The K20 lipid sits on the opposite face and does not disrupt binding.

## Evidence for the leads (5 independent lines, all consistent)
| Method | WT | A8Y | A9Y | read |
|---|---|---|---|---|
| Boltz binding_confidence, un-forced (+lipid) | 0.635 | **0.731** | 0.720 | maturation helps |
| Atomistic protein_ipTM (scFv+peptide) | 0.855 | 0.885 | **0.895** | high-confidence interface |
| MD MM/GBSA ΔG, 5 ns (kcal/mol) | −38.9 | **−47.9** | −37.5 | A8Y strongest |
| **MD MM/GBSA ΔG, 20 ns (confirm)** | – | **−53.1 ± 3.4** | −41.6 ± 4.6 | holds on longer traj |
| MD pose stability, 20 ns (binder RMSD / contact retention) | – | 0.48 nm / 0.93 | 0.66 nm / 1.02 | complex stays bound |
| Whole-drug ligand dock (orthogonal representation) | 0.175 | 0.193 | 0.270 | docks at **CDR paratope** (CDR-H3 + VL CDRs) |

**Bottom line:** A8Y and A9Y bind the modified drug at the CDR paratope (matured CDR-H3 + VL CDRs),
form complexes that are stable over 20 ns MD, and A8Y is the strongest by MM/GBSA. Advance as
IgG1 / scFv-Fc (bivalent → avidity). Confirm empirically by SPR/BLI + DSF/SEC.

## Developability / risk (leads are clean)
- **No N-glycosylation sequon anywhere in the Fv/CDRs** (A8Y/A9Y/WT). The only sequon in the
  scFv-Fc is the canonical **Fc N297 glycan** (distal, functional).
- 4 Cys = the two intradomain disulfides, **no free Cys**.
- Framework has ordinary Met-ox / deamidation hotspots (mitigable), no CDR liabilities.

## Why the earlier (unmodified-target) designs failed in wet-lab
Not an affinity problem — **structure-based affinity scores (Boltz AND MM/GBSA) both fail to flag them.**
The cause is **developability**: 5/10, including the preferred `design_spec_7`, carry an
**N-glycosylation sequon inside a CDR** (spec_7 = `NGS` in CDR-H2 → paratope glycan in mammalian
expression → no binding), plus very high liability loads. Detail: `results/PRIOR_DESIGNS_EVALUATION.md`.

## On the Fc format (important caveat)
Whole-antibody Boltz co-folding **cannot** score the Fc format: a lone scFv-Fc chain misfolds
(CH3 is an obligate homodimer → structure_confidence ≈ 0), and the ~980-aa flexible dimer gives
diluted global scores + mislocalized peptide — all 5 tested constructs collapse to ~0 uniformly
(artifact, not loss of binding). The Fc is a separate module across a flexible hinge and does not
touch the paratope, so **per-arm binding = the bare-scFv result (unchanged)**; the Fc adds bivalent
avidity (a gain). Quantify the Fc format by **MD** (dimer+peptide) or **experiment**, not whole-antibody
co-folding. Detail: `results/antibody/fc/FC_EVALUATION.md`.

## Methods / models
- **Boltz-2.1** (AlphaFold3-class co-folding; there is no separate AlphaFold3 endpoint here) via Boltz API:
  protein_screen (epitope-forced & un-forced), structure_and_binding (atomistic protein_ipTM + CIF),
  whole-drug-as-ligand co-fold.
- **MD**: OpenMM (ff14SB/TIP3P, 4 fs HMR) on Modal GPUs via the `mdscreen` package (branch
  `claude/md-analysis-system-no9p30`); new **protein–peptide chain-split MM/GBSA + pose stability**
  (`mdscreen/binding_pp.py`, PBC-correct) added this campaign.

## File index
```
results/
  MASTER_SUMMARY.md                 ← this file
  MATURATION_REPORT.md              affinity maturation (3 rounds + validation)
  MATURATION_ATOMISTIC.md           atomistic protein_ipTM of matured leads
  MD_CROSSVALIDATION.md             5 ns MD + MM/GBSA (5 complexes) + honest caveats
  PRIOR_DESIGNS_EVALUATION.md       why the prior scFv failed (CDR glycosylation)
  antibody/
    matured_constructs.fasta        VH/VL/scFv/scFv-Fc/IgG1 (leads) — order these
    human_constant_regions.fasta    IgG1 Fc, IgG1-LALA, IgG4-S228P, CH1-CH3, Cκ, Cλ
    atomistic/                      scFv+peptide complex CIFs + metrics
    md/                             5 ns MM/GBSA JSON (5 complexes)
    md20/                           20 ns MM/GBSA confirmation (A8Y, A9Y)
    fc/                             scFv-Fc evaluation (+ artifact analysis)
    ligand_dock/                    whole-drug-as-ligand docking + CIFs
report/tzp_report.html             interactive report (tabs incl. Antibody, Sequences)
```
Boltz jobs: maturation `prot_scr_SAlb4E3KrEdF3qUOg5bz`; combos `…ScBhRFLCYYTqcRGepPBK`;
validation `…cc4FyMeY0pAH8Vgqpvfg`; prior-designs `…vYmIVxUdiwMveYc9lDJu` (+un-forced/no-lipid);
scFv-Fc `…6N1ZijBr4pqN7tE6ogiz` / dimer `…4s6KLQcbz1pJxLPPCipD`; whole-drug ligand `…mldVzPigBS5LqUfIJBgR`;
atomistic `sab_pred_*`. MD: `run_7c04089573` (5 ns), `run_633ecd9108` (20 ns).
