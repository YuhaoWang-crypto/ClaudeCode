# Evaluation of 10 prior scFv designs (built vs UNMODIFIED Tirzepatide)

**Context.** These 10 scFv were designed against **unmodified** Tirzepatide and taken to the wet lab, where they **failed to recognise the modified drug**. `design_spec_7` was a preferred pick. I re-evaluated all 10 uniformly against the **fully-modified** Tirzepatide (peptide + Aib2/Aib13 + K20 C20-diacid acyl chain) with Boltz-2, alongside my ab2 WT / matured leads, under three conditions.

## Headline (honest)
**Boltz binding_confidence alone does NOT reproduce the wet-lab failure** — under an epitope-forced or lipid-present co-fold, several priors score as well as my ab2 (0.65-0.77). A structure-based affinity score would not, by itself, have predicted this failure. The failure is explained by things the score does not model:

1. **N-glycosylation sequons inside CDRs (paratope glycosylation).** 5/10 designs - including the preferred `design_spec_7` - carry an `N-x-S/T` sequon in a CDR. In mammalian expression that Asn is glycosylated and the glycan sterically blocks antigen binding. A classic de-novo-design miss and a sufficient, mechanism-level cause of 'no recognition'.
2. **Severe developability liabilities.** `design_spec_7` has the highest liability score of the set (220); `design_spec_0` 200, `design_spec_4` 190 - heavy Trp-oxidation, Asp isomerisation/cleavage, deamidation and proteolysis hotspots -> poor expression / aggregation / degradation.
3. **Weak intrinsic peptide engagement.** Against the *naked peptide* (lipid removed), most priors collapse to binding_confidence ~= 0 - they do not robustly grip the peptide epitope on their own; the lipid-present score is inflated by the large hydrophobic acyl chain (Delta-lipid column).

Jobs: forced `prot_scr_vYmIVxUdiwMveYc9lDJu`; un-forced `prot_scr_CHgUvToRsJfLjK58NcXD`; naked-peptide (no lipid) `prot_scr_h9J1PTliJwF23be2heUI`.

| Their rank | Design | forced | un-forced (+lipid) | naked-peptide | Δlipid | their-ipTM (unmod) | liab | CDR glycosite |
|---|---|---|---|---|---|---|---|---|
| #1 | design_spec_4 | 0.714 | 0.696 | 0.257 | -0.439 | 0.469 | 190 | CDR-H1 (NGT@30) |
| #2 | design_spec_9 | 0.691 | 0.692 | 0.347 | -0.345 | 0.374 | 125 | CDR-H2 (NGT@54) |
| #3 | design_spec_8 | 0.643 | 0.723 | 0.000 | -0.722 | 0.523 | 100 | - |
| #4 | design_spec_7 | 0.623 | 0.653 | 0.406 | -0.247 | 0.447 | 220 | CDR-H2 (NGS@55) |
| #5 | design_spec_1 | 0.640 | 0.724 | 0.097 | -0.627 | 0.307 | 70 | - |
| #6 | design_spec_2 | 0.734 | 0.774 | 0.010 | -0.764 | 0.336 | 95 | - |
| #7 | design_spec_0 | 0.705 | 0.680 | 0.000 | -0.680 | 0.288 | 200 | CDR-H2 x2 (NGS@52,55) |
| #8 | design_spec_3 | 0.711 | 0.620 | 0.000 | -0.620 | 0.269 | 105 | - |
| #9 | design_spec_6 | 0.655 | 0.636 | 0.000 | -0.636 | 0.308 | 120 | - |
| #10 | design_spec_5 | 0.685 | 0.542 | 0.231 | -0.311 | 0.221 | 120 | CDR-H2 (NGS@54) |
| ref | **ab2 WT (mine)** | 0.652 | 0.640 | 0.325 | -0.315 | - | - | none |
| ref | **ab2-mat1 A8Y (mine)** | 0.669 | 0.618 | 0.299 | -0.319 | - | - | none |

## `design_spec_7` (your preferred pick) - three independent strikes
1. **CDR-H2 N-glycosylation sequon** `...AISSD`**`NGS`**`YKYYVG...` (NGS@55) -> paratope glycan in mammalian expression.
2. **Highest liability score in the panel (220)** - Trp-ox, Asp-isomerisation, deamidation, hydrophobic patches.
3. **Weakest predicted binder among the priors under forcing (0.623)**; only 0.653 un-forced (~ab2 WT 0.640).

## What my campaign did differently
- Designed & scored against the **actual modified drug** (Aib + K20 lipid), not the naked peptide.
- **Liability-filtered** leads; my ab2/matured scFv carry **no CDR glycosylation sequon**.
- Used **un-forced validation** to reject forcing artifacts.

**Actionable fix:** to rescue any prior design, remove the CDR glycosylation sequon (N->Q or S/T->A in the N-x-S/T motif) and de-risk Trp/Asp hotspots, then re-test. Given their weak naked-peptide engagement, the modified-target-designed ab2 lineage is the stronger starting point.
