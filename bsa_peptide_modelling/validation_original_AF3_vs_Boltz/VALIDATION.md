# Validation — original peptide, Boltz-2.1 vs AlphaFold3

To check that the **Boltz-2.1** structure-prediction step used for Peptides A and B is
consistent with the **AlphaFold3** models used in the original order (CCB100725-YW01), the
original peptide `CFAGTPSILMLAGGGS` (Biotin @ res 1, NH2 @ res 16) was re-run through the
*identical* Boltz-2.1 + PRODIGY pipeline and compared with the original AF3-based results.

> Note: my PRODIGY step was already confirmed to reproduce the original report's numbers
> **byte-for-byte** on the supplied AF3 structures. This validation isolates the effect of
> the *structure predictor* (Boltz-2.1 vs AF3) on the same peptide.

## Summary

| Metric | AF3 (original CCB100725-YW01) | Boltz-2.1 (this run) | Agreement |
|---|---|---|---|
| **Mean ΔG** (5 models) | −10.02 kcal/mol | **−10.36 kcal/mol** | within 0.34 kcal/mol |
| **Best ΔG** | −12.4 (Model 5) | **−14.4 (Model 5)** | both Model 5, strongest |
| ΔG range | −8.1 … −12.4 | −8.6 … −14.4 | comparable spread |
| **Max contacts** | 102 (Model 5) | 142 (Model 5) | both densest in Model 5 |
| Interface character | 90–96% hydrophobic | 91–97% hydrophobic | essentially identical |
| Charged–charged contacts | 0 | 0 | identical |

The predicted **mean affinity agrees to ~0.3 kcal/mol**, both predictors independently rank
**Model 5** as the tightest and most-contacted pose, and both describe a strongly
hydrophobic, charge-free interface. See `AF3_vs_Boltz_comparison.png`.

## Binding site

Both methods place the strongest pose on the **same canonical BSA hydrophobic patch**:

- **AF3** best model top residues: ASN390, LEU397, ARG409, LEU386, ALA405, LEU452, MET547, GLN393 …
- **Boltz-2.1** best model top residues: LEU452, TRP213, ARG409, LEU386, ASN390, ALA405, LEU406 …

Shared anchor residues **LEU386 / ASN390 / ALA405 / ARG409 / LEU452** appear in both — i.e.
Boltz-2.1's best pose converges on the AF3 sub-domain-II hydrophobic pocket. Lower-ranked
models in *both* predictors also sample a secondary groove around residues ~205–350, so the
site heterogeneity is a shared property of the peptide/target, not an artefact of one method.

## Interpretation

Boltz-2.1 and AlphaFold3 give **quantitatively consistent** BSA-binding predictions for the
original peptide (mean ΔG within ~0.3 kcal/mol, same top-ranked pose, same hydrophobic
anchor). This supports using the Boltz-2.1 + PRODIGY pipeline for Peptides A and B as a
faithful continuation of the original AF3-based workflow.

The main quantitative difference is that Boltz-2.1's single best model is ~2 kcal/mol stronger
(−14.4 vs −12.4) with more contacts (142 vs 102) — expected model-to-model / seed-to-seed
variation, not a systematic offset (the 5-model means are nearly identical). As always these
are predictions and absolute ΔG/Kd should be read as *relative* guidance pending experiment.

## Files

```
AF3_vs_Boltz_comparison.png     Per-model ΔG and contact-count bar charts (AF3 vs Boltz)
structures/model_{1..5}.cif     Boltz-2.1 complexes (raw)
structures/model_{1..5}.pdb     Boltz-2.1 complexes (PDB; chain A = BSA, B = peptide)
figures/model_{1..5}.png        Ray-traced PyMOL cartoons
prodigy/model_{1..5}/*.csv      PRODIGY output.csv + ic.csv per model
analysis.json                   Parsed rankings / consensus
boltz_confidence.json           Boltz pLDDT / pTM / ipTM per sample
```
(The AF3 reference structures/numbers are the ones supplied in the original DATA_2 bundle.)
