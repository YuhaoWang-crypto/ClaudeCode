# Boltz-2.1 co-folding workflow (aptamer × protein)

Uses the `Boltz_API` MCP server. Estimates are free; runs are paid (~$0.025/complex).

## 1. Estimate (free) — always first
`boltz_estimate_structure_and_binding` with the same `input` you plan to run. Report the
cost to the user and get approval before any paid submission.

## 2. Input shape
Two entities: the cropped protein target (chain A) and the aptamer (chain B).
```json
{
  "entities": [
    { "type": "protein", "chain_ids": ["A"], "value": "<TARGET_CROP_AA>" },
    { "type": "dna",     "chain_ids": ["B"], "value": "<APTAMER_DNA>" }
  ]
}
```
- Use `"type": "rna"` (A/C/G/U) for RNA aptamers, `"type": "dna"` (A/C/G/T) for DNA.
- **Omit `msa`** on the protein → automatic MSA (needed for correct folding of a
  disulfide-rich domain). Do not add a `binding` block — that head is for
  small-molecule/protein binders and returns 0 for nucleic acids.
- `num_samples` default 1 is fine for screening; raise to 3–5 to reduce ipTM noise on finalists.

## 3. Submit (paid)
`boltz_start_structure_and_binding` with a stable, unique `idempotency_key` per
(candidate, target) pair — e.g. `"<target>-<aptamerID>-v1"`. Submit all candidates in
parallel (independent calls). Runs typically finish in ~0.5–3 min each.

## 4. Poll
`boltz_get_structure_and_binding_prediction(prediction_id)` until `status == "succeeded"`.
Prefer waiting via a single background timer over rapid re-polling. Structure download
URLs are signed and expire (~30 min) — fetch/save promptly if you need the .cif files.

## 5. Metrics to record (per best_sample.metrics)
| Field | Meaning | Use |
|---|---|---|
| `iptm` | interface predicted TM | **primary** interface-confidence rank signal |
| `ptm` | global predicted TM | overall fold plausibility |
| `complex_iplddt` | interface pLDDT | local interface confidence |
| `complex_ipde` | interface PDE (lower better) | pose distance error |
| `structure_confidence` | composite | secondary |
| `ligand_iptm` / `protein_iptm` | 0 for nucleic acids — ignore | — |

## 6. Controls & comparability
- **Always run a scrambled-sequence decoy** to see the ipTM baseline for a non-designed
  ligand. A candidate is only interesting if it clears the decoy by a clear margin.
- ipTM is only comparable **within the same target crop + same MSA**. If you mix
  campaigns, note it; treat differences < ~0.03 as ties.

## 7. Specificity screen
Re-run finalists against paralog crops (same length/region) and compare on-target vs
off-target ipTM. Record `offtarget_iptm` for the ranker.
