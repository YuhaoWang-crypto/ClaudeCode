# Boltz-2.1 co-fold workflow (exact MCP sequence)

The structure/affinity engine for this skill. The Boltz MCP tools are invoked by the model (not from
a script); the bundled scripts handle the deterministic post-processing (CIF download, contacts,
liabilities).

## Tools
- `boltz_get_account_context` — call once; confirms the org. Omit `organization_id` afterward if one
  org is already selected.
- `boltz_estimate_structure_and_binding` — validate input + get $ cost **before** paying.
- `boltz_start_structure_and_binding` — start a paid job; needs a stable `idempotency_key`.
- `boltz_get_structure_and_binding_prediction` — poll by `prediction_id`; returns metrics + presigned
  CIF/archive URLs (URLs expire ~30 min — download promptly).

## Input shape (antibody Fab + peptide antigen)
```json
{
  "entities": [
    {"type":"protein","chain_ids":["H"],"value":"<Fd = VH-CH1>","msa":{"type":"empty"}},
    {"type":"protein","chain_ids":["L"],"value":"<light = VL-CL>","msa":{"type":"empty"}},
    {"type":"protein","chain_ids":["P"],"value":"<antigen>","msa":{"type":"empty"},
     "modifications":[{"residue_index":1,"type":"ccd","value":"AIB"},
                      {"residue_index":12,"type":"ccd","value":"AIB"}]}
  ],
  "binding": {"type":"protein_protein_binding","binder_chain_ids":["H","L"]},
  "num_samples": 5
}
```
Conventions that matter:
- **Antibodies use single-sequence mode** (`"msa":{"type":"empty"}`). Auto-MSA and custom/empty MSA
  cannot be mixed — if any chain is empty, set all protein chains empty.
- `residue_index` in `modifications` is **0-based**; non-standard residues go in as CCD codes
  (`AIB`, `SEP`, `MSE`, …). For a bulky PTM (lipid/glycan), add it as a ligand entity
  (`ligand_smiles`/`ligand_ccd`) plus a `bonds` constraint to the anchor atom.
- Co-fold the **Fab**, not the whole IgG (see methodology § 4).
- `num_samples` 5 is a good robustness/cost balance; cost ≈ $0.033/sample (~$0.17 for a ~470-residue
  Fab+peptide complex, 5 samples). Estimate first.

## Metrics returned (per sample + best)
- `protein_iptm` — atomistic antibody↔antigen interface confidence (0–1; >0.8 confident). **Average
  over samples**; don't cherry-pick the best.
- `binding_confidence` (in `binding_metrics`) — Boltz binder-head affinity proxy (higher = stronger).
- `structure_confidence`, `complex_plddt`, `complex_pde` — fold quality.

## Loop
1. `estimate` → check cost.
2. `start` with a descriptive `idempotency_key` (e.g. `tirz-fab-ab2mat1-v1`); reuse it verbatim on
   retries so you never double-pay.
3. Poll `get_...prediction` until `status:"succeeded"` (a ~470-res complex finishes in ~2–3 min).
4. `curl` the `best_sample.structure.url` (and/or `all_sample_results[*]`) to a local `.cif`
   **immediately** (URLs expire). If a download returns a tiny XML `SignatureDoesNotMatch`, the URL
   was truncated/mis-copied — re-fetch a fresh URL by polling again.
5. Post-process with `scripts/contacts.py` (epitope + occlusion) and rank per the methodology rubric.

## Cross-validation
Run the same complex through a second model (Chai-1) when available and compare **structural**
outputs (burial, hotspots, epitope), not absolute ipTM — the two disagree on magnitude by design.
