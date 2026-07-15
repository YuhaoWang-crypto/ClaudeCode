# In-silico validation with Boltz

How the `biosensor_pipeline` maps biosensor-design questions onto Boltz-2.1
structure+binding predictions, and how to read the results honestly.

## What structure prediction can and cannot tell you

Boltz predicts a **static structure + a binding metric**. It does **not** predict
an allosteric **dynamic range** (a kinetic, entropic quantity). So the in-silico
layer answers only the *necessary-condition* questions that gate bench work:

| Design question | Boltz signal | Rigor |
|---|---|---|
| Does the chimera fold at all? | complex pLDDT, pTM | ✅ number from model |
| Does the receptor still bind its analyte after CP + insertion? | ligand_iptm, binding_confidence (holo vs native control) | ✅ number |
| Is the reporter active site left intact? | Cα–Cα constellation of catalytic residues vs native enzyme | ✅ geometric measurement |
| Any hint of ligand-coupled ordering? | apo→holo pLDDT / pTM change | ⚠️ weak, interpretive |
| Which library member to test first? | switch proxy (combines the above) | ⚠️ illustrative ranking |

## The three predictions per construct

For each candidate chimera (built by `screen.build_library`):

1. **holo** — chimera protein (chain A) + analyte ligand (chain L), with
   `binding: {type: ligand_protein_binding, binder_chain_id: "L"}`.
   → fold + does the receptor bind the analyte inside the chimera?
2. **apo** — chimera protein only. → OFF-state scaffold fold, apo/holo contrast.
3. **control** — native receptor + analyte. → positive control that the
   binder/ligand pair is modellable at all; the denominator for
   "binding retention vs native".

`boltz_io.py` builds these payloads. Runs use **single-sequence mode**
(`msa: {type: empty}`) because de-novo/engineered chains have no meaningful MSA.

## Running it (Boltz remote MCP)

Predictions are ~$0.05 each and finish in well under a minute at these sizes.

1. `boltz_estimate_structure_and_binding` — free cost/validity check.
2. `boltz_start_structure_and_binding` — pass a stable `idempotency_key`
   (e.g. `biosensor-<system>-<site>-<holo|apo>-r1`) so retries don't double-bill.
3. `boltz_get_structure_and_binding_prediction` — poll; on `succeeded`, read
   `output.best_sample.metrics` and `output.binding_metrics`, and download the
   `structure.url` CIF (presigned, ~30-min TTL) for the geometric check.
4. Save the compact metric table to `biosensor_out/boltz_results.json` and run
   `python3 -m biosensor_pipeline.analyze_boltz`.

## Reading the numbers

- **binding_retention_vs_native** = `chimera ligand_iptm / native ligand_iptm`.
  ≈ 1 means circular permutation + insertion did not damage the receptor. This
  is the single most informative in-silico check of the recipe.
- **catalytic constellation Δmax vs native** — a few Å is fine (models wobble);
  a large expansion means the insertion disrupted the active site → deprioritize.
- **switch proxy** — a transparent weighted combination (see `scoring.py`).
  Treat it as a **triage rank only**. A high proxy is *permission to test at the
  bench*, never a claimed working biosensor.

## Honest failure modes to expect

- Low-MSA de-novo binders (e.g. DIG10.3) get **modest absolute confidence**
  (~0.66). That is a property of single-sequence modelling, not evidence the
  design fails — compare chimera **to its own native control**, not to an
  absolute threshold.
- The chimera **pTM** often drops vs the native receptor because the two domains'
  relative orientation is genuinely uncertain (they are only loosely coupled).
  Per-domain **pLDDT** and **ligand_iptm** staying high is the meaningful signal.
- A high apo pLDDT does **not** disprove switchability — the switch is entropic
  and largely invisible to a single static prediction. Absence of an apo→holo
  shift is *uninformative*, not negative.
