# Drug repurposing — reuse an existing drug for a new indication

Discovering a target asks *"what gene reverses the disease?"*. Repurposing asks
a tighter question the Connectivity Map was built for: *"which **already-approved
drug** reverses the disease **and is not already used for it**?"* The skill
answers it with two ingredients on top of the ordinary reverser ranking
(`perturbomics/repurpose.py`).

## The two extra ingredients

1. **Existing-drug filter.** A repurposing candidate must already be a real drug.
   Keep only *launched* / late-clinical compounds. Source of truth: the **Broad
   Drug Repurposing Hub** annotation file (`repurposing_drugs_*.txt`), column
   `clinical_phase` ∈ {Launched, Phase 3/2/1, Preclinical, Withdrawn}.

2. **Indication novelty.** Compare the drug's **current** indication /
   disease-area (Hub columns `indication`, `disease_area`) to the target disease:
   - already indicated for it → **not** repurposing (excluded),
   - same disease area, different indication → adjacent (novelty 0.5),
   - different disease area → true indication shift (novelty 1.0).

```
repurpose_score = |WTCS| (reversal) × phase_weight (de-risking) × indication_novelty
```

`phase_weight` rewards de-risking: Launched 1.0 → Phase 3 0.85 → … → Preclinical
0.15. So the ranking favours drugs that both **reverse the disease** and are
**already approved for something else**.

## Data

| File | Where | Columns used |
|---|---|---|
| `repurposing_drugs_20200324.txt` | `data.clue.io/repurposing/downloads/` | `pert_iname, clinical_phase, moa, target, disease_area, indication` |

`load_repurposing_hub(path)` parses it (skips the `!`-comment header) into
`{drug → annotation}`. Drug names from L1000 match the Hub `pert_iname`
(salt-suffixes normalised).

## Two modes of repurposing (surface both, don't collapse)

Cross the candidate against **ClinicalTrials** ("is it already tried for THIS
disease?"):

- **De-risked / confirmatory** — reverses the signature **and** already in trials
  for the disease. The method recovered a real repurposing bet → high confidence.
- **Novel white-space** — reverses the signature, approved elsewhere, **no**
  trials for the disease yet → a genuinely new hypothesis to test.

## Worked real run (IPF) — `examples/repurpose_ipf.py`

Real L1000 IPF reversers × the real Hub → launched drugs reversing IPF from a
non-respiratory indication:

| score | drug | now approved for | mode (live ClinicalTrials) |
|---|---|---|---|
| 0.79 | **neratinib** (EGFR inhibitor) | breast cancer | **novel** — 0 IPF trials; aligns with the EGFR/canertinib finding |
| 0.75 | **dasatinib** (BCR-ABL/SRC TKI) | CML / ALL | **de-risked** — NCT02874989 (Dasatinib+Quercetin, IPF Ph1), NCT00764309 (scleroderma lung fibrosis) |
| 0.73 | **vorinostat** (HDAC inhibitor) | cutaneous T-cell lymphoma | novel — echoes the trichostatin/HDAC antifibrotic angle |
| 0.81 | betamethasone / budesonide (steroids) | dermatoses / Crohn's | anti-inflammatory signature match (⚠️ steroids are not effective in IPF specifically) |

The two anchor hits are exactly what a good repurposing screen should do:
**dasatinib** — which the screen nominated from transcriptomics alone — is
*genuinely* in IPF trials (independent confirmation); **neratinib/EGFR** is an
untried, mechanistically-coherent new hypothesis.

## Chaining with the rest of the skill

Repurposing is a specialisation of the integration funnel
(`reference/integration.md`): the Hub gives the **clinical/engageability**
evidence (launched + target + indication), the network-biomarker skill adds
**network control** (is the target on-pathway?), and `integrated_leads` fuses
them. Use `screen_repurposing` to get the launched-drug shortlist, then
`integrated_leads` to rank the shortlist across all axes.

## Rigor
- ✅ **rigorous**: the |WTCS| reversal, the launched status and current
  indication (read from the Hub), the trial status (ClinicalTrials), the score.
- ⚠️ **hypothesis**: *"reverses the signature ⇒ treats the new indication"*. A
  repurposing hypothesis to validate experimentally — **never** a clinical or
  treatment recommendation. (The Hub itself states: do not use it to make
  clinical treatment decisions.)
