# Calibrating our Boltz score against real wet-lab labels

We spent this whole project scoring designs with Boltz `binding_confidence` / `ligand_iptm`
without knowing what those numbers buy in the lab. The **Anthropic
`claude-protein-binder-design`** release (1,320 de-novo miniprotein designs, 15 targets,
two-CRO binding measurements, CC BY 4.0) lets us calibrate. Reproduced here with the
`binder-benchmark-audit` skill's kernel — every number below is computed from the release,
not quoted.

## Reproduced headline (sanity check — matches the paper)
- analysis set: **1,320 designs / 15 targets**, overall hit rate **26.8 %** (paper 26.8 %).
- 3-predictor campaign score, within-target mean average precision **0.508** (paper 0.51),
  beats chance on **12/13** targets.

## What our predictor (Boltz-2) actually buys  ✅ real
`ipsae_min_boltz2` (the Boltz-2 analogue of the interface score we use), within-target:

| keep top N% per target (rank by Boltz-2) | n | hit rate |
|---|---|---|
| 100 % (no filter) | 1320 | 26.8 % |
| 50 % | 660 | 35.3 % |
| 30 % | 396 | 38.9 % |
| 10 % | 132 | 40.2 % |
| 5 % | 61 | 42.6 % |

- Boltz-2 within-target **mAP 0.462**, beats chance on **13/13 targets** — a robust triage
  signal (~1.5× enrichment), slightly weaker than the 3-predictor ensemble but more uniform.

## The caveat that changes how we read every Boltz number  ⚠️
Among actual binders, **Boltz-2 score vs log KD Spearman ρ = −0.05** (paper 0.16 for the
ensemble). **The score predicts *whether* a design binds, not *how tightly*.** So our
`binding_confidence` is a *feasibility / triage* signal — never read a higher Boltz score as
"tighter binder" or quote it as an affinity proxy. This is exactly the honesty label the
biosensor skill already uses; now it's quantified.

## Actionable rules for this skill's Boltz scoring
1. **Triage by within-target percentile, not an absolute cutoff.** Absolute thresholds
   silently re-weight toward whichever targets co-fold confidently (in the release, score
   ≥0.90 *drops* hit rate vs ≥0.80 because one bad target dominates the tail).
2. **Rank to decide *which to test*, not to estimate affinity or dynamic range.** Kd /
   ΔF·F₀ still require the wet-lab titration.
3. Expect **~1.5× enrichment** from Boltz triage on a fresh target — real but uneven per
   target (leave-one-target-out on the release: top-K beats base rate on ~9/15 targets, loses
   on the rest). Not a guarantee.
4. **A score gate is not an epitope/pose gate.** Interface confidence says it docks
   *somewhere*, not *where* (release: ipTM vs on-target-face fraction ρ≈0). For our
   specificity work this reinforces the whole thread's finding — co-folding confidence can't
   substitute for a structurally-defined discriminating contact.

## Reproduce
```
python3 binder_audit/run_audit.py      # fetches the release tables, writes audit_boltz2_calibration.json
```
Artifacts: `binder_audit/audit_boltz2_calibration.json` (numbers above).
Source of truth: `huggingface.co/datasets/Anthropic/claude-protein-binder-design`.
