# Reproduction of "Autonomous de novo protein binder design with Claude"

Claude Science & A. Shanehsazzadeh, Anthropic, 18 August 2026.
Data release: <https://huggingface.co/datasets/Anthropic/claude-protein-binder-design> (CC-BY-4.0).

## What this reproduces

`reproduce_claims.py` re-derives **97 numerical claims** from the paper directly out of the
released per-design tables — hit rates, rank calibration, co-folding average precision,
affinity distributions, cross-species reactivity, method attribution, backbone-level
de-duplication and fold diversity. **All 97 match.**

```bash
pip install pandas pyarrow scipy
curl -L -o tables.zip \
  https://huggingface.co/datasets/Anthropic/claude-protein-binder-design/resolve/main/protein_binder_design_data_release_docs_tables.zip
unzip tables.zip
python reproduce_claims.py --release-dir protein_binder_design_data_release
```

Outputs `reproduction_report.{json,csv}` (one row per claim) and `per_target_hit_rates.csv`.
`run_output.txt` is the recorded run.

## Selected results

| Claim | Paper | Recomputed |
|---|---|---|
| Designs delivered / with interpretable data | 1,440 / 1,320 | 1,440 / 1,320 |
| Binders (integrated call) | 354 (26.8%) | 354 (26.8%) |
| Targets with ≥1 binder | 14 of 15 | 14 of 15 |
| Opus 4.8 multi-target | 88/390 (22.6%) | 88/390 (22.6%) |
| Mythos Preview multi-target | 104/390 (26.7%) | 104/390 (26.7%) |
| Mythos Preview single-target | 158/450 (35.1%) | 158/450 (35.1%) |
| Single vs multi, same model & 13 targets | 143 vs 104, p = 0.003 | 143 vs 104, p = 0.003 |
| Hit rate top-1 / top-5 / top-10 / all 30 | 49 / 44 / 39 / 28 % | 49 / 44 / 39 / 28 % |
| AP of Claude's delivered rank vs permutation null | 0.48 vs 0.35, p < 0.001 | 0.50 vs 0.355, p < 0.001 |
| AP of the re-scored co-folding score | 0.52 (chance 0.31) | 0.51 (chance 0.305) |
| AP of the 7-predictor held-out ensemble | 0.57 | 0.57 |
| Binders below 100 / 10 / 1 nM | 194 / 90 / 42 | 194 / 90 / 42 |
| RBX1 top design, Adaptyv K_D | 3.9 nM | 3.9 nM |
| RBX1 top design, K_D of record | 7.0 nM | 7.0 nM |
| TNFα binders, all from Opus 4.8 | 12/150 (8.0%) | 12/150 (8.0%) |
| Mouse / cyno cross-reactivity among binders | 130/233, 154/179 | 130/233, 154/179 |
| Backbone-level hit rate | 200/809 (24.7%) | 200/809 (24.7%) |
| Designs not all-α | 126 of 1,320 | 126 of 1,320 |
| Tested designs per generator (PXDesign … ProteinHunter) | 358, 267, 185, 135, 134, 118, 100, 14, 2, 2 | identical |

Two aggregation details the paper does not fully specify, resolved here:

* **The paper's "tested" denominator is 1,315, not 1,320** — five designs reached neither CRO
  (`vendor_agreement == "not_tested_either"`). Method-level counts and hit rates only match the
  paper on the 1,315 subset (e.g. FreeBindCraft 58/135 = 43%, matching the stated 22–43% range).
* **Rank average precision.** The paper reports 0.48 for AP "computed within each target and
  averaged over the 13 targets" without saying how the several campaigns per target are pooled.
  Pooling the three main campaigns' designs into one list per target gives 0.496; including the
  two Opus 4.8 single-target campaigns gives 0.466. Both bracket 0.48, and both sit far above the
  0.355 permutation null (p < 0.001) and just under the co-folding score's 0.51–0.52 — the paper's
  actual conclusion (Claude's ranking is calibrated, but does not beat the score it was based on).

Two further caveats worth stating:

* Because the delivered designs were already filtered on the same three predictors, the
  re-scored AP (0.51) understates pre-selection enrichment, as the paper itself notes.
* Every structure in the paper is a prediction; no design was solved experimentally. The
  reproduction inherits that limit.

## What this does *not* reproduce

The paper has three layers; this repository closes only the third.

1. **The campaigns** — Claude Opus 4.8 / Mythos Preview running 24–48 h autonomously from the
   protocol prompt, building each open-source tool from source on Modal GPUs.
   Cost as run: USD $10,000 per single-target campaign, $50,000 per multi-target campaign.
   The prompts are released (`prompts/` in the data release), so this is *runnable* by anyone
   with the GPU budget — but a rerun samples fresh designs, so it can be checked only in
   distribution (hit rate), never design-by-design.
2. **The wet lab** — every design synthesized and measured by SPR/BLI at Adaptyv Bio and Twist
   Bioscience. Months and a five-to-six-figure budget; not reproducible computationally by
   construction. This is where all binding evidence comes from.
3. **The analysis** — everything in this repository. Fully reproducible, and reproduced.
