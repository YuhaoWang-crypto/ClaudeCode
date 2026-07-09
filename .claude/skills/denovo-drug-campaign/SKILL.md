---
name: denovo-drug-campaign
description: >-
  Run an end-to-end de novo multi-modal drug-discovery campaign for ANY disease,
  from market-based indication scoring through multi-level target-axis discovery,
  biomarker mining, a target x modality feasibility matrix, GPU-executed designs
  (antibody / mini-binder / gene / small-molecule / PROTAC), consolidated
  safety-developability ranking, targeted delivery formulation, and a
  network-biomarker irreducibility + critical-transition analysis. Use when the
  task is "pick an indication and design drugs from scratch", multi-modal target
  triage, biomarker discovery for enrollment/monitoring, or reproducing this
  pipeline on a new disease by changing one config. Config-driven: retarget by
  editing the disease + candidate list. Enforces grounded data (Open Targets,
  cBioPortal, CIViC, ClinicalTrials) and rigorous vs hypothesis labeling.
---

# De novo multi-modal drug-discovery campaign

A reusable, config-driven pipeline that takes a **disease** as input and produces
a ranked, multi-modal drug-design dossier — every number computed from grounded
data or a real GPU design run, never asserted. Validated end-to-end on NSCLC
(KEAP1/NRF2 redox-metabolic axis).

Helper functions ship in `kernel.py` (auto-loaded). MCP data-access patterns and
GPU-dispatch recipes are in `reference/`.

## The one input that changes everything

Edit `config` (see `reference/config-example.md`):
- `disease` (name + MONDO/EFO id), `candidate_indications` (for Stage 1 scoring),
- optional `seed_targets` / `pathway_edges` if you already know the axis.

Everything downstream keys off this. Nothing else needs editing to retarget.

## The 11 stages (run in order; each ends with save_artifacts + a figure)

| Stage | What | Tools | Key helper |
|---|---|---|---|
| 1 | Indication scoring | Open Targets assoc counts + ClinicalTrials recruiting counts (MCP) | `composite_score`, `normalize_log` |
| 2 | Multi-level target axis (gene->protein->metabolism) | Open Targets top targets + cBioPortal mut freq (MCP) | — |
| 2b | **Network-biomarker layer** | fibration -> CRNT deficiency -> critical slowing | `input_fibers`, `crnt_deficiency`, `schlogl_states`, `critical_slowing`, `langevin_earlywarning` |
| 3 | Biomarker mining + prioritization | CIViC evidence (MCP) + mechanistic supplementation | — |
| 4 | Synthetic + genomic biomarker | Proto cassette + Evo2 naturalness (GPU) | — |
| 5 | target x modality feasibility matrix | PDB/AlphaFold availability (MCP) | — |
| 6 | Protein / antibody / gene GPU designs | RFdiffusion+ProteinMPNN+RF2, Boltz-2, oligo rules | — |
| 7 | Small-molecule / PROTAC GPU designs | Boltz-2 co-fold + RDKit ADMET | — |
| 8 | Consolidated safety/developability ranking | weighted composite | `rank_leads` |
| 9 | Targeted delivery formulation | Mihaila 5-ODE siRNA-LNP (lnp-delivery-kinetics skill) | `lnp_knockdown_from_escape` |
| 10 | Master dossier | HTML embedding all figures via artifact markers | — |

Stage 2b is the network-biomarker cross-check: it independently validates the
target axis (irreducible core) and yields early-warning biomarkers. See
`reference/network-biomarker.md`.

## Cost control (respect the user budget)

- Antibody campaigns: **cap at 100 designs** (25 backbones x 4 seqs) per campaign.
- Small-molecule / mini-binder co-folds: batch **4 per Boltz job** on one A100.
- Execute a **representative real subset** of the modality matrix, not all 25 cells
  — pick the top cell per modality so all modalities are covered with real GPU
  evidence, then rank. Document which cells were executed vs matrix-only.

## Non-negotiable discipline (learned the hard way)

1. **Ground before you claim.** Every market/target/biomarker number comes from a
   live MCP pull (Open Targets, cBioPortal, CIViC, ClinicalTrials), cited by id.
   Never assert prevalence/counts from memory.
2. **View before you assert legibility.** After saving any figure, call
   `read_file(version_id=...)` and actually look BEFORE writing "clean/legible".
   A geometric overlap count is a pre-check, not a substitute for viewing.
3. **Scope your metrics honestly.** If a GPU run recovers a metric for n<N designs
   (e.g. pLDDT for 2 of 100 from a log tail), say "n=2 recoverable", never imply
   the full set. Saved artifacts must carry the same caveat as the prose.
4. **Label rigorous vs hypothesis.** Deterministic computation / measured statistic
   = rigorous; illustrative parameterization or untested clinical mapping =
   hypothesis. Put both in the figure caption and the dossier.
5. **ipTM/pLDDT validate pose, not affinity.** Reference chemotypes co-folded at
   high ipTM are pose validations of known binders, not novel-scaffold potency.

## Helper quick reference (kernel.py)

```python
composite_score(row, weights=None)          # Stage 1 indication composite
normalize_log(values)                        # counts -> 0-5
input_fibers(edges)                          # 2b: nodes -> irreducible fibers
crnt_deficiency(n_complexes, l, s)           # 2b: delta = n - l - s
schlogl_states(k1,k2,k3,p)                   # 2b: (roots, is_stable)
critical_slowing(k1,k2,k3,p, branch="lo")    # 2b: lam/tau/variance/ar1 analytic
langevin_earlywarning(k1,k2,k3,p, seed=0)    # 2b: MEASURED var/ar1 on SDE
rank_leads(leads, weights=None)              # Stage 8 weighted ranking
```

See `reference/config-example.md` to retarget, `reference/network-biomarker.md` for
the Stage-2b math and gotchas, `reference/gpu-recipes.md` for Modal dispatch patterns.
