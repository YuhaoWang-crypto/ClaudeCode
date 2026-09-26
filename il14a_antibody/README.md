# IL-14α / TXLNA humanised blocking antibody design

Reproducible pipeline for selecting an epitope on IL-14α (α-taxilin, UniProt
P40222), triaging existing antibody parents, and humanising them into
deliverable constructs.

**Read [`REPORT.md`](REPORT.md) first** — it states what was and was not
delivered, and why the affinity axis is deliberately gated.

## Quick start

```bash
pip install numpy biopython
python3 src/step1_epitope_scan.py
python3 src/step1b_specificity_residues.py
python3 src/step2_epitope_geometry_and_parents.py
python3 src/step3_epitope_selectivity.py      # ~2 min
python3 src/step4_humanization.py
python3 src/step5_constructs.py
python3 src/step6_boltz_configs.py
```

## What each step does

| step | purpose | key output |
|---|---|---|
| `step1_epitope_scan` | scan epitopes; audit the four legacy G034 epitopes | `results/step1_epitope_scan.json` |
| `step1b_specificity_residues` | locate exposed TXLNA-unique residues vs β/γ-taxilin | `results/step1b_sdr_map.json` |
| `step2_epitope_geometry_and_parents` | measure the epitope's 3D footprint; grade existing parents | `results/step2_parents.json` |
| `step3_epitope_selectivity` | 4 directed rounds per epitope; best achievable paralog-selectivity margin | `results/step3_epitope_selectivity.json` |
| `step4_humanization` | 4 gated rounds of liability repair + germlining | `results/step4_humanization.json` |
| `step5_constructs` | assemble VHH-Fc and full IgG; QC | `results/IL14A_humanised_candidates.fasta` |
| `step6_boltz_configs` | epitope-constrained Boltz-2 payloads (not submitted) | `results/boltz_configs/` |

## Honesty conventions

Following the convention already used in this repository, every claim is
labelled by what backs it:

- **Computed from primary data** — UniProt sequences, the AlphaFold model,
  germline repertoires. Reproducible by rerunning the step.
- **Heuristic** — the `pair_score` interaction model in `step3`. Explicitly
  stated chemistry, not a trained potential. Only the *specificity margin*
  (a difference between three near-identical surfaces) is interpretable;
  absolute values are not affinities.
- **Gated** — anything needing an antigen–paratope interface, MHC-II
  prediction, or wet-lab measurement. Not estimated, not guessed.

No Kd, IC50, ΔΔG or "N-fold improvement" appears anywhere in this project.
