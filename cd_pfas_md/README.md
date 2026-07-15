# cd-pfas-md — β-cyclodextrin / PFAS host-guest binding free energy pipeline

Molecular-dynamics tooling to support synthesis decisions for **β-cyclodextrin
(β-CD) based PFAS sensors**. The assay is a host-guest displacement system (no
proteins): a sulfonated reporter **dye** and an anionic **PFAS** compete for the
modified β-CD cavity; PFAS displaces the dye, and the signal quantifies PFAS.

This package computes the two numbers that drive design:

| Deliverable | What it answers | Module(s) |
|---|---|---|
| **Heuristic prefilter** | Cheap CPU triage: which modifications are even worth FEP, and is displacement favorable? | `prefilter` |
| **APR absolute ΔG_bind** | Does our force field + protocol reproduce the *measured* dye·CD affinity? (calibration) | `parameterize → build_apr → run_apr → analyze_apr` |
| **FEP/TI ΔΔG screen** | Which CD modification binds PFAS (or the dye) tighter than plain β-CD? (design triage) | `fep_ti_ddg` |

β-CD host-guest binding is a **standard benchmark** for these free-energy methods
(the pAPRika β-CD tutorials, the SAMPL host-guest challenges), so the approach is
well-trodden — the work is getting *your* dye, *your* PFAS, and *your*
modifications parameterized and calibrated.

## Install

```bash
conda env create -f cd_pfas_md/environment.yml
conda activate cd-pfas-md
```

A **GPU** (or HPC allocation) is needed for the production MD legs. Parameterization,
setup, analysis, and the unit tests run on CPU.

## 0) Cheap triage first (no GPU): the heuristic prefilter

Before spending FEP cycles, rank the whole modification library on CPU:

```bash
cd cd_pfas_md
python -m src.prefilter          # -> work/prefilter/prefilter_report.json + a table
```

This computes `ΔG_bind ≈ cavity + electrostatics + fluorophilic` (all constants
exposed in `config/system.yaml → prefilter:`) and the displacement figure of
merit `ΔG_bind(PFAS) − ΔG_bind(dye)`. It is **order-of-magnitude, not converged
MD** — its job is to tell you which designs deserve the real run. A worked demo
(TNS dye + a cationic vs fluorophilic modification comparison) with the actual
output and interpretation is in [`RESULTS_demo.md`](RESULTS_demo.md); the headline
is that a symmetric rim charge tunes *affinity* but not PFAS-vs-dye *selectivity*
(both are −1 anions) — the fluorous-lined cavity is the selectivity lever.

## 1) Calibrate on your dye·CD system (APR)

```bash
# put your dye structure/Ka in config/system.yaml first
cd cd_pfas_md
scripts/run_calibration.sh dye
# -> work/apr_dye/binding_result.json  (ΔG_bind, calc vs experimental Ka, verdict)
```

The report prints `ΔG_bind ± SEM`, the back-computed `Ka`, and a
GOOD/FAIR/POOR verdict vs your measured value. **Only trust the ΔΔG screen once
this calibration is GOOD (<1 kcal/mol) or at least FAIR.**

## 2) Screen CD modifications (FEP/TI ΔΔG)

Edit `config/modifications.yaml` (a starter library of cationic / anionic /
methylated rims is provided), then:

```bash
scripts/run_ddg_screen.sh pfoa     # rank modifications by ΔΔG vs PFOA
scripts/run_ddg_screen.sh dye      # repeat vs the dye to check selectivity
# -> work/ddg/ddg_screen_<guest>.json  + a ranked table on stdout
```

**Interpretation:** `ΔΔG_bind = ΔG_complex(WT→MOD) − ΔG_apo(WT→MOD)`.
A **negative** ΔΔG means the modification binds that guest **more tightly** than
plain β-CD. For a sensitive, robust assay you generally want modifications that
are **strongly negative for PFAS** and **near-zero/positive for the dye**, so the
displacement is thermodynamically favorable and the dynamic range is wide.

## Design rationale baked in

PFAS is anionic, so the modification library leans on **cationic rim groups**
(amino, trimethylammonium) as the sensitivity lever, with a **sulfobutylether**
anionic control (expected to *repel* PFAS) and a neutral **methylated** cavity
tweak. The ΔΔG screen quantifies these *before* you commit to synthesis.

## Honesty / scope notes

- The MD legs (`run_apr.py`, `fep_ti_ddg.run_alchemical_leg`) are real OpenMM /
  openmmtools drivers but are configured with **short smoke-test lengths**. Scale
  `production_ps` / `per_window_ps` to your accuracy target before publishing.
- Two clearly-marked **integration points** need your input:
  1. APR anchor atom selections (`apr_manifest.json`) must match your host's real
     Amber atom names.
  2. Non-trivial CD grafts should be supplied as explicit SMILES / pre-built mol2
     (the auto-graft raises a helpful error rather than guessing connectivity).
- **Fluorine parameters** are the largest accuracy risk for PFAS; `parameterize.py`
  flags any GAFF2 parameters it had to guess for review.

## Test

```bash
pytest cd_pfas_md/tests -q     # config schema + thermodynamics + ranking logic (no MD engine needed)
```

## Layout

```
cd_pfas_md/
  environment.yml
  config/
    system.yaml          # host, dye, PFAS, force field, APR + sim protocol
    modifications.yaml   # CD modification library for the ΔΔG screen
  src/
    parameterize.py      # GAFF2 + AM1-BCC/RESP; PFAS fluorine review
    build_apr.py         # pAPRika attach-pull-release window setup + tleap solvation
    run_apr.py           # OpenMM MD per window
    analyze_apr.py       # pAPRika fe_calc -> ΔG_bind + calibration vs experiment
    fep_ti_ddg.py        # alchemical WT->MOD host edit, complex+apo legs -> ΔΔG, batch
    utils.py             # config, units (Ka<->ΔG), platform selection
  scripts/
    run_calibration.sh   # APR end-to-end
    run_ddg_screen.sh    # ΔΔG batch
  data/hosts, data/guests   # drop your structures here
  tests/
```
