# Running on Modal (modal.com) GPUs

`modal_app.py` runs the whole pipeline on Modal: CPU containers for
parameterization/build/analysis, one **GPU container per APR window and per FEP
λ-leg** (fanned out with `.starmap`, so wall-clock ≈ one window, not N), and a
persistent **Volume** (`cd-pfas-work`) holding all outputs under `/work`.

## One-time setup

```bash
pip install modal
modal setup                       # authenticates to your Modal workspace
```

## Commands

```bash
# 0) instant CPU triage (no GPU, no structures needed) — sanity check the graph
modal run cd_pfas_md/modal_app.py::triage

# 1) SMOKE test the full APR path end-to-end (tiny sims, ~minutes on a T4).
#    Proves the container image + CUDA + restraints + analysis all wire up.
modal run cd_pfas_md/modal_app.py --guest dye --mode smoke

# 2) PRODUCTION APR calibration on the dye·β-CD system (real sim lengths).
modal run cd_pfas_md/modal_app.py --guest dye --mode prod

# 3) FEP/TI ΔΔG screen of the modification library vs a guest.
modal run cd_pfas_md/modal_app.py::screen --guest pfoa --mode prod

# 4) pull results to your laptop
modal volume get cd-pfas-work /apr_dye ./pulled
modal volume get cd-pfas-work /ddg    ./pulled
```

## GPU sizing & rough cost

Set `GPU_TYPE` at the top of `modal_app.py`.

| System | GPU | ~ per-window prod (5 ns) | APR (18 windows) | Notes |
|---|---|---|---|---|
| β-CD + small guest (~5–8k atoms) | **T4** | ~20–40 min | parallel ⇒ ~40 min wall | cheapest; fine for this size |
| same, faster | **A10G** (default) | ~8–15 min | ~15 min wall | best price/perf here |
| large modified host / dimeric CD | **A100** | — | — | only if the box is big |

Because windows run **in parallel**, wall-clock is ~one window; the *cost* is the
sum (≈18 window-hours for a full APR, a few $ on T4/A10G at current rates). The
FEP screen cost scales with `#modifications × 2 legs × lambda_windows`.

Start every new system with `--mode smoke` (seconds/window) to confirm it runs
before committing prod GPU-hours.

## What runs turn-key vs. needs your input

| Stage | Status on Modal |
|---|---|
| `triage` (prefilter) | ✅ runs now, no inputs |
| Parameterize **TNS dye** (from SMILES) | ✅ automated (RDKit → antechamber) |
| Parameterize **PFOA/PFOS** | ✅ automated, with fluorine parameter review |
| Parameterize **β-CD host** | ⚠️ needs `data/hosts/bcd.mol2` — drop a validated β-CD structure (pAPRika β-CD tutorial, or a GLYCAM build). RDKit can't reliably 3D-embed the macrocycle. |
| APR anchors | ✅ auto-resolved by `src/anchors.py` (principal-axis rim + head charge) — **inspect `apr_manifest.json` before prod** |
| APR windows + analysis | ✅ once the host structure is in |
| FEP screen: cationic/methylated mods | ⚠️ each modified host needs an explicit structure (`smiles:`/mol2). The prefilter already says the symmetric-charge mods don't buy selectivity, so prioritize the fluorophilic ones. |

**Bottom line:** the only blocker to a real GPU APR number for the **TNS·β-CD**
calibration is dropping one validated `bcd.mol2` into `data/hosts/`. Everything
downstream (parameterize → build → 18 parallel GPU windows → calibrated ΔG vs the
−4.9 kcal/mol anchor) then runs unattended.

## Troubleshooting

- **CUDA/OpenMM platform error**: the conda-forge `openmm` build must match
  Modal's driver. If OpenMM can't find CUDA, pin a CUDA version in the image
  (`.micromamba_install(..., "cuda-version=12.4")`) or fall back to `gpu="T4"`.
- **Volume not updating**: functions call `vol.commit()` after writing and
  `vol.reload()` before reading; if you add stages, keep that pattern.
- **First run is slow**: Modal builds the conda image once, then caches it.
