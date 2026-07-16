---
name: cd-pfas-md
description: >-
  Compute β-cyclodextrin / PFAS host-guest binding free energies to guide
  cyclodextrin-sensor synthesis decisions — NO proteins, pure small-molecule
  host-guest. Three tiers: a cheap CPU displacement prefilter, an APR
  attach-pull-release absolute ΔG calibration against a measured dye·CD Ka, and
  an FEP/TI ΔΔG screen of CD modifications, all runnable on Modal GPUs. Use when
  designing/ranking modified cyclodextrins (cationic, fluorophilic, methylated
  rims) for a dye-displacement PFAS assay, calibrating a force field against a
  known Ka, parameterizing anionic/fluorinated guests, or extending the
  cd_pfas_md package with a new dye/PFAS/modification. Enforces ✅-validated vs
  ⚠️-needs-review labeling on every number.
---

# β-CD / PFAS host-guest binding free energy pipeline

A reusable methodology (and the working `cd_pfas_md/` package) for the design
question behind a cyclodextrin PFAS sensor: **which modified β-CD binds PFAS
tightly enough to displace the reporter dye, giving a sensitive, robust assay?**

The assay is host-guest displacement (no proteins): a sulfonated **dye** and an
anionic **PFAS** compete for the β-CD cavity; PFAS displaces the dye and the
signal quantifies PFAS. The physics that matters is **ΔG_bind(PFAS) vs
ΔG_bind(dye)** for each candidate host.

## Three tiers — pick by how much confidence you need

| Tier | Question | Cost | Module | Rigor |
|---|---|---|---|---|
| **Prefilter** | which modifications are even worth simulating? is displacement favorable? | seconds, CPU | `prefilter` | ⚠️ heuristic (transparent physics, all constants in config) |
| **APR** | does our FF+protocol reproduce the *measured* dye·CD Ka? (calibration) | GPU-hours | `parameterize→build_apr→run_apr→analyze_apr` | ✅ TI+MBAR, benchmark-grade for CD |
| **FEP/TI ΔΔG** | which modification binds PFAS/dye tighter than plain β-CD? | GPU-hours | `fep_ti_ddg` | ✅ relative FE; ⚠️ graft-decouple approximation (see below) |

Always run the **prefilter first** to triage, then **calibrate with APR** on a
system with a known Ka, then **screen with FEP** — never trust the screen before
the APR calibration is GOOD (<1 kcal/mol) or FAIR (<2).

## Execution status — what has actually run (be honest about this)

✅ **Verified here on CPU — reproducible right now, backed by the test suite:**
```bash
python -m cd_pfas_md.src.prefilter            # ranked displacement table (8 hosts)
python -m pytest cd_pfas_md/tests -q          # 16 passed
python -m cd_pfas_md.src.build_modified_host  # 3 modified hosts, exact formulas/charges
```
- Prefilter produces the displacement ranking (fluorophilic hosts FAVORABLE,
  symmetric-charge hosts UNFAVORABLE — the selectivity insight).
- The APR restraint force + unit conversion is validated against the real OpenMM
  API to machine precision (1e-16); the graft stoichiometry against RDKit.
- Host = real PDB CCD β-CD (147 atoms); 3 modified hosts built + shipped.

⏳ **Wired + smoke-testable, first GPU run on Modal still to be confirmed:**
`::check`, `--mode smoke`, `--mode prod`, `::screen`. Everything CPU-side of the
MD engine (config, parameterization plumbing, window/λ build, anchor resolution,
restraint construction, analysis aggregation) is validated; the GPU legs are the
one thing that needs a real Modal run to sign off. Paste the first `check`+`smoke`
logs to close this out.

**Rule for this skill: never label an APR/FEP number as trustworthy until (a) the
Modal smoke run is green and (b) the APR calibration vs the measured Ka is
GOOD/FAIR.** The CPU tier and all structure/parameter prep are the proven core.

## The key design insight the prefilter encodes

PFAS and a sulfonated dye are *both* −1 anions, so a symmetric cationic rim
charge stabilizes them almost equally — it tunes **overall affinity/robustness,
not PFAS-vs-dye selectivity**. Selectivity comes from what PFAS has that the dye
doesn't: the **fluorophilic tail**. So a fluorous-lined cavity is the selectivity
lever; charge is the affinity/detection-limit lever. Combine both.

## How to run (Modal GPUs)

```bash
pip install modal && modal setup
modal run cd_pfas_md/modal_app.py::triage                     # CPU prefilter, instant
modal run cd_pfas_md/modal_app.py::check                      # preflight: CUDA + toolchain
modal run cd_pfas_md/modal_app.py --guest dye --mode smoke    # validate APR path (T4)
modal run cd_pfas_md/modal_app.py --guest dye --mode prod     # real APR calibration
modal run cd_pfas_md/modal_app.py::screen --guest pfoa --mode prod   # FEP ΔΔG screen
modal volume get cd-pfas-work /apr_dye ./pulled               # fetch results
```

`--mode smoke` shrinks every sim to seconds — always smoke a new system before
prod. Windows/λ-legs fan out one-GPU-container-each, so wall-clock ≈ one window.
Full run guide + cost table: `cd_pfas_md/MODAL.md`.

## Adding a new dye / PFAS / modification

- **Dye or PFAS**: add an entry under `guests:` in `config/system.yaml` with a
  SMILES (deprotonated anion — assert the formal charge), `cavity_dG_kcal`, and
  `fluorinated:`. Anchor the dye to its measured Ka via `experimental.Ka_M_inv`.
- **CD modification**: add to `config/modifications.yaml`. For the three built-in
  designs (mono_6_trimethylammonium, fluorous_tagged, fluorous_cationic) the
  structure is auto-grafted onto `data/hosts/bcd.sdf` by `build_modified_host.py`.
  For others, supply an explicit `smiles:`/`mol2:` or add a builder.

## What is validated vs what needs your review

✅ **Validated (CPU, in the test suite — `pytest cd_pfas_md/tests`)**
- APR restraint force + kcal→kJ / Å→nm conversion, checked against real OpenMM to
  machine precision.
- Modified-host graft stoichiometry + formal charge (C45H78NO34⁺, C54H73F21O35,
  C57H81F21NO34⁺ — exactly the expected chemistry).
- Ka↔ΔG conversions, config schema, ΔΔG ranking, attach/pull k+r0 schedule.
- The host structure is the real PDB CCD β-CD (BCD, C42H70O35), not an RDKit guess.

⚠️ **Needs your review before trusting production numbers**
1. **PFAS fluorine parameters** — GAFF2 can misrepresent C–F; `parameterize.py`
   flags any guessed torsions. Prefer RESP for publication ΔG.
2. **APR anchors** — auto-resolved (principal-axis rim + head charge); inspect
   `apr_manifest.json` before prod.
3. **Grafted geometry** — unminimized starting points; tleap+min relaxes them.
4. **FEP topology** — the graft is treated as an alchemically decoupled group (a
   fast RANKING scheme), not a rigorous single-topology O6H↔graft morph. Confirm
   the winning design with a dual-topology run (perses/openfe) before synthesis.

## Package map

```
cd_pfas_md/
  config/system.yaml         # host (bcd.sdf), TNS dye, PFOA/PFOS, FF, APR + prefilter params
  config/modifications.yaml  # CD modification library
  src/prefilter.py           # tier-1 heuristic displacement triage
  src/parameterize.py        # GAFF2 + AM1-BCC/RESP; sdf/pdb/smiles; fluorine review
  src/build_apr.py + anchors.py + run_apr.py + analyze_apr.py   # APR tier
  src/build_modified_host.py # graft modified hosts onto bcd.sdf
  src/fep_ti_ddg.py          # tier-3 ΔΔG screen
  modal_app.py               # Modal GPU orchestration (triage/check/main/screen)
  data/hosts/bcd.sdf         # real 3D β-CD; data/hosts/modified/*.sdf  built modifications
  tests/                     # CPU validation (16 tests)
  MODAL.md, RESULTS_demo.md
```
