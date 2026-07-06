---
name: mof-drug-modeling
description: >
  Run validated materials & drug computational-chemistry predictions on Modal
  GPUs using Meta FAIR's UMA machine-learned potential (fairchem) and RASPA
  GCMC. Covers MOF adsorption energy, gas adsorption isotherms, MOF water
  uptake vs relative humidity (吸水率/吸湿), drug polymorph stability ranking
  (结晶/多晶型), MOF drug-loading affinity (载药), and inorganic material
  relaxation. Use when the user wants to predict MOF adsorption or water
  uptake, rank drug crystal polymorphs, estimate drug loading in a MOF, or
  relax/score a crystal with UMA. Runs on the user's Modal account.
---

# MOF & Drug Modeling (UMA + RASPA on Modal)

Validated pipelines living in `fairchem-modal/`. Compute runs on the user's
**Modal** account (GPU); this environment is CPU-only and only orchestrates.

## Prerequisites (check first)

1. `MODAL_TOKEN_ID` and `MODAL_TOKEN_SECRET` must be in the environment.
   Verify: `python3 -c "import os;print(bool(os.getenv('MODAL_TOKEN_ID')))"`
2. `pip install "modal[api-proxy-support]"` (the `python-socks` extra is
   required behind the agent proxy, else `modal` errors on connect).
3. A Modal secret named `huggingface` holding an `HF_TOKEN` with access to the
   gated `facebook/UMA` (+ OMol25/ODAC25) repos. Create once:
   `modal secret create huggingface HF_TOKEN=hf_xxx`
4. **Always `cd fairchem-modal` before any `modal run`/`modal deploy`** — the
   entrypoints reference local files by relative path.

## Two ways to run

- **Deployed endpoints (persistent, preferred for one-offs):** apps
  `fairchem-uma` and `raspa-gcmc` are deployed and scale to zero. Call from any
  machine with Modal creds via `fairchem-modal/client.py` or
  `modal.Function.from_name("fairchem-uma","_mof_adsorption").remote(...)`.
  Re-deploy after code edits: `modal deploy modal_app.py && modal deploy gcmc_raspa.py`.
- **`modal run` (ephemeral CLI, good for iterating):** commands below.

## Pipelines & commands (run from `fairchem-modal/`)

```bash
# 0) smoke test — cheapest end-to-end check (H2O energy on GPU)
modal run modal_app.py::smoke_test

# 1) MOF single-site adsorption energy (UMA odac head)
modal run modal_app.py::mof_adsorption --cif inputs/MOF-5.cif --adsorbate H2O

# 2) gas adsorption isotherm (RASPA GCMC)
modal run gcmc_raspa.py::isotherm --cif inputs/MOF-5.cif --molecule CO2 \
    --temperature 298 --pressures 1e3,1e4,5e4,1e5 --cycles 2000 --init-cycles 1000

# 2b) WATER uptake vs relative humidity (吸水率) — TIP4P + Ewald + framework charges
modal run gcmc_raspa.py::water_isotherm --cif inputs/HKUST-1.cif \
    --rh 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 --cycles 2500 --init-cycles 1500

# 3) drug polymorph stability ranking (UMA omc head)
modal run modal_app.py::polymorph_rank --cif-dir drugs/paracetamol_polymorphs \
    --atoms-per-molecule 20            # atoms per molecule (paracetamol C8H9NO2 = 20)

# 4) MOF drug-loading host–guest affinity (UMA odac head)
modal run modal_app.py::drug_loading --mof-cif inputs/ZIF-8.cif --drug-dir drugs/loading

# inorganic material relaxation + energy (UMA omat head)
modal run modal_app.py::relax_material --cif your_crystal.cif

# warm the model cache (UMA / OMol25 / ODAC25 into the hf-cache volume)
modal run setup_models.py::warm
```

## Getting input structures

- **MOF CIFs**: iRASPA/RASPA2 repo, e.g.
  `https://raw.githubusercontent.com/iRASPA/RASPA2/master/structures/mofs/cif/<NAME>.cif`
  (have: MOF-5=IRMOF-1, HKUST-1=Cu-BTC, ZIF-8). For **water GCMC the CIF must
  carry `_atom_site_charge`** — HKUST-1/ZIF-8 from that repo do.
- **Drug 3D molecules**: PubChem 3D SDF —
  `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/<NAME>/record/SDF?record_type=3d`
- **Drug/organic crystals (polymorphs)**: Crystallography Open Database (COD) —
  search `https://www.crystallography.net/cod/result?text=<name>&format=lst`,
  fetch `https://www.crystallography.net/cod/<ID>.cif`. Distinguish polymorphs
  by space group / cell angles.

## Key facts baked in (don't rediscover)

- **UMA model tag is `uma-s-1p2`** (not `uma-s-1`). Others: `uma-m-1p1`,
  `esen-sm-conserving-all-omol` (OMol25), `esen-sm-full-odac25` (ODAC25).
  OMC25/OMat24 have no standalone tag — use UMA's `omc`/`omat` task heads.
- **UMA task heads**: `omol` (isolated molecule), `omc` (molecular crystal /
  polymorph), `odac` (MOF host + guest: adsorption & drug loading), `omat`
  (inorganic material), `oc20` (catalysis). Use ONE head across all terms of an
  energy difference so references stay consistent.
- **RASPA data lives at `/opt/conda/share/raspa`** (dir name `raspa`, not
  `raspa2`); set `RASPA_DIR=/opt/conda`. Forcefield `ExampleMOFsForceField`
  has Cu/Zn/etc.; molecule set `ExampleDefinitions` (CO2 yes, **no water**).
- **Water GCMC** needs a hand-authored TIP4P model + pseudo-atoms/LJ injected
  into a runtime-built RASPA_DIR + `UseChargesFromCIFFile yes` + Ewald — all
  already implemented in `gcmc_raspa.py::_water_isotherm`.
- `polymorph_rank` infers Z from atom count → pass `--atoms-per-molecule`.

## Interpreting results (state these caveats)

- **Adsorption/drug-loading energy** (eV or kJ/mol): negative = favorable.
  Single centred pose + fixed framework → screening-grade; sample multiple
  poses for reliable numbers. Drug loading = thermodynamic binding only, no
  diffusion/kinetics (e.g. ZIF-8's small windows gate real uptake).
- **GCMC isotherms**: generic force fields capture the *shape* (e.g. HKUST-1's
  water step near RH 40–50%) but underestimate absolute uptake; tune the force
  field + lengthen sampling for quantitative values. Pressure↔RH via
  P_sat(298 K)=3169 Pa.
- **Polymorph ranking**: gaps <1–2 kJ/mol are within MLIP error — report
  "near-degenerate" rather than over-reading the sign; confirm ordering with
  dispersion-corrected DFT. Large unrelaxed→relaxed drops are X-ray H-position
  corrections (expected).

## Reference files

`fairchem-modal/`: `modal_app.py` (UMA), `gcmc_raspa.py` (RASPA GCMC),
`client.py` (endpoint client), `setup_models.py` (cache warmer),
`REPORT.md` (validated results), `DEPLOY.md` (endpoint usage).
