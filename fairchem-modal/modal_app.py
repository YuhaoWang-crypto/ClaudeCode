"""
FAIRChem (UMA) on Modal GPUs
============================

Runs Meta FAIR's UMA universal interatomic potential (from HuggingFace
facebook/UMA) on Modal's GPUs for three workloads:

  1. MOF adsorption     -> UMA task head "odac"  (Open Direct Air Capture:
                            MOFs + H2O/CO2 adsorption)
  2. Drug polymorphs    -> UMA task head "omc"   (organic molecular crystals)
  3. Inorganic materials-> UMA task head "omat"  (bulk crystals)

Why Modal: the Claude Code container is CPU-only, so all heavy inference is
offloaded here. Auth uses MODAL_TOKEN_ID / MODAL_TOKEN_SECRET (env vars).

Prereqs (see README.md):
  - `modal secret create huggingface HF_TOKEN=hf_xxx`  (UMA is gated; your
    token must have accepted the facebook/UMA license)

Run:
  modal run modal_app.py::smoke_test
  modal run modal_app.py::mof_adsorption --cif path/to/mof.cif
  modal run modal_app.py::polymorph_rank --cif-dir path/to/polymorphs/
  modal run modal_app.py::relax_material --cif path/to/crystal.cif

NOTE: This is a scaffold. The UMA model tag ("uma-s-1") and task-head names
are set to the current fairchem-core v2 API; confirm against your access on
first run (smoke_test is the cheapest way to validate the whole chain).
"""

import modal

# --- Image: fairchem-core (ships UMA support) + ASE/pymatgen for I/O ---------
image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        "torch",
        "fairchem-core>=2.0.0",
        "ase",
        "pymatgen",
        "huggingface_hub",
    )
)

app = modal.App("fairchem-uma", image=image)

# HF token needed to pull the gated facebook/UMA weights. Create with:
#   modal secret create huggingface HF_TOKEN=hf_xxxxx
hf_secret = modal.Secret.from_name("huggingface")

# Cache HF downloads so the ~GB of weights aren't re-pulled every run.
model_cache = modal.Volume.from_name("hf-cache", create_if_missing=True)
CACHE_DIR = "/root/.cache/huggingface"

GPU = "A10G"          # cheap; bump to "A100" for large MOF supercells
MODEL_TAG = "uma-s-1" # confirm the exact tag you have access to


def _load_calc(task_name: str):
    """Load a UMA calculator bound to a given task head."""
    from fairchem.core import pretrained_mlip, FAIRChemCalculator

    predictor = pretrained_mlip.get_predict_unit(MODEL_TAG, device="cuda")
    return FAIRChemCalculator(predictor, task_name=task_name)


# ---------------------------------------------------------------------------
# 0. Smoke test: cheapest end-to-end validation of GPU + gated model access
# ---------------------------------------------------------------------------
@app.function(gpu=GPU, secrets=[hf_secret], volumes={CACHE_DIR: model_cache}, timeout=1800)
def _smoke():
    from ase.build import molecule

    calc = _load_calc("omol")  # isolated-molecule head
    atoms = molecule("H2O")
    atoms.calc = calc
    e = atoms.get_potential_energy()
    model_cache.commit()
    return {"system": "H2O", "energy_eV": float(e), "gpu": GPU, "model": MODEL_TAG}


@app.local_entrypoint()
def smoke_test():
    print("Loading UMA on Modal GPU and computing a single-point energy...")
    out = _smoke.remote()
    print(out)
    print("OK -- chain works: Modal GPU + gated UMA weights + fairchem inference.")


# ---------------------------------------------------------------------------
# 1. MOF adsorption energy (task head: odac)
#    E_ads = E(MOF+adsorbate) - E(MOF) - E(adsorbate)
#    A negative value means favorable adsorption.
# ---------------------------------------------------------------------------
@app.function(gpu=GPU, secrets=[hf_secret], volumes={CACHE_DIR: model_cache}, timeout=3600)
def _mof_adsorption(cif_text: str, adsorbate: str = "H2O", fmax: float = 0.05):
    import io
    from ase.io import read
    from ase.build import molecule
    from ase.optimize import BFGS
    from ase.constraints import FixAtoms

    calc = _load_calc("odac")

    # --- relaxed MOF framework ---
    mof = read(io.StringIO(cif_text), format="cif")
    mof.calc = calc
    BFGS(mof, logfile=None).run(fmax=fmax, steps=200)
    e_mof = mof.get_potential_energy()

    # --- isolated adsorbate (H2O / CO2) ---
    ads = molecule(adsorbate)
    ads.calc = calc
    BFGS(ads, logfile=None).run(fmax=fmax, steps=200)
    e_ads = ads.get_potential_energy()

    # --- adsorbate placed inside the pore, framework fixed during relax ---
    combo = mof.copy()
    ads_placed = ads.copy()
    ads_placed.translate(mof.get_center_of_mass() - ads_placed.get_center_of_mass())
    combo += ads_placed
    combo.calc = calc
    combo.set_constraint(FixAtoms(indices=list(range(len(mof)))))
    BFGS(combo, logfile=None).run(fmax=fmax, steps=300)
    e_combo = combo.get_potential_energy()

    model_cache.commit()
    e_binding = e_combo - e_mof - e_ads
    return {
        "adsorbate": adsorbate,
        "E_MOF_eV": float(e_mof),
        "E_adsorbate_eV": float(e_ads),
        "E_complex_eV": float(e_combo),
        "E_adsorption_eV": float(e_binding),
        "note": "Naive center placement + fixed framework. For screening/isotherms "
                "sample many sites (GCMC/RASPA) -- see README.",
    }


@app.local_entrypoint()
def mof_adsorption(cif: str, adsorbate: str = "H2O"):
    with open(cif) as f:
        cif_text = f.read()
    print(f"Computing {adsorbate} adsorption energy on {cif} via UMA/odac ...")
    print(_mof_adsorption.remote(cif_text, adsorbate))


# ---------------------------------------------------------------------------
# 2. Drug polymorph energy ranking (task head: omc)
#    Relax each candidate crystal, rank by lattice energy per molecule.
#    CAVEAT: polymorph energy gaps are often <1 kJ/mol; treat as a screen,
#    confirm the top candidates with dispersion-corrected DFT.
# ---------------------------------------------------------------------------
@app.function(gpu=GPU, secrets=[hf_secret], volumes={CACHE_DIR: model_cache}, timeout=3600)
def _polymorph_rank(cifs: dict, z_per_cell: dict | None = None, fmax: float = 0.03):
    import io
    from ase.io import read
    from ase.optimize import BFGS
    from ase.filters import FrechetCellFilter

    calc = _load_calc("omc")
    z_per_cell = z_per_cell or {}
    results = []
    for name, cif_text in cifs.items():
        atoms = read(io.StringIO(cif_text), format="cif")
        atoms.calc = calc
        # relax atoms AND cell (variable-cell) for molecular crystals
        BFGS(FrechetCellFilter(atoms), logfile=None).run(fmax=fmax, steps=400)
        e = atoms.get_potential_energy()
        z = z_per_cell.get(name, 1)
        results.append({"name": name, "E_total_eV": float(e),
                        "Z": z, "E_per_molecule_eV": float(e) / z})
    model_cache.commit()
    results.sort(key=lambda r: r["E_per_molecule_eV"])
    lowest = results[0]["E_per_molecule_eV"]
    for r in results:
        r["dE_meV_per_mol"] = (r["E_per_molecule_eV"] - lowest) * 1000.0
    return {"ranking": results,
            "warning": "MLIP polymorph gaps are approximate; verify top hits with DFT-D."}


@app.local_entrypoint()
def polymorph_rank(cif_dir: str):
    import os, json
    cifs = {}
    for fn in sorted(os.listdir(cif_dir)):
        if fn.lower().endswith(".cif"):
            with open(os.path.join(cif_dir, fn)) as f:
                cifs[fn] = f.read()
    if not cifs:
        print(f"No .cif files found in {cif_dir}")
        return
    print(f"Ranking {len(cifs)} polymorph candidates via UMA/omc ...")
    print(json.dumps(_polymorph_rank.remote(cifs), indent=2))


# ---------------------------------------------------------------------------
# 3. Inorganic material relaxation + energy (task head: omat)
# ---------------------------------------------------------------------------
@app.function(gpu=GPU, secrets=[hf_secret], volumes={CACHE_DIR: model_cache}, timeout=3600)
def _relax_material(cif_text: str, fmax: float = 0.02):
    import io
    from ase.io import read
    from ase.optimize import BFGS
    from ase.filters import FrechetCellFilter

    calc = _load_calc("omat")
    atoms = read(io.StringIO(cif_text), format="cif")
    atoms.calc = calc
    e0 = atoms.get_potential_energy()
    BFGS(FrechetCellFilter(atoms), logfile=None).run(fmax=fmax, steps=400)
    model_cache.commit()
    return {
        "n_atoms": len(atoms),
        "E_initial_eV": float(e0),
        "E_relaxed_eV": float(atoms.get_potential_energy()),
        "E_per_atom_eV": float(atoms.get_potential_energy()) / len(atoms),
    }


@app.local_entrypoint()
def relax_material(cif: str):
    with open(cif) as f:
        cif_text = f.read()
    print(f"Relaxing {cif} via UMA/omat ...")
    print(_relax_material.remote(cif_text))
