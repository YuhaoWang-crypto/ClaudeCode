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
MODEL_TAG = "uma-s-1p2"  # latest UMA-small; other options: uma-m-1p1, esen-sm-full-odac25


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


def _parse_sdf(sdf_text):
    """Minimal SDF V2000 -> (symbols, positions). Robust to PubChem 3D records."""
    lines = sdf_text.splitlines()
    natoms = int(lines[3][:3])
    symbols, pos = [], []
    for ln in lines[4:4 + natoms]:
        t = ln.split()
        pos.append([float(t[0]), float(t[1]), float(t[2])])
        symbols.append(t[3])
    return symbols, pos


def _largest_void(atoms, n=12):
    """Coarse grid search for the point in the cell farthest from any framework
    atom (minimum-image) -- a cheap proxy for the largest pore centre."""
    import numpy as np
    cell = np.array(atoms.cell)
    frac_atoms = atoms.get_scaled_positions()
    best, best_d = None, -1.0
    grid = [(i + 0.5) / n for i in range(n)]
    for fx in grid:
        for fy in grid:
            for fz in grid:
                fp = np.array([fx, fy, fz])
                df = frac_atoms - fp
                df -= np.round(df)                    # minimum image
                d = np.linalg.norm(df @ cell, axis=1).min()
                if d > best_d:
                    best_d, best = d, fp
    return best @ cell, best_d


# ---------------------------------------------------------------------------
# 1b. Drug loading in a MOF: host-guest interaction energy
#     E_int = E(MOF+drug) - E(MOF) - E(drug).  More negative = stronger
#     framework affinity -> higher loading / slower release.
# ---------------------------------------------------------------------------
@app.function(gpu=GPU, secrets=[hf_secret], volumes={CACHE_DIR: model_cache}, timeout=5400)
def _drug_loading(mof_cif: str, drugs: dict, fmax: float = 0.05):
    import io
    import numpy as np
    from ase import Atoms
    from ase.io import read
    from ase.optimize import BFGS
    from ase.constraints import FixAtoms

    calc = _load_calc("odac")  # MOF host + molecular guest; one head = consistent refs

    mof = read(io.StringIO(mof_cif), format="cif")
    mof.calc = calc
    BFGS(mof, logfile=None).run(fmax=fmax, steps=200)   # relax framework (cell fixed)
    e_mof = mof.get_potential_energy()
    void_cart, void_r = _largest_void(mof)

    results = []
    for name, sdf in drugs.items():
        symbols, positions = _parse_sdf(sdf)
        drug = Atoms(symbols=symbols, positions=positions)

        # isolated drug reference (same cell, vacuum), relaxed
        iso = drug.copy(); iso.set_cell(mof.cell); iso.set_pbc(True); iso.center()
        iso.calc = calc
        BFGS(iso, logfile=None).run(fmax=fmax, steps=300)
        e_drug = iso.get_potential_energy()

        # place drug COM at the largest void, relax guest with framework fixed
        g = drug.copy(); g.translate(void_cart - g.get_positions().mean(axis=0))
        combo = mof.copy(); combo += g
        combo.set_pbc(True); combo.calc = calc
        combo.set_constraint(FixAtoms(indices=list(range(len(mof)))))
        BFGS(combo, logfile=None).run(fmax=fmax, steps=400)
        e_combo = combo.get_potential_energy()

        e_int = e_combo - e_mof - e_drug
        results.append({"drug": name, "n_atoms": len(drug),
                        "E_MOF_eV": float(e_mof), "E_drug_eV": float(e_drug),
                        "E_complex_eV": float(e_combo),
                        "E_interaction_eV": float(e_int),
                        "E_interaction_kJ_per_mol": float(e_int) * 96.485})
    model_cache.commit()
    results.sort(key=lambda r: r["E_interaction_eV"])  # most favorable first
    return {"void_radius_A": float(void_r), "ranking": results,
            "note": "E_int<0 = favorable loading. Single centred pose + fixed framework; "
                    "for capacity/selectivity sample multiple poses and see GCMC."}


@app.local_entrypoint()
def drug_loading(mof_cif: str, drug_dir: str):
    import os, json
    with open(mof_cif) as f:
        cif_text = f.read()
    drugs = {}
    for fn in sorted(os.listdir(drug_dir)):
        if fn.lower().endswith(".sdf"):
            with open(os.path.join(drug_dir, fn)) as f:
                drugs[os.path.splitext(fn)[0]] = f.read()
    if not drugs:
        print(f"No .sdf drug files in {drug_dir}"); return
    print(f"Drug loading in {os.path.basename(mof_cif)} for {list(drugs)} via UMA/odac ...")
    print(json.dumps(_drug_loading.remote(cif_text, drugs), indent=2))


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
def _polymorph_rank(cifs: dict, atoms_per_molecule: int, fmax: float = 0.03):
    import io
    from ase.io import read
    from ase.optimize import BFGS
    from ase.filters import FrechetCellFilter

    calc = _load_calc("omc")
    results = []
    for name, cif_text in cifs.items():
        atoms = read(io.StringIO(cif_text), format="cif")
        atoms.calc = calc
        e0 = atoms.get_potential_energy()  # as-published (unrelaxed) energy
        # relax atoms AND cell (variable-cell) for molecular crystals
        BFGS(FrechetCellFilter(atoms), logfile=None).run(fmax=fmax, steps=400)
        e = atoms.get_potential_energy()
        n = len(atoms)
        # infer molecules-per-cell (Z) from the actual atom count -> robust to
        # whether the CIF reader expanded symmetry or not.
        z = max(1, round(n / atoms_per_molecule))
        results.append({"name": name, "n_atoms": n, "Z": z,
                        "E_unrelaxed_eV": float(e0),
                        "E_relaxed_eV": float(e),
                        "E_per_molecule_eV": float(e) / z})
    model_cache.commit()
    results.sort(key=lambda r: r["E_per_molecule_eV"])
    lowest = results[0]["E_per_molecule_eV"]
    for r in results:
        r["dE_kJ_per_mol"] = (r["E_per_molecule_eV"] - lowest) * 96.485  # eV -> kJ/mol
    return {"ranking": results,
            "warning": "MLIP polymorph gaps are approximate (often <2 kJ/mol, near MLIP "
                       "accuracy limit); verify ordering with dispersion-corrected DFT."}


@app.local_entrypoint()
def polymorph_rank(cif_dir: str, atoms_per_molecule: int = 0):
    import os, json
    if atoms_per_molecule <= 0:
        raise SystemExit("Pass --atoms-per-molecule (e.g. paracetamol C8H9NO2 = 20)")
    cifs = {}
    for fn in sorted(os.listdir(cif_dir)):
        if fn.lower().endswith(".cif"):
            with open(os.path.join(cif_dir, fn)) as f:
                cifs[fn] = f.read()
    if not cifs:
        print(f"No .cif files found in {cif_dir}")
        return
    print(f"Ranking {len(cifs)} polymorph candidates via UMA/omc ...")
    print(json.dumps(_polymorph_rank.remote(cifs, atoms_per_molecule), indent=2))


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
