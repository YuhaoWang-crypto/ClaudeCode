"""
Client for the deployed FAIRChem/RASPA Modal endpoints
======================================================

The compute functions are DEPLOYED (persistent) on Modal apps `fairchem-uma`
and `raspa-gcmc`. They scale to zero (no idle cost) and spin up a GPU on call.

Call them from ANY machine that has your Modal creds
(MODAL_TOKEN_ID / MODAL_TOKEN_SECRET) -- no need to clone the compute code:

    from client import drug_loading, water_isotherm
    print(drug_loading("inputs/ZIF-8.cif", ["drugs/loading/ibuprofen.sdf"]))
    print(water_isotherm("inputs/HKUST-1.cif", rh=(0.2, 0.5, 0.8)))

Or as a quick check:  python client.py smoke
"""

import os
import sys
import modal

UMA_APP = "fairchem-uma"
RASPA_APP = "raspa-gcmc"


def _fn(app, name):
    return modal.Function.from_name(app, name)


def _read(path):
    with open(path) as f:
        return f.read()


# --- UMA (fairchem-uma) -----------------------------------------------------
def smoke():
    """Cheapest end-to-end check: single-point energy of H2O on a GPU."""
    return _fn(UMA_APP, "_smoke").remote()


def mof_adsorption(cif_path, adsorbate="H2O"):
    """Single-site adsorption energy of a small molecule in a MOF."""
    return _fn(UMA_APP, "_mof_adsorption").remote(_read(cif_path), adsorbate)


def drug_loading(mof_cif_path, drug_sdf_paths):
    """Host-guest interaction energy of drug(s) inside a MOF pore."""
    drugs = {os.path.splitext(os.path.basename(p))[0]: _read(p) for p in drug_sdf_paths}
    return _fn(UMA_APP, "_drug_loading").remote(_read(mof_cif_path), drugs)


def polymorph_rank(cif_paths, atoms_per_molecule):
    """Rank molecular-crystal polymorphs by lattice energy per molecule."""
    cifs = {os.path.basename(p): _read(p) for p in cif_paths}
    return _fn(UMA_APP, "_polymorph_rank").remote(cifs, atoms_per_molecule)


def relax_material(cif_path):
    """Variable-cell relaxation + energy of an inorganic crystal."""
    return _fn(UMA_APP, "_relax_material").remote(_read(cif_path))


# --- RASPA GCMC (raspa-gcmc) ------------------------------------------------
def gas_isotherm(cif_path, molecule="CO2", temperature=298.0,
                 pressures=(1e3, 1e4, 5e4, 1e5), cycles=2000, init_cycles=1000,
                 forcefield="ExampleMOFsForceField", unit_cells="1 1 1"):
    """GCMC adsorption isotherm (loading vs pressure) of a gas in a MOF."""
    fw = os.path.splitext(os.path.basename(cif_path))[0]
    return _fn(RASPA_APP, "_isotherm").remote(
        _read(cif_path), fw, molecule, temperature, list(pressures),
        cycles, init_cycles, forcefield, unit_cells)


def water_isotherm(cif_path, temperature=298.0, rh=(0.1, 0.3, 0.5, 0.7, 0.9),
                   cycles=2500, init_cycles=1500, unit_cells="1 1 1"):
    """Water uptake vs relative humidity (TIP4P + Ewald) in a MOF."""
    fw = os.path.splitext(os.path.basename(cif_path))[0]
    return _fn(RASPA_APP, "_water_isotherm").remote(
        _read(cif_path), fw, temperature, list(rh), cycles, init_cycles, unit_cells)


if __name__ == "__main__":
    import json
    cmd = sys.argv[1] if len(sys.argv) > 1 else "smoke"
    fns = {
        "smoke": lambda: smoke(),
        "mof_adsorption": lambda: mof_adsorption("inputs/MOF-5.cif", "H2O"),
        "drug_loading": lambda: drug_loading("inputs/ZIF-8.cif",
            ["drugs/loading/ibuprofen.sdf", "drugs/loading/5-fluorouracil.sdf"]),
        "water_isotherm": lambda: water_isotherm("inputs/HKUST-1.cif", rh=(0.3, 0.6, 0.9)),
    }
    if cmd not in fns:
        print(f"usage: python client.py [{'|'.join(fns)}]"); sys.exit(1)
    print(json.dumps(fns[cmd](), indent=2))
