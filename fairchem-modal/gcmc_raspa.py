"""
GCMC water/CO2 isotherms in MOFs via RASPA2 on Modal
====================================================

This is the "macroscopic uptake" layer that complements modal_app.py:

  modal_app.py (UMA)  -> accurate per-site adsorption ENERGY + relaxation
  gcmc_raspa.py (RASPA)-> macroscopic LOADING / isotherm (the "吸水率" curve)

RASPA runs classical grand-canonical Monte Carlo with standard force fields
(framework: UFF/DREIDING; adsorbates: TIP water models, TraPPE CO2). It is the
mature, fast route to uptake isotherms. UMA is NOT used inside the MC loop here
(that would be research-grade MLIP-driven MC -- one UMA call per MC step).

Run:
  # first, verify RASPA installs and locate its data dir
  modal run gcmc_raspa.py::check_raspa
  # then an isotherm (pressures in Pa)
  modal run gcmc_raspa.py::isotherm --cif inputs/HKUST-1.cif --molecule tip5p \
      --temperature 298 --pressures 1e3,5e3,1e4,5e4,1e5

NOTE: force-field / molecule names below use RASPA's shipped examples. For
publication work you'd supply vetted force fields; treat defaults as a
working starting point, not a validated parameter set.
"""

import modal

# RASPA2 from conda-forge (ships the `simulate` binary + example FF/molecules).
image = (
    modal.Image.micromamba()
    .micromamba_install("raspa2", "python=3.11", channels=["conda-forge"])
    .pip_install("ase")
)

app = modal.App("raspa-gcmc", image=image)


def _raspa_dir():
    """Locate the RASPA data root (holds share/raspa2/{forcefield,molecules,...})."""
    import os, shutil
    sim = shutil.which("simulate")
    if not sim:
        raise RuntimeError("RASPA `simulate` binary not found on PATH")
    # binary at <prefix>/bin/simulate ; data at <prefix>/share/raspa2
    prefix = os.path.dirname(os.path.dirname(sim))
    return prefix, sim


# ---------------------------------------------------------------------------
# 0. Install / environment check -- run this FIRST
# ---------------------------------------------------------------------------
@app.function(timeout=600)
def _check():
    import os, glob
    prefix, sim = _raspa_dir()
    data = os.path.join(prefix, "share", "raspa2")
    ff = sorted(os.path.basename(p) for p in glob.glob(os.path.join(data, "forcefield", "*")))
    mols = sorted(os.path.basename(p) for p in glob.glob(os.path.join(data, "molecules", "*")))
    return {
        "simulate_bin": sim,
        "raspa_data": data,
        "data_exists": os.path.isdir(data),
        "forcefields": ff[:20],
        "molecule_sets": mols[:20],
    }


@app.local_entrypoint()
def check_raspa():
    import json
    print(json.dumps(_check.remote(), indent=2))


# ---------------------------------------------------------------------------
# 1. GCMC isotherm: loading vs pressure at fixed T
# ---------------------------------------------------------------------------
@app.function(timeout=7200)
def _isotherm(cif_text: str, framework_name: str, molecule: str,
              temperature: float, pressures: list, cycles: int,
              init_cycles: int, forcefield: str, unit_cells: str):
    import os, re, subprocess, glob, tempfile

    prefix, sim = _raspa_dir()
    data = os.path.join(prefix, "share", "raspa2")

    workroot = tempfile.mkdtemp(prefix="gcmc_")
    # RASPA looks up the framework CIF relative to RASPA_DIR/share/raspa2/structures/cif
    # or the working dir; we place it in the working dir and point RASPA_DIR at data.
    cif_path = os.path.join(workroot, f"{framework_name}.cif")
    with open(cif_path, "w") as f:
        f.write(cif_text)

    env = dict(os.environ, RASPA_DIR=prefix)

    def run_point(pressure):
        d = os.path.join(workroot, f"P_{pressure:.0f}")
        os.makedirs(d, exist_ok=True)
        # symlink CIF so RASPA finds it in cwd
        os.symlink(cif_path, os.path.join(d, f"{framework_name}.cif"))
        sim_input = f"""SimulationType                MonteCarlo
NumberOfCycles                {cycles}
NumberOfInitializationCycles  {init_cycles}
PrintEvery                    {max(1, cycles // 10)}
RestartFile                   no

Forcefield                    {forcefield}
CutOff                        12.0
ChargeMethod                  Ewald
EwaldPrecision                1e-6

Framework 0
FrameworkName {framework_name}
UnitCells {unit_cells}
ExternalTemperature {temperature}
ExternalPressure {pressure}

Component 0 MoleculeName             {molecule}
            MoleculeDefinition       TraPPE
            TranslationProbability   1.0
            RotationProbability      1.0
            ReinsertionProbability   1.0
            SwapProbability          1.0
            CreateNumberOfMolecules  0
"""
        with open(os.path.join(d, "simulation.input"), "w") as f:
            f.write(sim_input)
        proc = subprocess.run([sim], cwd=d, env=env, capture_output=True, text=True, timeout=3600)
        # parse the Output data file for absolute loading
        outs = glob.glob(os.path.join(d, "Output", "System_0", "*.data"))
        loading_molkg = loading_uc = None
        raw_tail = proc.stdout[-500:] + proc.stderr[-500:]
        if outs:
            txt = open(outs[0]).read()
            m1 = re.search(r"Average loading absolute \[mol/kg framework\]\s+([-\d.eE]+)", txt)
            m2 = re.search(r"Average loading absolute \[molecules/unit cell\]\s+([-\d.eE]+)", txt)
            if m1: loading_molkg = float(m1.group(1))
            if m2: loading_uc = float(m2.group(1))
            raw_tail = "parsed OK"
        return {"pressure_Pa": pressure, "loading_mol_per_kg": loading_molkg,
                "loading_molec_per_cell": loading_uc,
                "returncode": proc.returncode, "note": raw_tail}

    results = [run_point(p) for p in pressures]
    return {"framework": framework_name, "molecule": molecule,
            "temperature_K": temperature, "forcefield": forcefield,
            "isotherm": results}


@app.local_entrypoint()
def isotherm(cif: str, molecule: str = "CO2", temperature: float = 298.0,
             pressures: str = "1e3,1e4,1e5", cycles: int = 5000,
             init_cycles: int = 2000, forcefield: str = "GenericMOFs",
             unit_cells: str = "1 1 1"):
    import os, json
    with open(cif) as f:
        cif_text = f.read()
    fw = os.path.splitext(os.path.basename(cif))[0]
    plist = [float(x) for x in pressures.split(",")]
    print(f"GCMC: {molecule} in {fw} at {temperature} K, {len(plist)} pressure points ...")
    out = _isotherm.remote(cif_text, fw, molecule, temperature, plist,
                           cycles, init_cycles, forcefield, unit_cells)
    print(json.dumps(out, indent=2))
