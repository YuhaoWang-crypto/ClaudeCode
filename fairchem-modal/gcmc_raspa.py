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
    import os, glob, subprocess
    prefix, sim = _raspa_dir()
    # RASPA data is marked by pseudo_atoms.def / force_field*.def and a molecules dir.
    # conda-forge may not place it at <prefix>/share/raspa2 -- search for it.
    hits = subprocess.run(
        ["find", "/opt/conda", "-name", "pseudo_atoms.def", "-o", "-name", "force_field_mixing_rules.def"],
        capture_output=True, text=True, timeout=120).stdout.strip().splitlines()
    # directories that contain example molecule defs (CO2.def, tip5p.def, ...)
    moldirs = subprocess.run(
        ["find", "/opt/conda", "-type", "d", "-name", "molecules"],
        capture_output=True, text=True, timeout=120).stdout.strip().splitlines()
    ffdirs = subprocess.run(
        ["find", "/opt/conda", "-type", "d", "-name", "forcefield"],
        capture_output=True, text=True, timeout=120).stdout.strip().splitlines()
    # guess RASPA_DIR = parent of share/raspa2
    guess = None
    for d in ffdirs:
        # d = <root>/share/raspa2/forcefield/<name>  OR  <root>/share/raspa2/forcefield
        parts = d.split(os.sep)
        if "share" in parts:
            guess = os.sep.join(parts[: parts.index("share")])
            break
    DATA = "/opt/conda/share/raspa"
    ff_names = sorted(os.path.basename(p) for p in glob.glob(os.path.join(DATA, "forcefield", "*")))
    mol_sets = sorted(os.path.basename(p) for p in glob.glob(os.path.join(DATA, "molecules", "*")))
    # locate CO2 / water model defs and report which molecule-set holds them
    wanted = {}
    for target in ("CO2.def", "Tip5p.def", "tip5p.def", "Tip4p.def", "tip4p.def", "Water.def", "H2O.def"):
        found = subprocess.run(["find", os.path.join(DATA, "molecules"), "-name", target],
                               capture_output=True, text=True, timeout=60).stdout.strip().splitlines()
        if found:
            wanted[target] = [p.replace(DATA + "/molecules/", "") for p in found]
    return {
        "simulate_bin": sim,
        "RASPA_DIR": "/opt/conda",
        "DATA": DATA,
        "forcefields": ff_names,
        "molecule_sets": mol_sets,
        "adsorbate_defs_found": wanted,
    }


@app.local_entrypoint()
def check_raspa():
    import json
    print(json.dumps(_check.remote(), indent=2))


@app.function(timeout=300)
def _inspect(forcefield: str = "ExampleMOFsForceField"):
    import os, glob
    DATA = "/opt/conda/share/raspa"
    md = os.path.join(DATA, "molecules", "ExampleDefinitions")
    mols = sorted(os.path.basename(p)[:-4] for p in glob.glob(os.path.join(md, "*.def")))
    ffd = os.path.join(DATA, "forcefield", forcefield)
    out = {"molecule_defs": mols, "forcefield": forcefield}
    for fn in ("pseudo_atoms.def", "force_field_mixing_rules.def"):
        p = os.path.join(ffd, fn)
        out[fn] = open(p).read() if os.path.exists(p) else "MISSING"
    return out


@app.local_entrypoint()
def inspect(forcefield: str = "ExampleMOFsForceField"):
    out = _inspect.remote(forcefield)
    print("MOLECULE DEFS:", ", ".join(out["molecule_defs"]))
    print("\n===== pseudo_atoms.def =====\n" + out["pseudo_atoms.def"][:4000])
    print("\n===== force_field_mixing_rules.def =====\n" + out["force_field_mixing_rules.def"][:4000])


# ---------------------------------------------------------------------------
# 1. GCMC isotherm: loading vs pressure at fixed T
# ---------------------------------------------------------------------------
@app.function(timeout=7200)
def _isotherm(cif_text: str, framework_name: str, molecule: str,
              temperature: float, pressures: list, cycles: int,
              init_cycles: int, forcefield: str, unit_cells: str):
    import os, re, subprocess, glob, tempfile

    prefix, sim = _raspa_dir()  # prefix = /opt/conda; RASPA reads share/raspa/ under it

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
            MoleculeDefinition       ExampleDefinitions
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


# ---------------------------------------------------------------------------
# 2b. WATER isotherm: loading vs relative humidity (the "吸水率" curve)
#     Ships a hand-authored TIP4P water model + matching pseudo-atoms/LJ,
#     reads framework partial charges from the CIF, uses Ewald electrostatics.
# ---------------------------------------------------------------------------

# TIP4P (Jorgensen 1983): O + 2H + massless M-site on the H-O-H bisector.
_TIP4P_DEF = """# critical constants: Temperature [T], Pressure [Pa], and Acentric factor [-]
647.0
22064000.0
0.3449
#Number Of Atoms
 4
# Number of groups
1
# water-group
rigid
# number of atoms
4
# atomic positions
0 Ow   0.0        0.0        0.0
1 Hw   0.756950   0.0        0.585882
2 Hw  -0.756950   0.0        0.585882
3 Mw   0.0        0.0        0.15
# Chiral Bond BondDipoles Bend UB InvBend Torsion ImpTorsion Bond/Bond Stretch/Bend Bend/Bend Stretch/Torsion Bend/Torsion IntraVDW IntraCoulomb
      0    3           0    0  0       0       0          0         0            0         0               0            0        0            0
# Bond stretch: atom n1-n2, type, parameters
0 1 RIGID_BOND
0 2 RIGID_BOND
0 3 RIGID_BOND
# Number of config moves
0
"""

# pseudo-atoms for TIP4P water (Ow neutral+LJ, Hw +0.52, Mw massless -1.04)
_WATER_PSEUDO = [
    "Ow         yes     O     O     0           15.9994     0.0      0.0          1.0      0.7    0            0           relative           0",
    "Hw         yes     H     H     0           1.00794     0.52     0.0          1.0      0.7    0            0           relative           0",
    "Mw         no      O     -     0           0.0        -1.04     0.0          1.0      0.7    0            0           relative           0",
]
_WATER_MIX = [
    "Ow             lennard-jones    78.0      3.15365      // TIP4P, Jorgensen et al. J.Chem.Phys. 1983",
    "Hw             none",
    "Mw             none",
]
P_SAT_298 = 3169.0  # water saturation vapour pressure at 298.15 K [Pa]


def _build_water_raspa_dir():
    """Clone the conda RASPA data into a writable tree with TIP4P water added."""
    import os, shutil, tempfile
    src = "/opt/conda/share/raspa"
    root = tempfile.mkdtemp(prefix="raspa_water_")
    raspa = os.path.join(root, "share", "raspa")
    os.makedirs(raspa)
    for name in os.listdir(src):
        if name not in ("forcefield", "molecules"):
            os.symlink(os.path.join(src, name), os.path.join(raspa, name))

    # forcefield: symlink all, then a real modified copy "WaterMOF"
    ffroot = os.path.join(raspa, "forcefield"); os.makedirs(ffroot)
    for ff in os.listdir(os.path.join(src, "forcefield")):
        os.symlink(os.path.join(src, "forcefield", ff), os.path.join(ffroot, ff))
    os.remove(os.path.join(ffroot, "WaterMOF")) if os.path.islink(os.path.join(ffroot, "WaterMOF")) else None
    dst = os.path.join(ffroot, "WaterMOF")
    shutil.copytree(os.path.join(src, "forcefield", "ExampleMOFsForceField"), dst)

    # append water pseudo-atoms (bump the count)
    pa = os.path.join(dst, "pseudo_atoms.def")
    lines = open(pa).read().splitlines()
    for i, ln in enumerate(lines):
        if ln.strip().isdigit():
            lines[i] = str(int(ln.strip()) + len(_WATER_PSEUDO)); break
    open(pa, "w").write("\n".join(lines + _WATER_PSEUDO) + "\n")

    # append water LJ before the trailing "general mixing rule" block (bump count)
    mr = os.path.join(dst, "force_field_mixing_rules.def")
    lines = open(mr).read().splitlines()
    for i, ln in enumerate(lines):
        if ln.strip().isdigit():
            lines[i] = str(int(ln.strip()) + len(_WATER_MIX)); break
    ins = next(i for i, ln in enumerate(lines) if ln.startswith("# general mixing rule"))
    lines[ins:ins] = _WATER_MIX
    open(mr, "w").write("\n".join(lines) + "\n")

    # molecules: symlink existing sets, add a WaterDefs set with Tip4p.def
    mroot = os.path.join(raspa, "molecules"); os.makedirs(mroot)
    for ms in os.listdir(os.path.join(src, "molecules")):
        os.symlink(os.path.join(src, "molecules", ms), os.path.join(mroot, ms))
    wd = os.path.join(mroot, "WaterDefs"); os.makedirs(wd)
    open(os.path.join(wd, "Tip4p.def"), "w").write(_TIP4P_DEF)
    return root


@app.function(timeout=10800)
def _water_isotherm(cif_text: str, framework_name: str, temperature: float,
                    rh_points: list, cycles: int, init_cycles: int, unit_cells: str):
    import os, re, subprocess, glob, tempfile

    root = _build_water_raspa_dir()
    _, sim = _raspa_dir()
    env = dict(os.environ, RASPA_DIR=root)

    workroot = tempfile.mkdtemp(prefix="wateriso_")
    cif_path = os.path.join(workroot, f"{framework_name}.cif")
    open(cif_path, "w").write(cif_text)
    p_sat = P_SAT_298  # 298 K reference; fine for the demo temperature

    def run_point(rh):
        pressure = rh * p_sat
        d = os.path.join(workroot, f"RH_{int(rh*100)}"); os.makedirs(d, exist_ok=True)
        os.symlink(cif_path, os.path.join(d, f"{framework_name}.cif"))
        sim_input = f"""SimulationType                MonteCarlo
NumberOfCycles                {cycles}
NumberOfInitializationCycles  {init_cycles}
PrintEvery                    {max(1, cycles // 10)}
RestartFile                   no

Forcefield                    WaterMOF
CutOff                        12.0
ChargeMethod                  Ewald
EwaldPrecision                1e-6
UseChargesFromCIFFile         yes

Framework 0
FrameworkName {framework_name}
UnitCells {unit_cells}
ExternalTemperature {temperature}
ExternalPressure {pressure}

Component 0 MoleculeName             Tip4p
            MoleculeDefinition       WaterDefs
            TranslationProbability   1.0
            RotationProbability      1.0
            ReinsertionProbability   1.0
            SwapProbability          1.0
            CreateNumberOfMolecules  0
"""
        open(os.path.join(d, "simulation.input"), "w").write(sim_input)
        proc = subprocess.run([sim], cwd=d, env=env, capture_output=True, text=True, timeout=7200)
        outs = glob.glob(os.path.join(d, "Output", "System_0", "*.data"))
        molkg = uc = None
        note = (proc.stdout[-600:] + proc.stderr[-600:]) or "no output"
        if outs:
            txt = open(outs[0]).read()
            m1 = re.search(r"Average loading absolute \[mol/kg framework\]\s+([-\d.eE]+)", txt)
            m2 = re.search(r"Average loading absolute \[molecules/unit cell\]\s+([-\d.eE]+)", txt)
            if m1: molkg = float(m1.group(1))
            if m2: uc = float(m2.group(1))
            note = "parsed OK" if m1 else note
        return {"RH": rh, "pressure_Pa": round(pressure, 1),
                "loading_mol_per_kg": molkg, "loading_molec_per_cell": uc,
                "returncode": proc.returncode, "note": note}

    return {"framework": framework_name, "molecule": "TIP4P water",
            "temperature_K": temperature, "P_sat_Pa": p_sat,
            "isotherm": [run_point(rh) for rh in rh_points]}


@app.local_entrypoint()
def water_isotherm(cif: str, temperature: float = 298.0,
                   rh: str = "0.1,0.2,0.3,0.5,0.7,0.9",
                   cycles: int = 3000, init_cycles: int = 2000,
                   unit_cells: str = "1 1 1"):
    import os, json
    with open(cif) as f:
        cif_text = f.read()
    fw = os.path.splitext(os.path.basename(cif))[0]
    rh_points = [float(x) for x in rh.split(",")]
    print(f"WATER GCMC: {fw} at {temperature} K, RH points {rh_points} (TIP4P) ...")
    out = _water_isotherm.remote(cif_text, fw, temperature, rh_points,
                                 cycles, init_cycles, unit_cells)
    print(json.dumps(out, indent=2))


@app.local_entrypoint()
def isotherm(cif: str, molecule: str = "CO2", temperature: float = 298.0,
             pressures: str = "1e3,1e4,1e5", cycles: int = 5000,
             init_cycles: int = 2000, forcefield: str = "ExampleMOFsForceField",
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
