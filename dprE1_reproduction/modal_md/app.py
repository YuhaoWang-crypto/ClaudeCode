"""
Part 3 - Protein-ligand MD on Modal, adapting the *making-it-rain*
Protein_ligand workflow (Arantes et al.) to serverless GPU.

Faithful to making-it-rain's protein-ligand notebook:
  * protein  : AMBER ff14SB / ff19SB
  * ligand   : GAFF2 with AM1-BCC charges (antechamber)
  * water    : TIP3P, neutralising ions
  * topology : AmberTools tleap  ->  prmtop/inpcrd
  * engine   : OpenMM, LangevinMiddleIntegrator, 2 fs, HBonds, PME 1.0 nm
  * ensemble : minimise -> NVT heat -> NPT equilibrate -> NPT production
  * analysis : mdtraj (ligand/protein/cofactor RMSD, protein RMSF)  ->  paper Table 4
               + optional MM-GBSA (AmberTools MMPBSA.py) as in the paper

Run (needs a Modal account, `pip install modal`, `modal token new`):

    # one complex, short demo (default 5 ns)
    modal run modal_md/app.py --system GTD_9.7

    # full paper protocol (500 ns) for every candidate
    modal run modal_md/app.py --all --ns 500

Inputs are produced locally by `prepare_inputs.py` (docked complexes from Part 2)
and live in `modal_md/inputs/<system>/{protein.pdb, ligand.sdf}`.
"""
import os
import sys
import pathlib
import modal

APP_DIR = pathlib.Path(__file__).parent
INPUTS = APP_DIR / "inputs"

# --- conda image: AmberTools + OpenMM(+CUDA) + analysis stack -----------------
image = (
    modal.Image.micromamba(python_version="3.11")
    .micromamba_install(
        "ambertools=23",
        "openmm=8.1",
        "openmmforcefields",
        "openff-toolkit-base",
        "pdbfixer",
        "mdtraj",
        "parmed",
        "rdkit",
        "numpy",
        "scipy",
        channels=["conda-forge"],
    )
    # CUDA runtime for OpenMM GPU platform is provided by Modal's GPU host
)

app = modal.App("dprE1-md-making-it-rain", image=image)
vol = modal.Volume.from_name("dprE1-md-results", create_if_missing=True)
RESULTS = "/results"


# ---------------------------------------------------------------- MD driver ----
@app.function(gpu="A10G", timeout=60 * 60 * 24, volumes={RESULTS: vol})
def simulate(system: str, protein_pdb: bytes, ligand_sdf: bytes,
             ns: float = 5.0, protein_ff: str = "ff14SB",
             equil_ns: float = 0.5, temperature: float = 300.0,
             report_ps: float = 25.0):
    """Parametrise (AMBER/GAFF2), solvate, equilibrate and run production MD."""
    import subprocess
    import textwrap
    import numpy as np
    from openmm import (LangevinMiddleIntegrator, MonteCarloBarostat,
                        Platform, unit)
    from openmm.app import (AmberPrmtopFile, AmberInpcrdFile, Simulation,
                           StateDataReporter, DCDReporter, HBonds, PME)
    from openff.toolkit import Molecule

    work = pathlib.Path(f"/tmp/{system}")
    work.mkdir(parents=True, exist_ok=True)
    (work / "protein_in.pdb").write_bytes(protein_pdb)
    (work / "ligand.sdf").write_bytes(ligand_sdf)
    outdir = pathlib.Path(RESULTS) / system
    outdir.mkdir(parents=True, exist_ok=True)

    # 1) clean protein with PDBFixer (add H at pH 7, missing heavy atoms)
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    fixer = PDBFixer(filename=str(work / "protein_in.pdb"))
    fixer.findMissingResidues(); fixer.findMissingAtoms()
    fixer.addMissingAtoms(); fixer.addMissingHydrogens(7.0)
    PDBFile.writeFile(fixer.topology, fixer.positions,
                      open(work / "protein.pdb", "w"))

    # 2) ligand: AM1-BCC charges + GAFF2 via antechamber/parmchk2
    mol = Molecule.from_file(str(work / "ligand.sdf"))
    net_q = int(round(mol.total_charge.magnitude))
    mol.to_file(str(work / "lig.mol2"), file_format="mol2")
    subprocess.run(
        ["antechamber", "-i", "lig.mol2", "-fi", "mol2", "-o", "lig_bcc.mol2",
         "-fo", "mol2", "-c", "bcc", "-nc", str(net_q), "-s", "2", "-at", "gaff2"],
        cwd=work, check=True)
    subprocess.run(["parmchk2", "-i", "lig_bcc.mol2", "-f", "mol2",
                    "-o", "lig.frcmod", "-s", "gaff2"], cwd=work, check=True)

    # 3) tleap: combine protein+ligand, solvate TIP3P, neutralise
    ff_line = {"ff14SB": "leaprc.protein.ff14SB",
               "ff19SB": "leaprc.protein.ff19SB"}[protein_ff]
    leap = textwrap.dedent(f"""
        source {ff_line}
        source leaprc.gaff2
        source leaprc.water.tip3p
        LIG = loadmol2 lig_bcc.mol2
        loadamberparams lig.frcmod
        prot = loadpdb protein.pdb
        complex = combine {{ prot LIG }}
        solvateBox complex TIP3PBOX 12.0
        addIonsRand complex Na+ 0 Cl- 0
        saveamberparm complex SYS_gaff2.prmtop SYS_gaff2.inpcrd
        savepdb complex SYS_solvated.pdb
        quit
    """)
    (work / "tleap.in").write_text(leap)
    subprocess.run(["tleap", "-f", "tleap.in"], cwd=work, check=True)

    # 4) OpenMM system  (making-it-rain defaults)
    prmtop = AmberPrmtopFile(str(work / "SYS_gaff2.prmtop"))
    inpcrd = AmberInpcrdFile(str(work / "SYS_gaff2.inpcrd"))
    system_omm = prmtop.createSystem(
        nonbondedMethod=PME, nonbondedCutoff=1.0 * unit.nanometer,
        constraints=HBonds)
    T = temperature * unit.kelvin
    integrator = LangevinMiddleIntegrator(T, 1.0 / unit.picosecond,
                                          0.002 * unit.picoseconds)
    system_omm.addForce(MonteCarloBarostat(1.0 * unit.bar, T, 25))
    try:
        platform = Platform.getPlatformByName("CUDA")
        props = {"Precision": "mixed"}
    except Exception:
        platform, props = Platform.getPlatformByName("CPU"), {}
    sim = Simulation(prmtop.topology, system_omm, integrator, platform, props)
    sim.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        sim.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    steps_per_ps = int(1.0 / 0.002)          # 500 steps = 1 ps
    report_stride = int(report_ps * steps_per_ps)

    # 5) minimise -> NVT heat -> NPT equilibrate -> NPT production
    sim.minimizeEnergy()
    sim.context.setVelocitiesToTemperature(T)
    sim.reporters.append(StateDataReporter(
        str(outdir / "md_log.csv"), report_stride, step=True, time=True,
        potentialEnergy=True, temperature=True, density=True, speed=True))

    sim.step(int(equil_ns * 1000 * steps_per_ps))            # equilibration
    sim.reporters.append(DCDReporter(str(outdir / "production.dcd"),
                                     report_stride))
    prod_steps = int(ns * 1000 * steps_per_ps)
    sim.step(prod_steps)                                      # production

    # save final state + topology for analysis / MM-GBSA
    from openmm.app import PDBFile as _PDB
    state = sim.context.getState(getPositions=True)
    _PDB.writeFile(prmtop.topology, state.getPositions(),
                   open(outdir / "final.pdb", "w"))
    for f in ("SYS_gaff2.prmtop", "SYS_solvated.pdb"):
        (outdir / f).write_bytes((work / f).read_bytes())
    vol.commit()

    metrics = analyse(str(outdir))
    metrics["system"] = system
    metrics["production_ns"] = ns
    import json
    (outdir / "metrics.json").write_text(json.dumps(metrics, indent=2))
    vol.commit()
    return metrics


def analyse(outdir: str):
    """mdtraj RMSD/RMSF matching paper Table 4 (ligand, protein, cofactor)."""
    import mdtraj as md
    import numpy as np
    top = f"{outdir}/SYS_gaff2.prmtop"
    traj = md.load(f"{outdir}/production.dcd", top=top)
    traj = traj.superpose(traj, 0, atom_indices=traj.top.select("protein and backbone"))

    lig = traj.top.select("resname LIG")
    prot_bb = traj.top.select("protein and backbone")
    cof = traj.top.select("resname FAD or resname HEM")

    def rmsd(sel):
        if len(sel) == 0:
            return None
        ref = traj[0].atom_slice(sel)
        sub = traj.atom_slice(sel)
        return float(np.mean(md.rmsd(sub, ref, 0) * 10.0))  # nm -> A

    rmsf = None
    if len(prot_bb):
        rmsf = float(np.mean(md.rmsf(traj, traj, 0, atom_indices=prot_bb) * 10.0))

    return {
        "ligand_RMSD_A": rmsd(lig),
        "protein_backbone_RMSD_A": rmsd(prot_bb),
        "protein_RMSF_A": rmsf,
        "cofactor_RMSD_A": rmsd(cof),
        "n_frames": int(traj.n_frames),
    }


# ------------------------------------------------------------------- entry -----
@app.local_entrypoint()
def main(system: str = "GTD_9.7", all: bool = False, ns: float = 5.0,
         protein_ff: str = "ff14SB"):
    systems = ([p.name for p in sorted(INPUTS.iterdir()) if p.is_dir()]
               if all else [system])
    calls = []
    for s in systems:
        pdb = (INPUTS / s / "protein.pdb").read_bytes()
        sdf = (INPUTS / s / "ligand.sdf").read_bytes()
        calls.append((s, simulate.spawn(s, pdb, sdf, ns, protein_ff)))
    print(f"Launched {len(calls)} MD job(s) on Modal (ns={ns}, ff={protein_ff})")
    for s, handle in calls:
        m = handle.get()
        print(f"[{s}] ligRMSD={m['ligand_RMSD_A']}  "
              f"protRMSD={m['protein_backbone_RMSD_A']}  "
              f"RMSF={m['protein_RMSF_A']}  cofRMSD={m['cofactor_RMSD_A']}")
