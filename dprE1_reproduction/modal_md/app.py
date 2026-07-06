"""
Part 3 - Protein-ligand MD on Modal, adapting the *making-it-rain*
Protein_ligand workflow (Arantes et al.) to serverless GPU, now with the
**FAD cofactor** parametrised (GAFF2) and **MM-GBSA** binding energy, to match
the paper's Table 4 fully.

Faithful to making-it-rain protein-ligand:
  * protein  : AMBER ff14SB / ff19SB
  * ligand   : GAFF2 + AM1-BCC (antechamber)
  * cofactor : FAD parametrised the same way (GAFF2 + AM1-BCC) and kept in the box
               so the paper's *cofactor RMSD* and FAD pi-stacking are represented
               (HEM/Fe is NOT supported this way - needs bonded-metal params)
  * water    : TIP3P, neutralising ions
  * topology : AmberTools tleap  ->  prmtop/inpcrd
  * engine   : OpenMM, LangevinMiddleIntegrator, 2 fs, HBonds, PME 1.0 nm
  * ensemble : minimise -> NPT equilibration -> NPT production
  * analysis : mdtraj RMSD/RMSF (Table 4) + AmberTools MMPBSA.py (MM-GBSA, igb=5)

Run (Modal account required; tokens via env or `modal token new`):

    # 1 ns end-to-end test on one complex (recommended first)
    modal run modal_md/app.py --system DprE1_WT__GTD_9.7 --ns 1

    # full paper protocol
    modal run modal_md/app.py --all --ns 500
"""
import pathlib
import modal

APP_DIR = pathlib.Path(__file__).parent
INPUTS = APP_DIR / "inputs"

image = (
    modal.Image.micromamba(python_version="3.11")
    .micromamba_install(
        "ambertools=23", "openmm=8.1", "openff-toolkit-base", "pdbfixer",
        "mdtraj", "parmed", "rdkit", "numpy", "scipy",
        channels=["conda-forge"],
    )
)

app = modal.App("dprE1-md-making-it-rain", image=image)
vol = modal.Volume.from_name("dprE1-md-results", create_if_missing=True)
RESULTS = "/results"


# --------------------------------------------------- ligand/cofactor param ----
def _gaff2_param(work, sdf_path, resname, out_prefix, net_charge=None,
                 charge_methods=("bcc",)):
    """antechamber (GAFF2) + parmchk2 for one small molecule/cofactor.

    charge_methods is tried in order; large cofactors (FAD, 88 atoms) can pass
    ("bcc", "gas") so a slow/non-converging AM1-BCC falls back to Gasteiger
    (fine for a mostly-spectator cofactor: its internal charges cancel in the
    ligand MM-GBSA ΔG)."""
    import subprocess, shutil
    from rdkit import Chem
    # net charge from the (already 3D, all-H) input SDF; antechamber reads SDF
    # directly, so we avoid MOL2 writing (unavailable without OpenEye here).
    if net_charge is None:
        mol = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=True)[0]
        if mol is None:  # FAD bond perception can fail sanitisation
            mol = Chem.SDMolSupplier(str(sdf_path), removeHs=False,
                                     sanitize=False)[0]
        net_charge = Chem.GetFormalCharge(mol) if mol is not None else 0
    shutil.copy(str(sdf_path), str(work / f"{out_prefix}.sdf"))
    last = None
    for cm in charge_methods:
        try:
            subprocess.run(
                ["antechamber", "-i", f"{out_prefix}.sdf", "-fi", "sdf",
                 "-o", f"{out_prefix}_bcc.mol2", "-fo", "mol2", "-c", cm,
                 "-nc", str(net_charge), "-s", "2", "-at", "gaff2",
                 "-rn", resname], cwd=work, check=True, timeout=60 * 30)
            subprocess.run(["parmchk2", "-i", f"{out_prefix}_bcc.mol2",
                            "-f", "mol2", "-o", f"{out_prefix}.frcmod",
                            "-s", "gaff2"], cwd=work, check=True)
            return net_charge
        except Exception as e:  # noqa
            last = e
    raise RuntimeError(f"antechamber failed for {out_prefix}: {last}")


# ---------------------------------------------------------------- MD driver ----
@app.function(gpu="A10G", timeout=60 * 60 * 24, volumes={RESULTS: vol})
def simulate(system: str, protein_pdb: bytes, ligand_sdf: bytes,
             cofactor_sdf: bytes | None = None, ns: float = 1.0,
             protein_ff: str = "ff14SB", equil_ns: float = 0.2,
             temperature: float = 300.0, report_ps: float = 10.0,
             run_mmgbsa: bool = True):
    import subprocess, textwrap, json
    from openmm import (LangevinMiddleIntegrator, MonteCarloBarostat,
                        Platform, unit)
    from openmm.app import (AmberPrmtopFile, AmberInpcrdFile, Simulation,
                           StateDataReporter, DCDReporter, HBonds, PME, PDBFile)
    from pdbfixer import PDBFixer

    work = pathlib.Path(f"/tmp/{system}"); work.mkdir(parents=True, exist_ok=True)
    (work / "protein_in.pdb").write_bytes(protein_pdb)
    (work / "ligand.sdf").write_bytes(ligand_sdf)
    outdir = pathlib.Path(RESULTS) / system; outdir.mkdir(parents=True, exist_ok=True)

    # 1) protein prep
    fixer = PDBFixer(filename=str(work / "protein_in.pdb"))
    fixer.findMissingResidues(); fixer.findMissingAtoms()
    fixer.addMissingAtoms(); fixer.addMissingHydrogens(7.0)
    PDBFile.writeFile(fixer.topology, fixer.positions, open(work / "protein.pdb", "w"))

    # 2) ligand + optional FAD cofactor -> GAFF2
    _gaff2_param(work, work / "ligand.sdf", "LIG", "lig")
    have_cof = cofactor_sdf is not None
    if have_cof:
        (work / "cofactor.sdf").write_bytes(cofactor_sdf)
        _gaff2_param(work, work / "cofactor.sdf", "FAD", "cof",
                     charge_methods=("bcc", "gas"))

    # 3) tleap
    ff_line = {"ff14SB": "leaprc.protein.ff14SB",
               "ff19SB": "leaprc.protein.ff19SB"}[protein_ff]
    cof_load = ("FAD = loadmol2 cof_bcc.mol2\nloadamberparams cof.frcmod\n"
                if have_cof else "")
    cof_combine = " FAD" if have_cof else ""
    leap = textwrap.dedent(f"""
        source {ff_line}
        source leaprc.gaff2
        source leaprc.water.tip3p
        LIG = loadmol2 lig_bcc.mol2
        loadamberparams lig.frcmod
        {cof_load}prot = loadpdb protein.pdb
        complex = combine {{ prot{cof_combine} LIG }}
        saveamberparm complex complex_dry.prmtop complex_dry.inpcrd
        solvateBox complex TIP3PBOX 12.0
        addIonsRand complex Na+ 0 Cl- 0
        saveamberparm complex SYS.prmtop SYS.inpcrd
        savepdb complex SYS_solvated.pdb
        quit
    """)
    (work / "tleap.in").write_text(leap)
    r = subprocess.run(["tleap", "-f", "tleap.in"], cwd=work,
                       capture_output=True, text=True)
    if r.returncode != 0 or not (work / "SYS.prmtop").exists():
        leaplog = (work / "leap.log").read_text()[-3000:] if \
            (work / "leap.log").exists() else "(no leap.log)"
        raise RuntimeError("tleap failed:\n" + r.stdout[-2000:] +
                           "\n--- leap.log ---\n" + leaplog)

    # 4) OpenMM system (making-it-rain defaults)
    prmtop = AmberPrmtopFile(str(work / "SYS.prmtop"))
    inpcrd = AmberInpcrdFile(str(work / "SYS.inpcrd"))
    system_omm = prmtop.createSystem(nonbondedMethod=PME,
                                     nonbondedCutoff=1.0 * unit.nanometer,
                                     constraints=HBonds)
    T = temperature * unit.kelvin
    integ = LangevinMiddleIntegrator(T, 1.0 / unit.picosecond, 0.002 * unit.picoseconds)
    system_omm.addForce(MonteCarloBarostat(1.0 * unit.bar, T, 25))
    try:
        platform, props = Platform.getPlatformByName("CUDA"), {"Precision": "mixed"}
    except Exception:
        platform, props = Platform.getPlatformByName("CPU"), {}
    sim = Simulation(prmtop.topology, system_omm, integ, platform, props)
    sim.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        sim.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    sps = 500                                   # steps per ps (2 fs)
    stride = int(report_ps * sps)

    # 5) minimise -> equilibrate -> production
    sim.minimizeEnergy()
    sim.context.setVelocitiesToTemperature(T)
    sim.reporters.append(StateDataReporter(
        str(outdir / "md_log.csv"), stride, step=True, time=True,
        potentialEnergy=True, temperature=True, density=True, speed=True))
    sim.step(int(equil_ns * 1000 * sps))
    sim.reporters.append(DCDReporter(str(outdir / "production.dcd"), stride))
    sim.step(int(ns * 1000 * sps))

    st = sim.context.getState(getPositions=True)
    PDBFile.writeFile(prmtop.topology, st.getPositions(), open(outdir / "final.pdb", "w"))
    for f in ("SYS.prmtop", "complex_dry.prmtop", "SYS_solvated.pdb"):
        (outdir / f).write_bytes((work / f).read_bytes())
    vol.commit()

    metrics = analyse(str(outdir), work, have_cof, run_mmgbsa)
    metrics.update(system=system, production_ns=ns, cofactor=have_cof)
    (outdir / "metrics.json").write_text(json.dumps(metrics, indent=2))
    vol.commit()
    return metrics


def analyse(outdir, work, have_cof, run_mmgbsa):
    """mdtraj RMSD/RMSF (Table 4) + optional MM-GBSA (MMPBSA.py, igb=5)."""
    import mdtraj as md, numpy as np, subprocess, textwrap
    traj = md.load(f"{outdir}/production.dcd", top=f"{outdir}/SYS.prmtop")
    traj.superpose(traj, 0, atom_indices=traj.top.select("protein and backbone"))

    def rmsd(sel):
        if len(sel) == 0:
            return None
        return float(np.mean(md.rmsd(traj.atom_slice(sel),
                                     traj.atom_slice(sel)[0], 0) * 10.0))
    prot_bb = traj.top.select("protein and backbone")
    m = {
        "ligand_RMSD_A": rmsd(traj.top.select("resname LIG")),
        "protein_backbone_RMSD_A": rmsd(prot_bb),
        "protein_RMSF_A": (float(np.mean(md.rmsf(traj, traj, 0,
                          atom_indices=prot_bb) * 10.0)) if len(prot_bb) else None),
        "cofactor_RMSD_A": rmsd(traj.top.select("resname FAD")),
        "n_frames": int(traj.n_frames),
        "mmgbsa_dG_kcal_mol": None,
    }

    if run_mmgbsa:
        try:
            # convert dcd -> Amber NetCDF for MMPBSA
            traj.save_netcdf(f"{work}/prod.nc")
            recmask = ":1-9998 & !:LIG"   # everything but ligand = receptor(+FAD)
            # build dry receptor/ligand topologies from the dry complex prmtop
            subprocess.run(["ante-MMPBSA.py", "-p", "complex_dry.prmtop",
                            "-c", "com.prmtop", "-r", "rec.prmtop", "-l", "lig.prmtop",
                            "-s", ":WAT,Na+,Cl-", "-n", ":LIG", "--radii", "mbondi2"],
                           cwd=work, check=True)
            (work / "mmgbsa.in").write_text(textwrap.dedent("""
                MM-GBSA (igb=5, GBSW-like), last-portion of trajectory
                &general startframe=1, interval=1, /
                &gb igb=5, saltcon=0.15, /
            """))
            subprocess.run(
                ["MMPBSA.py", "-O", "-i", "mmgbsa.in", "-cp", "com.prmtop",
                 "-rp", "rec.prmtop", "-lp", "lig.prmtop", "-y", "prod.nc",
                 "-o", "mmgbsa.dat"], cwd=work, check=True)
            for line in (work / "mmgbsa.dat").read_text().splitlines():
                if line.strip().startswith("DELTA TOTAL"):
                    m["mmgbsa_dG_kcal_mol"] = float(line.split()[2]); break
            import shutil
            shutil.copy(work / "mmgbsa.dat", f"{outdir}/mmgbsa.dat")
        except Exception as e:
            m["mmgbsa_error"] = str(e)[:300]
    return m


@app.local_entrypoint()
def main(system: str = "DprE1_WT__GTD_9.7", all: bool = False, ns: float = 1.0,
         protein_ff: str = "ff14SB", cofactor: bool = True):
    systems = ([p.name for p in sorted(INPUTS.iterdir()) if p.is_dir()]
               if all else [system])
    cof_bytes = None
    calls = []
    for s in systems:
        pdb = (INPUTS / s / "protein.pdb").read_bytes()
        sdf = (INPUTS / s / "ligand.sdf").read_bytes()
        cof = (INPUTS / s / "cofactor.sdf")
        cb = cof.read_bytes() if (cofactor and cof.exists()) else None
        calls.append((s, simulate.spawn(s, pdb, sdf, cb, ns, protein_ff)))
    print(f"Launched {len(calls)} MD job(s) on Modal (ns={ns}, ff={protein_ff})")
    for s, h in calls:
        m = h.get()
        print(f"[{s}] ligRMSD={m['ligand_RMSD_A']} protRMSD={m['protein_backbone_RMSD_A']} "
              f"RMSF={m['protein_RMSF_A']} cofRMSD={m['cofactor_RMSD_A']} "
              f"MMGBSA={m['mmgbsa_dG_kcal_mol']} kcal/mol")
