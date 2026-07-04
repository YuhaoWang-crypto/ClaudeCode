import modal, json

app = modal.App("is621-md")
vol = modal.Volume.from_name("is621-md-vol", create_if_missing=True)
img = (modal.Image.micromamba(python_version="3.11")
       .micromamba_install("ambertools", "openmm", "cudatoolkit=11.8", "mdtraj",
                           "parmed", "pdbfixer", "numpy",
                           channels=["conda-forge"]))

DDGB_NOTE = "single-trajectory MM-GBSA (OBC2), initial estimate"


@app.function(image=img, gpu="A10G", timeout=7200, volumes={"/vol": vol})
def run_variant(complex_pdb: str, variant: str, prod_ns: float = 2.0):
    import os, subprocess, numpy as np
    from pdbfixer import PDBFixer
    import openmm as mm
    import openmm.app as app_mm
    from openmm import unit
    wd = f"/work/{variant}"; os.makedirs(wd, exist_ok=True); os.chdir(wd)
    log = {"variant": variant}

    # ---- 1. PDBFixer: mutate (if needed), complete heavy atoms, add H ----
    open("in.pdb", "w").write(complex_pdb)
    fixer = PDBFixer(filename="in.pdb")
    if variant == "A248K":
        fixer.applyMutations(["ALA-248-LYS"], "A")
    fixer.findMissingResidues(); fixer.missingResidues = {}   # don't model absent loops
    fixer.findMissingAtoms(); fixer.addMissingAtoms()
    # NOTE: no addMissingHydrogens -> let tleap protonate (avoids HIS H-name clashes)
    with open("fixed.pdb", "w") as fh:
        app_mm.PDBFile.writeFile(fixer.topology, fixer.positions, fh)

    # ---- 2. clean for tleap: drop H, strip 5'-terminal DNA phosphates, per-chain TER ----
    DNA = {"DA", "DT", "DG", "DC"}
    first_resi = {}   # first residue number seen per chain (= 5' end for DNA)
    lines, last = [], None
    for l in open("fixed.pdb"):
        if not l.startswith("ATOM"):
            continue
        elem = l[76:78].strip()
        if elem == "H" or l[12:16].strip().startswith("H"):
            continue  # let tleap add hydrogens
        ch = l[21]; resn = l[17:20].strip(); resi = int(l[22:26]); an = l[12:16].strip()
        first_resi.setdefault(ch, resi)
        if resn in DNA and resi == first_resi[ch] and an in ("P", "OP1", "OP2"):
            continue  # remove 5'-terminal phosphate (tleap DX5 template has none)
        if last is not None and ch != last:
            lines.append("TER\n")
        lines.append(l); last = ch
    lines.append("TER\nEND\n")
    open("clean.pdb", "w").writelines(lines)

    leapin = f"""source leaprc.protein.ff19SB
source leaprc.DNA.OL21
source leaprc.water.tip3p
set default PBRadii mbondi2
m = loadpdb clean.pdb
savepdb m dry.pdb
saveamberparm m dry.prmtop dry.inpcrd
solvatebox m TIP3PBOX 10.0
addionsrand m Na+ 0
addionsrand m Cl- 0
saveamberparm m sol.prmtop sol.inpcrd
quit
"""
    open("leap.in", "w").write(leapin)
    r = subprocess.run(["tleap", "-f", "leap.in"], capture_output=True, text=True)
    log["tleap_tail"] = r.stdout[-1500:]
    if not os.path.exists("sol.prmtop"):
        log["error"] = "tleap failed to build solvated system"
        return log

    # ---- 3. OpenMM MD (explicit solvent, PME, HMR 4 fs) ----
    prm = app_mm.AmberPrmtopFile("sol.prmtop")
    inp = app_mm.AmberInpcrdFile("sol.inpcrd")
    system = prm.createSystem(nonbondedMethod=app_mm.PME, nonbondedCutoff=1.0*unit.nanometer,
                              constraints=app_mm.HBonds, hydrogenMass=3.0*unit.amu)
    system.addForce(mm.MonteCarloBarostat(1*unit.bar, 300*unit.kelvin))
    integ = mm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.004*unit.picoseconds)
    sim = app_mm.Simulation(prm.topology, system, integ)
    sim.context.setPositions(inp.positions)
    if inp.boxVectors is not None:
        sim.context.setPeriodicBoxVectors(*inp.boxVectors)
    sim.minimizeEnergy(maxIterations=5000)
    sim.context.setVelocitiesToTemperature(300*unit.kelvin)
    sim.step(50000)   # 200 ps equilibration
    n_frames = 100
    steps_per = int(prod_ns * 1e6 / 4 / n_frames)  # 4 fs steps
    sim.reporters.append(app_mm.DCDReporter("prod.dcd", steps_per))
    sim.reporters.append(app_mm.StateDataReporter("md.csv", steps_per*10, step=True,
                         potentialEnergy=True, temperature=True))
    sim.step(steps_per * n_frames)
    log["prod_ns"] = prod_ns; log["frames"] = n_frames

    # ---- 4. MM-GBSA (single trajectory, OBC2) on stripped snapshots ----
    import parmed as pmd, mdtraj as md
    dry = pmd.load_file("dry.prmtop")
    # receptor = protein (strip nucleic), ligand = DNA (strip protein)
    prot_mask = "@%" ; # build masks via residue names
    rec = pmd.load_file("dry.prmtop"); rec.strip(":DA,DT,DG,DC")
    lig = pmd.load_file("dry.prmtop"); lig.strip("!(:DA,DT,DG,DC)")
    def mk(struct):
        s = struct.createSystem(implicitSolvent=app_mm.OBC2, nonbondedMethod=app_mm.NoCutoff)
        integ = mm.VerletIntegrator(1*unit.femtosecond)
        c = mm.Context(s, integ, mm.Platform.getPlatformByName("CPU"))
        return c
    c_cpx, c_rec, c_lig = mk(dry), mk(rec), mk(lig)
    n_prot = sum(1 for a in dry.atoms if a.residue.name not in ("DA","DT","DG","DC"))
    # align frames: dry complex atom order = protein then DNA (as loaded); mdtraj trajectory
    traj = md.load("prod.dcd", top="sol.prmtop")
    # keep only solute (protein + DNA) to match dry complex atom order
    sel = traj.top.select("protein or nucleic")
    traj = traj.atom_slice(sel)
    assert traj.n_atoms == len(dry.atoms), (traj.n_atoms, len(dry.atoms))
    def energy(c, xyz_nm):
        c.setPositions(xyz_nm)
        return c.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
    dg = []
    for i in range(0, traj.n_frames, max(1, traj.n_frames//60)):
        xyz = traj.xyz[i]  # nm
        E_c = energy(c_cpx, xyz)
        E_r = energy(c_rec, xyz[:n_prot])
        E_l = energy(c_lig, xyz[n_prot:])
        dg.append(E_c - E_r - E_l)
    dg = np.array(dg)
    log["dG_bind_kcal"] = float(dg.mean())
    log["dG_sem_kcal"] = float(dg.std()/np.sqrt(len(dg)))
    log["dg_per_frame"] = dg.tolist()
    log["n_snapshots"] = len(dg)
    log["method"] = DDGB_NOTE
    # persist to volume so results survive local-client disconnects
    small = {k: v for k, v in log.items() if k != "tleap_tail"}
    with open(f"/vol/{variant}.json", "w") as f:
        json.dump(small, f)
    vol.commit()
    return log


@app.local_entrypoint()
def main():
    pdb = open("mdqm/complex_wt.pdb").read()
    # run both variants concurrently; results are written to the volume by each call
    h_wt = run_variant.spawn(pdb, "WT", 1.5)
    h_mut = run_variant.spawn(pdb, "A248K", 1.5)
    print("spawned WT:", h_wt.object_id, " A248K:", h_mut.object_id)
    print("results will be written to volume is621-md-vol as WT.json / A248K.json")
