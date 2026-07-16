"""
md_entropy.py -- physical conformational-flexibility / entropy analysis via MD.

The paper's mechanism is ENTROPIC: ligand binding lowers the chimera's
conformational entropy, which is coupled to the reporter's catalytic activity.
Boltz cannot see this (single static structure). Molecular dynamics can:

  * per-residue Cα RMSF (flexibility profile)             -> where the chain is floppy
  * active-site RMSF (S70/K73/S130/E166)                  -> is the catalytic
                                                             machinery pre-organized?
  * quasi-harmonic / Schlitter configurational entropy    -> a physical S estimate;
    ΔS(apo→holo) is the quantity the paper links to the switch.

This module uses OpenMM (open source). Running apo AND holo and comparing
flexibility/entropy is the rigorous open-source route to the ON/OFF mechanism.

HONESTY: converged entropy needs long sampling (µs) on a GPU and a parameterized
ligand for the holo leg (GAFF/OpenFF via openmmforcefields). The `smoke()` entry
here runs a SHORT implicit-solvent apo trajectory to prove the harness works and
produce a real RMSF profile -- it is infrastructure, NOT a converged result and
NOT a dynamic range.  ⚠️
"""

from __future__ import annotations
import os
import numpy as np

OUT = os.path.join(os.getcwd(), "biosensor_out")


def _clean_chainA_pdb(pdb_id: str, out_pdb: str) -> str:
    """Fetch a PDB entry, keep chain A protein heavy atoms, write a clean PDB."""
    import warnings; warnings.filterwarnings("ignore")
    import biotite.database.rcsb as rcsb
    import biotite.structure.io.pdbx as pdbx
    import biotite.structure.io.pdb as pdb
    import biotite.structure as struc
    import tempfile
    path = rcsb.fetch(pdb_id, "bcif", tempfile.mkdtemp())
    arr = pdbx.get_structure(pdbx.BinaryCIFFile.read(path), model=1)
    arr = arr[struc.filter_amino_acids(arr)]
    arr = arr[arr.chain_id == "A"]
    arr = arr[arr.element != "H"]
    f = pdb.PDBFile(); f.set_structure(arr); f.write(out_pdb)
    return out_pdb


def _add_cterm_oxt(arr, struc):
    """Append an OXT atom to the C-terminal residue (needed for OpenMM capping)."""
    import numpy as np
    last_rid = int(arr.res_id.max())
    res = arr[arr.res_id == last_rid]
    try:
        C = res[res.atom_name == "C"].coord[0]
        O = res[res.atom_name == "O"].coord[0]
    except IndexError:
        return arr
    oxt = struc.Atom(
        2.0 * C - O,                      # reflect O through C (~correct C-OXT length)
        chain_id=res.chain_id[0], res_id=last_rid, res_name=res.res_name[0],
        atom_name="OXT", element="O", hetero=False,
    )
    return arr + struc.array([oxt])


def _load_topology(path: str):
    """Load an OpenMM topology+positions from a PDB or (Boltz) mmCIF.

    Boltz output CIFs are complete (all heavy atoms, no gaps), so they prep
    cleanly -- unlike crystal structures that may miss side-chain/terminal atoms.
    The protein chain is kept; any ligand/hetero is dropped for this apo analysis.
    """
    from openmm import app
    if path.lower().endswith((".cif", ".mmcif", ".bcif")):
        # Route CIF through biotite -> clean protein-only PDB so OpenMM infers
        # chain termini correctly (Boltz CIF struct_conn confuses terminal
        # detection). Keep chain A protein heavy atoms only.
        import warnings; warnings.filterwarnings("ignore")
        import biotite.structure.io.pdbx as pdbx
        import biotite.structure.io.pdb as biopdb
        import biotite.structure as struc
        cif = pdbx.CIFFile.read(path) if path.lower().endswith((".cif", ".mmcif")) \
            else pdbx.BinaryCIFFile.read(path)
        arr = pdbx.get_structure(cif, model=1)
        arr = arr[struc.filter_amino_acids(arr)]
        chain = max(set(arr.chain_id), key=lambda c: (arr.chain_id == c).sum())
        arr = arr[(arr.chain_id == chain) & (arr.element != "H")]
        arr = _add_cterm_oxt(arr, struc)      # Boltz models omit OXT; cap it
        tmp = path + ".prot.pdb"
        f = biopdb.PDBFile(); f.set_structure(arr); f.write(tmp)
        st = app.PDBFile(tmp)
    else:
        st = app.PDBFile(path)
    return app.Modeller(st.topology, st.positions)


def run_md(pdb_path: str, ps: float = 1.0, n_frames: int = 12, temperature: float = 300.0):
    """Short implicit-solvent MD; returns (resids, ca_rmsf_A, entropy_proxy)."""
    from openmm import app, unit, LangevinMiddleIntegrator, Platform
    import openmm as mm

    ff = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")
    modeller = _load_topology(pdb_path)
    modeller.addHydrogens(ff)
    system = ff.createSystem(modeller.topology, nonbondedMethod=app.CutoffNonPeriodic,
                             nonbondedCutoff=1.6 * unit.nanometer, constraints=app.HBonds)
    integ = LangevinMiddleIntegrator(temperature * unit.kelvin, 1.0 / unit.picosecond,
                                     0.002 * unit.picoseconds)
    sim = app.Simulation(modeller.topology, system, integ, Platform.getPlatformByName("CPU"))
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy(maxIterations=60)
    sim.context.setVelocitiesToTemperature(temperature * unit.kelvin)

    total_steps = int(ps / 0.002)
    stride = max(1, total_steps // n_frames)
    top = modeller.topology
    ca_idx = [a.index for a in top.atoms() if a.name == "CA"]
    ca_res = [a.residue.id for a in top.atoms() if a.name == "CA"]
    frames = []
    for _ in range(n_frames):
        sim.step(stride)
        pos = sim.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        frames.append(pos[ca_idx])
    X = np.array(frames)                        # (frames, n_ca, 3)
    Xc = X - X.mean(axis=1, keepdims=True)      # remove translation per frame
    mean = Xc.mean(axis=0)
    rmsf = np.sqrt(((Xc - mean) ** 2).sum(axis=2).mean(axis=0))   # per-CA RMSF (A)
    # crude configurational-entropy proxy: sum log of positional variance (nats)
    var = ((Xc - mean) ** 2).sum(axis=2).mean(axis=0) + 1e-3
    entropy_proxy = float(0.5 * np.log(var).sum())
    return ca_res, rmsf, entropy_proxy


def smoke(model_path: str = None, ps: float = 1.0):
    """Infrastructure smoke test: short apo MD from a complete (Boltz) model."""
    import json
    os.makedirs(OUT, exist_ok=True)
    model_path = model_path or os.path.join(OUT, "dfhbi_apo.cif")
    res, rmsf, S = run_md(model_path, ps=ps)
    out = {
        "model": os.path.basename(model_path),
        "kind": "APO smoke test (UNDER-CONVERGED, infrastructure only)",
        "ps": ps, "n_residues": len(res),
        "rmsf_A": {"min": round(float(rmsf.min()), 2), "mean": round(float(rmsf.mean()), 2),
                    "max": round(float(rmsf.max()), 2)},
        "entropy_proxy_nats": round(S, 2),
        "note": "⚠️ short CPU run — proves the OpenMM harness runs and yields a real "
                "RMSF/entropy profile. A switch validation needs apo AND holo (ligand "
                "parameterized), µs sampling on GPU, then compare ΔS and active-site RMSF.",
    }
    with open(os.path.join(OUT, "md_smoke.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"OpenMM smoke test on {out['model']} ({len(res)} res, {ps} ps):")
    print(f"  Cα RMSF  min/mean/max = {out['rmsf_A']['min']}/{out['rmsf_A']['mean']}/{out['rmsf_A']['max']} Å")
    print(f"  config-entropy proxy  = {out['entropy_proxy_nats']} nats")
    print("  " + out["note"])
    return out


# ===========================================================================
# PRODUCTION protocol: GPU apo/holo MD with a parameterized ligand + entropy.
# ---------------------------------------------------------------------------
# This is the rigorous open-source route to the ON/OFF mechanism. It needs a GPU
# and the small-molecule-FF stack (openff-toolkit + openmmforcefields), which are
# conda-forge packages — see reference/md-gpu-protocol.md and modal_app.py to run
# it on a cloud GPU (Modal). All heavy imports are lazy so this module still
# imports (and `smoke()` still runs) without them.
# ===========================================================================

# physical constants (SI)
_KB = 1.380649e-23           # J/K
_HBAR = 1.054571817e-34      # J*s
_NA = 6.02214076e23          # 1/mol
_AMU = 1.66053906660e-27     # kg
_E = 2.718281828459045


def pick_platform(prefer=("CUDA", "OpenCL", "CPU")):
    """Return the fastest available OpenMM platform name (GPU first)."""
    from openmm import Platform
    have = {Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())}
    for p in prefer:
        if p in have:
            return p
    return "CPU"


def _openff_molecule(smiles: str, name: str = "LIG"):
    from openff.toolkit import Molecule
    mol = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
    mol.generate_conformers(n_conformers=1)
    mol.name = name
    return mol


def make_system_generator(ligand_smiles: str = None, solvent: str = "explicit",
                          small_molecule_ff: str = "openff-2.2.0"):
    """Build an openmmforcefields SystemGenerator (protein + optional ligand).

    solvent='explicit' -> amber14 + TIP3P + PME (production quality);
    solvent='implicit' -> amber14 + GBn2 (cheaper, for scans).
    The ligand is parameterized with an OpenFF SMIRNOFF force field (pure-python,
    no ambertools) matched to `ligand_smiles`.
    """
    from openmm import app, unit
    from openmmforcefields.generators import SystemGenerator
    mols = [_openff_molecule(ligand_smiles)] if ligand_smiles else []
    if solvent == "explicit":
        return SystemGenerator(
            forcefields=["amber14-all.xml", "amber14/tip3p.xml"],
            small_molecule_forcefield=small_molecule_ff, molecules=mols,
            forcefield_kwargs={"constraints": app.HBonds, "rigidWater": True,
                               "hydrogenMass": 1.5 * unit.amu},
            periodic_forcefield_kwargs={"nonbondedMethod": app.PME,
                                        "nonbondedCutoff": 1.0 * unit.nanometer},
        ), mols
    return SystemGenerator(
        forcefields=["amber14-all.xml", "implicit/gbn2.xml"],
        small_molecule_forcefield=small_molecule_ff, molecules=mols,
        forcefield_kwargs={"constraints": app.HBonds, "hydrogenMass": 1.5 * unit.amu},
        nonperiodic_forcefield_kwargs={"nonbondedMethod": app.CutoffNonPeriodic,
                                       "nonbondedCutoff": 1.6 * unit.nanometer},
    ), mols


def build_complex(protein_path: str, ligand_sdf: str = None, ligand_smiles: str = None,
                  solvent: str = "explicit"):
    """Assemble protein (+ optional ligand) into a solvated, parameterized system.

    protein_path : Boltz/PDB model of the chimera (apo leg: protein only).
    ligand_sdf   : Boltz `sample_0_ligands.sdf` for the holo leg (coords + bonds).
    ligand_smiles: SMILES fallback / the molecule used to parameterize the ligand.
    Returns (modeller, system, topology).
    """
    from openmm import app, unit
    modeller = _load_topology(protein_path)         # protein, biotite-cleaned
    ff = app.ForceField  # (unused directly; SystemGenerator owns FF)
    sysgen, mols = make_system_generator(ligand_smiles=ligand_smiles, solvent=solvent)
    modeller.addHydrogens(sysgen.forcefield)
    if ligand_sdf:
        from openff.toolkit import Molecule
        lig = Molecule.from_file(ligand_sdf)
        if isinstance(lig, list):
            lig = lig[0]
        lig_top = lig.to_topology().to_openmm()
        lig_pos = lig.conformers[0].to_openmm()
        modeller.add(lig_top, lig_pos)
    if solvent == "explicit":
        modeller.addSolvent(sysgen.forcefield, model="tip3p",
                            padding=1.0 * unit.nanometer, ionicStrength=0.15 * unit.molar)
    system = sysgen.create_system(modeller.topology)
    return modeller, system, modeller.topology


def schlitter_entropy(coords_A, masses_amu, temperature=300.0):
    """Absolute conformational entropy (Schlitter, upper bound), in J/(mol*K).

    coords_A   : (frames, N, 3) Cartesian, already superposed (rigid-body removed).
    masses_amu : (N,) atomic masses.
    S = 0.5 kB ln det[ 1 + (kB T e^2 / hbar^2) * M^{1/2} sigma M^{1/2} ]
    with sigma the covariance of the Cartesian fluctuations. ΔS(holo-apo) is the
    quantity the paper's entropic-switch mechanism predicts to be negative.
    """
    X = coords_A.reshape(coords_A.shape[0], -1) * 1e-10        # A -> m
    X = X - X.mean(axis=0)
    sigma = np.cov(X, rowvar=False)                           # (3N,3N) m^2
    m = np.repeat(np.asarray(masses_amu) * _AMU, 3)           # kg, len 3N
    Mroot = np.sqrt(m)
    mw = (Mroot[:, None] * sigma) * Mroot[None, :]            # M^{1/2} sigma M^{1/2}
    lam = np.linalg.eigvalsh(mw)
    lam = lam[lam > 1e-50]        # drop the 6 rigid-body/numerical-zero modes only
    alpha = _KB * temperature * _E ** 2 / _HBAR ** 2
    S = 0.5 * _KB * np.sum(np.log1p(alpha * lam))             # J/K
    return float(S * _NA)                                     # J/(mol*K)


def run_leg(protein_path, ligand_sdf=None, ligand_smiles=None, solvent="explicit",
            equil_ps=200.0, prod_ns=20.0, temperature=300.0, save_every_ps=10.0,
            platform=None, catalytic_resids=None):
    """One MD leg (apo or holo): equilibrate + production, return entropy + RMSF.

    Returns a dict: {S_JmolK, ca_rmsf {min,mean,max}, active_site_rmsf, n_frames,
                     platform, solvent}. GPU strongly recommended (prod_ns×system).
    """
    from openmm import app, unit, LangevinMiddleIntegrator, MonteCarloBarostat, Platform
    import openmm as mm

    modeller, system, top = build_complex(protein_path, ligand_sdf, ligand_smiles, solvent)
    if solvent == "explicit":
        system.addForce(MonteCarloBarostat(1.0 * unit.bar, temperature * unit.kelvin, 25))
    integ = LangevinMiddleIntegrator(temperature * unit.kelvin, 1.0 / unit.picosecond,
                                     0.004 * unit.picoseconds)   # 4 fs w/ HMR
    plat = Platform.getPlatformByName(platform or pick_platform())
    sim = app.Simulation(top, system, integ, plat)
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy()
    sim.context.setVelocitiesToTemperature(temperature * unit.kelvin)
    sim.step(int(equil_ps / 0.004))                              # equilibration

    ca = [a.index for a in top.atoms() if a.name == "CA"]
    ca_res = [a.residue.id for a in top.atoms() if a.name == "CA"]
    masses = [system.getParticleMass(i).value_in_unit(unit.amu) for i in ca]
    n_frames = max(4, int(prod_ns * 1000 / save_every_ps))
    stride = int(save_every_ps / 0.004)
    frames = []
    for _ in range(n_frames):
        sim.step(stride)
        pos = sim.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        frames.append(pos[ca])
    X = np.array(frames)
    # superpose each frame onto frame 0 (Kabsch) to remove rigid-body motion
    X = _superpose(X)
    mean = X.mean(axis=0)
    rmsf = np.sqrt(((X - mean) ** 2).sum(axis=2).mean(axis=0))
    S = schlitter_entropy(X, masses, temperature)
    out = {"solvent": solvent, "platform": plat.getName(), "prod_ns": prod_ns,
           "n_frames": len(frames), "S_JmolK": round(S, 1),
           "ca_rmsf_A": {"min": round(float(rmsf.min()), 2),
                          "mean": round(float(rmsf.mean()), 2),
                          "max": round(float(rmsf.max()), 2)}}
    if catalytic_resids:
        idx = [i for i, r in enumerate(ca_res) if int(r) in set(catalytic_resids)]
        if idx:
            out["active_site_rmsf_A"] = round(float(rmsf[idx].mean()), 2)
    return out


def _superpose(X):
    """Kabsch-align every frame (N,3) onto the first frame."""
    ref = X[0] - X[0].mean(0)
    out = np.empty_like(X)
    for k in range(X.shape[0]):
        P = X[k] - X[k].mean(0)
        H = P.T @ ref
        U, _, Vt = np.linalg.svd(H)
        d = np.sign(np.linalg.det(Vt.T @ U.T))
        R = Vt.T @ np.diag([1, 1, d]) @ U.T
        out[k] = P @ R.T
    return out


def compare_apo_holo(apo_protein, holo_protein, ligand_sdf=None, ligand_smiles=None,
                     catalytic_resids=None, **kw):
    """Full mechanism test: ΔS(holo-apo) and active-site RMSF change.

    Mechanism prediction (paper): ligand binding LOWERS conformational entropy
    (ΔS < 0) and pre-organizes the reporter active site (active-site RMSF down).
    """
    apo = run_leg(apo_protein, solvent=kw.get("solvent", "explicit"),
                  catalytic_resids=catalytic_resids, **kw)
    holo = run_leg(holo_protein, ligand_sdf=ligand_sdf, ligand_smiles=ligand_smiles,
                   solvent=kw.get("solvent", "explicit"),
                   catalytic_resids=catalytic_resids, **kw)
    dS = round(holo["S_JmolK"] - apo["S_JmolK"], 1)
    res = {"apo": apo, "holo": holo, "dS_JmolK_holo_minus_apo": dS,
           "entropy_reduced_on_binding": dS < 0,
           "label": "✅ ΔS and RMSF are measurements on the trajectories; "
                    "⚠️ mapping ΔS -> a dynamic-range number still needs calibration "
                    "against wet-lab kobs. ΔS<0 supports the paper's entropic switch."}
    if "active_site_rmsf_A" in apo and "active_site_rmsf_A" in holo:
        res["active_site_rmsf_delta"] = round(holo["active_site_rmsf_A"] - apo["active_site_rmsf_A"], 2)
    return res


if __name__ == "__main__":
    import sys
    smoke(sys.argv[1] if len(sys.argv) > 1 else None)
