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


if __name__ == "__main__":
    import sys
    smoke(sys.argv[1] if len(sys.argv) > 1 else None)
