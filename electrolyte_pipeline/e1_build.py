"""
E1a — from a recipe to an explicit particle set and box   [digest p4, panel A]

    counts (E0)  -> per-species conformers + NAGL charges
                 -> packmol random packing (no severe overlaps)
                 -> OpenFF Sage parametrisation (Interchange)
                 -> OpenMM System XML + PDB + system_record.json

Ion charge scaling is applied to the *molecule* charges before parametrisation
so that intramolecular 1-4 exceptions stay consistent with the scaled charges.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

from electrolyte_pipeline.e0_systems import (SERIES, SERIES_BY_LABEL, SPECIES, WORK,
                                             ForceFieldSpec, DynamicsSpec,
                                             system_record, write_record)


def _species_molecule(name: str, ff: ForceFieldSpec):
    from openff.toolkit import Molecule
    from openff.units import unit
    spec = SPECIES[name]
    mol = Molecule.from_smiles(spec["smiles"], allow_undefined_stereo=True)
    mol.name = name
    if mol.n_atoms == 1:                       # monatomic ion: formal charge
        q = float(spec["charge"])
        mol.partial_charges = np.array([q]) * unit.elementary_charge
    else:
        mol.generate_conformers(n_conformers=1)
        mol.assign_partial_charges(ff.charge_model)
    if spec["charge"] != 0 and ff.ion_charge_scale != 1.0:
        mol.partial_charges = mol.partial_charges * ff.ion_charge_scale
    return mol


def build(label: str, ff: ForceFieldSpec | None = None, dyn: DynamicsSpec | None = None,
          outdir: str | None = None, seed: int = 2024) -> str:
    """Build one composition; returns the output directory."""
    from openff.toolkit import ForceField
    from openff.interchange import Interchange
    from openff.interchange.components._packmol import pack_box
    from openff.units import unit
    import openmm
    from openmm import app

    ff = ff or ForceFieldSpec()
    dyn = dyn or DynamicsSpec()
    comp = SERIES_BY_LABEL[label]
    outdir = outdir or os.path.join(WORK, label)
    os.makedirs(outdir, exist_ok=True)

    names = list(comp.counts)
    mols = [_species_molecule(n, ff) for n in names]
    counts = [comp.counts[n] for n in names]

    # --- packmol -----------------------------------------------------------
    L = comp.box_nm * unit.nanometer
    box = np.eye(3) * comp.box_nm * unit.nanometer
    topology = pack_box(molecules=mols, number_of_copies=counts,
                        box_vectors=box, tolerance=2.0 * unit.angstrom)
    topology.box_vectors = box

    # --- parametrise ---------------------------------------------------------
    sage = ForceField(ff.name)
    ic = Interchange.from_smirnoff(sage, topology, charge_from_molecules=mols,
                                   allow_nonintegral_charges=(ff.ion_charge_scale != 1.0))
    system = ic.to_openmm(combine_nonbonded_forces=True)
    omm_top = ic.to_openmm_topology()
    positions = ic.positions.to_openmm()

    # cutoff / switching / PME as recorded in the spec
    for f in system.getForces():
        if isinstance(f, openmm.NonbondedForce):
            f.setNonbondedMethod(openmm.NonbondedForce.PME if ff.pme else openmm.NonbondedForce.CutoffPeriodic)
            f.setCutoffDistance(ff.cutoff_nm * openmm.unit.nanometer)
            f.setUseSwitchingFunction(True)
            f.setSwitchingDistance(ff.switch_nm * openmm.unit.nanometer)
            f.setUseDispersionCorrection(True)
            nb = f
    # net charge sanity
    qtot = sum(nb.getParticleParameters(i)[0].value_in_unit(openmm.unit.elementary_charge)
               for i in range(nb.getNumParticles()))

    with open(os.path.join(outdir, "system.xml"), "w") as fh:
        fh.write(openmm.XmlSerializer.serialize(system))
    with open(os.path.join(outdir, "initial.pdb"), "w") as fh:
        app.PDBFile.writeFile(omm_top, positions, fh)

    # atom-level bookkeeping used by every analysis module
    atoms = []
    for res in omm_top.residues():
        for a in res.atoms():
            atoms.append({"index": a.index, "element": a.element.symbol,
                          "species": res.name, "molecule": res.index})
    # OpenFF residue names: use molecule order to map back to species labels
    mol_index = 0
    for name, n in zip(names, counts):
        for k in range(n):
            pass
    # residue index -> species label (Interchange keeps molecule order)
    res_species = []
    for name, n in zip(names, counts):
        res_species += [name] * n
    for a in atoms:
        a["species"] = res_species[a["molecule"]]
    with open(os.path.join(outdir, "atoms.json"), "w") as fh:
        json.dump(atoms, fh)

    rec = system_record(comp, ff, dyn, extra={
        "build": {"packmol_tolerance_A": 2.0, "seed": seed,
                  "n_atoms": system.getNumParticles(),
                  "net_charge_e": round(qtot, 6),
                  "species_order": names, "counts_order": counts,
                  "partial_charges_note": "NAGL AM1-BCC (GNN) for molecules; "
                                          f"ions formal charge x {ff.ion_charge_scale}"}
    })
    write_record(os.path.join(outdir, "system_record.json"), rec)
    print(f"[E1 build] {label}: {system.getNumParticles()} atoms, box {comp.box_nm} nm, "
          f"net charge {qtot:+.3f} e -> {outdir}")
    return outdir


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("labels", nargs="*", default=[c.label for c in SERIES])
    p.add_argument("--scale", type=float, default=None, help="override ion charge scale")
    p.add_argument("--suffix", default="")
    a = p.parse_args(argv)
    ff = ForceFieldSpec()
    if a.scale is not None:
        ff.ion_charge_scale = a.scale
    for lab in a.labels:
        build(lab, ff=ff, outdir=os.path.join(WORK, lab + a.suffix))


if __name__ == "__main__":
    main()
