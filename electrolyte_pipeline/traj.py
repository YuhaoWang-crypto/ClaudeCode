"""Shared trajectory access: MDAnalysis Universe + species bookkeeping."""
from __future__ import annotations

import json
import os

import numpy as np

from electrolyte_pipeline.e0_systems import WORK


class Traj:
    def __init__(self, label: str, workdir: str = WORK, start_frac: float = 0.0):
        import MDAnalysis as mda
        self.label = label
        self.dir = os.path.join(workdir, label)
        self.u = mda.Universe(os.path.join(self.dir, "initial.pdb"), os.path.join(self.dir, "prod.dcd"))
        with open(os.path.join(self.dir, "atoms.json")) as fh:
            atoms = json.load(fh)
        self.species = np.array([a["species"] for a in atoms])
        self.element = np.array([a["element"] for a in atoms])
        self.molid = np.array([a["molecule"] for a in atoms])
        with open(os.path.join(self.dir, "system_record.json")) as fh:
            self.record = json.load(fh)
        rr = os.path.join(self.dir, "run_record.json")
        self.run_record = json.load(open(rr)) if os.path.exists(rr) else {}
        self.frame_ps = self.run_record.get("frame_ps", 2.0)
        self.start = int(len(self.u.trajectory) * start_frac)
        self.cation = "Na" if "Na" in self.record["composition"]["counts"] else "Li"
        # masses from OpenMM system for COM / dielectric
        self.mass = np.array([a.mass for a in self.u.atoms])
        # partial charges are not in the PDB: pull them from system.xml
        self.charge = self._charges_from_xml()

    def _charges_from_xml(self):
        import openmm
        with open(os.path.join(self.dir, "system.xml")) as fh:
            s = openmm.XmlSerializer.deserialize(fh.read())
        nb = [f for f in s.getForces() if isinstance(f, openmm.NonbondedForce)][0]
        return np.array([nb.getParticleParameters(i)[0].value_in_unit(openmm.unit.elementary_charge)
                         for i in range(nb.getNumParticles())])

    def sel(self, species: str | None = None, element: str | None = None) -> np.ndarray:
        m = np.ones(len(self.species), bool)
        if species:
            m &= self.species == species
        if element:
            m &= self.element == element
        return np.where(m)[0]

    @property
    def n_frames(self):
        return len(self.u.trajectory) - self.start

    def frames(self, stride: int = 1):
        for ts in self.u.trajectory[self.start::stride]:
            yield ts

    def box(self, ts):
        return ts.dimensions[:3].astype(float)   # Å, orthorhombic

    def counts(self):
        return self.record["composition"]["counts"]
