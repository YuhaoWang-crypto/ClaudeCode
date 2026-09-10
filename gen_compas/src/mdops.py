"""MD operations used by Gen-COMPAS: targeted MD and committor shooting."""

import os
import sys

import numpy as np
import openmm
from openmm import unit
from openmm.app import Simulation

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402


class Engine:
    """One OpenMM context reused for every short simulation.

    Carries a switchable targeted-MD restraint.  Schlitter's original TMD uses
    an RMSD restraint to a reference structure; OpenMM cannot update the
    reference of an RMSDForce that lives inside a CustomCVForce, so the
    restraint is written as a harmonic pull of the heavy atoms toward
    reference positions that are re-superposed onto the instantaneous
    structure at every steering window.  This is equivalent to a best-fit RMSD
    restraint, and unlike a distance-matrix restraint it distinguishes a
    structure from its mirror image -- which matters here, because flipping
    the backbone dihedral phi is close to a reflection of the heavy-atom
    skeleton.  The force constant is zero for unbiased runs.
    """

    def __init__(self, platform="CPU", threads=1, seed=0):
        top, system, pos, masses = common.get_system()
        self.masses = masses
        self.heavy = common.heavy_atoms()

        # deep-copy so that the restraint is not added to the shared system
        system = openmm.XmlSerializer.deserialize(
            openmm.XmlSerializer.serialize(system)
        )
        ext = openmm.CustomExternalForce(
            "0.5*k_tmd*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        ext.addGlobalParameter("k_tmd", 0.0)
        for p in ("x0", "y0", "z0"):
            ext.addPerParticleParameter(p)
        for i in self.heavy:
            ext.addParticle(int(i), [0.0, 0.0, 0.0])
        system.addForce(ext)
        self._ext = ext

        integ = openmm.LangevinMiddleIntegrator(
            common.TEMPERATURE, common.FRICTION, common.TIMESTEP
        )
        integ.setRandomNumberSeed(seed)
        plat = openmm.Platform.getPlatformByName(platform)
        props = {"Threads": str(threads)} if platform == "CPU" else {}
        self.sim = Simulation(top, system, integ, plat, props)
        self.sim.context.setPositions(pos)
        self.n_atoms = system.getNumParticles()
        self.steps_used = 0

    # -- helpers ------------------------------------------------------------
    def positions(self):
        return self.sim.context.getState(getPositions=True).getPositions(
            asNumpy=True).value_in_unit(unit.nanometer)

    def set_positions(self, x):
        self.sim.context.setPositions(np.asarray(x) * unit.nanometer)

    def randomize_velocities(self, seed=None):
        if seed is None:
            self.sim.context.setVelocitiesToTemperature(common.TEMPERATURE)
        else:
            self.sim.context.setVelocitiesToTemperature(common.TEMPERATURE, seed)

    def set_reference(self, ref_heavy_xyz):
        """Set the restraint reference positions (nm, lab frame)."""
        r = np.asarray(ref_heavy_xyz, dtype=np.float64)
        for a, i in enumerate(self.heavy):
            self._ext.setParticleParameters(a, int(i), [float(r[a, 0]),
                                                        float(r[a, 1]),
                                                        float(r[a, 2])])
        self._ext.updateParametersInContext(self.sim.context)

    def set_k(self, k):
        self.sim.context.setParameter("k_tmd", float(k))

    def clear_restraint(self):
        self.sim.context.setParameter("k_tmd", 0.0)

    def potential_energy(self, x=None):
        """Unrestrained potential energy in kcal/mol."""
        if x is not None:
            self.clear_restraint()
            self.set_positions(x)
        return self.sim.context.getState(getEnergy=True).getPotentialEnergy(
        ).value_in_unit(unit.kilocalorie_per_mole)

    def step(self, n):
        self.sim.step(n)
        self.steps_used += n

    def ns_used(self):
        return self.steps_used * 0.002 / 1000.0

    def recover(self):
        """Reset the context after a blow-up (NaN coordinates)."""
        _, _, pos, _ = common.get_system()
        self.clear_restraint()
        self.sim.context.setPositions(pos * unit.nanometer)
        self.sim.context.setVelocitiesToTemperature(common.TEMPERATURE)
        self.sim.minimizeEnergy(maxIterations=200)


# ---------------------------------------------------------------------------
# Targeted MD
# ---------------------------------------------------------------------------


_GEOM = None


def _geometry_reference():
    """Heavy-atom covalent bonds and the loosest thresholds compatible with
    ordinary thermal fluctuation (calibrated on unbiased MD frames)."""
    global _GEOM
    if _GEOM is None:
        heavy = common.heavy_atoms()
        top, _, _, _ = common.get_system()
        hset = {int(i): k for k, i in enumerate(heavy)}
        bonds = [(hset[a.index], hset[b.index]) for a, b in top.bonds()
                 if a.index in hset and b.index in hset]
        nonbonded = [(i, j) for i in range(len(heavy))
                     for j in range(i + 1, len(heavy))
                     if (i, j) not in bonds and (j, i) not in bonds]
        _GEOM = (bonds, nonbonded)
    return _GEOM


def target_is_physical(target_heavy_xyz, bond_range=(0.09, 0.20),
                       min_sep=0.16):
    """Reject generated geometries that could not be a molecule at all.

    Only catastrophes are rejected: a heavy-atom covalent bond outside
    0.9-2.0 A, or two non-bonded heavy atoms closer than 1.6 A.  In 5 ns of
    unbiased MD the observed ranges are 1.16-1.66 A and >2.13 A, so the test
    never rejects a real conformation.  Moderately distorted targets are kept
    on purpose: the targeted-MD step is what projects them back onto the
    physical manifold, exactly as in the paper.
    """
    bonds, nonbonded = _geometry_reference()
    t = np.asarray(target_heavy_xyz, dtype=np.float64)
    d = np.linalg.norm(t[:, None] - t[None, :], axis=-1)
    for i, j in bonds:
        if not (bond_range[0] < d[i, j] < bond_range[1]):
            return False
    for i, j in nonbonded:
        if d[i, j] < min_sep:
            return False
    return True


def target_rmsd(engine, positions, target_heavy_xyz):
    """Heavy-atom RMSD (nm) to the target after optimal superposition."""
    x = np.asarray(positions)[engine.heavy]
    y = np.asarray(target_heavy_xyz)
    xa = common.kabsch_align(x, y - y.mean(0))
    return float(np.sqrt(((xa - (y - y.mean(0))) ** 2).sum(-1).mean()))


def targeted_md(engine, start_positions, target_heavy_xyz, n_steps=6000,
                k=200000.0, relax_steps=100, n_path=10):
    """Pull `start_positions` toward a generated heavy-atom target.

    `target_heavy_xyz` is a (n_heavy, 3) array produced by the diffusion model.
    The restraint force constant is ramped from 0 to `k`, driving the system
    onto the target geometry, and is then switched off for a short unbiased
    relaxation so that the resulting structure is physical.

    Returns dict(positions, rmsd, path) or None if the target was rejected or
    the run blew up.  `path` holds structures sampled along the steering, which
    is what makes the separatrix findable: the steering carries the molecule
    all the way across the barrier, so the crossing is bracketed by consecutive
    path points even when the committor is still poorly calibrated.
    """
    if not target_is_physical(target_heavy_xyz):
        return None
    try:
        return _targeted_md(engine, start_positions, target_heavy_xyz,
                            n_steps, k, relax_steps, n_path)
    except openmm.OpenMMException:
        engine.recover()
        return None


def _targeted_md(engine, start_positions, target_heavy_xyz, n_steps, k,
                 relax_steps, n_path):
    engine.set_positions(start_positions)
    engine.randomize_velocities()

    tgt = np.asarray(target_heavy_xyz, dtype=np.float64)
    tgt = tgt - tgt.mean(0)

    # Schlitter-style steering: at each window the target is best-fit onto the
    # instantaneous structure, and the restraint reference is placed a fraction
    # s of the way from the current heavy-atom positions to that fitted target.
    n_windows = 30
    per = max(1, n_steps // n_windows)
    keep_every = max(1, n_windows // n_path)
    engine.set_k(k)
    path = []
    for i in range(n_windows):
        s = (i + 1) / n_windows
        cur = engine.positions()[engine.heavy]
        cen = cur.mean(0)
        fitted = common.kabsch_align(tgt, cur - cen) + cen
        engine.set_reference((1.0 - s) * cur + s * fitted)
        engine.step(per)
        if i % keep_every == 0 or i == n_windows - 1:
            path.append(engine.positions())

    engine.clear_restraint()
    x_pulled = engine.positions()
    # Only a brief relaxation: the restraint artefacts disappear in tens of
    # femtoseconds, whereas a structure sitting on the barrier falls into a
    # basin within about a picosecond, which would defeat the whole point.
    engine.step(relax_steps)
    x = engine.positions()
    if not np.all(np.isfinite(x)):
        return None
    e = engine.sim.context.getState(getEnergy=True).getPotentialEnergy()
    if e.value_in_unit(unit.kilocalorie_per_mole) > 200.0:
        return None
    path.append(x)
    return dict(positions=x, rmsd=target_rmsd(engine, x_pulled,
                                              target_heavy_xyz),
                path=np.asarray(path, dtype=np.float32))


# ---------------------------------------------------------------------------
# Committor shooting
# ---------------------------------------------------------------------------


def shoot(engine, x0, total_ps=20.0, save_ps=0.2, seed=None):
    """Fixed-length unbiased trajectory from x0.

    The run is *not* stopped when the trajectory commits: a stopping rule tied
    to the states would censor the data and bias any Markov model built from
    it later.  The first core reached is recorded as the commitment outcome
    (0 = A, 1 = B, -1 = never committed within total_ps).
    """
    engine.clear_restraint()
    engine.set_positions(x0)
    engine.randomize_velocities(seed)
    try:
        return _shoot(engine, total_ps, save_ps)
    except openmm.OpenMMException:
        engine.recover()
        return None


def _shoot(engine, total_ps, save_ps):
    save_every = max(1, int(save_ps / 0.002))
    n_saves = int(total_ps / save_ps)

    coords = []
    for _ in range(n_saves):
        engine.step(save_every)
        x = engine.positions()
        if not np.all(np.isfinite(x)):
            break
        coords.append(x)
    if not coords:
        return None
    coords = np.asarray(coords, dtype=np.float32)
    phi, psi = common.phi_psi(coords)
    core = common.which_core(phi)
    hit = np.where(core >= 0)[0]
    outcome = int(core[hit[0]]) if len(hit) else -1
    commit_ps = float((hit[0] + 1) * save_ps) if len(hit) else float("nan")

    # Per-frame committor label: the state this frame's own future reaches
    # first.  Labelling every frame with the trajectory's first commitment
    # would be wrong for frames that come after it -- a run that touches A at
    # 5 ps and B at 18 ps has frames in between whose future is B, not A.
    labels = np.full(len(coords), -1, dtype=np.int8)
    nxt = -1
    for i in range(len(coords) - 1, -1, -1):
        if core[i] >= 0:
            nxt = int(core[i])
        labels[i] = nxt
    return dict(coords=coords, phi=phi, psi=psi, outcome=outcome,
                commit_ps=commit_ps, labels=labels)
