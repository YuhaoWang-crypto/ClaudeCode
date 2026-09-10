"""Shared utilities: system construction, MD driver, CVs, state definitions."""

import os

import numpy as np
import openmm
from openmm import unit
from openmm.app import PDBFile, ForceField, Simulation, NoCutoff

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")
RESULTS = os.path.join(HERE, "results")
FIGURES = os.path.join(HERE, "figures")
for _d in (DATA, RESULTS, FIGURES):
    os.makedirs(_d, exist_ok=True)

TEMPERATURE = 300.0 * unit.kelvin
FRICTION = 1.0 / unit.picosecond
TIMESTEP = 0.002 * unit.picoseconds
KT_KCAL = 0.0019872041 * 300.0  # kcal/mol at 300 K

# ---------------------------------------------------------------------------
# System
# ---------------------------------------------------------------------------

_TOPOLOGY = None
_SYSTEM = None
_MASSES = None


def get_system():
    """Return (topology, system, reference_positions_nm, masses)."""
    global _TOPOLOGY, _SYSTEM, _MASSES
    pdb = PDBFile(os.path.join(DATA, "nanma.pdb"))
    if _SYSTEM is None:
        ff = ForceField("amber14-all.xml")
        _SYSTEM = ff.createSystem(
            pdb.topology,
            nonbondedMethod=NoCutoff,
            constraints=openmm.app.HBonds,
            rigidWater=False,
        )
        _TOPOLOGY = pdb.topology
        _MASSES = np.array(
            [
                _SYSTEM.getParticleMass(i).value_in_unit(unit.dalton)
                for i in range(_SYSTEM.getNumParticles())
            ]
        )
    pos = np.array(pdb.positions.value_in_unit(unit.nanometer))
    return _TOPOLOGY, _SYSTEM, pos, _MASSES


def atom_index_map():
    top, _, _, _ = get_system()
    idx = {}
    for atom in top.atoms():
        idx[f"{atom.name}_{atom.residue.index}"] = atom.index
    return idx


def make_simulation(platform_name="CPU", threads=None):
    top, system, pos, _ = get_system()
    integ = openmm.LangevinMiddleIntegrator(TEMPERATURE, FRICTION, TIMESTEP)
    plat = openmm.Platform.getPlatformByName(platform_name)
    props = {}
    if platform_name == "CPU" and threads:
        props["Threads"] = str(threads)
    sim = Simulation(top, system, integ, plat, props)
    sim.context.setPositions(pos)
    return sim


# ---------------------------------------------------------------------------
# Collective variables (used only for analysis / state definition, never to bias)
# ---------------------------------------------------------------------------


def _dihedral(p0, p1, p2, p3):
    """Vectorised dihedral in degrees.  Inputs (..., 3)."""
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
    b1n = b1 / np.linalg.norm(b1, axis=-1, keepdims=True)
    v = b0 - (b0 * b1n).sum(-1, keepdims=True) * b1n
    w = b2 - (b2 * b1n).sum(-1, keepdims=True) * b1n
    x = (v * w).sum(-1)
    y = (np.cross(b1n, v) * w).sum(-1)
    return np.degrees(np.arctan2(y, x))


_IDX = None


def phi_psi(coords):
    """coords: (..., natoms, 3) -> (phi, psi) in degrees."""
    global _IDX
    if _IDX is None:
        m = atom_index_map()
        _IDX = (
            m["C_0"], m["N_1"], m["CA_1"], m["C_1"], m["N_2"],
        )
    cace, n, ca, c, nnme = _IDX
    coords = np.asarray(coords)
    phi = _dihedral(coords[..., cace, :], coords[..., n, :],
                    coords[..., ca, :], coords[..., c, :])
    psi = _dihedral(coords[..., n, :], coords[..., ca, :],
                    coords[..., c, :], coords[..., nnme, :])
    return phi, psi


# Metastable states of NANMA in vacuum.
#
# The slow degree of freedom is the backbone dihedral phi.  The reference
# metadynamics free energy along phi shows two basins: one centred near
# phi = -57 (C7eq / alpha_R) and one near phi = +72 (C7ax), the latter
# merging into a broad extended region that continues past phi = 180.  The
# basins are separated by a 7.4 kcal/mol saddle at phi ~ 0 and, going the
# other way round, by a 13 kcal/mol saddle near phi = -127, so essentially
# all of the flux goes through the phi ~ 0 channel.  The core boundaries
# below lie inside each basin, short of both saddles.
#
# These definitions enter only as the *reactant* and *product* the method is
# asked to connect, and as the labels of the commitment outcome.  They are
# never used to bias the dynamics or to parameterise a pathway.
#
# State B wraps across phi = 180, so the test is periodic: B is everything
# from +25 round through 180 to -150.
# Two nested definitions are needed, and conflating them is a mistake that
# quietly ruins the committor.
#
# BASIN_* are the thermodynamic macrostates: everything on one side of the
# barriers.  They are what the free energies are integrated over.
#
# CORE_* are strict subsets used to decide that a trajectory has *committed*.
# They must sit well away from the barrier: with a boundary only 25 degrees
# from the saddle, ordinary libration carries a structure that is still
# genuinely on the barrier across it within one saved frame, and every
# measured committor collapses to 0 or 1.
BASIN_A = dict(name="C7eq / alpha_R", phi_lo=-105.0, phi_hi=-25.0)
BASIN_B = dict(name="C7ax / extended", phi_lo=25.0, phi_hi=210.0)
CORE_A = dict(name="C7eq / alpha_R", phi_lo=-95.0, phi_hi=-45.0)
CORE_B = dict(name="C7ax / extended", phi_lo=55.0, phi_hi=210.0)

# representative structures used to start the seed runs
START_A = (-82.0, 73.0)
START_B = (71.0, -62.0)


def _ang_dist(a, b):
    return (a - b + 180.0) % 360.0 - 180.0


def _in_range(phi, lo, hi):
    """Periodic interval test on phi (degrees), hi may exceed 180."""
    return ((phi - lo) % 360.0) < (hi - lo)


def _assign(phi, a, b):
    phi = np.asarray(phi, dtype=float)
    out = np.full(phi.shape, -1, dtype=int)
    out = np.where(_in_range(phi, a["phi_lo"], a["phi_hi"]), 0, out)
    out = np.where(_in_range(phi, b["phi_lo"], b["phi_hi"]), 1, out)
    return out


def which_core(phi, psi=None):
    """Strict commitment test: 0 = committed to A, 1 = to B, -1 = neither."""
    return _assign(phi, CORE_A, CORE_B)


def which_basin(phi, psi=None):
    """Thermodynamic macrostate: 0 = A, 1 = B, -1 = on a barrier."""
    return _assign(phi, BASIN_A, BASIN_B)


# ---------------------------------------------------------------------------
# Featurisation for the neural networks
# ---------------------------------------------------------------------------
# Both the diffusion model and the committor network act on a
# translation/rotation-invariant Cartesian representation: the heavy-atom
# coordinates after optimal (mass-weighted) superposition on a reference.
# No dihedral, no hand-picked CV enters the learning.

_HEAVY = None


def heavy_atoms():
    global _HEAVY
    if _HEAVY is None:
        top, _, _, _ = get_system()
        _HEAVY = np.array(
            [a.index for a in top.atoms() if a.element.symbol != "H"]
        )
    return _HEAVY


def kabsch_align(coords, ref, weights=None):
    """Superimpose (..., n, 3) onto ref (n, 3).  Returns aligned coords."""
    coords = np.asarray(coords, dtype=np.float64)
    single = coords.ndim == 2
    if single:
        coords = coords[None]
    w = np.ones(coords.shape[1]) if weights is None else np.asarray(weights)
    w = w / w.sum()
    refc = ref - (ref * w[:, None]).sum(0)
    out = np.empty_like(coords)
    for i, x in enumerate(coords):
        xc = x - (x * w[:, None]).sum(0)
        h = (xc * w[:, None]).T @ refc
        u, s, vt = np.linalg.svd(h)
        d = np.sign(np.linalg.det(vt.T @ u.T))
        r = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
        out[i] = xc @ r.T
    return out[0] if single else out


_REF_HEAVY = None


def featurize(coords):
    """(..., natoms, 3) nm -> (..., n_heavy*3) aligned heavy-atom coordinates."""
    global _REF_HEAVY
    heavy = heavy_atoms()
    _, _, ref, masses = get_system()
    if _REF_HEAVY is None:
        _REF_HEAVY = ref[heavy]
    x = np.asarray(coords)[..., heavy, :]
    shp = x.shape
    x = x.reshape(-1, shp[-2], 3)
    x = kabsch_align(x, _REF_HEAVY, weights=masses[heavy])
    return x.reshape(*shp[:-2], shp[-2] * 3)


def defeaturize_into(full_coords, feat):
    """Write heavy-atom features back into a full coordinate set."""
    heavy = heavy_atoms()
    out = np.array(full_coords, dtype=np.float64)
    out[..., heavy, :] = np.asarray(feat).reshape(*np.shape(feat)[:-1], len(heavy), 3)
    return out
