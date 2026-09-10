"""Build the NANMA (N-acetyl-N'-methylalaninamide, alanine dipeptide) system.

The paper (arXiv:2510.24979v1) uses NANMA in vacuum as its proof-of-concept
system for Gen-COMPAS.  OpenMM ships no test-system library on PyPI, so we
build ACE-ALA-NME from an internal-coordinate (Z-matrix) definition with the
NeRF algorithm, then relax it with the AMBER14 force field.

Outputs:
    data/nanma.pdb          minimized structure
"""

import os

import numpy as np

BOHR = None  # unused; kept for clarity that everything below is in angstrom

# ---------------------------------------------------------------------------
# NeRF placement
# ---------------------------------------------------------------------------


def place_atom(a, b, c, bond, angle_deg, dihedral_deg):
    """Place atom D given three placed atoms.

    D is at distance `bond` from `a`, angle D-a-b = `angle_deg`, and dihedral
    D-a-b-c = `dihedral_deg` (natural extension reference frame).
    """
    ang = np.deg2rad(angle_deg)
    tor = np.deg2rad(dihedral_deg)

    bc = a - b
    bc /= np.linalg.norm(bc)
    n = np.cross(b - c, bc)
    n /= np.linalg.norm(n)
    m = np.cross(n, bc)

    d2 = np.array(
        [
            -bond * np.cos(ang),
            bond * np.sin(ang) * np.cos(tor),
            bond * np.sin(ang) * np.sin(tor),
        ]
    )
    return a + d2[0] * bc + d2[1] * m + d2[2] * n


# (name, residue_name, residue_index, parent, angle_ref, dihedral_ref,
#  bond, angle, dihedral)
# Backbone is built first, then the hydrogens.
PHI0, PSI0 = -80.0, 80.0  # start in the C7eq basin


def zmat(PHI0=PHI0, PSI0=PSI0):
    return [
        # name        res  ridx  parent      angref      dihref      b      ang    dih
        ("CH3", "ACE", 0, None, None, None, 0.0, 0.0, 0.0),
        ("C", "ACE", 0, "CH3_0", None, None, 1.522, 0.0, 0.0),
        ("O", "ACE", 0, "C_0", "CH3_0", None, 1.229, 120.5, 0.0),
        ("N", "ALA", 1, "C_0", "O_0", "CH3_0", 1.335, 122.9, 180.0),
        ("CA", "ALA", 1, "N_1", "C_0", "O_0", 1.449, 121.9, 0.0),
        ("C", "ALA", 1, "CA_1", "N_1", "C_0", 1.522, 110.1, PHI0),
        ("N", "NME", 2, "C_1", "CA_1", "N_1", 1.335, 116.6, PSI0),
        ("O", "ALA", 1, "C_1", "CA_1", "N_1", 1.229, 120.5, PSI0 + 180.0),
        ("CH3", "NME", 2, "N_2", "C_1", "CA_1", 1.449, 121.9, 180.0),
        # hydrogens
        ("HH31", "ACE", 0, "CH3_0", "C_0", "O_0", 1.090, 109.5, 0.0),
        ("HH32", "ACE", 0, "CH3_0", "C_0", "O_0", 1.090, 109.5, 120.0),
        ("HH33", "ACE", 0, "CH3_0", "C_0", "O_0", 1.090, 109.5, 240.0),
        ("H", "ALA", 1, "N_1", "C_0", "O_0", 1.010, 119.0, 180.0),
        ("CB", "ALA", 1, "CA_1", "N_1", "C_0", 1.525, 110.5, PHI0 + 122.0),
        ("HA", "ALA", 1, "CA_1", "N_1", "C_0", 1.090, 109.5, PHI0 - 120.0),
        ("HB1", "ALA", 1, "CB_1", "CA_1", "N_1", 1.090, 109.5, 60.0),
        ("HB2", "ALA", 1, "CB_1", "CA_1", "N_1", 1.090, 109.5, 180.0),
        ("HB3", "ALA", 1, "CB_1", "CA_1", "N_1", 1.090, 109.5, 300.0),
        ("H", "NME", 2, "N_2", "C_1", "CA_1", 1.010, 119.0, 180.0),
        ("HH31", "NME", 2, "CH3_2", "N_2", "C_1", 1.090, 109.5, 0.0),
        ("HH32", "NME", 2, "CH3_2", "N_2", "C_1", 1.090, 109.5, 120.0),
        ("HH33", "NME", 2, "CH3_2", "N_2", "C_1", 1.090, 109.5, 240.0),
    ]


ZMAT = zmat()

# order the atoms as AMBER/PDB expects them within each residue
PDB_ORDER = [
    ("HH31", 0), ("CH3", 0), ("HH32", 0), ("HH33", 0), ("C", 0), ("O", 0),
    ("N", 1), ("H", 1), ("CA", 1), ("HA", 1), ("CB", 1), ("HB1", 1),
    ("HB2", 1), ("HB3", 1), ("C", 1), ("O", 1),
    ("N", 2), ("H", 2), ("CH3", 2), ("HH31", 2), ("HH32", 2), ("HH33", 2),
]
RESNAMES = {0: "ACE", 1: "ALA", 2: "NME"}


def build_coordinates(phi=PHI0, psi=PSI0):
    coords = {}
    for i, (name, res, ridx, par, aref, dref, b, ang, dih) in enumerate(zmat(phi, psi)):
        key = f"{name}_{ridx}"
        if i == 0:
            coords[key] = np.zeros(3)
        elif i == 1:
            coords[key] = coords[par] + np.array([b, 0.0, 0.0])
        elif i == 2:
            a = coords[par]
            bb = coords[aref]
            v = a - bb
            v /= np.linalg.norm(v)
            perp = np.array([0.0, 1.0, 0.0])
            perp = perp - np.dot(perp, v) * v
            perp /= np.linalg.norm(perp)
            th = np.deg2rad(ang)
            coords[key] = a + b * (-np.cos(th) * v + np.sin(th) * perp)
        else:
            coords[key] = place_atom(
                coords[par], coords[aref], coords[dref], b, ang, dih
            )
    return coords


def chirality_ok(coords):
    """L-amino acid: the improper N-C-CB-HA seen from CA must be right-handed."""
    n = coords["N_1"]
    c = coords["C_1"]
    cb = coords["CB_1"]
    ca = coords["CA_1"]
    # signed volume of (N-CA, C-CA, CB-CA); L-Ala is negative in this convention
    v = np.dot(np.cross(n - ca, c - ca), cb - ca)
    return v < 0


def write_pdb(path, coords):
    lines = []
    serial = 1
    for name, ridx in PDB_ORDER:
        key = f"{name}_{ridx}"
        x, y, z = coords[key]
        element = "H" if name.startswith("H") else name[0]
        lines.append(
            f"ATOM  {serial:5d} {name:<4s} {RESNAMES[ridx]:>3s} A{ridx + 1:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}"
        )
        serial += 1
    lines.append("TER")
    lines.append("END")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def main():
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data = os.path.join(here, "data")
    os.makedirs(data, exist_ok=True)

    coords = build_coordinates()
    if not chirality_ok(coords):
        # swap the HA / CB torsion offsets to obtain the L enantiomer
        cb, ha = coords["CB_1"].copy(), coords["HA_1"].copy()
        ca = coords["CA_1"]
        coords["CB_1"] = ca + (cb - ca) * 0 + (ha - ca) / np.linalg.norm(ha - ca) * 1.525
        coords["HA_1"] = ca + (cb - ca) / np.linalg.norm(cb - ca) * 1.090
        print("[build] flipped CA substituents to obtain L-alanine")

    raw = os.path.join(data, "nanma_raw.pdb")
    write_pdb(raw, coords)

    from openmm import LangevinMiddleIntegrator, unit
    from openmm.app import PDBFile, ForceField, Simulation, NoCutoff

    pdb = PDBFile(raw)
    ff = ForceField("amber14-all.xml")
    system = ff.createSystem(
        pdb.topology, nonbondedMethod=NoCutoff, constraints=None, rigidWater=False
    )
    integ = LangevinMiddleIntegrator(
        300 * unit.kelvin, 1.0 / unit.picosecond, 0.001 * unit.picoseconds
    )
    sim = Simulation(pdb.topology, system, integ)
    sim.context.setPositions(pdb.positions)
    e0 = sim.context.getState(getEnergy=True).getPotentialEnergy()
    sim.minimizeEnergy(maxIterations=5000)
    st = sim.context.getState(getEnergy=True, getPositions=True)
    print(f"[build] atoms={system.getNumParticles()}")
    print(f"[build] E(initial) = {e0.value_in_unit(unit.kilocalorie_per_mole):10.2f} kcal/mol")
    print(f"[build] E(minimized) = "
          f"{st.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole):10.2f} kcal/mol")

    out = os.path.join(data, "nanma.pdb")
    with open(out, "w") as fh:
        PDBFile.writeFile(sim.topology, st.getPositions(), fh)
    print(f"[build] wrote {out}")


if __name__ == "__main__":
    main()
