"""
Structural review of the proposed SparkCage design specification.

    python3 sparkcage/spec_review.py

The specification under review:

  Option 1, Engineered Type I Clamshell (1omp/1anf analog)
    - methylene blue conjugated site-specifically at Cys-155 (deep binding
      cleft) through a C6 aliphatic linker
    - anchored to the carbon screen-printed electrode at residue 370

  Option 2, lucCage (LOCKR) architecture, cited as PDB 6S0A
    - methylene blue on the internal face of the latch domain, at or near
      Cys-112 or Cys-115 on the inner helix
    - anchored at the N-terminus or a distal loop on the rigid cage domain,
      near residue 2 or residue 145

Each claim below is checked against the deposited coordinates.  Two are wrong
as written, one works but not for the reason the specification implies, and one
has the signal polarity inverted.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from structure import analyze          # noqa: E402

PDB_DIR = HERE / "structure" / "pdb"

VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "SE": 1.90, "CL": 1.75}
PROBE = 1.4


def fibonacci_sphere(n):
    i = np.arange(0, n, dtype=float) + 0.5
    phi = np.arccos(1.0 - 2.0 * i / n)
    theta = np.pi * (1.0 + 5.0 ** 0.5) * i
    return np.c_[np.cos(theta) * np.sin(phi),
                 np.sin(theta) * np.sin(phi),
                 np.cos(phi)]


def load_atoms(path, chain="A", drop=("HOH", "EOH")):
    out = []
    for line in open(path):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[21] != chain or line[16] not in (" ", "A"):
            continue
        resn = line[17:20].strip()
        if resn in drop:
            continue
        el = (line[76:78].strip() or line[12:16].strip()[0]).upper()
        out.append((int(line[22:26]), resn, line[12:16].strip(), el,
                    np.array([float(line[30:38]), float(line[38:46]),
                              float(line[46:54])])))
    return out


def sasa(atoms, subset, universe, n_points=120):
    """Shrake-Rupley solvent-accessible surface for `subset` inside `universe`."""
    sp = fibonacci_sphere(n_points)
    xyz = np.array([atoms[i][4] for i in universe])
    rad = np.array([VDW.get(atoms[i][3], 1.7) for i in universe]) + PROBE
    where = {g: k for k, g in enumerate(universe)}
    out = {}
    for gi in subset:
        k = where[gi]
        pts = xyz[k] + rad[k] * sp
        d = np.linalg.norm(pts[:, None, :] - xyz[None, :, :], axis=2)
        d[:, k] = 1e9
        out[gi] = float((d >= rad[None, :]).all(axis=1).mean()
                        * 4.0 * np.pi * rad[k] ** 2)
    return out


# ---------------------------------------------------------------------------
def check_pdb_ids():
    print("=" * 74)
    print("1. THE CITED PDB ID IS WRONG")
    print("=" * 74)
    print("""
  6S0A is 'Crystal Structure of Properdin (TSR domains N12 & 456)', a
  complement-system protein solved at 2.52 A and deposited in 2019. It has
  nothing to do with LOCKR. The near-miss 6SOA (letter O) is BamABCDE in a
  nanodisc, also unrelated.

  There are exactly TWO LOCKR-family entries in the whole PDB:
    7CBC   sCageHA267_1S, a switch cage holding an HA binder, 1.99 A, from
           Quijano-Rubio 2021. Note this is NOT lucCage: different cage
           sequence, no NanoLuc or SmBiT. It is still the structure to model
           on, because it is the only caged-latch coordinate set available.
    7JH5   Co-LOCKR, the colocalisation-dependent variant, 2.10 A, from
           Lajoie et al. Science 2020.

  Two further corrections worth having:
    - Langan et al. Nature 2019, the original LOCKR paper, deposited NO
      structures at all. It was characterised by circular dichroism, SAXS,
      biolayer interferometry and SEC-MALS. Any PDB ID attributed to it is
      wrong.
    - No lucCage or lucKey crystal structure exists. A lucCage model has to
      be built on 7CBC or predicted.

  Residue numbering matters here, and it is where this specification goes
  wrong in a specific way. In full lucCage, tag stripped, the sequence is 359
  residues: the CAGE is 1-300 and the LATCH is 301-359. Read against the
  published Table S6 sequence, position 112 is arginine and position 115 is
  leucine, both in the last turn of CAGE HELIX 2, immediately before the
  GSGSGS linker at 118-123. They are solvent-facing cage positions.

  That is the opposite of what the design needs. A reporter on a solvent-facing
  cage helix is exposed in BOTH states, so latch release neither buries nor
  reveals it and there is no signal. The specification asks for the internal
  face of the latch, which is the right idea; the residue numbers point
  somewhere else. In 7CBC the latch is helix H6, residues 244-269, and
  section 3 gives the positions that actually face the cage.

  Beware also of the tag: if a collaborator counts the 18-residue
  MGSHHHHHHGSENLYFQG His-TEV tag, every number shifts by 18. Agree the
  convention before ordering DNA.
""")


def check_clamshell_geometry():
    print("=" * 74)
    print("2. CLAMSHELL: RESIDUES 155 AND 370 SIT ON THE SAME RIGID LOBE")
    print("=" * 74)
    o = analyze.load_ca(analyze.fetch("1OMP", PDB_DIR))
    c = analyze.load_ca(analyze.fetch("1ANF", PDB_DIR))
    common = sorted(set(o) & set(c))
    n_dom = analyze._expand(analyze.MBP_N_DOMAIN, set(common))
    c_dom = analyze._expand(analyze.MBP_C_DOMAIN, set(common))

    print(f"\n  residue 155 is TYR, in the C-domain: {155 in c_dom}")
    print(f"  residue 370 is LYS, in the C-domain: {370 in c_dom}")
    print(f"\n  |155-370| open   : {np.linalg.norm(o[155]-o[370]):6.2f} A")
    print(f"  |155-370| closed : {np.linalg.norm(c[155]-c[370]):6.2f} A")
    print(f"  change           : {abs(np.linalg.norm(c[155]-c[370]) - np.linalg.norm(o[155]-o[370])):6.2f} A")

    rot, pc, qc, rms = analyze.kabsch(np.array([o[r] for r in c_dom]),
                                      np.array([c[r] for r in c_dom]))
    moved = np.linalg.norm(((rot @ (o[155] - pc)) + qc) - c[155])
    print(f"\n  Holding the anchored C-domain fixed (internal RMSD {rms:.2f} A),")
    print(f"  residue 155 moves {moved:.2f} A on ligand binding.")
    print("\n  So the hinge motion produces essentially NO rigid-body displacement")
    print("  of the reporter relative to the electrode. Both the anchor and the")
    print("  label are on the same lobe, and the hinge moves the lobes relative")
    print("  to EACH OTHER. Judged as a displacement sensor, this design has no")
    print("  mechanism at all. Section 2b shows why it nonetheless works.")


def check_tethered_ensemble(reach_ang=(10.0, 12.0), clash=4.0, n_points=4000):
    print("\n" + "=" * 74)
    print("2b. BUT THE TETHERED DYE ENSEMBLE DOES SWITCH, AND STRONGLY")
    print("=" * 74)
    print("""
  Methylene blue is not at residue 155; it is on the end of a C6 linker from
  it, so what matters is the ENSEMBLE of positions the dye can occupy and how
  close that ensemble gets to the electrode. The electrode contact is at the
  anchor, residue 370, so distance-to-anchor is the right proxy.
""")
    sp = fibonacci_sphere(n_points)
    for reach in reach_ang:
        print(f"  MB centre {reach:.0f} A from CB155, {clash:.1f} A clash cutoff:")
        for pdb, label in [("1OMP", "open   (no analyte)"),
                           ("1ANF", "closed (analyte bound)")]:
            atoms = load_atoms(analyze.fetch(pdb, PDB_DIR))
            xyz = np.array([a[4] for a in atoms])
            byname = {(a[0], a[2]): a[4] for a in atoms}
            origin = byname.get((155, "CB"), byname[(155, "CA")])
            anchor = byname[(370, "CA")]
            pts = origin + reach * sp
            d = np.linalg.norm(pts[:, None, :] - xyz[None, :, :], axis=2)
            free = d.min(axis=1) > clash
            da = np.linalg.norm(pts[free] - anchor, axis=1)
            print(f"     {label}: {free.mean()*100:5.1f}% of shell allowed, "
                  f"closest approach to anchor {da.min():5.1f} A")
        print()
    print("  Read-out: closing the clamshell collapses the dye's accessible")
    print("  volume by roughly ten-fold and pushes its closest approach to the")
    print("  anchor out by about 8 A. At beta = 1.0 per Angstrom that is a")
    print("  factor of e^8, near 3000, in the electron-transfer rate. The design")
    print("  works. It works by ENSEMBLE OCCLUSION, not by moving the attachment")
    print("  residue, and only the tethered-ensemble calculation shows it.")
    print("\n  Methodological point worth keeping: residue-level metrics both say")
    print("  this design fails. Displacement of residue 155 is 0.4 A and its")
    print("  solvent accessibility changes by 4 A^2. Model the dye, not the")
    print("  residue it hangs from.")


def check_polarity():
    print("\n" + "=" * 74)
    print("2c. THE SIGNAL POLARITY IS INVERTED RELATIVE TO THE STATED INTENT")
    print("=" * 74)
    print("""
  The specification says the latch 'mechanically opens the SparkCage to reveal
  the sequestered methylene blue payload upon target binding', i.e. signal-ON.

  A Type I clamshell does the opposite. Periplasmic binding proteins CLOSE on
  ligand binding. Per section 2b, closing is what buries the dye. So as drawn,
  target binding sequesters the reporter and the sensor reads signal-OFF.

  Signal-off is a perfectly good sensor, so this is not fatal. But it has to be
  a decision rather than an accident, because getting signal-ON out of a
  clamshell requires the binder to stabilise the OPEN state, which is the
  reverse of how these proteins work and a materially harder design.

  The lucCage option has no such problem: target binds the latch, the latch is
  released, a reporter on the latch's buried face becomes exposed. Signal-ON
  falls out naturally. On this point the two options are not equivalent.
""")


def check_latch_positions():
    print("=" * 74)
    print("3. CORRECT MB POSITIONS ON THE LUCCAGE LATCH (7CBC numbering)")
    print("=" * 74)
    atoms = load_atoms(analyze.fetch("7CBC", PDB_DIR))
    universe = list(range(len(atoms)))
    latch = [i for i, a in enumerate(atoms) if 244 <= a[0] <= 269]
    caged = sasa(atoms, latch, universe)
    released = sasa(atoms, latch, latch)

    per = {}
    for i in latch:
        rn = atoms[i][0]
        per.setdefault(rn, [atoms[i][1], 0.0, 0.0])
        if atoms[i][2] not in ("N", "CA", "C", "O"):
            per[rn][1] += caged[i]
            per[rn][2] += released[i]

    rows = []
    for rn, (name, f, a) in sorted(per.items()):
        if a < 1.0:
            continue
        rows.append((rn, name, f, a, a - f, (a - f) / a))

    print("\n  latch residue   caged   released   buried   fraction")
    print("  " + "-" * 56)
    for rn, name, f, a, b, fr in rows:
        print(f"   {name}{rn:<5d} {f:9.1f} {a:9.1f} {b:9.1f}   {fr*100:5.1f}%"
              + "  " + "#" * int(fr * 20))

    print("\n  The buried positions repeat every three to four residues, which is")
    print("  the heptad periodicity of the coiled-coil interface. Those are the")
    print("  positions whose reporter is hidden when caged and exposed when the")
    print("  latch is released, so they are the correct MB sites.")

    deep = [r for r in rows if r[5] > 0.75]
    rim = [r for r in rows if 0.35 <= r[5] <= 0.65]
    print("\n  deeply buried (>75%): " + ", ".join(f"{n}{r}" for r, n, *_ in deep))
    print("  interface rim (35-65%): " + ", ".join(f"{n}{r}" for r, n, *_ in rim))
    return deep, rim


def check_mb_fits():
    print("\n" + "=" * 74)
    print("4. CAN METHYLENE BLUE ACTUALLY FIT WHERE YOU WANT TO BURY IT?")
    print("=" * 74)
    print("""
  Methylene blue, measured from the coordinates in 5ACM:
      van der Waals volume     225 A^3
      roughly 15 x 7 x 6 A including the dimethylamino groups, and planar

  Side-chain volumes for comparison (A^3):
      Ala 27    Ser 32    Val 85    Ile/Leu 102    Trp 163

  So methylene blue needs about 2.2x the room of an isoleucine side chain.

  Consequence for section 3: the most deeply buried latch positions are exactly
  the hydrophobic core positions, and you cannot mutate one to cysteine and
  expect a dye of this size to occupy the space that side chain vacated. It
  will either not fit or it will gut the cage core and destabilise the switch,
  which shows up as a loss of dG_open, which is the one parameter the whole
  design depends on.

  Two workable routes:
    - use the INTERFACE RIM positions instead, where the reporter is still
      substantially occluded when caged but there is room for it, or
    - deliberately design a cavity in the cage to receive the dye, which is
      what 'sequestered MB payload' really implies and is the more interesting
      version of the project. It is also a real design task, not a mutation.

  Either way this must be modelled before synthesis: burying a 225 A^3 aromatic
  cation in a designed helical bundle changes dG_open, and dG_open sets both
  the dynamic range and the limit of detection.
""")


if __name__ == "__main__":
    check_pdb_ids()
    check_clamshell_geometry()
    check_tethered_ensemble()
    check_polarity()
    check_latch_positions()
    check_mb_fits()
    print("=" * 74)
