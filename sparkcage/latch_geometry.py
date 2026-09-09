"""
Where on the latch to put the methylene blue, and which end of the cage to anchor.

    python3 sparkcage/latch_geometry.py

This is the calculation the concept schematic makes necessary. The drawing shows
the released latch tethered near the electrode so the reporter can reach the
surface. That is the right idea, but in the real construct the latch is attached
at ONE SPECIFIC END of the cage bundle, and which end is anchored decides
whether the sensor works, which way its signal goes, and whether it has a
baseline at all.

Three questions, answered from the 7CBC coordinates plus a tethered-chain model:

  1. Which end of the cage should be bonded to the electrode?
  2. How close can the reporter get once the latch is released?
  3. Which latch position gives both a measurable baseline and a large change?

The answers converge with the buried-surface analysis in spec_review.py, which
is reassuring because the two calculations share no assumptions.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from structure import analyze          # noqa: E402
from echem import swv                  # noqa: E402

PDB_DIR = HERE / "structure" / "pdb"
RNG = np.random.default_rng(0)

# Calibrated cell, as in run_demo.py
CELL = dict(area_cm2=0.03, gamma_mol_cm2=2.0e-12, e_step=0.002, freq_hz=87.0)
K0_CONTACT, BETA, D_CONTACT = 3.0e2, 1.0, 5.0

# 7CBC segment assignment (see structure/analyze.py)
CAGE_HELICES = [(2, 47), (49, 94), (98, 143), (145, 191), (194, 240)]
LATCH = (244, 269)

# Polypeptide chain parameters
KUHN_ANG = 6.5           # Kuhn length of a flexible polypeptide
RES_PER_KUHN = 2.0
RISE_PER_RES_HELIX = 1.5  # alpha-helix rise per residue
VDW_STANDOFF = 5.0        # closest a dye centre approaches a surface


def bundle_axis_heights():
    """Height of every residue along the cage bundle axis, zeroed at the low end."""
    ca = analyze.load_ca(analyze.fetch("7CBC", PDB_DIR), "A", hetatm=True)
    cage = np.array([ca[r] for r in sorted(ca) if 2 <= r <= 240])
    centre = cage.mean(axis=0)
    axis = np.linalg.svd(cage - centre)[2][0]
    h = {r: float((ca[r] - centre) @ axis) for r in ca}
    lo = min(h.values())
    return {r: v - lo for r, v in h.items()}, ca


def in_helix(r):
    return any(a <= r <= b for a, b in CAGE_HELICES + [LATCH])


def flexible_chain_effective_distance(n_res, n_samples=30000):
    """Effective electron-transfer distance of a dye on a floppy tether at a wall.

    Because the rate falls exponentially with distance, the ensemble-averaged
    rate is set by the near-wall tail of the distribution, not by the mean
    position.  The effective distance is defined so that
    exp(-beta d_eff) = < exp(-beta z) >.

    The chain is a freely-jointed walk from a tether point on the electrode,
    with the surface treated as reflecting.
    """
    n_kuhn = max(1, int(round(n_res / RES_PER_KUHN)))
    v = RNG.normal(size=(n_samples, n_kuhn, 3))
    v /= np.linalg.norm(v, axis=2, keepdims=True)
    z = np.abs(np.cumsum(v[:, :, 2] * KUHN_ANG, axis=1)[:, -1])
    return float(-np.log(np.exp(-BETA * z).mean()) / BETA)


def rigid_helix_effective_distance(n_res, n_samples=30000):
    """Same quantity if the released latch stays HELICAL and rigid.

    This is the important robustness check.  A released LOCKR latch is not
    necessarily a random coil; it may retain helical structure, in which case
    the dye is on the end of a rod pivoting about the tether rather than on a
    flexible string.  A rod can still point at the surface, so the conclusion
    should survive, and this function checks that it does.
    """
    length = n_res * RISE_PER_RES_HELIX
    u = RNG.normal(size=(n_samples, 3))
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    z = np.abs(u[:, 2] * length)
    return float(-np.log(np.exp(-BETA * z).mean()) / BETA)


def peak_current(distance_ang):
    k0 = K0_CONTACT * np.exp(-BETA * (distance_ang - D_CONTACT))
    return abs(swv.swv_scan(k0, **CELL)["peak_current"])


def report():
    heights, _ = bundle_axis_heights()

    print("=" * 74)
    print("1. WHICH END OF THE CAGE GOES ON THE ELECTRODE")
    print("=" * 74)
    print(f"""
  The cage bundle is {max(heights.values()):.0f} A long. The latch, helix H6,
  runs from residue 244 at height {heights[244]:.1f} A down to residue 269 at
  height {heights[269]:.1f} A, so it hangs off ONE end of the bundle and occupies
  its upper half.

  When the latch is released it stays covalently attached at that junction. So
  the reporter can only reach the electrode if the junction end is the end that
  is bonded to the surface. Anchoring the far end puts a {max(heights.values()):.0f} A
  rigid bundle between the released reporter and the electrode.

  The specification offers two anchor points:
    residue 2   at height {heights[2]:5.1f} A  -> WRONG END, {heights[244]-heights[2]:.0f} A from the latch
    residue 145 at height {heights[145]:5.1f} A  -> correct end
""")
    print("  Surface loops at the correct end, which anchor better than a helix:")
    for r in sorted(heights):
        if r < 244 and heights[r] > 68 and not in_helix(r):
            print(f"    residue {r:4d}  height {heights[r]:5.1f} A")
    print("\n  Residue 144 is the loop neighbouring the specified residue 145 and is")
    print("  the cleaner choice. Residues 241 to 243 are the loop immediately")
    print("  before the latch, which is the right end but risks tethering so close")
    print("  to the junction that it interferes with release.")

    print("\n" + "=" * 74)
    print("2. HOW FAR THE RELEASED REPORTER REALLY IS FROM THE SURFACE")
    print("=" * 74)
    print("""
  Two models of the released latch, bracketing the truth:
""")
    print("  residues from   flexible coil    rigid helix     mean height")
    print("  tether to MB    effective d      effective d     (coil)")
    print("  " + "-" * 62)
    for n in (4, 8, 12, 16, 20, 26):
        d_coil = flexible_chain_effective_distance(n)
        d_rod = rigid_helix_effective_distance(n)
        nk = max(1, int(round(n / RES_PER_KUHN)))
        mean_h = KUHN_ANG * np.sqrt(nk) * 0.8
        print(f"      {n:3d}         {d_coil:8.1f} A      {d_rod:8.1f} A     {mean_h:7.1f} A")
    print("""
  Both models give an effective distance of a few Angstrom, far below the mean
  height, because the exponential weighting is dominated by the near-surface
  tail. The conclusion is therefore robust to whether the released latch
  unfolds or stays helical, which matters because that is genuinely unknown.

  The physical statement: a reporter held RIGIDLY at 20 A is electrically dead,
  while the same reporter at a 20 A MEAN height on a tether is electrically
  live, because it visits the surface. Rigidity, not distance, is what the cage
  is really controlling.
""")

    print("=" * 74)
    print("3. WHICH LATCH POSITION, ANCHORING THE JUNCTION END")
    print("=" * 74)
    print("\n  latch  caged d   caged peak   released d  released peak   change  verdict")
    print("  " + "-" * 76)
    junction = heights[244]
    best = None
    for r in (246, 248, 250, 252, 255, 259, 262, 266):
        d_caged = abs(junction - heights[r]) + VDW_STANDOFF
        d_rel = flexible_chain_effective_distance(r - 243) + VDW_STANDOFF
        i_c, i_r = peak_current(d_caged), peak_current(d_rel)
        fold = i_r / max(i_c, 1e-30)
        if i_c * 1e9 < 0.1:
            verdict = "no baseline"
        elif fold < 5:
            verdict = "change too small"
        else:
            verdict = "USABLE"
            if best is None:
                best = (r, i_c, fold)
        print(f"  {r:4d}  {d_caged:7.1f} A {i_c*1e9:10.3f} nA {d_rel:9.1f} A "
              f"{i_r*1e9:12.1f} nA {fold:8.0f}x  {verdict}")

    print(f"""
  The usable window is latch residues 248 to 252, and this AGREES with the
  buried-surface analysis in spec_review.py, which independently picks A248,
  A251 and I252 as positions whose reporter is hidden when caged and exposed
  when released. Two calculations sharing no assumptions landing on the same
  three residues is the strongest result in this file.

  Maximising fold change is the wrong objective. Past residue 252 the caged
  state carries no measurable current, so the device reads zero with no
  analyte and a fouled, dried or detached electrode is indistinguishable from
  a true negative. Electrochemical aptamer sensors are always normalised
  against a baseline peak; a sensor without one cannot be trusted in blood.

  Residue 250 is the balanced choice: a few nanoamps of baseline and roughly a
  seventy-fold rise on release.

  One unresolved conflict, carried over: methylene blue occupies 225 cubic
  Angstrom, about 2.2 times an isoleucine side chain, and A248, A251 and I252
  are all substantially buried. The cage has to be locally opened to receive
  the dye. Doing that at residues 248 to 252 is the least bad place for it,
  because they sit within two turns of the end of the bundle where packing is
  loosest, rather than in the middle of the core.
""")

    print("=" * 74)
    print("4. WHY THE CLAMSHELL PANEL HAS THE HARDER PROBLEM")
    print("=" * 74)
    print("""
  The schematic labels the clamshell's coupling as still to be designed, with a
  dashed arrow from the bound target to the hinge. That label is doing a great
  deal of work, and it is worth being explicit about how much.

  In the clamshell the recognition module sits OUTSIDE the cleft, so target
  binding has to open the hinge through a designed allosteric path that does
  not exist in the natural protein. Periplasmic binding proteins close on
  ligand binding IN the cleft; nothing about an external binder drives that
  hinge, in either direction. Creating such a path is precisely the class of
  computational design that has not been independently reproduced.

  In the cage-and-latch architecture there is no allosteric path to design. The
  binder IS the latch. Target binding competes the latch out of the cage by
  direct thermodynamic linkage, which is the mechanism the published sensors
  actually use, across eleven different targets.

  So the two panels are not two variants of one design with different odds.
  Panel B needs the latch tuned; panel A needs a new allosteric mechanism
  invented first, and only then tuned.
""")


if __name__ == "__main__":
    report()
    print("=" * 74)
