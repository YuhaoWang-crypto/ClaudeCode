"""
Methylene blue placement: a screening library that decides the mechanism.

    python3 sparkcage/mb_placement_design.py

SCOPE. Monod Bio supplies a latch that forms a working switch with the cage.
Our deliverable is where to put the reporter and how to prove the choice was
right. This file produces a small variant library and the decision rules that
make each construct informative rather than merely tested.

THREE CONSTRAINTS, AND WHY THEY BARELY INTERSECT.

  1. Occlusion. The reporter must be substantially hidden when caged and
     exposed when released, or there is no signal. That wants HIGH burial.
  2. Volume. Methylene blue is 225 cubic Angstrom, about 2.2 times an
     isoleucine side chain. The buried fraction of it has to fit where the
     native side chain was. That wants a LARGE native residue and MODERATE
     burial.
  3. Distance. The caged reporter must be close enough to the electrode to
     give a measurable baseline peak, because an electrochemical sensor with
     no baseline cannot be distinguished from a fouled or detached one. That
     wants a position within roughly 10-15 Angstrom of the anchored end.

Constraint 1 pulls against 2, and 3 restricts both to the first two turns of
the latch. On the 7CBC scaffold exactly one position satisfies all three, and
it does so at the cost of a salt bridge. That is not a comfortable margin, and
it is the reason this has to be screened rather than predicted.

WHY VOLUME CANNOT BE BOUGHT. The obvious fix for a buried position is to
truncate neighbouring core residues and make a cavity. A single large-to-small
substitution in a packed protein core costs roughly 4-5 kcal/mol, from the
classic cavity-creation studies on T4 lysozyme. The entire dG_open budget for
a picomolar-sensitivity sensor is 1-2 kcal/mol (see clinical_targets.py). One
truncation therefore costs several times the whole switching free energy, and
the cage stops being a cage. The reporter has to fit where it is put.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from structure import analyze     # noqa: E402
from echem import swv             # noqa: E402

PDB_DIR = HERE / "structure" / "pdb"
MB_VOLUME_A3 = 225.0

# Tsai/Zamyatnin side-chain volumes, cubic Angstrom
SIDECHAIN_VOL = {
    "GLY": 0, "ALA": 27, "SER": 32, "CYS": 55, "THR": 62, "PRO": 72,
    "VAL": 85, "ASN": 73, "ASP": 70, "ILE": 102, "LEU": 102, "MET": 105,
    "GLU": 96, "GLN": 99, "LYS": 119, "ARG": 148, "HIS": 98, "PHE": 135,
    "TYR": 141, "TRP": 163,
}

# Side-chain burial fraction of each latch residue when caged, measured from
# 7CBC in spec_review.py (SASA inside the full protein vs the latch alone).
BURIAL = {
    244: 0.082, 245: 0.432, 246: 0.000, 247: 0.000, 248: 0.950, 249: 0.589,
    250: 0.000, 251: 0.759, 252: 0.905, 253: 0.000, 254: 0.127, 255: 0.785,
    256: 0.448, 257: 0.000, 258: 0.201, 259: 0.983, 260: 0.356, 261: 0.000,
    262: 0.908, 263: 0.843, 264: 0.000, 265: 0.058, 266: 0.918, 267: 0.129,
}

CELL = dict(area_cm2=0.03, gamma_mol_cm2=2.0e-12, e_step=0.002, freq_hz=87.0)
K0_CONTACT, BETA, D_CONTACT = 3.0e2, 1.0, 5.0
BASELINE_FLOOR_NA = 0.1     # smallest caged peak that can serve as a reference


def load_scaffold():
    ca = analyze.load_ca(analyze.fetch("7CBC", PDB_DIR), "A", hetatm=True)
    identities = {}
    for line in open(PDB_DIR / "7cbc.pdb"):
        if line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == "CA" \
                and line[21] == "A":
            identities[int(line[22:26])] = line[17:20].strip()
    cage = np.array([ca[r] for r in sorted(ca) if 2 <= r <= 240])
    centre = cage.mean(axis=0)
    axis = np.linalg.svd(cage - centre)[2][0]
    height = {r: float((ca[r] - centre) @ axis) for r in ca}
    return identities, height


def caged_distance(height, residue, standoff=5.0):
    """Distance from the electrode with the latch-junction end anchored."""
    return abs(height[244] - height[residue]) + standoff


def peak_na(distance_ang):
    k0 = K0_CONTACT * np.exp(-BETA * (distance_ang - D_CONTACT))
    return abs(swv.swv_scan(k0, **CELL)["peak_current"]) * 1e9


def score_positions():
    ident, height = load_scaffold()
    rows = []
    for r, burial in sorted(BURIAL.items()):
        name = ident.get(r, "UNK")
        native_vol = SIDECHAIN_VOL.get(name, 80)
        need = burial * MB_VOLUME_A3
        deficit = need - native_vol
        d = caged_distance(height, r)
        rows.append(dict(res=r, name=name, burial=burial, native_vol=native_vol,
                         deficit=deficit, distance=d, baseline_na=peak_na(d)))
    return rows


def report():
    rows = score_positions()

    print("=" * 78)
    print("1. THE THREE CONSTRAINTS SCORED ON EVERY LATCH POSITION")
    print("=" * 78)
    print("\n  pos  aa   burial   native   MB deficit   caged d   baseline   passes")
    print("  " + "-" * 74)
    for row in rows:
        if row["burial"] < 0.30:
            continue
        p = []
        p.append("occl" if row["burial"] >= 0.35 else "----")
        p.append("vol" if row["deficit"] <= 0 else "---")
        p.append("dist" if row["baseline_na"] >= BASELINE_FLOOR_NA else "----")
        n_pass = sum(1 for x in p if not x.startswith("-"))
        mark = "  <== ALL THREE" if n_pass == 3 else ""
        print(f"  {row['res']:3d}  {row['name']}  {row['burial']*100:5.1f}%  "
              f"{row['native_vol']:5.0f} A3  {row['deficit']:+8.0f} A3  "
              f"{row['distance']:6.1f} A  {row['baseline_na']:8.2f} nA   "
              f"{'/'.join(p)}{mark}")

    print("""
  Read the deficit column as: how much extra room methylene blue needs beyond
  what the native side chain vacates. Positive means a cavity has to be made,
  which the introduction explains you cannot afford.

  Only three positions have room for the dye at all, and they are the three
  with a large native side chain and partial burial. Of those, only one is
  close enough to the electrode to leave a baseline.

  One position looks tempting and is not a candidate. PRO245 sits closest to
  the electrode with by far the largest baseline, but it is the proline of the
  TDP motif that caps the N terminus of the latch helix. Mutating it would
  unpick the start of the helix it is supposed to report on, and the dye still
  would not fit. Structural role beats the score sheet.
""")

    print("=" * 78)
    print("2. THE LEAD CANDIDATE, AND WHAT IT COSTS")
    print("=" * 78)
    print("""
  R249C is the only position on this scaffold that satisfies all three
  constraints: 59% buried when caged, an arginine side chain whose 148 cubic
  Angstrom is just enough for the 133 that must be buried, and 12.8 Angstrom
  from the anchored end, which leaves a few nanoamps of baseline.

  It is not free. ARG249 makes a salt bridge to ASP246 at 3.33 Angstrom, an
  i-to-i-minus-3 pair within the latch helix, and its guanidinium also reaches
  GLU46 on cage helix 1. Mutating it removes both. An intrahelical salt bridge
  is worth roughly half to one and a half kcal/mol, which is a real fraction of
  a dG_open budget that is only one to two kcal/mol in total.

  That cost is not necessarily bad. It weakens the latch helix and therefore
  lowers dG_open, and section 3 of clinical_targets.py shows that a picomolar
  target wants a LOW dG_open anyway. It may push the design in the right
  direction. But it changes the switch thermodynamics at the same time as the
  reporter placement, so the two effects must be separated by measurement, not
  assumed. Include the R249C-without-dye control for exactly this reason.

  LYS256 and LYS260 also have room for the dye and make no salt bridge, so they
  are cleaner mutations. They sit 23 and 30 Angstrom out, which on this
  scaffold means no baseline. They become viable if the anchor chemistry ends
  up placing the bundle differently than modelled, which is why they stay in
  the library as the distance arm.
""")

    print("=" * 78)
    print("3. THE SCREENING LIBRARY")
    print("=" * 78)
    print("""
  Twelve constructs, chosen so that each one changes exactly one thing and the
  result of the set identifies the mechanism rather than just ranking guesses.

  ARM A - distance, at constant chemistry (all Cys, all C6 linker)
    A1  A248C    11.6 A   deeply buried, dye does not fit; expect little signal
    A2  R249C    12.8 A   the lead
    A3  I252C    17.7 A   buried, dye does not fit, 5 A further out
    A4  K256C    23.4 A   dye fits, expect no baseline
    A5  K260C    29.5 A   dye fits, expect no baseline at all
        -> If signal tracks DISTANCE, A2 wins and A4/A5 are silent.
        -> If it tracks BURIAL, A1 and A3 win despite the volume problem.
        Those two outcomes are distinguishable in one experiment.

  ARM B - linker length at the lead position
    B1  R249C + C2 linker      dye held close to the backbone
    B2  R249C + C6 linker      = A2, the reference
    B3  R249C + C11 linker     dye can swing further
        -> A long linker helps the released state reach the electrode and hurts
           the caged state's occlusion. The optimum is a real experiment; the
           model cannot call it because it depends on how the linker packs.

  ARM C - the geometry control, the single most informative construct
    C1  R249C, anchored at the FAR end of the cage instead of the junction end
        -> The geometric model says the released reporter then has a 76 A
           bundle between it and the electrode and the sensor must be dead.
           If C1 signals anyway, the model is wrong and we learn it for the
           price of one construct instead of a design cycle.

  ARM D - the controls that make the rest interpretable
    D1  R249C, no dye conjugated        electrode background, fouling baseline
    D2  wild-type latch + dye on a surface cysteine elsewhere on the cage
                                        reports non-specific dye signal
    D3  R249C with the latch deleted    dye on a free cage, no switching
    D4  R249C plus a non-binding target does the signal need real binding

  Twelve constructs is one synthesis batch and one plate. The information comes
  from the pattern across arms, not from any single number.
""")

    print("=" * 78)
    print("4. WHAT TO MEASURE, AND THE DECISION RULES")
    print("=" * 78)
    print("""
  Measure in this order. Each step gates the next, so a failure stops the run
  rather than propagating into an uninterpretable data set.

  STEP 1  Coverage. Integrate the square-wave peak with no analyte and convert
          to surface coverage. Compare with an independent measure of how much
          protein is on the surface. If the electroactive coverage is far below
          the protein coverage, most of the monolayer is dark and the whole
          signal analysis is about a small subpopulation.

  STEP 2  Frequency sweep, 5 to 1000 Hz, before any analyte. Find the peak in
          current versus frequency. That maximum locates k0/f near unity and
          gives the interrogation frequency. Working on the wrong side of it
          makes a released reporter read as MORE current instead of less, so
          this must be done before any dose-response is interpreted.

  STEP 3  Voltammogram SHAPE, still with no analyte. Fit with a distribution of
          rate constants rather than a single one. A single narrow population
          means one well-defined reporter environment. A broad or bimodal
          distribution means a mixture of responsive and non-responsive probes,
          and it is the difference between the reporter barely moving and most
          of the monolayer being dead. Those need opposite fixes, so establish
          which one before redesigning anything.

  STEP 4  Dose-response against the real target, with D4 alongside. Extract the
          apparent dissociation constant and the fractional signal change.

  STEP 5  Only if steps 1 to 4 are clean, characterise in plasma.

  Decision rules, stated in advance so the result cannot be rationalised after
  the fact:

    - Baseline peak below about 0.1 nA: reject the position regardless of fold
      change. Without a reference peak the device cannot tell a true negative
      from a dead electrode.
    - Fold change below 3 at the frequency chosen in step 2: reject.
    - C1 gives signal: stop and revisit the geometric model before building
      anything else.
    - Step 3 shows a broad rate distribution: fix immobilisation chemistry and
      surface density first. Moving the dye will not help.
    - Apparent dissociation constant more than tenfold above the clinical
      decision threshold: the problem is the binder or dG_open, not placement,
      and it goes back to the thermodynamic design rather than to another
      round of conjugation.
""")

    print("=" * 78)
    print("5. THE TRANSFERABLE RULE")
    print("=" * 78)
    print("""
  The residue numbers above belong to 7CBC. The latch that arrives from Monod
  Bio will have a different sequence, so what has to transfer is the rule, not
  the numbers. Applied to any latch:

    Choose the position that is 40 to 65 percent buried at the cage-latch
    interface, carries a large native side chain (arginine, lysine, methionine,
    glutamine or a long aromatic), and sits 10 to 15 Angstrom from the end of
    the bundle that is bonded to the electrode.

    Reject positions above 75 percent burial: they switch well but the dye does
    not fit and buying the space costs more free energy than the switch has.
    Reject positions below 35 percent burial: the dye is exposed in both states.
    Reject anything beyond about 20 Angstrom from the anchored end: no baseline.

  Every quantity in that rule is computable from a model of the delivered latch
  in the cage, before any DNA is ordered, using the code in this repository.
  The screening library then exists to test the rule rather than to search
  blindly, which is what keeps it to twelve constructs instead of fifty.
""")


if __name__ == "__main__":
    report()
    print("=" * 78)
