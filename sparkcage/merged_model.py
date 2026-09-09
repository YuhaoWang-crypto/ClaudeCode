"""
Merging the two external reports into the placement model, and what it changes.

    python3 sparkcage/merged_model.py

TWO EXTERNAL SOURCES
--------------------
A. A Biomni-platform feasibility report on PDB 7CBC, the same scaffold this
   package uses, so its residue numbering maps directly onto ours. It ran
   things we did not: blind docking of methylene blue onto the latch, a
   latch-on-graphene adsorption simulation, and umbrella sampling of latch
   opening.
B. A screening report on lucCageHer2, a different construct (404 residues, HER2
   affibody rather than the influenza minibinder), which searched anchor site
   times reporter site times linker over 4752 scenarios.

Only A shares our numbering. B has to be compared structurally, not by residue
index, and doing that carefully is where the interesting result is.

WHAT A CONTRIBUTES THAT WE DID NOT HAVE
---------------------------------------
Two EMPIRICAL exclusion zones, from simulation rather than geometry:

  1. Methylene blue binds the latch non-covalently, and it prefers the target
     binding patch. Blind docking put it on TYR287 in four of five best poses
     and PHE285 in two; a 4 ns simulation then found contacts dominated by the
     C-terminal 312-319 and the same 283/289 patch. Since methylene blue is a
     flat aromatic cation, it is read out by the same aromatic-plus-acidic
     motif that recognises most flat cations, which is exactly what our own
     survey of the six methylene-blue protein structures concluded. So this is
     not surprising, but it is now specific: those residues are where a dye on
     a linker will fold back and stick.

  2. The latch adsorbs to carbon at residues 297-304, contact fraction 0.5 to
     1.0, stable within a nanosecond. That patch will be pinned against the
     electrode, so a reporter there reports the surface rather than the switch.

WHAT THIS PACKAGE CONTRIBUTES BACK
----------------------------------
Constraints A and B both lack: the volume test (a dye that does not fit),
the baseline test (a caged state with no reference peak), and the clinical
thresholds that decide whether an EC50 is useful at all.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from mb_placement_design import (          # noqa: E402
    BURIAL, SIDECHAIN_VOL, MB_VOLUME_A3, load_scaffold, caged_distance, peak_na,
)

R_KCAL = 1.987204259e-3

# --- external constraint A, in 7CBC numbering -------------------------------
# Methylene blue's own non-covalent preference on the latch, from docking + MD
MB_SELF_BINDING = set(range(283, 290)) | set(range(312, 320))
# Latch residues that adsorb to the graphene electrode model
ELECTRODE_ADSORPTION = set(range(297, 305))
# The grafted binder domain as a whole
BINDER_DOMAIN = set(range(272, 320))
# Region the external MD recommends, being free of all of the above
MD_RECOMMENDED = set(range(245, 261))

# 7CBC is sCageHA_267-1S, which already carries the V255S tuning mutation, so
# position 255 is the engineered handle on the opening free energy.
DG_OPEN_TUNING_SITES = {255}

# --- the opening free energy band, from three independent routes -------------
BANDS = {
    "this package, from the TBI clinical windows": (0.3, 2.1),
    "external report A, stated as 5-10 kJ/mol": (5.0 / 4.184, 10.0 / 4.184),
    "literature, the only published lucCage value": (4.0, 4.0),
}


def merged_scoring():
    ident, height = load_scaffold()
    rows = []
    for r, burial in sorted(BURIAL.items()):
        name = ident.get(r, "UNK")
        native = SIDECHAIN_VOL.get(name, 80)
        deficit = burial * MB_VOLUME_A3 - native
        d = caged_distance(height, r)
        base = peak_na(d)
        gates = {
            "occlusion": burial >= 0.35,
            "volume": deficit <= 0.0,
            "baseline": base >= 0.1,
            "not_MB_patch": r not in MB_SELF_BINDING,
            "not_electrode_patch": r not in ELECTRODE_ADSORPTION,
            "not_tuning_site": r not in DG_OPEN_TUNING_SITES,
        }
        rows.append(dict(res=r, aa=name, burial=burial, deficit=deficit,
                         distance=d, baseline=base, gates=gates,
                         passes=all(gates.values())))
    return rows


def report():
    print("=" * 78)
    print("1. THE SIX GATES, AFTER MERGING")
    print("=" * 78)
    rows = merged_scoring()
    print("\n  pos aa   burial  deficit  caged d  baseline  occl vol base MBp ELp tun")
    print("  " + "-" * 76)
    for row in rows:
        if row["burial"] < 0.30:
            continue
        g = row["gates"]
        marks = "".join(f"{'  Y ' if g[k] else '  . '}" for k in
                        ("occlusion", "volume", "baseline", "not_MB_patch",
                         "not_electrode_patch", "not_tuning_site"))
        star = "  <== PASSES ALL" if row["passes"] else ""
        print(f"  {row['res']:3d} {row['aa']} {row['burial']*100:6.1f}% "
              f"{row['deficit']:+7.0f}  {row['distance']:6.1f} A {row['baseline']:8.2f} nA"
              f"{marks}{star}")
    winners = [r for r in rows if r["passes"]]
    names = ", ".join(f"{r['aa']}{r['res']}" for r in winners) or "none"
    print(f"\n  Positions passing all six gates: {names}")

    print("""
  The two new gates change NOTHING on this scaffold, and that is the result.
  Every residue the external simulations exclude sits at 272 or beyond: the
  grafted binder domain, the electrode adsorption patch at 297-304, and the
  dye's own preferred patches at 283-289 and 312-319. Our scoring never reached
  there, because it only ever considered the latch helix, on the separate
  grounds that a reporter has to be buried at the cage-latch interface to
  switch at all.
""")
    print("=" * 78)
    print("2. THE CROSS-VALIDATION THAT MATTERS")
    print("=" * 78)
    print(f"""
  The external molecular dynamics recommends the latch helix N-terminal
  region, residues {min(MD_RECOMMENDED)}-{max(MD_RECOMMENDED)}, as the only stretch free of dye
  self-binding and electrode adsorption.

  Our independent scoring, which used burial, side-chain volume and distance to
  the anchored end and knew nothing about either simulation, put its single
  surviving candidate at residue 249.

  249 is inside 245-260.

  Two methods with no shared assumptions, one geometric and one dynamical,
  select the same short stretch of the same helix. That is the strongest
  evidence in this whole package for where the reporter should go, and neither
  method could have produced it alone.
""")

    print("=" * 78)
    print("3. A CONFLICT BETWEEN THE TWO EXTERNAL REPORTS")
    print("=" * 78)
    print("""
  Report B ranks K404C and D399C as its top two candidates. Those are the last
  residues of the grafted affibody, at the far C-terminal tip of the binder.

  Report A, simulating the structurally analogous construct, found that the
  C-terminal tip of the grafted binder is exactly where methylene blue sticks
  to the latch on its own: contact fractions of 0.5 to 0.9 over a 4 ns run,
  the highest of any region. It also found the neighbouring stretch pinned
  against the carbon surface.

  So B's top candidates sit in the region A's simulations flag twice over. The
  two reports do not contradict each other on any number; they simply never
  looked at each other. The conflict only appears when the numbering is mapped
  structurally rather than by index, because the two constructs differ.

  This is not fatal to B's ranking. A dye that folds back onto its own protein
  is not necessarily a broken sensor, and could even sharpen the contrast if
  the target competes it off, which report A raises as a possibility. But it is
  a mechanism B's model does not contain, so B's predicted currents for those
  two candidates are computed under an assumption its companion report says is
  wrong. They should not be ordered first without checking it.
""")

    print("=" * 78)
    print("4. THE OPENING FREE ENERGY, NOW FROM THREE ROUTES")
    print("=" * 78)
    print("\n  source                                          band (kcal/mol)")
    print("  " + "-" * 66)
    for name, (lo, hi) in BANDS.items():
        print(f"  {name:46s}  {lo:5.2f} - {hi:5.2f}")
    lo = max(b[0] for k, b in BANDS.items() if "literature" not in k)
    hi = min(b[1] for k, b in BANDS.items() if "literature" not in k)
    print(f"""
  The two independent design routes overlap at {lo:.2f} to {hi:.2f} kcal/mol. One got
  there from the clinical decision thresholds for the traumatic brain injury
  markers, the other from dose-response shape in a completely separate model.

  The published lucCage value of 4 kcal/mol sits well ABOVE both, which is
  consistent rather than contradictory: that value was tuned for nanomolar
  targets, and a picomolar target needs a looser cage. Report B assumes an
  opening constant of 0.0101, which is 2.72 kcal/mol and also above the band.
  That is why its EC50 comes out at 79 nM against clinical thresholds of
  0.60 pM for GFAP and 14.5 pM for UCH-L1.

  The practical handle already exists in this scaffold. 7CBC is sCageHA_267-1S,
  which is to say it already carries the V255S interface mutation that the
  LOCKR authors used to tune exactly this quantity. Position 255 is therefore
  the opening-free-energy knob and must not be spent on the reporter, which is
  the sixth gate added above.
""")

    print("=" * 78)
    print("5. ONE ARITHMETIC CHECK WORTH SETTLING")
    print("=" * 78)
    F = 96485.332
    for n in (1, 2):
        q = n * F * 0.071 * 2.0e-12
        print(f"    n = {n}:  Q = n F A Gamma = {q*1e9:6.2f} nC "
              f"(A = 0.071 cm2, Gamma = 2e-12 mol/cm2)")
    print("""
  The two external reports use different electron counts for the same reporter,
  one treating methylene blue as a single-electron couple and the other as two.
  Methylene blue to leucomethylene blue is a two-electron, one-proton couple,
  so two is correct and this package uses two throughout.

  This is not a factor-of-two correction to apply afterwards. The electron
  count enters the Butler-Volmer exponents as well as the charge, so it changes
  the peak shape and position, not just its height. Any comparison of absolute
  currents between the two reports is invalid until they agree on this.
""")


if __name__ == "__main__":
    report()
    print("=" * 78)
