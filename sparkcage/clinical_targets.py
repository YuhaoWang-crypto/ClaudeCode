"""
What the actual clinical numbers for GFAP, UCH-L1 and NfL demand of the design.

    python3 sparkcage/clinical_targets.py

All cutoffs below are from primary FDA labeling, not from review articles.
Everything is converted to molar, because that is the unit the physics cares
about and it is where the decisive asymmetry between the targets shows up.

The short version: the three named targets are not one problem. UCH-L1 is a
tractable electrochemical target, GFAP at its mild-TBI threshold is about
twenty times harder for reasons no protein design can touch, and NfL is not an
acute marker at all.
"""

from __future__ import annotations

import numpy as np

R_KCAL = 1.987204259e-3

# UniProt-verified masses
MASS_DA = {"GFAP": 49880, "UCH-L1": 24824, "NfL": 61517}

# FDA-cleared decision thresholds, pg/mL.
# Logic in every cleared test: positive if EITHER marker is at or above cutoff.
CUTOFFS = {
    "Banyan BTI (DEN170045, 2018, serum/plasma)": {"GFAP": 22, "UCH-L1": 327},
    "i-STAT TBI Plasma (2021, EDTA plasma)": {"GFAP": 30, "UCH-L1": 360},
    "Alinity i / ARCHITECT TBI (2023)": {"GFAP": 35, "UCH-L1": 400},
    "i-STAT TBI Cartridge (2024, whole blood)": {"GFAP": 65, "UCH-L1": 360},
}

# Reportable ranges, i-STAT plasma cartridge, pg/mL
REPORTABLE = {"GFAP": (30, 10000), "UCH-L1": (200, 3200)}

# Representative measured concentrations, pg/mL (medians)
SCENARIOS = [
    ("healthy control", "GFAP", 15),
    ("mTBI, CT-negative", "GFAP", 78),
    ("mTBI, CT-positive", "GFAP", 739),
    ("moderate-severe TBI", "GFAP", 764),
    ("healthy control", "UCH-L1", 71),
    ("mTBI, CT-negative", "UCH-L1", 164),
    ("mTBI, CT-positive", "UCH-L1", 211),
    ("moderate-severe TBI", "UCH-L1", 705),
]

# Research-platform limits of detection, pg/mL
SIMOA_LOD = {"GFAP": 0.221, "UCH-L1": 1.74, "NfL": 0.104}


def molar(pg_ml, marker):
    return pg_ml * 1e-12 / MASS_DA[marker] * 1000.0


def time_to_occupancy(conc_m, k_on, occupancy=0.10):
    """Seconds for a tight binder to reach `occupancy`, kinetic (non-equilibrium).

    For Kd well below the analyte concentration the off-rate is negligible over
    the measurement, so phi(t) = 1 - exp(-k_on c t).  A rate-based readout does
    not have to wait for equilibrium, which is the only reason a picomolar
    measurement in minutes is conceivable at all.

    Note the direction of the affinity trade: a TIGHTER binder has a smaller
    off-rate and therefore equilibrates more slowly, not faster.  Affinity buys
    sensitivity, never speed.
    """
    return -np.log(1.0 - occupancy) / (k_on * conc_m)


def report():
    print("=" * 74)
    print("1. THE CUTOFFS IN MOLAR TERMS, WHERE THE ASYMMETRY APPEARS")
    print("=" * 74)
    print("\n  cleared test                                 GFAP        UCH-L1")
    print("  " + "-" * 68)
    for name, cut in CUTOFFS.items():
        g, u = molar(cut["GFAP"], "GFAP"), molar(cut["UCH-L1"], "UCH-L1")
        print(f"  {name:42s} {g*1e12:6.2f} pM   {u*1e12:6.2f} pM")
    ratio = molar(360, "UCH-L1") / molar(30, "GFAP")
    print(f"\n  The UCH-L1 threshold is {ratio:.0f}x higher in molar terms than GFAP's.")
    print("  Since a cleared test calls a positive if EITHER marker is elevated,")
    print("  a device that does UCH-L1 well is already clinically meaningful.")

    print("\n  Measured concentrations for reference:\n")
    print("  scenario                marker      pg/mL        molar")
    print("  " + "-" * 60)
    for scen, marker, pg in SCENARIOS:
        print(f"  {scen:22s} {marker:8s} {pg:7d}    {molar(pg, marker)*1e12:7.2f} pM")

    print("\n" + "=" * 74)
    print("2. THE CLOCK: CAN IT ANSWER IN THE 15 MINUTES A POINT-OF-CARE TEST HAS?")
    print("=" * 74)
    print("""
  A cleared i-STAT cartridge returns a result in about 15 minutes from 20 uL of
  whole blood. That is the bar. Below, the time for a tight binder to reach 10%
  occupancy at a diffusion-limited on-rate, which is the physical best case.
""")
    print("  scenario                          conc     k_on=1e6    k_on=1e7")
    print("  " + "-" * 68)
    rows = [("GFAP, mTBI cutoff", "GFAP", 30),
            ("GFAP, CT-positive mTBI", "GFAP", 739),
            ("GFAP, moderate-severe", "GFAP", 764),
            ("UCH-L1, mTBI cutoff", "UCH-L1", 360),
            ("UCH-L1, moderate-severe", "UCH-L1", 705)]
    for label, marker, pg in rows:
        c = molar(pg, marker)
        t6 = time_to_occupancy(c, 1e6) / 60.0
        t7 = time_to_occupancy(c, 1e7) / 60.0
        flag = "  <- fits" if t7 <= 15 else ""
        print(f"  {label:32s} {c*1e12:6.2f} pM {t6:8.0f} min {t7:8.1f} min{flag}")

    c_gfap = molar(30, "GFAP")
    print(f"\n  At the GFAP mild-TBI threshold the best case is "
          f"{time_to_occupancy(c_gfap, 1e7)/3600:.1f} hours, and with a")
    print("  typical rather than diffusion-limited on-rate it is "
          f"{time_to_occupancy(c_gfap, 1e6)/3600:.0f} hours.")
    print("""
  This is not an engineering defect to be optimised away. It is what binding at
  sub-picomolar concentration costs in time, and affinity makes it worse rather
  than better because a tighter binder releases more slowly.

  Everything at roughly 15 pM and above fits the window: the UCH-L1 threshold,
  UCH-L1 in moderate-severe injury, and GFAP in CT-positive or moderate-severe
  injury. What does not fit is GFAP at the mild-TBI decision threshold, which
  is the case the cleared tests are built around.
""")

    print("=" * 74)
    print("3. WHAT THIS DEMANDS OF THE BINDER AND THE CAGE")
    print("=" * 74)
    print("""
  The sensor responds where its apparent dissociation constant sits, and the
  cage shifts that upward:

      Kd_apparent = Kd_binder x exp(dG_open / RT)

  To respond at a clinical threshold you need Kd_apparent near that threshold.
  Since exp(dG_open/RT) >= 1 always, the binder alone must already be at or
  below the threshold concentration, and every kcal/mol of cage stability makes
  the requirement tighter.
""")
    print("  marker    threshold    required binder Kd (dG_open = 0)   ... at 2 kcal/mol")
    print("  " + "-" * 74)
    rt = R_KCAL * 298.15
    for marker, pg in [("GFAP", 30), ("UCH-L1", 360)]:
        c = molar(pg, marker)
        print(f"  {marker:8s} {c*1e12:6.2f} pM      <= {c*1e12:6.2f} pM"
              f"                     <= {c/np.exp(2.0/rt)*1e12:6.3f} pM")
    print("""
  For scale, the tightest binder in the LucCage paper was LCB1 against the
  SARS-CoV-2 receptor binding domain at 500 pM. Meeting the UCH-L1 threshold
  needs roughly 35x better than that; meeting the GFAP mild-TBI threshold needs
  roughly 800x better. This is the single most important number to agree with
  Monod Bio before design starts, because it decides which targets are in scope.

  The corollary for the cage is that dG_open must be SMALL. Optimising the
  usable signal swing across each reportable range gives:
""")
    from thermo import switch_model as sw
    for marker in ("GFAP", "UCH-L1"):
        lo, hi = [molar(p, marker) for p in REPORTABLE[marker]]
        print(f"  {marker} reportable range {lo*1e12:.2f} to {hi*1e12:.0f} pM:")
        for kd in (1e-12, 1e-11, 1e-10, 1e-9):
            dgo, swing = sw.optimal_dg_open((lo, hi), sw.dg_from_kd(kd, 298.15), 298.15)
            print(f"    binder Kd {kd*1e12:7.1f} pM -> optimal dG_open {dgo:4.2f} kcal/mol,"
                  f" swing {swing*100:5.1f}%")
        print()
    print("  Note how far this is from the LucCage regime. That paper's only")
    print("  published free energy is dG_open = 4 kcal/mol, tuned for nanomolar")
    print("  targets. A picomolar target forces a much weaker cage, which means")
    print("  higher background and less dynamic range. The design is being pushed")
    print("  into the corner of its own trade-off, and that is a consequence of")
    print("  the target concentrations rather than of any modelling choice.")

    print("\n" + "=" * 74)
    print("4. NfL DOES NOT BELONG IN AN ACUTE TRIAGE DEVICE")
    print("=" * 74)
    print("""
  NfL has no FDA-cleared traumatic brain injury cutoff. The Siemens threshold
  of 12.9 pg/mL that circulates in the literature is a multiple-sclerosis
  disease-activity cut-point and must not be quoted as a TBI threshold.

  More decisively, its kinetics are wrong for triage. NfL peaks between ten
  days and six weeks after injury and has a blood half-life around three weeks,
  so it is still abnormal a year later. UCH-L1 peaks at about eight hours with
  a seven to nine hour half-life, and GFAP peaks at twenty to twenty-four hours
  with a twenty-four to forty-eight hour half-life. The twelve-hour window in
  the cleared labeling is set by UCH-L1's short half-life.

  NfL is a subacute and chronic marker of axonal injury. It is a good target
  for a monitoring device and a poor one for an emergency-department decision.
  Keeping all three targets in one acute cartridge conflates two different
  products.
""")

    print("=" * 74)
    print("5. WHERE THIS TECHNOLOGY ACTUALLY WINS")
    print("=" * 74)
    print("""
  Section 2 says a reagentless equilibrium sensor cannot beat a cleared
  cartridge at fifteen-minute triage on the hardest marker. That is worth
  saying plainly rather than discovering in year two.

  But the comparison is the wrong one. A sandwich immunoassay is a single
  end-point measurement: it consumes its sample and reports one number. A
  reagentless conformational-switch sensor is reversible and can be interrogated
  continuously for hours or days, which no immunoassay can do at any price.

  For traumatic brain injury the clinically unmet question is often not the
  single triage number but the TRAJECTORY: is this patient's GFAP still rising
  at hour six, is the secondary injury progressing, does the curve bend after
  intervention. The slow binding kinetics that disqualify this technology for
  triage are irrelevant for a sensor that is already equilibrated and tracking.

  That reframing also relaxes the two hardest constraints at once. Continuous
  interrogation means integration over hours, which lifts the charge bound, and
  a flowing or indwelling format removes the depletion problem that a static
  drop on a screen-printed electrode cannot escape.
""")


if __name__ == "__main__":
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    report()
    print("=" * 74)
