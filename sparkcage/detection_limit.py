"""
What sets the detection limit for GFAP, UCH-L1 and NfL on a carbon electrode.

    python3 sparkcage/detection_limit.py

The three named targets are traumatic-brain-injury biomarkers measured in
plasma at picogram-per-millilitre concentrations.  That single fact imposes
constraints that no amount of protein design can relax, and they are worth
settling before any money is spent on cages and latches.

The chain of reasoning:

1. A switched sensor gives up a fixed amount of charge (2 electrons for the
   methylene blue / leucomethylene blue couple).
2. No sensor can switch more molecules than the sample contains.
3. Therefore the maximum modulated charge is set by the SAMPLE, not the sensor:
       dQ_max = 2 e c V N_A
   independent of electrode area, surface coverage, binder affinity, cage
   stability, or interrogation frequency.
4. Separately, if the immobilised sensor population outnumbers the analyte, the
   fractional occupancy collapses even for a perfect binder, because the sensor
   itself depletes the sample.

Together these make sample volume and electrode size first-class design
parameters, on a par with anything in the protein engineering.
"""

from __future__ import annotations

import numpy as np

N_A = 6.02214076e23
E_CHARGE = 1.602176634e-19
N_ELECTRONS = 2          # methylene blue / leucomethylene blue

# UniProt masses, verified: P14136, P09936, P07196
TARGETS = {
    "GFAP": dict(uniprot="P14136", mass_da=49880, length=432),
    "UCH-L1": dict(uniprot="P09936", mass_da=24824, length=223),
    "NfL": dict(uniprot="P07196", mass_da=61518, length=543),
}

# Stokes-Einstein at 37 C with the standard globular-protein radius scaling
# Rh(nm) = 0.0515 * M(Da)^0.392 .  NfL is an intermediate-filament protein and
# is elongated, so its true D is lower than this and the value is an upper bound.
KB = 1.380649e-23
T_K = 310.15
ETA = 0.7e-3


def hydrodynamic_radius_m(mass_da):
    return 0.0515 * mass_da ** 0.392 * 1e-9


def diffusion_cm2_s(mass_da):
    return KB * T_K / (6.0 * np.pi * ETA * hydrodynamic_radius_m(mass_da)) * 1e4


def molar_from_pg_ml(pg_ml, mass_da):
    return (pg_ml * 1e-12) / mass_da * 1000.0


def molecules(pg_ml, mass_da, volume_l):
    return molar_from_pg_ml(pg_ml, mass_da) * volume_l * N_A


def max_modulated_charge_c(pg_ml, mass_da, volume_l):
    """Upper bound on the faradaic charge the analyte can modulate, in coulomb."""
    return N_ELECTRONS * E_CHARGE * molecules(pg_ml, mass_da, volume_l)


def concentration_for_charge(charge_c, mass_da, volume_l):
    """Concentration (pg/mL) whose charge bound equals `charge_c`."""
    c_molar = charge_c / (N_ELECTRONS * E_CHARGE * volume_l * N_A)
    return c_molar * mass_da / 1000.0 * 1e12


def sensor_count(area_cm2, gamma_mol_cm2):
    return gamma_mol_cm2 * area_cm2 * N_A


def max_area_before_depletion(pg_ml, mass_da, volume_l, gamma_mol_cm2,
                              sensor_to_analyte=0.1):
    """Largest electrode that keeps the sensor from depleting the sample."""
    n_analyte = molecules(pg_ml, mass_da, volume_l)
    return sensor_to_analyte * n_analyte / (gamma_mol_cm2 * N_A)


def report():
    print("=" * 74)
    print("1. TARGET PROPERTIES")
    print("=" * 74)
    print("\n  target    UniProt   length    mass      R_h      D (cm^2/s)")
    print("  " + "-" * 60)
    for name, t in TARGETS.items():
        rh = hydrodynamic_radius_m(t["mass_da"]) * 1e9
        print(f"  {name:8s}  {t['uniprot']}   {t['length']:5d}  "
              f"{t['mass_da']/1000:6.1f} kDa  {rh:4.2f} nm  "
              f"{diffusion_cm2_s(t['mass_da']):.2e}")
    print("\n  NfL is an intermediate-filament protein, not globular, so its true")
    print("  diffusion coefficient is lower than this estimate.")

    print("\n" + "=" * 74)
    print("2. THE CHARGE BOUND: A CEILING NO SENSOR DESIGN CAN LIFT")
    print("=" * 74)
    print("\n  Maximum modulated charge for GFAP, by sample volume:\n")
    print("   concentration      10 uL          100 uL           1 mL")
    print("  " + "-" * 62)
    for pg in [1, 10, 100, 1000, 10000]:
        row = [max_modulated_charge_c(pg, TARGETS["GFAP"]["mass_da"], v)
               for v in (10e-6, 100e-6, 1e-3)]
        print(f"  {pg:7d} pg/mL " + "  ".join(f"{q*1e12:11.3f} pC" for q in row))

    print("\n  For scale, a typical E-AB square-wave peak carries 1 to 300 nC,")
    print("  and a good potentiostat resolves roughly 1 pC with averaging.")
    print("\n  Concentration at which the bound equals 1 pC:\n")
    print("   target        10 uL         100 uL          1 mL")
    print("  " + "-" * 58)
    for name, t in TARGETS.items():
        row = [concentration_for_charge(1e-12, t["mass_da"], v)
               for v in (10e-6, 100e-6, 1e-3)]
        print(f"  {name:8s} " + "  ".join(f"{v:10.2f} pg/mL" for v in row))
    print("\n  A 100 uL sample puts the ceiling near 1 to 3 pg/mL, and that")
    print("  assumes every analyte molecule is captured, every captured molecule")
    print("  fully switches a sensor, and the measurement is noise-free. Capture")
    print("  efficiency, modulation depth and noise each cost about a decade, so")
    print("  a realistic floor is tens to hundreds of pg/mL at 100 uL.")

    print("\n" + "=" * 74)
    print("3. THE DEPLETION CONSTRAINT: THE ELECTRODE CAN BE TOO BIG")
    print("=" * 74)
    gamma = 5e-12
    area = 0.03
    n_sensor = sensor_count(area, gamma)
    print(f"\n  A {area} cm^2 electrode at Gamma = {gamma:.0e} mol/cm^2 carries")
    print(f"  {n_sensor:.2e} sensor molecules.\n")
    print("  GFAP molecules in 100 uL, and the occupancy they could reach even")
    print("  if binding were irreversible and complete:\n")
    print("   concentration    molecules      max occupancy")
    print("  " + "-" * 48)
    for pg in [1, 10, 100, 1000, 10000]:
        n = molecules(pg, TARGETS["GFAP"]["mass_da"], 100e-6)
        print(f"  {pg:7d} pg/mL   {n:.2e}      {n/n_sensor*100:8.4f}%")
    print("\n  At 100 pg/mL the sensors outnumber the analyte about 750 to 1, so")
    print("  occupancy cannot exceed roughly 0.1% and the signal change cannot")
    print("  exceed 0.1% of full modulation. The E-AB noise floor is about 1%.")
    print("  No binder and no cage fixes this; only fewer sensors or more sample.")

    print("\n  Largest electrode that keeps sensors at or below 10% of the analyte")
    print("  count, 100 uL of GFAP, full monolayer:\n")
    print("   concentration      max area        equivalent disc")
    print("  " + "-" * 52)
    for pg in [1, 10, 100, 1000, 10000]:
        a = max_area_before_depletion(pg, TARGETS["GFAP"]["mass_da"], 100e-6, gamma)
        d = 2.0 * np.sqrt(a / np.pi) * 1e4
        print(f"  {pg:7d} pg/mL   {a:.3e} cm^2   {d:9.2f} um")

    print("\n" + "=" * 74)
    print("4. WHAT THIS MEANS FOR THE PROPOSED FORMAT")
    print("=" * 74)
    print("""
  The two constraints pull in opposite directions. The electrode must be small
  enough that its sensor population does not swamp the analyte, and large
  enough that the switched fraction carries measurable charge. For picogram
  targets those two requirements meet near the clinical range, which is the
  real reason this is hard.

  A screen-printed macroelectrode is on the wrong side of the depletion
  constraint by two to three orders of magnitude at clinically relevant GFAP,
  UCH-L1 and NfL concentrations. Three ways out, in order of how much they buy:

    - increase sample volume; 10 uL to 1 mL is two decades, and it is free
    - shrink the working electrode toward tens of microns, at the cost of
      absolute current and a harder electronics problem
    - pre-concentrate or capture-and-release the analyte before interrogation

  None of these is a protein-design question, and all of them should be settled
  before the cage and latch work starts, because they determine what dG_open
  and what binder affinity are worth designing for.
""")


if __name__ == "__main__":
    report()
    print("=" * 74)
