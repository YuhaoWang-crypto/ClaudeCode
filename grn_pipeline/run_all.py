"""
Run the full pipeline end to end and emit a consolidated summary.

Pipeline:
  M1  symmetry reduction  (automorphism group -> quotient -> irreducible core)
  M2  CRNT deficiency     (topology gates switching capacity / bistability)
  M3  elementary flux modes (irreducible generators of the flux cone)
  M4  DNB / critical slowing down / Lyapunov exponent (tipping biomarker)
"""
from grn_pipeline import (m1_symmetry, m2_crnt, m3_efm, m4_dnb_lyapunov,
                          m5_kras_real, m6_integrate)


def main():
    print("\n" + "#" * 68)
    print("#  GENE-REGULATORY / METABOLIC NETWORK — MATHEMATICAL PIPELINE")
    print("#" * 68 + "\n")

    r1 = m1_symmetry.report();          print()
    m2_crnt.report();                   print()
    r3 = m3_efm.report();               print()
    r4 = m4_dnb_lyapunov.report();      print()
    r5 = m5_kras_real.report();         print()
    r6 = m6_integrate.report();         print()

    print("=" * 68)
    print("CONSOLIDATED SUMMARY")
    print("=" * 68)
    print(f"M1  symmetry : |Aut(G)|={r1['group_order']}, "
          f"network {r1['n_nodes']}->{r1['n_core_nodes']} nodes "
          f"(3 RAS paralogues -> 1 core node)")
    print(f"M2  CRNT     : deficiency delta gates switching; delta=0 weakly "
          f"reversible => monostable, delta>=1 => bistable switch")
    print(f"M3  EFM      : {len(r3['efms'])} irreducible flux generators "
          f"span all steady-state flux")
    print(f"M4  biomarker: LLE code validated on Rossler = "
          f"{r4['lle_rossler']:+.4f}; leading eigenvalue -> 0 with rising "
          f"SD/autocorr/DNB at the tipping point")
    print(f"M5  real KRAS: covalent G12C drug breaks paralog symmetry "
          f"S_3(6)->S_2({r5['drugged']['order']}); sotorasib 1217x "
          f"G12C-selective (ChEMBL)")
    print(f"M6  integrate: ChEMBL+Boltz binding -> engagement -> mu -> "
          f"network-stability biomarker (adagrasib slightly closer to edge)")
    print("\nAbstract, measurable biomarker candidates produced:")
    print("  * irreducible-core node identity        (M1 quotient)")
    print("  * deficiency delta / distance-to-bistability (M2)")
    print("  * EFM usage / rate-limiting generator    (M3)")
    print("  * leading Lyapunov exponent & DNB index  (M4)")
    print("  * paralog-symmetry order as drug-selectivity readout (M5)")
    print("  * binding-driven network-stability score (M6)")


if __name__ == "__main__":
    main()
