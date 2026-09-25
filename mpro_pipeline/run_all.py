"""
Run the Mpro enzymatic pipeline end to end.

  M1  label QC          pool + quarantine the discordant tail; preincub tiers
  M2  null-model gate   the descriptor bar every score must clear
  M3  Boltz calibration pre-registered criteria; which metric to trust
  M4  enzymatic QSAR     scaffold-split, null-gated, tier-reported

Scope: the purified-enzyme axis only (a BPS #79955-type FRET assay). Cellular
antiviral activity is a different axis and is deliberately out of scope --
enzymatic pIC50 explains only ~14% of the variance in cellular potency
(rho = +0.43 over 615 paired compounds), and the missing part is
permeability / efflux / glutathione / esterase, not anything this pipeline
can see.
"""
from mpro_pipeline import m1_labels, m2_null_model, m3_boltz_calib, m4_qsar


def main():
    print("\n" + "#" * 70)
    print("#  SARS-CoV-2 Mpro (3CLpro) — ENZYMATIC QSAR PIPELINE")
    print("#" * 70 + "\n")

    r1 = m1_labels.report();      print()
    r2 = m2_null_model.report();  print()
    r3 = m3_boltz_calib.report(); print()
    r4 = m4_qsar.report();        print()

    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"trainable compounds        {len(r1['trainable'])}")
    print(f"null-model bar (calib set) |rho| = {r2['best_null']:.3f}")
    print(f"Boltz optimization_score   rho = {r3['rho']:+.3f} "
          f"CI [{r3['ci'][0]:+.3f}, {r3['ci'][1]:+.3f}]  -> trust for ranking")
    best = max(r4, key=lambda k: r4[k]["rho"])
    print(f"QSAR ({best})       rho = {r4[best]['rho']:+.3f}  "
          f"RMSE = {r4[best]['rmse']:.2f} log  margin over null "
          f"{r4[best]['margin']:+.3f}")
    print("\nnot established by this pipeline:")
    print("  - dimer > monomer (direction consistent, CI on delta spans 0)")
    print("  - stereo-SAR (Boltz scores the nirmatrelvir epimer within 0.008)")
    print("  - anything cellular (different axis, see module docstrings)")


if __name__ == "__main__":
    main()
