"""
Run the proteomic-aging-clock reproduction end to end.

  M1  synthetic Olink-shaped cohorts with a planted ground truth
  M2  NPX preprocessing, QC, panel harmonisation
  M3  the six published clocks (real weights if supplied, else surrogates)
  M4  age acceleration on a common scale
  M5  trial statistics + negative control  (M5b: power sweep)
  M6  differential abundance + aging enrichment
  M7  pathway analysis

Usage:
    python3 -m protclock_pipeline.run_all              # full run
    python3 -m protclock_pipeline.run_all --quick      # skip the power sweep
    python3 -m protclock_pipeline.run_all --figures    # also write figures/

Every number this prints comes from simulated data. See REPRODUCTION.md for
what would be needed to make it a real reproduction.
"""
import argparse
import sys

import numpy as np

from protclock_pipeline import (m1_cohort, m2_preprocess, m3_clocks,
                                m4_ageaccel, m5_trialstats, m6_enrichment,
                                m7_pathways)

N_REFERENCE = 2500


def banner(text):
    print("\n" + "=" * 68)
    print(f"  {text}")
    print("=" * 68)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true",
                    help="skip the power sweep (the slow part)")
    ap.add_argument("--figures", action="store_true",
                    help="write figures to figures/protclock/")
    ap.add_argument("--n-reference", type=int, default=N_REFERENCE)
    ap.add_argument("--weight-dir", default=None,
                    help="directory of published clock coefficient CSVs")
    args = ap.parse_args(argv)

    print("#" * 68)
    print("#  PROTEOMIC AGING CLOCKS IN A PHASE 2a TRIAL")
    print("#  reproduction scaffold for doi:10.1038/s41587-026-03286-y")
    print("#  DATA: SYNTHETIC. The trial proteome is controlled-access.")
    print("#" * 68)

    banner("M1  cohorts")
    r1 = m1_cohort.report()

    banner("M2  preprocessing")
    pre = m2_preprocess.run(n_reference=args.n_reference)
    q = pre.qc
    print(f"  shared analysis panel : {q['n_shared']} proteins "
          f"(paper: {m1_cohort.N_SHARED})")
    print(f"  trial matrix          : {pre.trial_npx.shape}")
    print(f"  reference matrix      : {pre.ref_npx.shape}")
    print(f"  imputed below LOD     : {q['trial_missing_imputed']} values")

    banner("M3  clocks")
    clocks = m3_clocks.fit_all(pre, weight_dir=args.weight_dir)
    hold = m3_clocks.evaluate_holdout(pre)
    for name, d in hold.items():
        mae = "  -  " if np.isnan(d["mae_years"]) else f"{d['mae_years']:.2f}"
        print(f"    {name:<19} held-out r = {d['r']:+.3f}   MAE = {mae} yr")

    banner("M4  age acceleration")
    scored = m4_ageaccel.score_trial(clocks, pre, verbose=True)
    val = m4_ageaccel.validate_against_truth(scored, pre)
    print("\n  recovery of the planted biological-age offset (baseline):")
    for name, d in val.items():
        print(f"    {name:<19} r = {d['r']:+.3f}")

    banner("M5  trial statistics")
    r5 = m5_trialstats.report(scored=scored, run_null=True)

    if not args.quick:
        banner("M5b power sweep")
        sweep = m5_trialstats.report_power(n_rep=6, n_reference=1500)
    else:
        sweep = None
        print("\n  (power sweep skipped: --quick)")

    banner("M6  differential abundance + enrichment")
    r6 = m6_enrichment.report(pre=pre)

    banner("M7  pathways")
    da = m6_enrichment.differential_abundance(pre, contrast="paired_treated")
    ag = m6_enrichment.aging_associated(pre)
    r7 = m7_pathways.report(pre=pre, differential=da, aging=ag)

    banner("SUMMARY")
    _summary(r1, pre, scored, r5, r6, r7, sweep, clocks, val)

    if args.figures:
        from protclock_pipeline import figures
        figures.write_all(pre, scored, r5, r6, r7, sweep)

    return 0


def _summary(r1, pre, scored, r5, r6, r7, sweep, clocks, val):
    truth = r1["truth"]
    conc = r5["concordance"]
    best = conc.sort_values("mean_diff_z").iloc[0]
    nc = r5.get("null_control", {})
    enr_vp = r6["enrichment"]
    enr_pt = r6["paired_treated"]["enrichment"]

    n_track = sum(1 for d in val.values() if d["r"] > 0.3)
    print("  WHAT THE PIPELINE RECOVERED FROM SIMULATED DATA")
    print(f"    clocks scored                  : {scored['clock'].nunique()}/6")
    print(f"    tracking the planted bio-age")
    print(f"      offset (r > 0.3 at baseline) : {n_track}/{len(val)}")
    print(f"    strongest arm                  : {best['arm']} "
          f"({best['mean_diff_z']:+.2f} SD, {best['n_agree']}/6 clocks "
          f"{best['direction']})")
    print(f"    planted strongest arm          : rento_30mg_BID")
    if nc:
        print(f"    null false-positive rate       : "
              f"{nc['false_positive_rate']:.3f} (nominal 0.05)")
    print(f"    enrichment OR, vs placebo      : "
          f"{enr_vp['odds_ratio']:.2f} "
          f"[{enr_vp['ci'][0]:.2f}, {enr_vp['ci'][1]:.2f}]")
    print(f"    enrichment OR, paired treated  : "
          f"{enr_pt['odds_ratio']:.2f} "
          f"[{enr_pt['ci'][0]:.2f}, {enr_pt['ci'][1]:.2f}]")
    print(f"    planted OR                     : {truth['odds_ratio']:.2f}")
    print(f"    paper's reported OR            : 1.74")
    rej = r7[(r7["aligned_q"] < 0.05) & (r7["aligned_npx"] < 0)]["pathway"]
    print(f"    pathways moved toward younger  : {', '.join(rej) or 'none'}")

    print("\n  WHAT THIS DOES NOT SHOW")
    print("    Nothing here is evidence about rentosertib. The drug effect,")
    print("    the dose ordering and the aging enrichment were written into")
    print("    the simulation by M1 and then measured back out. That tests")
    print("    the code, not the biology.")
    n_sur = sum(1 for c in clocks.values() if c.is_surrogate)
    print(f"\n    {n_sur} of the {len(clocks)} clocks are surrogates fitted on")
    print("    the synthetic reference, not the published models, and the")
    print("    trial proteome is controlled access. See REPRODUCTION.md.")


if __name__ == "__main__":
    sys.exit(main())
