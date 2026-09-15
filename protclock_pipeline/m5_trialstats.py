"""
Module 5 - Trial statistics on age acceleration.

Three analyses, matching what the paper's design supports:

  1. WITHIN-PATIENT CHANGE. Each patient is their own control:
         delta_i = accel(week 12) - accel(week 0)
     then each active arm is compared with placebo by Welch t-test. This is
     the primary readout. It has to be within-patient because M4 showed the
     arms are not balanced at baseline and every patient carries a large
     fixed biological-age offset; a cross-sectional comparison at week 12
     would be measuring that offset, not the drug.

  2. MIXED-EFFECTS MODEL over all four visits:
         accel ~ week + week:arm,  random intercept per patient
     The week:arm coefficient is the per-week divergence of an arm from
     placebo, and it uses the week 2 and week 4 visits the delta analysis
     throws away.

  3. CROSS-CLOCK CONCORDANCE. The paper's central claim is that six
     independently built clocks agree in direction. Sign agreement is counted
     per arm and given a binomial p-value.

CAVEAT on that binomial p-value: it assumes the six clocks are independent.
They are not. They are trained on overlapping proteins from overlapping
cohorts, so the true probability of accidental agreement is higher than the
nominal figure and the p-value is anticonservative. It is reported as a
descriptive statistic, not as evidence.

NEGATIVE CONTROL: `null_control()` reruns the whole pipeline on data built
with no drug effect at all. If that run produces significant results, the
pipeline is manufacturing signal and nothing else it says can be trusted.
"""
import warnings

import numpy as np
import pandas as pd
from scipy import stats

ACTIVE_ARMS = ("rento_30mg_QD", "rento_30mg_BID", "rento_60mg_QD")
PLACEBO = "placebo"
PRIMARY_WEEK = 12


def within_patient_delta(scored, week=PRIMARY_WEEK):
    """
    Per patient and clock, age acceleration at `week` minus at baseline.

    Both scales are carried through: `delta_z` for inference and
    `delta_years` for interpretation (see M4 on why years is unstable for
    the mortality clocks).
    """
    base = scored[scored["week"] == 0].set_index(["patient_id", "clock"])
    late = scored[scored["week"] == week].set_index(["patient_id", "clock"])
    common = base.index.intersection(late.index)

    out = pd.DataFrame({
        "delta_years": late.loc[common, "accel_years"]
                       - base.loc[common, "accel_years"],
        "delta_z": late.loc[common, "accel_z"] - base.loc[common, "accel_z"],
        "arm": late.loc[common, "arm"],
    }).reset_index()
    return out


def delta_vs_placebo(delta, value="delta_z"):
    """Welch t-test of each active arm against placebo, per clock."""
    rows = []
    for clock, grp in delta.groupby("clock"):
        ctrl = grp.loc[grp["arm"] == PLACEBO, value].to_numpy()
        yrs = grp.groupby("arm")["delta_years"].mean()
        for arm in ACTIVE_ARMS:
            act = grp.loc[grp["arm"] == arm, value].to_numpy()
            t, p = stats.ttest_ind(act, ctrl, equal_var=False)
            diff = act.mean() - ctrl.mean()
            # Welch-Satterthwaite CI on the difference of means.
            se = np.sqrt(act.var(ddof=1) / len(act) + ctrl.var(ddof=1) / len(ctrl))
            df = se ** 4 / (
                (act.var(ddof=1) / len(act)) ** 2 / (len(act) - 1)
                + (ctrl.var(ddof=1) / len(ctrl)) ** 2 / (len(ctrl) - 1))
            crit = stats.t.ppf(0.975, df)
            rows.append({"clock": clock, "arm": arm,
                         "mean_delta_arm": act.mean(),
                         "mean_delta_placebo": ctrl.mean(),
                         "diff_z": diff,
                         "diff_years": float(yrs[arm] - yrs[PLACEBO]),
                         "ci_low": diff - crit * se,
                         "ci_high": diff + crit * se,
                         "t": t, "p": p,
                         "n_arm": len(act), "n_placebo": len(ctrl)})
    return pd.DataFrame(rows)


def mixed_model(scored):
    """
    accel ~ week + week:arm with a per-patient random intercept.

    The arm main effect is deliberately omitted: arm is a between-patient
    factor, so it is absorbed by the random intercept and is not identifiable
    alongside it. The estimand is the week:arm interaction.
    """
    import statsmodels.formula.api as smf

    rows = []
    for clock, grp in scored.groupby("clock"):
        d = grp.copy()
        d["arm"] = pd.Categorical(d["arm"],
                                  categories=[PLACEBO] + list(ACTIVE_ARMS))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                fit = smf.mixedlm("accel_z ~ week + week:arm", d,
                                  groups=d["patient_id"]).fit(reml=True)
            except Exception as exc:                      # pragma: no cover
                rows.append({"clock": clock, "arm": None, "coef": np.nan,
                             "p": np.nan, "note": f"failed: {exc}"})
                continue
        for arm in ACTIVE_ARMS:
            key = [k for k in fit.params.index
                   if k.startswith("week:arm") and arm in k]
            if not key:
                continue
            rows.append({"clock": clock, "arm": arm,
                         "coef_z_per_week": float(fit.params[key[0]]),
                         "p": float(fit.pvalues[key[0]]),
                         "note": ""})
    return pd.DataFrame(rows)


def concordance(stats_tab):
    """
    Per arm, how many of the six clocks point the same way as the majority.

    Returns the count, the direction, and a two-sided binomial p-value under
    the (false, see module docstring) assumption of independent clocks.
    """
    rows = []
    for arm, grp in stats_tab.groupby("arm"):
        signs = np.sign(grp["diff_z"].to_numpy())
        n = len(signs)
        n_neg = int((signs < 0).sum())
        n_agree = max(n_neg, n - n_neg)
        direction = "younger" if n_neg >= n - n_neg else "older"
        p = stats.binomtest(n_agree, n, 0.5, alternative="two-sided").pvalue
        rows.append({"arm": arm, "n_clocks": n, "n_agree": n_agree,
                     "direction": direction,
                     "mean_diff_z": grp["diff_z"].mean(),
                     "mean_diff_years": grp["diff_years"].mean(),
                     "binom_p_ANTICONSERVATIVE": p})
    return pd.DataFrame(rows)


def run(scored):
    delta = within_patient_delta(scored)
    tab = delta_vs_placebo(delta)
    conc = concordance(tab)
    mixed = mixed_model(scored)
    return {"delta": delta, "vs_placebo": tab, "concordance": conc,
            "mixed": mixed}


def null_control(n_reference=1500, seed=0, n_rep=5):
    """
    Rerun everything with effect_scale=0, i.e. no drug effect whatsoever.

    Averaged over `n_rep` independent redraws. A single null replicate is a
    coin toss: 3 of 18 tests reaching p<0.05 by chance is unremarkable, so a
    one-draw negative control cannot tell a calibrated pipeline from a leaky
    one. The false-positive RATE across replicates can.
    """
    from protclock_pipeline import m2_preprocess, m3_clocks, m4_ageaccel

    rates, shifts, tests = [], [], 0
    for rep in range(n_rep):
        pre = m2_preprocess.run(n_reference=n_reference, effect_scale=0.0,
                                seed=seed + 101 + 613 * rep)
        clocks = m3_clocks.fit_all(pre, verbose=False)
        scored = m4_ageaccel.score_trial(clocks, pre)
        res = run(scored)
        tab = res["vs_placebo"]
        rates.append(float((tab["p"] < 0.05).mean()))
        shifts.append(float(res["concordance"]["mean_diff_z"].abs().max()))
        tests = len(tab)
    return {"n_tests": tests, "n_rep": n_rep,
            "false_positive_rate": float(np.mean(rates)),
            "max_abs_shift_z": float(np.mean(shifts))}


def power_sweep(scales=(0.0, 0.5, 1.0, 2.0, 4.0), n_rep=3,
                n_reference=1500, seed=0):
    """
    How large must the effect be before this 42-patient design detects it?

    For each effect scale the trial is regenerated `n_rep` times with fresh
    noise and the pipeline is rerun end to end. Reported per scale:
      detect_rate  fraction of arm x clock tests with p < 0.05
      concordant   fraction of replicates where >= 5 of 6 clocks agree in
                   direction on the strongest arm
      mean_shift_z mean standardised shift of the 30 mg BID arm

    Scale 0.0 is the null. This is the honest answer to "did the pipeline
    work": not a single significant p-value, but a curve showing the effect
    size at which detection becomes reliable.
    """
    from protclock_pipeline import m2_preprocess, m3_clocks, m4_ageaccel

    rows = []
    for scale in scales:
        det, conc_hits, shifts, best_hits = [], [], [], []
        for rep in range(n_rep):
            pre = m2_preprocess.run(n_reference=n_reference,
                                    effect_scale=scale,
                                    seed=seed + 1000 * rep + 7)
            clocks = m3_clocks.fit_all(pre, verbose=False)
            scored = m4_ageaccel.score_trial(clocks, pre)
            res = run(scored)
            tab, conc = res["vs_placebo"], res["concordance"]
            det.append(float((tab["p"] < 0.05).mean()))
            bid = conc[conc["arm"] == "rento_30mg_BID"].iloc[0]
            conc_hits.append(bid["n_agree"] >= 5 and bid["direction"] == "younger")
            shifts.append(float(bid["mean_diff_z"]))
            best_hits.append(conc.sort_values("mean_diff_z").iloc[0]["arm"]
                             == "rento_30mg_BID")
        rows.append({"effect_scale": scale,
                     "planted_years": scale * 30.0,
                     "detect_rate": float(np.mean(det)),
                     "concordant_frac": float(np.mean(conc_hits)),
                     "best_arm_correct": float(np.mean(best_hits)),
                     "mean_shift_z_BID": float(np.mean(shifts)),
                     "n_rep": n_rep})
    return pd.DataFrame(rows)


def report(scored=None, pre=None, clocks=None, run_null=True):
    from protclock_pipeline import m2_preprocess, m3_clocks, m4_ageaccel
    if scored is None:
        if pre is None:
            pre = m2_preprocess.run(n_reference=2500)
        if clocks is None:
            clocks = m3_clocks.fit_all(pre, verbose=False)
        scored = m4_ageaccel.score_trial(clocks, pre)

    res = run(scored)

    print("M5  TRIAL STATISTICS ON AGE ACCELERATION")
    print("-" * 68)
    print("  (1) within-patient change, baseline -> week 12, vs placebo")
    print("      standardised (reference SD units); negative = YOUNGER\n")
    piv = res["vs_placebo"].pivot(index="clock", columns="arm",
                                  values="diff_z")
    print(piv.round(3).to_string())
    print("\n      same contrast in years (unstable for mortality clocks):")
    print(res["vs_placebo"].pivot(index="clock", columns="arm",
                                  values="diff_years").round(2).to_string())

    print("\n  (2) cross-clock concordance per arm")
    print(res["concordance"].round(4).to_string(index=False))

    print("\n  (3) mixed model  accel_z ~ week + week:arm + (1|patient)")
    mpiv = res["mixed"].pivot(index="clock", columns="arm",
                              values="coef_z_per_week")
    print(mpiv.round(4).to_string())

    best = res["concordance"].sort_values("mean_diff_z").iloc[0]
    print(f"\n  strongest arm by mean shift: {best['arm']}  "
          f"({best['mean_diff_z']:+.3f} SD, {best['n_agree']}/"
          f"{best['n_clocks']} clocks {best['direction']})")
    print("  PLANTED strongest arm was rento_30mg_BID (dose weight 1.00).")
    print("  Do NOT read the arm ranking off one replicate: with 10-11")
    print("  patients per arm the ranking is unstable even when the overall")
    print("  effect is detected. m5b's power sweep measures how often it is")
    print("  recovered.")

    if run_null:
        print("\n  NEGATIVE CONTROL - same pipeline, drug effect set to zero:")
        nc = null_control()
        print(f"    false-positive rate over {nc['n_rep']} redraws x "
              f"{nc['n_tests']} tests: {nc['false_positive_rate']:.3f}  "
              f"(nominal 0.05)")
        print(f"    mean of the largest |arm shift| per replicate: "
              f"{nc['max_abs_shift_z']:.3f} SD")
        res["null_control"] = nc

    return res


def report_power(**kw):
    print("M5b POWER SWEEP - what can 42 patients actually detect?")
    print("-" * 68)
    sweep = power_sweep(**kw)
    print(sweep.round(3).to_string(index=False))
    print("\n  detect_rate is the fraction of the 18 arm x clock tests"
          " reaching p<0.05.")
    print("  At effect_scale=0 it should sit near 0.05; anything much"
          " higher means")
    print("  the pipeline invents signal.")
    return sweep


if __name__ == "__main__":
    report()
