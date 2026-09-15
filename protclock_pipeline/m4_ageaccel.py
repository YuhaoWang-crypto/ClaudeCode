"""
Module 4 - Age acceleration, on one scale for all six clocks.

A clock's raw output is not comparable across clocks: ProtAge emits years,
OrganAge_mortality emits a log hazard, PAC emits a mortality-informed score.
The paper's whole design rests on reading six clocks side by side, so every
clock has to be projected onto a common axis first.

Calibration, done once in the reference cohort and then frozen:

  1. Regress the clock's raw output on chronological age and sex:
         pred = a + b*age + c*sex + residual
  2. AGE ACCELERATION is that residual. It is the part of the proteomic
     signal not explained by how old the subject actually is.
  3. Divide by b, the clock's output per year of age, to express the residual
     in YEARS. A mortality clock's b is in logit units per year, so the
     quotient is a years-equivalent and the six clocks become comparable.

The reference model is fit on reference subjects ONLY and then applied
unchanged to the trial. Re-fitting it on trial samples would absorb the
treatment effect into the intercept and erase what the trial is measuring.

Validation - because M1 planted a true biological-age offset per subject,
age acceleration can be checked against it directly. A clock whose
acceleration does not track the planted offset is not measuring anything.
"""
import numpy as np
import pandas as pd


class AgeAccelCalibration:
    """Frozen reference-cohort calibration for one clock."""

    def __init__(self, name, coef, intercept, slope_per_year, resid_sd):
        self.name = name
        self.coef = coef                    # [age, sex]
        self.intercept = intercept
        self.slope_per_year = slope_per_year
        self.resid_sd = resid_sd            # SD of the reference residual

    def residual(self, pred, age, sex):
        expected = self.intercept + self.coef[0] * age + self.coef[1] * sex
        return pred - expected

    def years(self, pred, age, sex):
        """
        Age acceleration in years.

        Interpretable, but UNSTABLE for mortality-trained clocks: their output
        is a log hazard whose slope against age is small (order 0.01 logit per
        year), and dividing by a small slope inflates the noise by the same
        factor. In the synthetic run this makes the years-scale SD of
        OrganAge_mortality roughly 20 years against 2-5 for the age clocks.
        Use `z` for anything inferential and `years` for interpretation.
        """
        return self.residual(pred, age, sex) / self.slope_per_year

    def z(self, pred, age, sex):
        """
        Age acceleration in reference-cohort SD units.

        Standardising by the residual SD rather than the age slope puts all
        six clocks on one axis without the small-slope amplification above,
        so this is the scale the statistics in M5 run on.
        """
        return self.residual(pred, age, sex) / self.resid_sd


def calibrate(clock, pre):
    """Fit the age+sex reference model for one clock and freeze it."""
    pred = clock.predict(pre.ref_npx)
    age = pre.ref_meta["age"].to_numpy()
    sex = pre.ref_meta["sex"].to_numpy()

    x = np.column_stack([np.ones_like(age), age, sex])
    beta, *_ = np.linalg.lstsq(x, pred, rcond=None)
    slope = beta[1]
    resid_sd = float(np.std(pred - x @ beta, ddof=3))

    if abs(slope) < 1e-9:
        raise RuntimeError(
            f"{clock.spec.name}: output does not vary with age in the "
            f"reference cohort, so it cannot be converted to years")
    if resid_sd < 1e-9:
        raise RuntimeError(
            f"{clock.spec.name}: reference residual has no variance")

    return AgeAccelCalibration(clock.spec.name, coef=beta[1:],
                               intercept=beta[0], slope_per_year=slope,
                               resid_sd=resid_sd)


def score_trial(clocks, pre, verbose=False):
    """
    Return a tidy frame: one row per (sample, clock) with age acceleration.

    Columns: sample_id, patient_id, arm, week, age, sex, clock,
             pred_raw, accel_years, accel_z.
    """
    meta = pre.trial_meta
    rows = []
    for name, clk in clocks.items():
        cal = calibrate(clk, pre)
        pred = clk.predict(pre.trial_npx)
        age_v, sex_v = meta["age"].to_numpy(), meta["sex"].to_numpy()
        block = meta.copy()
        block["clock"] = name
        block["pred_raw"] = pred
        block["accel_years"] = cal.years(pred, age_v, sex_v)
        block["accel_z"] = cal.z(pred, age_v, sex_v)
        rows.append(block)
        if verbose:
            print(f"    {name:<19} slope = {cal.slope_per_year:+8.4f} /yr   "
                  f"resid SD = {cal.resid_sd:7.3f}   "
                  f"1 SD = {cal.resid_sd / abs(cal.slope_per_year):6.2f} yr")
    return pd.concat(rows, ignore_index=True)


def validate_against_truth(scored, pre):
    """
    Does age acceleration track the biological-age offset M1 planted?

    Uses baseline (week 0) samples only, because from week 2 onward the drug
    term is also moving the proteome and would contaminate the comparison.
    """
    base = scored[scored["week"] == 0]
    out = {}
    for name, grp in base.groupby("clock"):
        r = float(np.corrcoef(grp["accel_years"],
                              grp["bio_offset_true"])[0, 1])
        slope = float(np.polyfit(grp["bio_offset_true"],
                                 grp["accel_years"], 1)[0])
        out[name] = {"r": r, "slope": slope, "n": len(grp)}
    return out


def report(pre=None, clocks=None):
    from protclock_pipeline import m2_preprocess, m3_clocks
    if pre is None:
        pre = m2_preprocess.run(n_reference=2500)
    if clocks is None:
        clocks = m3_clocks.fit_all(pre, verbose=False)

    print("M4  AGE ACCELERATION (all six clocks, in years)")
    print("-" * 68)
    print("  reference calibration  pred ~ a + b*age + c*sex, frozen:")
    scored = score_trial(clocks, pre, verbose=True)

    print("\n  recovery of the PLANTED biological-age offset "
          "(baseline visits only):")
    val = validate_against_truth(scored, pre)
    for name, d in val.items():
        print(f"    {name:<19} r = {d['r']:+.3f}   slope = {d['slope']:+.3f} "
              f"yr/yr   n = {d['n']}")

    ok = [n for n, d in val.items() if d["r"] > 0.3]
    print(f"\n  clocks whose acceleration tracks the planted offset: "
          f"{len(ok)}/{len(val)}")
    print("  baseline age acceleration by arm (should be BALANCED, "
          "randomisation check):")
    base = scored[scored["week"] == 0]
    tab = base.pivot_table(index="arm", columns="clock",
                           values="accel_years", aggfunc="mean")
    print(tab.round(2).to_string())
    return {"scored": scored, "validation": val}


if __name__ == "__main__":
    report()
