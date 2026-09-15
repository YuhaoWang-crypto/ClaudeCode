"""
Module 6 - Differential abundance and aging enrichment.

The paper's second claim is that the proteins rentosertib moves are enriched
for aging-associated proteins: 326 proteins altered, and those proteins are
1.74 times more likely than background to be aging-linked.

This module reproduces that logic without ever consulting M1's ground truth:

  1. DRUG-ALTERED SET. Per protein, take each patient's within-patient change
     from baseline to week 12, then test active arms against placebo
     (Welch t-test). Benjamini-Hochberg FDR across all 2,832 proteins.

  2. AGING SET. Per protein, regress NPX on age and sex in the REFERENCE
     cohort and keep the proteins with a significant age term, again with
     BH FDR. This is discovered from data, not read off the simulation.

  3. ENRICHMENT. Fisher exact test on the 2x2 of altered x aging, giving the
     odds ratio the paper reports as 1.74.

Both sets are DISCOVERED, so the recovered odds ratio is a genuine
measurement of the synthetic cohort. It is compared against the odds ratio
M1 planted (about 1.73) to check the machinery, and against the paper's 1.74
only to show the scale matches. The paper's number cannot be reproduced
without the controlled trial data.
"""
import numpy as np
import pandas as pd
from scipy import stats

FDR = 0.05
ACTIVE_ARMS = ("rento_30mg_QD", "rento_30mg_BID", "rento_60mg_QD")


def benjamini_hochberg(p):
    """BH-adjusted p-values (q-values), same length and order as `p`."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    # Enforce monotonicity from the largest p downwards.
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.clip(ranked, 0, 1)
    return q


def patient_deltas(pre, week=12):
    """
    Per patient, NPX(week) - NPX(baseline), for every protein.

    Returns (delta matrix [n_patients x n_proteins], arm labels).
    """
    meta = pre.trial_meta
    base = meta.index[meta["week"] == 0]
    late = meta.index[meta["week"] == week]

    b = meta.loc[base].reset_index().set_index("patient_id")
    l = meta.loc[late].reset_index().set_index("patient_id")
    pids = b.index.intersection(l.index)

    delta = (pre.trial_npx[l.loc[pids, "index"].to_numpy(), :]
             - pre.trial_npx[b.loc[pids, "index"].to_numpy(), :])
    arms = l.loc[pids, "arm"].to_numpy()
    return delta, arms, np.asarray(pids)


def differential_abundance(pre, week=12, fdr=FDR, contrast="vs_placebo"):
    """
    Per-protein test of the baseline-to-week-`week` change, BH corrected.

    Two contrasts, and the choice matters more than anything else in this
    module:

      "vs_placebo"      treated change against placebo change (Welch).
                        Controls for anything that drifts over 12 weeks
                        regardless of treatment, but spends its power on an
                        11-patient control group.

      "paired_treated"  one-sample test that the treated change differs from
                        zero. Far more powerful at this sample size, but it
                        cannot separate a drug effect from disease
                        progression, regression to the mean, or assay drift,
                        because nothing holds those constant.

    The larger altered-protein counts come from the paired contrast. Any
    count quoted from it is an upper bound on what the drug did.
    """
    delta, arms, _ = patient_deltas(pre, week=week)
    treated = delta[np.isin(arms, ACTIVE_ARMS), :]
    control = delta[arms == "placebo", :]

    if contrast == "vs_placebo":
        t, p = stats.ttest_ind(treated, control, axis=0, equal_var=False)
        effect = treated.mean(0) - control.mean(0)
    elif contrast == "paired_treated":
        t, p = stats.ttest_1samp(treated, 0.0, axis=0)
        effect = treated.mean(0)
    else:
        raise ValueError(f"unknown contrast: {contrast}")

    q = benjamini_hochberg(p)
    return pd.DataFrame({
        "protein": pre.protein_ids,
        "mean_delta_treated": treated.mean(0),
        "mean_delta_placebo": control.mean(0),
        "effect_npx": effect,
        "t": t, "p": p, "q": q,
        "altered": q < fdr,
    })


def aging_associated(pre, fdr=FDR):
    """
    Proteins whose NPX tracks age in the reference cohort.

    Regresses each protein on [1, age, sex] and tests the age coefficient.
    Vectorised across all proteins in one least-squares solve.
    """
    age = pre.ref_meta["age"].to_numpy()
    sex = pre.ref_meta["sex"].to_numpy()
    x = np.column_stack([np.ones_like(age), age, sex])
    y = pre.ref_npx

    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    resid = y - x @ beta
    n, k = x.shape
    dof = n - k
    sigma2 = (resid ** 2).sum(0) / dof
    xtx_inv = np.linalg.inv(x.T @ x)
    se_age = np.sqrt(sigma2 * xtx_inv[1, 1])

    t = beta[1] / se_age
    p = 2 * stats.t.sf(np.abs(t), dof)
    q = benjamini_hochberg(p)

    return pd.DataFrame({
        "protein": pre.protein_ids,
        "age_slope": beta[1], "t": t, "p": p, "q": q,
        "is_aging": q < fdr,
    })


def enrichment(altered, aging):
    """Fisher exact test of altered x aging, plus the 2x2 it came from."""
    a = int(np.sum(altered & aging))
    b = int(np.sum(altered & ~aging))
    c = int(np.sum(~altered & aging))
    d = int(np.sum(~altered & ~aging))

    table = [[a, b], [c, d]]
    if min(a, b, c, d) == 0:
        return {"table": table, "odds_ratio": float("nan"),
                "p": float("nan"), "note": "a zero cell; OR undefined"}

    odds, p = stats.fisher_exact(table, alternative="two-sided")
    # Woolf log-OR standard error for a confidence interval.
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    lo, hi = np.exp(np.log(odds) - 1.96 * se), np.exp(np.log(odds) + 1.96 * se)
    return {"table": table, "odds_ratio": float(odds), "ci": (float(lo),
            float(hi)), "p": float(p), "n_altered": a + b, "n_aging": a + c,
            "note": ""}


def run(pre, week=12, fdr=FDR, contrast="vs_placebo"):
    da = differential_abundance(pre, week=week, fdr=fdr, contrast=contrast)
    ag = aging_associated(pre, fdr=fdr)
    enr = enrichment(da["altered"].to_numpy(), ag["is_aging"].to_numpy())
    return {"differential": da, "aging": ag, "enrichment": enr,
            "contrast": contrast}


def _recovery(discovered, truth):
    """Precision and recall of a discovered set against the planted one."""
    tp = int(np.sum(discovered & truth))
    fp = int(np.sum(discovered & ~truth))
    fn = int(np.sum(~discovered & truth))
    prec = tp / (tp + fp) if tp + fp else float("nan")
    rec = tp / (tp + fn) if tp + fn else float("nan")
    return {"tp": tp, "fp": fp, "fn": fn, "precision": prec, "recall": rec}


def report(pre=None):
    from protclock_pipeline import m1_cohort, m2_preprocess
    if pre is None:
        pre = m2_preprocess.run(n_reference=2500)

    truth = pre.panel.truth_table()
    print("M6  DIFFERENTIAL ABUNDANCE + AGING ENRICHMENT")
    print("-" * 68)

    both = {}
    for contrast in ("vs_placebo", "paired_treated"):
        res = run(pre, contrast=contrast)
        da, ag, enr = res["differential"], res["aging"], res["enrichment"]
        r_alt = _recovery(da["altered"].to_numpy(), pre.is_responsive_true)
        r_age = _recovery(ag["is_aging"].to_numpy(), pre.is_aging_true)
        res["recovery"] = {"altered": r_alt, "aging": r_age}
        both[contrast] = res

        print(f"\n  CONTRAST = {contrast}")
        print(f"    proteins tested            : {len(da)}")
        print(f"    drug-altered  (BH q<{FDR})  : {int(da['altered'].sum())}"
              f"   (paper: {m1_cohort.N_RESPONSIVE})")
        print(f"    aging-assoc   (BH q<{FDR})  : {int(ag['is_aging'].sum())}")
        (a, b), (c, d) = enr["table"]
        print(f"    2x2 [altered x aging]      : "
              f"[[{a}, {b}], [{c}, {d}]]")
        if np.isnan(enr["odds_ratio"]):
            print(f"    odds ratio                 : undefined ({enr['note']})")
        else:
            print(f"    RECOVERED odds ratio       : {enr['odds_ratio']:.3f}  "
                  f"95% CI [{enr['ci'][0]:.2f}, {enr['ci'][1]:.2f}]  "
                  f"p = {enr['p']:.2e}")
        print(f"    altered-set  precision={r_alt['precision']:.3f} "
              f"recall={r_alt['recall']:.3f}")
        print(f"    aging-set    precision={r_age['precision']:.3f} "
              f"recall={r_age['recall']:.3f}")

    print(f"\n  PLANTED odds ratio           : {truth['odds_ratio']:.3f}")
    print(f"  paper's reported value       : 1.74")
    print("\n  The two contrasts disagree on how many proteins the drug")
    print("  altered, by roughly an order of magnitude, from identical data.")
    print("  That gap is a property of a 42-patient design, not of the code.")

    out = both["vs_placebo"]
    out["paired_treated"] = both["paired_treated"]
    out["planted_or"] = truth["odds_ratio"]
    return out


if __name__ == "__main__":
    report()
