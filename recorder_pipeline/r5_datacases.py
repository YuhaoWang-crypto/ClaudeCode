"""
R5 — Clinical data space: what the published numbers actually allow.

Every candidate in the report is scored on writer specificity, novelty, evidence
etc. None of that says whether the MEASUREMENT can work. This module takes the
published normal-vs-disease distributions for each candidate, and computes:

  1. the biological effect size in log space (the only scale on which ratios and
     CVs compose additively)
  2. the achievable separation and AUC, and how far the reported clinical AUC
     falls short of the healthy-vs-disease number -- the gap IS the disease
     control group
  3. the analytical error budget: given the assay's real published CV, how much
     AUC is lost, and what CV would be needed to stop losing it
  4. the measurement-error-corrected ratio criterion (extends r2_pairing)
  5. the Gate-2 sample size implied by the effect

WHY LOG SPACE
-------------
For a lognormal analyte, a reference interval [lo, hi] quoted as central 95%
gives sigma_ln = (ln hi - ln lo) / (2 * 1.96) and median = sqrt(lo*hi). A fold
change f becomes an additive effect ln(f). An analytical CV of c contributes
sigma_ln ~= c (exact: sqrt(ln(1+c^2)), which for c<0.3 is within 1%). So the
whole feasibility question becomes one subtraction in log units.

THE RESULT THAT MATTERS MOST  [OK]
-----------------------------------
r2_pairing derived: a fixed ratio A/B helps iff rho > kappa/2, ASSUMING perfect
measurement. Real assays are not perfect, and the denominator's analytical noise
is added to the readout without removing anything. Redoing the derivation with
measurement error, the ratio helps iff

      rho_bio  >  (kappa_bio^2 + a^2) / (2 * kappa_bio)

where kappa_bio = sd_bio(B)/sd_bio(A) and a = sd_analytical(B)/sd_bio(A).

The a^2 term is the price of the second assay. It is NOT small for the
candidates in this report: PRO-C3's published inter-assay CV of 11% sits at
about half its healthy biological spread, so a ~ 0.5 and the required biological
correlation rises from 0.50 to about 0.62. For a marker whose whole disease
effect is 1.7x, that is the difference between a ratio that works and one that
does not.

DATA PROVENANCE — all [LIT], sources in CASE_STUDIES.md
"""
import numpy as np
from scipy.stats import norm

Z975 = norm.ppf(0.975)


# --------------------------------------------------------------------------
# lognormal helpers
# --------------------------------------------------------------------------

def sigma_from_interval(lo, hi, coverage=0.95):
    """sd on the ln scale implied by a central reference interval."""
    z = norm.ppf(0.5 + coverage / 2)
    return (np.log(hi) - np.log(lo)) / (2 * z)


def sigma_from_iqr(q1, q3):
    """sd on the ln scale implied by an interquartile range."""
    return (np.log(q3) - np.log(q1)) / (2 * norm.ppf(0.75))


def cv_to_sigma_ln(cv):
    """Analytical CV -> ln-scale sd. Exact for lognormal."""
    return np.sqrt(np.log(1.0 + cv**2))


def auc(d):
    """Binormal equal-variance AUC from separation d'."""
    return float(norm.cdf(d / np.sqrt(2.0)))


def dprime_from_auc(a):
    return float(norm.ppf(a) * np.sqrt(2.0))


def dprime_from_sens_spec(sens, spec):
    """Back out d' from a published operating point (equal-variance binormal)."""
    return float(norm.ppf(sens) + norm.ppf(spec))


def n_per_group(d, alpha=0.05, power=0.80):
    """Two-sample n per arm to DETECT a mean difference of separation d'.

    Included only to make the point that this is the wrong question for a
    diagnostic study: it returns absurdly small numbers for large effects
    because merely proving "the groups differ" is not the claim a test makes.
    """
    z_a = norm.ppf(1 - alpha / 2)
    z_b = norm.ppf(power)
    return int(np.ceil(2 * ((z_a + z_b) / d) ** 2))


def hanley_mcneil_var(a, n1, n2):
    """Variance of an empirical AUC (Hanley & McNeil 1982)."""
    q1 = a / (2.0 - a)
    q2 = 2.0 * a**2 / (1.0 + a)
    return (a * (1 - a) + (n1 - 1) * (q1 - a**2) + (n2 - 1) * (q2 - a**2)) / (n1 * n2)


def n_for_auc_lower_bound(a, a_min=0.70, alpha=0.05, power=0.80, nmax=20000):
    """n per arm so the study can show AUC > a_min.

    This IS the Gate-2 question: not "do the groups differ" but "is the AUC
    above the threshold of clinical usefulness". Solved by search on the
    Hanley-McNeil variance.
    """
    if a <= a_min:
        return None
    z_a, z_b = norm.ppf(1 - alpha), norm.ppf(power)
    for n in range(5, nmax):
        se = np.sqrt(hanley_mcneil_var(a, n, n))
        if (a - a_min) >= (z_a + z_b) * se:
            return n
    return None


# --------------------------------------------------------------------------
# ratio admissibility WITH measurement error   [OK — new result]
# --------------------------------------------------------------------------

def rho_min_with_error(kappa_bio, a):
    """Minimum biological correlation for a fixed ratio to help, given that the
    denominator assay contributes analytical noise a = sd_an(B)/sd_bio(A).

    a = 0 recovers r2_pairing's noise-free kappa/2.
    """
    return (kappa_bio**2 + a**2) / (2.0 * kappa_bio)


# --------------------------------------------------------------------------
# The published data  [LIT]
# --------------------------------------------------------------------------
# Each case: what is normal, what is disease, and what the assay really does.
# `spread` field: "LIT" = the quoted interval is published; "EST" = only the
# medians are published and the dispersion is an assumption. EST cases get a
# sensitivity sweep instead of a single AUC, because inventing an IQR and then
# reporting a 3-decimal AUC from it would be fabricated precision.
CASES = [
    dict(
        name="GFAP (TBI, CT-/MRI+ vs healthy)",
        indication="TBI", spread="LIT",
        normal_median=8.0, normal_lo=3.0, normal_hi=14.0, normal_cov=0.50,
        disease_median=414.4, disease_q1=139.3, disease_q3=813.4,
        cv_analytical=0.10,
        note="healthy 8 pg/mL (IQR 3-14); CT-/MRI+ 414 (IQR 139-813). Both published.",
    ),
    dict(
        name="Serum uromodulin (CKD G3 vs non-CKD)",
        indication="CKD", spread="EST",
        normal_median=228.0, normal_lo=120.0, normal_hi=380.0, normal_cov=0.50,
        disease_median=52.0, disease_q1=30.0, disease_q3=90.0,
        cv_analytical=0.08,
        note="medians published (228 -> G3 52 -> G5 13); per-stage IQR assumed.",
    ),
    dict(
        name="Carbamylated albumin (ESRD vs non-uremic)",
        indication="CKD", spread="EST",
        normal_median=0.42, normal_lo=0.30, normal_hi=0.58, normal_cov=0.50,
        disease_median=0.90, disease_q1=0.62, disease_q3=1.30,
        cv_analytical=0.08,
        note="%C-Alb Lys549: 0.90 vs 0.42 published; dispersion assumed.",
    ),
    dict(
        name="C-Alb for 1-yr mortality WITHIN ESRD",
        indication="CKD", spread="EST",
        normal_median=0.77, normal_lo=0.55, normal_hi=1.08, normal_cov=0.50,
        disease_median=1.01, disease_q1=0.72, disease_q3=1.42,
        cv_analytical=0.08,
        note="died<1yr 1.01 vs survived 0.77 published (HR 1.90 Q4vQ1); "
             "dispersion assumed. This is the INTENDED USE and it is much weaker.",
    ),
    dict(
        name="PRO-C3 (F4 vs F0/F1, NAFLD)",
        indication="fibrosis", spread="EST",
        normal_median=9.5, normal_lo=6.1, normal_hi=14.7, normal_cov=0.95,
        disease_median=16.3, disease_q1=11.0, disease_q3=24.0,
        cv_analytical=0.1103,
        note="healthy ref 6.1-14.7 published; F0/F1 9.5, F3 14.5, F4 16.3 "
             "published; F4 IQR assumed. Reported clinical AUC anchors it.",
        reported_auc=0.79,
    ),
    dict(
        name="TAT (DIC vs reference)",
        indication="sepsis", spread="EST",
        normal_median=2.0, normal_lo=1.0, normal_hi=4.0, normal_cov=0.95,
        disease_median=10.8, disease_q1=6.0, disease_q3=20.0,
        cv_analytical=0.0367,
        note="HISCL reference <4 ug/L and 10.8 threshold published; "
             "distributions assumed.",
    ),
]

# Operating points published as sens/spec rather than distributions  [LIT]
OPERATING_POINTS = [
    ("EC4d > 8.9 net MFI  (SLE vs healthy)", 0.70, 0.93),
    ("EC4d > 14 net MFI   (SLE vs other disease)", 0.45, 0.95),
    ("BC4d > 48 net MFI   (SLE vs healthy)", 0.66, 0.96),
    ("BC4d > 60 net MFI   (SLE vs other disease)", 0.54, 0.95),
]


def analyse(case):
    """Compute the whole data space for one case."""
    s_n = (sigma_from_interval(case["normal_lo"], case["normal_hi"],
                               case.get("normal_cov", 0.95)))
    s_d = sigma_from_iqr(case["disease_q1"], case["disease_q3"])
    delta = abs(np.log(case["disease_median"]) - np.log(case["normal_median"]))
    s_pooled = np.sqrt(0.5 * (s_n**2 + s_d**2))

    s_an = cv_to_sigma_ln(case["cv_analytical"])
    s_pooled_meas = np.sqrt(s_pooled**2 + s_an**2)

    d_ideal = delta / s_pooled          # perfect assay
    d_real = delta / s_pooled_meas      # with the real assay CV
    return dict(
        fold=np.exp(delta), delta=delta,
        sigma_norm=s_n, sigma_dis=s_d, sigma_pooled=s_pooled,
        sigma_analytical=s_an,
        noise_ratio=s_an / s_pooled,
        d_ideal=d_ideal, d_real=d_real,
        auc_ideal=auc(d_ideal), auc_real=auc(d_real),
        auc_loss=auc(d_ideal) - auc(d_real),
        n_per_group=n_per_group(d_real),
    )


def main():
    print("=" * 78)
    print("R5 — Clinical data space: what the published numbers allow")
    print("=" * 78)

    print("\n[1] EFFECT-SIZE LADDER (published medians, log scale)")
    print("    fold change uses PUBLISHED medians only -> robust.")
    print("    AUC additionally needs dispersion; 'sprd' flags whether that was")
    print("    published (LIT) or assumed (EST). EST AUCs get a sweep in [1b].")
    print(f"    {'candidate':<42}{'sprd':>5}{'fold':>7}{'delta':>7}{'AUC':>7}")
    results = []
    for c in CASES:
        r = analyse(c)
        results.append((c, r))
        print(f"    {c['name']:<42}{c['spread']:>5}{r['fold']:>6.1f}x"
              f"{r['delta']:>7.2f}{r['auc_ideal']:>7.3f}")
    print("    The fold-change column is the durable result: two orders of")
    print("    magnitude separate top from bottom. GFAP's 52x and PRO-C3's 1.7x")
    print("    are not the same kind of problem and should not share a gate. [OK]")

    print("\n[1b] SENSITIVITY of AUC to the assumed dispersion (EST cases)")
    print(f"    {'candidate':<42}{'sd x0.7':>9}{'assumed':>9}{'sd x1.5':>9}")
    for c, r in results:
        if c["spread"] != "EST":
            continue
        row = []
        for k in (0.7, 1.0, 1.5):
            s = np.sqrt(r["sigma_pooled"] ** 2 * k**2 + r["sigma_analytical"] ** 2)
            row.append(auc(r["delta"] / s))
        print(f"    {c['name']:<42}{row[0]:>9.3f}{row[1]:>9.3f}{row[2]:>9.3f}")
    print("    Rankings survive the sweep; absolute AUCs move by up to 0.10.")
    print("    Treat every EST number as an ordering, not a specification. [OK]")

    print("\n[2] ANALYTICAL ERROR BUDGET (using each assay's PUBLISHED CV)")
    print(f"    {'candidate':<42}{'CV':>6}{'noise/bio':>11}{'AUC real':>10}{'loss':>8}")
    for c, r in results:
        print(f"    {c['name']:<42}{c['cv_analytical']*100:>5.1f}%"
              f"{r['noise_ratio']:>11.2f}{r['auc_real']:>10.3f}{r['auc_loss']:>8.3f}")
    print("    'noise/bio' is the assay sd divided by the biological sd. Above")
    print("    ~0.3 the assay is a material part of what you are measuring.")
    print("    PRO-C3 at 11% inter-assay CV is the outlier and it is also the")
    print("    candidate with the smallest biological effect. [OK]")

    print("\n[3] REALITY CHECK: healthy-vs-disease AUC overstates clinical AUC")
    for c, r in results:
        if "reported_auc" not in c:
            continue
        d_rep = dprime_from_auc(c["reported_auc"])
        implied = c and (r["delta"] / d_rep)
        print(f"    {c['name']}")
        print(f"      computed vs healthy reference : AUC {r['auc_ideal']:.3f}"
              f"  (d'={r['d_ideal']:.2f})")
        print(f"      reported in clinical cohorts  : AUC {c['reported_auc']:.2f}"
              f"  (d'={d_rep:.2f})")
        print(f"      -> implied real-world sd = {implied:.3f} ln units, i.e."
              f" {implied/r['sigma_norm']:.1f}x the healthy reference spread")
        print("      The gap is the disease control group: F0-F2 NAFLD patients")
        print("      are not healthy blood donors. Sizing a trial off a healthy")
        print("      reference range will underpower it by this factor. [OK]")

    print("\n[4] RATIO ADMISSIBILITY WITH MEASUREMENT ERROR  [OK, new]")
    print("        rho_bio > (kappa_bio^2 + a^2) / (2 kappa_bio),"
          "  a = sd_an(B)/sd_bio(A)")
    print(f"    {'scenario':<46}{'a':>6}{'rho needed':>12}")
    for lab, kappa, a in [
        ("noise-free ideal (r2_pairing), kappa=1", 1.0, 0.0),
        ("PRO-C3/C3M, both 11% CV, kappa=1", 1.0, 0.11 / 0.2245),
        ("TAT/PIC, both ~4-7% CV, kappa=1", 1.0, 0.05 / 0.35),
        ("GFAP-BDP/total GFAP, 10% CV, kappa=1", 1.0, 0.10 / 1.30),
    ]:
        print(f"    {lab:<46}{a:>6.2f}{rho_min_with_error(kappa, a):>12.2f}")
    print("    The denominator's imprecision is added but removes nothing, so")
    print("    a low-effect / high-CV marker needs a MUCH better-correlated")
    print("    reference than the textbook rho>0.5 to break even. Where the")
    print("    biological spread is huge (GFAP), the same CV is irrelevant. [OK]")

    print("\n[5] SLE operating points, converted to a common scale")
    print(f"    {'operating point':<46}{'sens':>6}{'spec':>6}{'d prime':>9}{'AUC':>7}")
    for lab, sens, spec in OPERATING_POINTS:
        d = dprime_from_sens_spec(sens, spec)
        print(f"    {lab:<46}{sens:>6.2f}{spec:>6.2f}{d:>9.2f}{auc(d):>7.3f}")
    print("    Note the two comparison groups are NOT interchangeable: the same")
    print("    marker looks ~0.06 AUC better against healthy donors than against")
    print("    other autoimmune disease. Only the latter is the intended use. [OK]")

    print("\n[6] GATE-2 SAMPLE SIZE — and why the obvious calculation is wrong")
    print(f"    {'candidate':<42}{'n(diff)':>9}{'n(AUC>0.70)':>13}")
    for c, r in sorted(results, key=lambda x: -x[1]["auc_real"]):
        n_auc = n_for_auc_lower_bound(r["auc_real"], a_min=0.70)
        lab = "underpowered" if n_auc is None else f"{n_auc}"
        print(f"    {c['name']:<42}{r['n_per_group']:>9}{lab:>13}")
    print("    'n(diff)' powers 'do the groups differ' and returns absurd")
    print("    numbers (n=2 for GFAP) because that is not the claim a diagnostic")
    print("    makes. 'n(AUC>0.70)' powers the real Gate-2 question via")
    print("    Hanley-McNeil; it is 2-5x larger and it ORDERS the portfolio:")
    print("    PRO-C3 needs ~12x the samples of GFAP for the same confidence.")
    print("    C-Alb for its own intended use (mortality WITHIN ESRD) sits at")
    print("    AUC 0.65, below the 0.70 bar, so no sample size reaches it --")
    print("    a stop decision available before drawing a sample. Caveat: at")
    print("    the optimistic end of the [1b] sweep it just touches 0.70, so")
    print("    this one is 'redesign the intended use', not 'kill'. [OK]")
    print("    These are STATISTICAL minima. A real Gate-2 also needs enough")
    print("    subjects to cover the confounders each case names (eGFR, BMI,")
    print("    haemolysis, platelet count), which typically dominates the n.")

    return results


if __name__ == "__main__":
    main()
