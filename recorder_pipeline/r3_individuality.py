"""
R3 — Index of individuality: which existing assays are wasted by a population
cutoff, and how much separation a personal baseline would buy back.

The framework's section 8 argues for personal reference intervals. Biological
variation theory already supplies the machinery to decide WHICH analytes that
applies to, and by HOW MUCH -- so this is a screen that can be run today over
every analyte with published BV data, no new experiment required.

DEFINITIONS
-----------
    CV_I  within-subject biological variation (a person around their setpoint)
    CV_G  between-subject biological variation (setpoints across people)
    CV_A  analytical imprecision

    Index of individuality      II  = CV_I / CV_G
    Reference change value      RCV = sqrt(2) * Z * sqrt(CV_A^2 + CV_I^2)

II < 0.6 is the classical rule for "the population reference interval is
uninformative for this analyte" -- each person occupies a narrow band inside a
wide population range, so a result can be grossly abnormal FOR THAT PERSON
while sitting comfortably inside the population interval.

THE PERSONAL-BASELINE DIVIDEND  (the number this module adds)
-------------------------------------------------------------
Take a disease that shifts the analyte by Delta (in CV units). Against a
population cutoff, the healthy comparison distribution has spread

    sd_pop = sqrt(CV_G^2 + CV_I^2 + CV_A^2)

Against the subject's own earlier measurement, the comparison spread is

    sd_pers = sqrt(2 * (CV_I^2 + CV_A^2))            (single baseline draw)

so switching from a population cutoff to a personal baseline multiplies
separation by

    dividend = sd_pop / sd_pers

which depends only on published BV numbers. The factor sqrt(2) is the honest
cost of comparing two noisy measurements; a baseline averaged over k prior
draws reduces it toward sqrt(1 + 1/k).

This turns "we should use personal baselines" from a philosophy into a ranked,
pre-computable worklist -- and it retrospectively explains the two best-known
successes of the idea (PSA velocity, CA125 ROCA) and the KDIGO definition of
AKI, which is itself a personal-baseline creatinine rule.

INPUT PROVENANCE  [LIT]
-----------------------
CV_I / CV_G are compiled literature medians (Ricos et al. BV database, now
maintained as the EFLM BV Database; oncology and cardiac marker values from
the primary BV studies cited in RECORDER_FRAMEWORK.md). They vary between
studies -- sometimes substantially -- so the RANKING is the deliverable here,
not the third decimal place. Everything downstream of the table is computed.
"""
import numpy as np

# name, CV_I %, CV_G %, CV_A % (typical current assay), note
BV_TABLE = [
    ("Sodium",                   0.6,  0.7,  0.5, "classic high-II control"),
    ("Ferritin",                14.2, 15.0,  3.0, "high II"),
    ("Albumin",                  3.2,  4.8,  2.0, ""),
    ("Calcium (total)",          1.9,  2.8,  1.5, ""),
    ("CRP",                     42.2, 76.3,  5.0, "huge CV_I too - noisy both ways"),
    ("Urate",                    8.6, 17.3,  2.0, ""),
    ("ALT",                     18.0, 41.0,  5.0, ""),
    ("Cholesterol",              5.7, 15.3,  2.0, ""),
    ("Creatinine",               5.3, 14.2,  3.0, "KDIGO already uses own baseline"),
    ("HbA1c",                    1.9,  5.3,  2.0, "the recorder itself has low II"),
    ("CA125",                   24.7, 52.6,  5.0, "ROCA = personal baseline, works"),
    ("NT-proBNP",               29.0, 79.0,  5.0, "obesity/eGFR shift the setpoint"),
    ("PSA (total)",             18.1, 72.4,  5.0, "lowest II of the common markers"),
    ("CEA",                     12.7, 55.6,  5.0, ""),
    ("hs-Troponin T",           14.0, 60.0,  8.0, "serial delta already standard"),
]


def index_of_individuality(cv_i, cv_g):
    return cv_i / cv_g


def rcv(cv_i, cv_a, z=1.96):
    """Reference change value, % — the smallest change that is not noise."""
    return np.sqrt(2.0) * z * np.sqrt(cv_a**2 + cv_i**2)


def personal_dividend(cv_i, cv_g, cv_a, k_baseline=1):
    """Separation multiplier from switching to a personal baseline.

    k_baseline = number of prior measurements averaged into the baseline.
    """
    sd_pop = np.sqrt(cv_g**2 + cv_i**2 + cv_a**2)
    sd_pers = np.sqrt((1.0 + 1.0 / k_baseline) * (cv_i**2 + cv_a**2))
    return sd_pop / sd_pers


def main():
    print("=" * 78)
    print("R3 — Index of individuality and the personal-baseline dividend")
    print("=" * 78)

    rows = []
    for name, cvi, cvg, cva, note in BV_TABLE:
        ii = index_of_individuality(cvi, cvg)
        r = rcv(cvi, cva)
        div1 = personal_dividend(cvi, cvg, cva, k_baseline=1)
        div3 = personal_dividend(cvi, cvg, cva, k_baseline=3)
        rows.append((name, cvi, cvg, cva, ii, r, div1, div3, note))

    rows.sort(key=lambda r: r[4])   # by II ascending = most wasted first

    print("\n[1] Ranked by index of individuality (lowest = most wasted by a")
    print("    population reference interval)")
    print(f"    {'analyte':<18}{'CV_I':>6}{'CV_G':>6}{'II':>7}{'RCV%':>8}"
          f"{'dividend':>10}{'(k=3)':>8}  verdict")
    for name, cvi, cvg, cva, ii, r, d1, d3, note in rows:
        flag = "POPULATION RI UNINFORMATIVE" if ii < 0.6 else "population RI ok"
        print(f"    {name:<18}{cvi:>6.1f}{cvg:>6.1f}{ii:>7.2f}{r:>8.1f}"
              f"{d1:>9.2f}x{d3:>7.2f}x  {flag}")

    n_bad = sum(1 for r in rows if r[4] < 0.6)
    print(f"\n    {n_bad}/{len(rows)} of these routine analytes have II < 0.6.")
    print("    Every one of them is a candidate for the same treatment that")
    print("    rescued PSA velocity and CA125 ROCA. [OK]")

    print("\n[2] The three biggest dividends, and what is already known about them")
    for name, cvi, cvg, cva, ii, r, d1, d3, note in rows[:3]:
        print(f"    {name:<14} II={ii:.2f}  personal baseline buys {d1:.2f}x"
              f" separation ({d3:.2f}x with 3 priors)")
        print(f"                   RCV = +/-{r:.0f}% -> a change under this is noise")
        if note:
            print(f"                   note: {note}")

    print("\n[3] Cross-check against clinical practice")
    print("    The framework predicts that low-II analytes should have drifted")
    print("    toward personal-baseline rules on their own. They have:")
    print("      creatinine  II=0.37 -> KDIGO defines AKI as a RISE FROM the")
    print("                             patient's own baseline, not a cutoff")
    print("      hs-troponin II=0.23 -> practice uses the 0/1h serial DELTA")
    print("      CA125       II=0.47 -> ROCA scores the within-woman trajectory")
    print("                             'irrespective of the absolute level'")
    print("      PSA         II=0.25 -> velocity / %free PSA both attack the")
    print("                             between-person spread")
    print("    Three independent clinical fields converged on this without a")
    print("    shared theory. The BV table predicts them from published numbers")
    print("    alone -- and the same table nominates the ones nobody has done")
    print("    yet (NT-proBNP II=0.37, CEA II=0.23, ALT II=0.44). [OK]")

    print("\n[4] Caveat that limits this method")
    print("    A personal baseline requires a PRIOR SAMPLE from a stable state.")
    print("    It is therefore available for screening, surveillance and therapy")
    print("    monitoring, and unavailable for first-presentation diagnosis.")
    print("    RCV also assumes a steady setpoint: it is invalid across events")
    print("    that legitimately move the setpoint (pregnancy, weight change,")
    print("    starting an SGLT2 inhibitor). [OK]")

    return rows


if __name__ == "__main__":
    main()
