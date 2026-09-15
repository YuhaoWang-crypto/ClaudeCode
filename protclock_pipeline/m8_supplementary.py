"""
Module 8 - REAL data. The paper's supplementary tables.

Everything in M1 through M7 runs on simulation. This module does not. The
supplementary workbook is openly downloadable even though the article is
paywalled and the proteome is controlled, and it carries enough to reproduce
parts of the paper directly:

  S2   real per-sample predictions from all 43 clock variants, 168 rows
       (42 patients x 4 visits), keyed by a SHA256 of the NPX row
  S3   Wilcoxon paired statistics, baseline vs each visit, per arm and clock
  S4   Mann-Whitney U comparing treated arms against placebo
  S5   the headline tally of how many tests passed
  S8   enrichment of clock features in drug-modulated proteins
  S14  aggregated GSEA scores and leading-edge genes

WHAT S2 DOES AND DOES NOT ALLOW
-------------------------------
The sample identifier is a SHA256 hash of the sample's own NPX row, so the
predictions cannot be linked back to an arm or a visit without the
controlled OMIX008341 data. Per-arm statistics therefore cannot be
recomputed from S2 alone; that linkage is exactly what access control
protects.

What S2 does give, and what no amount of simulation could, is the JOINT
DISTRIBUTION of the six clocks over the real trial samples. That settles a
question M5 could only flag: the paper's central argument is that six
independently developed clocks agree, and a binomial p-value over six
agreeing clocks assumes they are independent. They are not, and S2 says by
how much.

The effective number of independent clocks is estimated from the eigenvalues
of their correlation matrix (Nyholt 2004, as corrected by Li and Ji 2005),
and the concordance p-value is then recomputed against that number instead
of six.
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

SUPP_URL = ("https://static-content.springer.com/esm/"
            "art%3A10.1038%2Fs41587-026-03286-y/MediaObjects/"
            "41587_2026_3286_MOESM3_ESM.xlsx")
DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
SUPP_PATH = os.path.join(DATA_DIR, "supplementary_tables.xlsx")

# The six headline clocks, as S2 spells them.
HEADLINE = {
    "ProtAge": "Argentieri_2024",
    "ipfP3GPT": "Galkin_2025(v6_standard)",
    "PAC": "Kuo_2024(none)",
    "OrganAge_chrono": "Goeminne_2025(chrono_conventional_none)",
    "OrganAge_mortality": "Goeminne_2025(mortality_conventional_none)",
    "PAOPAC": "Han_2025(Conventional)",
}


def fetch(path=SUPP_PATH, url=SUPP_URL, force=False):
    """Download the supplementary workbook if it is not already cached."""
    if os.path.exists(path) and not force:
        return path
    os.makedirs(os.path.dirname(path), exist_ok=True)
    import urllib.request
    urllib.request.urlretrieve(url, path)
    return path


def load_predictions(path=SUPP_PATH):
    """S2: real clock predictions, 168 samples x 43 clock variants."""
    d = pd.read_excel(path, sheet_name="S2", header=3)
    d = d.dropna(how="all")
    return d


def headline_matrix(pred):
    """The six headline clocks as a numeric frame, dropping incomplete rows."""
    missing = [k for k, v in HEADLINE.items() if v not in pred.columns]
    if missing:
        raise KeyError(f"S2 is missing expected clock columns: {missing}")
    m = pred[[HEADLINE[k] for k in HEADLINE]].astype(float)
    m.columns = list(HEADLINE)
    return m.dropna()


def effective_n_independent(corr):
    """
    Effective number of independent variables from a correlation matrix.

    Two standard estimators, both eigenvalue-based:

      Nyholt 2004   Meff = 1 + (M-1) * (1 - Var(lambda)/M)
      Li & Ji 2005  Meff = sum over eigenvalues of f(|lambda|), where
                    f(x) = I(x >= 1) + (x - floor(x))

    Li and Ji is the less anticonservative of the two and is preferred here.
    Both are reported because they disagree when correlations are strong, and
    the disagreement is itself informative.
    """
    lam = np.linalg.eigvalsh(np.asarray(corr, dtype=float))
    lam = np.clip(lam, 0, None)
    m = len(lam)

    nyholt = 1.0 + (m - 1) * (1.0 - np.var(lam, ddof=1) / m)
    li_ji = float(np.sum((np.abs(lam) >= 1).astype(float)
                         + (np.abs(lam) - np.floor(np.abs(lam)))))
    return {"n_clocks": m, "nyholt": float(nyholt),
            "li_ji": float(min(li_ji, m)),
            "eigenvalues": lam[::-1]}


def concordance_pvalue(n_agree, n_total, n_effective=None):
    """
    Two-sided binomial p for k of n clocks agreeing in direction.

    With `n_effective` the count is rescaled to the effective number of
    independent clocks before testing, which is the honest version when the
    clocks are correlated. Rescaling keeps the observed agreement FRACTION
    and reduces the sample size, so the p-value rises.
    """
    if n_effective is None:
        k, n = n_agree, n_total
    else:
        n = max(int(round(n_effective)), 1)
        k = int(round(n_agree / n_total * n))
    return float(stats.binomtest(k, n, 0.5, alternative="two-sided").pvalue)


def load_test_tally(path=SUPP_PATH):
    """S5: how many of the paper's tests actually passed."""
    raw = pd.read_excel(path, sheet_name="S5", header=1)
    raw = raw.dropna(how="all")
    raw.columns = [str(c) for c in raw.columns]
    return raw


def load_enrichment(path=SUPP_PATH):
    """S8: enrichment of clock features by response trajectory class."""
    raw = pd.read_excel(path, sheet_name="S8", header=1)
    raw = raw.dropna(how="all")
    raw["Arm"] = raw["Arm"].ffill()
    return raw


def report(path=None):
    path = fetch(path or SUPP_PATH)
    print("M8  THE PAPER'S SUPPLEMENTARY TABLES  (REAL DATA, not simulated)")
    print("-" * 68)
    print(f"  source: {SUPP_URL[:58]}...")

    pred = load_predictions(path)
    m = headline_matrix(pred)
    print(f"  S2 predictions          : {pred.shape[0]} samples x "
          f"{pred.shape[1] - 1} clock variants")
    print(f"  six headline clocks     : {m.shape[0]} complete samples")
    print(f"  sample key              : SHA256 of the NPX row, so arm and "
          f"visit\n                            cannot be recovered without "
          f"OMIX008341")

    corr = m.corr(method="pearson")
    print("\n  REAL cross-clock correlation on the trial samples:")
    print(corr.round(3).to_string())

    off = corr.to_numpy()[np.triu_indices(len(corr), k=1)]
    eff = effective_n_independent(corr)
    print(f"\n  mean pairwise r         : {off.mean():+.3f}")
    print(f"  range                   : {off.min():+.3f} to {off.max():+.3f}")
    print(f"  effective independent clocks:")
    print(f"    Nyholt 2004           : {eff['nyholt']:.2f} of 6")
    print(f"    Li & Ji 2005          : {eff['li_ji']:.2f} of 6")
    print(f"  eigenvalues             : "
          f"{np.round(eff['eigenvalues'], 3).tolist()}")

    p_naive = concordance_pvalue(6, 6)
    p_adj = concordance_pvalue(6, 6, eff["li_ji"])
    print(f"\n  what this does to the concordance argument:")
    print(f"    p for 6/6 agreeing, assuming independence : {p_naive:.4f}")
    print(f"    p using {eff['li_ji']:.2f} effective clocks           "
          f": {p_adj:.4f}")
    print("    The clocks are not six independent votes. M5 labelled its")
    print("    binomial column anticonservative; this is the measurement.")

    print("\n  S5, the paper's own tally of passed tests:")
    tally = load_test_tally(path)
    print(tally.head(4).to_string(index=False))

    print("\n  S8, enrichment of clock features by trajectory class:")
    enr = load_enrichment(path)
    print(enr.to_string(index=False))
    print("\n  Note the direction: enrichment is strong for SUSTAINED")
    print("  responders and BELOW 1 for delayed and transient ones, so a")
    print("  single headline odds ratio does not describe this table.")

    return {"predictions": pred, "headline": m, "corr": corr,
            "effective": eff, "tally": tally, "enrichment": enr}


if __name__ == "__main__":
    report()
