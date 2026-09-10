"""Translate the suite's reported accuracy metrics into screening-decision quantities.

A Spearman coefficient does not tell a lab how many plates to run. This converts each
reported (Spearman, RMSE) pair into: enrichment over random picking, hit rate in the
selected batch, and the misclassification risk of the VC-TRIAGE decision gates.
"""
import numpy as np, pandas as pd
from scipy.stats import norm, multivariate_normal


def pearson_from_spearman(rs):
    return 2 * np.sin(np.pi * rs / 6)


def enrichment(rs, top_frac, hit_frac):
    """P(true hit | selected in top `top_frac` by prediction) / hit_frac, Gaussian copula."""
    r = pearson_from_spearman(rs)
    a = norm.ppf(1 - top_frac)      # prediction cutoff
    b = norm.ppf(1 - hit_frac)      # truth cutoff (what counts as a hit)
    mv = multivariate_normal(mean=[0, 0], cov=[[1, r], [r, 1]])
    joint = mv.cdf([-a, -b], lower_limit=[a, b]) if hasattr(mv, "cdf") else None
    # P(T>b, P>a) via inclusion-exclusion on the bivariate normal CDF
    Phi = lambda x, y: multivariate_normal(mean=[0, 0], cov=[[1, r], [r, 1]]).cdf([x, y])
    joint = 1 - norm.cdf(a) - norm.cdf(b) + Phi(a, b)
    precision = joint / top_frac
    return precision / hit_frac, precision


# ---- reported metrics of every deployed model in the suite (from the report) ----
MODELS = [
    # name, endpoint, n_train, reported scaffold-CV Spearman, RMSE(log), page
    ("NA pIC50 (RF+qRASAR)",        "antiviral potency", 1143, 0.809, 1.070, "VC-AVI"),
    ("PA pIC50 (Ridge+qRASAR)",     "antiviral potency",  252, 0.734, 0.728, "VC-AVI"),
    ("CC50 (RF_msl10+qRASAR)",      "cytotoxicity",     33641, 0.704, 0.476, "VC-AVI"),
    ("MIC E. coli",                 "antibacterial",    31676, 0.719, 0.720, "VC-MIC"),
    ("MIC S. aureus",               "antibacterial",    48610, 0.714, 0.730, "VC-MIC"),
    ("MIC P. aeruginosa",           "antibacterial",    20022, 0.698, 0.680, "VC-MIC"),
    ("MIC K. pneumoniae",           "antibacterial",    12015, 0.754, 0.700, "VC-MIC"),
    ("MIC S. pyogenes",             "antibacterial",     6453, 0.824, 0.753, "VC-MIC2"),
    ("MIC M. tuberculosis",         "antibacterial",    19468, 0.643, 0.712, "VC-MIC2"),
    ("MIC E. cloacae",              "antibacterial",     2440, 0.785, 0.746, "VC-MIC2"),
    ("logS solubility",             "ADMET",            12930, 0.599, 1.898, "VC-ADMET"),
    ("hERG pIC50",                  "ADMET",             9716, 0.606, 0.704, "VC-ADMET"),
    ("CYP3A4 pIC50",                "ADMET",             5490, 0.573, 0.652, "VC-ADMET"),
    ("CYP2D6 pIC50",                "ADMET",             2993, 0.579, 0.695, "VC-ADMET"),
    ("HLM half-life",               "ADMET",              601, 0.449, 0.488, "VC-ADMET"),
]


def table():
    rows = []
    for name, ep, n, rs, rmse, page in MODELS:
        ef5, p5 = enrichment(rs, 0.05, 0.05)
        ef1, p1 = enrichment(rs, 0.01, 0.01)
        # how many compounds to test to find 10 true top-5% hits
        n_to_10 = 10 / p5 if p5 > 0 else np.inf
        rows.append(dict(model=name, page=page, endpoint=ep, n=n, spearman=rs, rmse_log=rmse,
                         EF_top5=round(ef5, 1), hit_rate_top5=round(p5 * 100, 1),
                         EF_top1=round(ef1, 1), hit_rate_top1=round(p1 * 100, 1),
                         tests_for_10_hits=int(round(n_to_10)),
                         fold_error_1sd=round(10 ** rmse, 1)))
    return pd.DataFrame(rows)


def gate_risk(true_margin, sigma, gate=0.0):
    """P(decision flips) for a compound whose true value sits `true_margin` logs
    above the gate, given combined model noise sigma."""
    return norm.cdf((gate - true_margin) / sigma)


if __name__ == "__main__":
    pd.set_option("display.width", 200)
    t = table()
    print("=== screening utility implied by each deployed model's reported accuracy ===")
    print(t.to_string(index=False))
    t.to_csv(__file__.rsplit("/", 1)[0] + "/results/utility_table.csv", index=False)

    print("\n=== error propagation through the suite's chained endpoints ===")
    chains = {
        "SI = CC50/EC50 (VC-AVI)": [("CC50", 0.476), ("NA pIC50", 1.070)],
        "antibacterial window = pMIC - pCC50 (VC-TRIAGE)": [("CC50", 0.476), ("MIC", 0.720)],
        "EC50 -> ODE peak reduction (VC-AVI/VC-VK)": [("NA pIC50", 1.070),
                                                      ("biochem IC50 -> cell EC50 offset", 1.0)],
    }
    for k, parts in chains.items():
        s = np.sqrt(sum(v ** 2 for _, v in parts))
        print(f"{k:52s} sigma = {s:.2f} log  ({10**s:.0f}-fold, 1 SD)")

    print("\n=== VC-TRIAGE gate reliability (window gate at 1.0 log) ===")
    sigma = np.sqrt(0.476 ** 2 + 0.720 ** 2)
    for margin in [0.0, 0.5, 1.0, 1.5, 2.0, 3.0]:
        print(f"  compound whose TRUE window exceeds the gate by {margin:.1f} log:"
              f" P(wrongly rejected) = {gate_risk(margin, sigma):.2f}")
