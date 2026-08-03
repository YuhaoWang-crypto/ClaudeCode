"""
T5 - is a tissue's evolutionary rate a property OF THE TISSUE, or just the
composition of its gene set?

This is the module the rest of the pipeline exists to support. The naive
finding - "testis genes have high omega, therefore testis evolves fast" - is
almost entirely reproducible from two gene-level covariates that have nothing
to do with testis:

  * expression LEVEL. The E-R anticorrelation (Drummond & Wilke 2008) is the
    single strongest known predictor of protein evolutionary rate: highly
    expressed proteins evolve slowly, plausibly because mistranslation-induced
    misfolding is costlier when there are more copies.
  * tissue SPECIFICITY tau. A gene used in one organ answers to fewer selective
    constraints than a gene used in thirty (pleiotropy).

Tissue-enriched genes are, by construction, high-tau; and organs differ wildly
in how many tissue-specific genes they carry. So "organ" and "gene-level
constraint" are collinear by design, and any analysis that does not break that
collinearity is measuring tau with extra steps.

Three complementary tests here:

  1. NESTED OLS on log10(omega): how much variance does organ identity add
     AFTER expression, tau and coding length are already in the model?
  2. PARTIAL SPEARMAN: rank-based, distribution-free version of the same
     question for the two continuous covariates.
  3. MATCHED SAMPLING (the decisive one): for each organ, draw control genes
     from the rest of the genome matched one-to-one on (tau, log expression),
     and ask whether the organ's genes are STILL faster. The raw effect splits
     into a COMPOSITION part (explained by matching) and an INTRINSIC part
     (survives matching). Only the intrinsic part is evidence about the tissue.
"""
import os

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree

from .common import DERIVED, FIGDIR, apply_figure_style, PALETTE, log
from . import t4_tissue_rates

MIN_ENRICHED = 40
RNG_SEED = 0
# "strict"     = HPA tissue-enriched only (>=4x the runner-up organ)
# "with_group" = also count HPA group-enriched sets containing the organ, which
#                is the only way brain/cerebellum/spinal cord get usable gene
#                sets at all (see t1_expression.group_enriched_sets)
SET_MODE = "with_group"


def organ_membership(d, organ, mode=SET_MODE):
    """Boolean mask: which genes count as this organ's gene set."""
    m = (d["enriched_in"] == organ)
    if mode == "with_group" and "group_enriched_in" in d.columns:
        g = d["group_enriched_in"].fillna("")
        m = m | g.str.split(",").apply(lambda parts: organ in parts)
    return m


def organ_list(d, mode=SET_MODE):
    organs = set(d["enriched_in"].dropna().unique())
    if mode == "with_group" and "group_enriched_in" in d.columns:
        for v in d["group_enriched_in"].dropna():
            organs.update(v.split(","))
    return sorted(organs)


def analysis_frame(sp_key=t4_tissue_rates.PRIMARY):
    df, mat, _ = t4_tissue_rates.build(sp_key)
    d = df[df["omega"].notna() & (df["omega"] > 0)].copy()
    if "group_enriched_in" not in d.columns:
        d["group_enriched_in"] = np.nan
    d["log_omega"] = np.log10(d["omega"])
    d = d[np.isfinite(d["log_omega"]) & d["tau"].notna()]
    return d


# --------------------------------------------------------------------------
# 1. nested OLS
# --------------------------------------------------------------------------
def nested_ols(d):
    import statsmodels.formula.api as smf

    sub = d[d["enriched_in"].notna() | d["group_enriched_in"].notna()].copy()
    sub["organ_set"] = sub["enriched_in"].fillna(sub["group_enriched_in"])
    models = [
        ("intercept only", "log_omega ~ 1", d),
        ("+ log peak expression", "log_omega ~ log_expr", d),
        ("+ tau", "log_omega ~ log_expr + tau", d),
        ("+ log coding length", "log_omega ~ log_expr + tau + log_codons", d),
    ]
    out = []
    for label, formula, data in models:
        m = smf.ols(formula, data=data).fit()
        out.append(dict(model=label, n=int(m.nobs), r2=float(m.rsquared)))

    # organ term is only defined for genes that ARE enriched somewhere, so the
    # last two rows are fit on that subset and must be compared to each other
    base = smf.ols("log_omega ~ log_expr + tau + log_codons", data=sub).fit()
    full = smf.ols("log_omega ~ log_expr + tau + log_codons + C(organ_set)",
                   data=sub).fit()
    ftest = full.compare_f_test(base)
    out.append(dict(model="covariates only (organ-assigned genes)", n=int(base.nobs),
                    r2=float(base.rsquared)))
    out.append(dict(model="+ organ identity", n=int(full.nobs),
                    r2=float(full.rsquared)))
    return pd.DataFrame(out), dict(f=float(ftest[0]), p=float(ftest[1]),
                                   df_diff=float(ftest[2]),
                                   delta_r2=float(full.rsquared - base.rsquared))


# --------------------------------------------------------------------------
# 2. partial Spearman
# --------------------------------------------------------------------------
def _partial_spearman(x, y, covars):
    """Spearman correlation of x,y after linearly removing covars from the ranks."""
    rx = stats.rankdata(x)
    ry = stats.rankdata(y)
    C = np.column_stack([stats.rankdata(c) for c in covars] + [np.ones(len(rx))])
    bx = np.linalg.lstsq(C, rx, rcond=None)[0]
    by = np.linalg.lstsq(C, ry, rcond=None)[0]
    ex, ey = rx - C @ bx, ry - C @ by
    r, p = stats.pearsonr(ex, ey)
    return float(r), float(p)


def partial_correlations(d):
    rows = []
    lo, ex, ta, lc = (d["log_omega"].to_numpy(), d["log_expr"].to_numpy(),
                      d["tau"].to_numpy(), d["log_codons"].to_numpy())
    exm = d["log_expr_mean"].to_numpy()
    r, p = stats.spearmanr(lo, ex)
    rows.append(dict(pair="omega ~ peak expression", control="none", rho=r, p=p))
    r, p = _partial_spearman(lo, ex, [ta, lc])
    rows.append(dict(pair="omega ~ peak expression", control="tau, length", rho=r, p=p))
    r, p = stats.spearmanr(lo, exm)
    rows.append(dict(pair="omega ~ mean expression", control="none", rho=r, p=p))
    r, p = stats.spearmanr(ta, exm)
    rows.append(dict(pair="tau ~ mean expression", control="none", rho=r, p=p))
    r, p = stats.spearmanr(lo, ta)
    rows.append(dict(pair="omega ~ tau", control="none", rho=r, p=p))
    r, p = _partial_spearman(lo, ta, [ex, lc])
    rows.append(dict(pair="omega ~ tau", control="peak expr, length", rho=r, p=p))
    r, p = stats.spearmanr(ta, ex)
    rows.append(dict(pair="tau ~ peak expression", control="none", rho=r, p=p))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# 3. matched sampling
# --------------------------------------------------------------------------
def matched_control(d, organ, seed=RNG_SEED, caliper=0.25):
    """Nearest-neighbour match on standardised (tau, log_expr), no replacement.

    `caliper` is in standard-deviation units; a case with no control inside the
    caliper is DROPPED rather than matched to a distant gene, and the number
    dropped is reported so a poorly-matchable organ cannot silently look clean.
    """
    member = organ_membership(d, organ)
    case = d[member]
    pool = d[~member]
    if len(case) < MIN_ENRICHED or len(pool) < len(case):
        return None
    feats = ["tau", "log_expr"]
    mu = d[feats].mean()
    sd = d[feats].std(ddof=0)
    C = ((case[feats] - mu) / sd).to_numpy()
    P = ((pool[feats] - mu) / sd).to_numpy()

    tree = cKDTree(P)
    used = set()
    idx_case, idx_ctrl = [], []
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(C))
    k = min(len(P), 200)
    dists, nbrs = tree.query(C, k=k)
    for i in order:
        for dist, j in zip(np.atleast_1d(dists[i]), np.atleast_1d(nbrs[i])):
            if j in used or dist > caliper:
                continue
            used.add(int(j))
            idx_case.append(i)
            idx_ctrl.append(int(j))
            break
    if len(idx_case) < MIN_ENRICHED:
        return None
    return case.iloc[idx_case], pool.iloc[idx_ctrl], len(case) - len(idx_case)


def matched_test(d, seed=RNG_SEED):
    genome_med = d["omega"].median()
    rows = []
    for organ in organ_list(d):
        m = matched_control(d, organ, seed=seed)
        if m is None:
            continue
        case, ctrl, n_dropped = m
        u, p = stats.mannwhitneyu(case["omega"], ctrl["omega"],
                                  alternative="two-sided")
        # rank-biserial effect size: +1 = every case gene faster than its control
        rb = 2.0 * u / (len(case) * len(ctrl)) - 1.0
        raw = np.log2(case["omega"].median() / genome_med)
        intrinsic = np.log2(case["omega"].median() / ctrl["omega"].median())
        rows.append(dict(
            organ=organ, n_matched=len(case), n_dropped=n_dropped,
            median_omega_organ=case["omega"].median(),
            median_omega_matched=ctrl["omega"].median(),
            median_tau_organ=case["tau"].median(),
            median_tau_matched=ctrl["tau"].median(),
            log2_raw_vs_genome=raw,
            log2_intrinsic_vs_matched=intrinsic,
            log2_composition=raw - intrinsic,
            rank_biserial=rb, p=p))
    res = pd.DataFrame(rows).set_index("organ")
    if len(res):
        res["q_BH"] = _bh(res["p"].to_numpy())
    return res.sort_values("log2_raw_vs_genome", ascending=False)


def _bh(p):
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank, i in enumerate(order[::-1]):
        k = n - rank
        prev = min(prev, p[i] * n / k)
        q[i] = prev
    return q


def make_figure(matched, sp_key):
    plt = apply_figure_style()
    d = matched.sort_values("log2_raw_vs_genome")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))

    y = np.arange(len(d))
    axes[0].barh(y - 0.2, d["log2_composition"], height=0.4,
                 color=PALETTE[7], label=r"composition ($\tau$ + expression)")
    axes[0].barh(y + 0.2, d["log2_intrinsic_vs_matched"], height=0.4,
                 color=PALETTE[1], label="intrinsic (survives matching)")
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(d.index)
    axes[0].axvline(0, color="k", lw=0.8)
    axes[0].set_xlabel(r"$\log_2$ fold-change in median $\omega$ vs genome")
    axes[0].set_title("decomposition of each organ's apparent rate")
    axes[0].legend(loc="lower right")
    axes[0].grid(axis="y", alpha=0)

    sig = d["q_BH"] < 0.05
    axes[1].scatter(d.loc[~sig, "log2_raw_vs_genome"],
                    d.loc[~sig, "log2_intrinsic_vs_matched"],
                    c="#999999", s=34, label="q >= 0.05")
    axes[1].scatter(d.loc[sig, "log2_raw_vs_genome"],
                    d.loc[sig, "log2_intrinsic_vs_matched"],
                    c=PALETTE[1], s=42, label="q < 0.05")
    lim = [min(d["log2_raw_vs_genome"].min(), d["log2_intrinsic_vs_matched"].min()) - 0.1,
           d["log2_raw_vs_genome"].max() + 0.1]
    axes[1].plot(lim, lim, "k--", lw=0.8, label="y = x (no confounding)")
    axes[1].axhline(0, color="k", lw=0.8)
    for organ in d.index:
        if abs(d.loc[organ, "log2_raw_vs_genome"]) > 0.25 or sig.get(organ, False):
            axes[1].annotate(organ, (d.loc[organ, "log2_raw_vs_genome"],
                                     d.loc[organ, "log2_intrinsic_vs_matched"]),
                             fontsize=6.5, xytext=(3, 3),
                             textcoords="offset points")
    axes[1].set_xlabel(r"raw $\log_2$ FC vs genome")
    axes[1].set_ylabel(r"matched $\log_2$ FC vs $(\tau,$ expression$)$-matched controls")
    axes[1].set_title("how much of the raw effect survives matching")
    axes[1].legend(loc="upper left")

    fig.suptitle(f"T5 - tissue rate: composition vs intrinsic (human-{sp_key})",
                 y=1.01)
    p = os.path.join(FIGDIR, f"t5_confounders_{sp_key}.png")
    fig.savefig(p)
    plt.close(fig)
    return p


def report(sp_key=t4_tissue_rates.PRIMARY):
    print("=" * 72)
    print(f"T5 - composition vs intrinsic tissue effect (human-{sp_key})")
    print("=" * 72)
    d = analysis_frame(sp_key)
    print(f"genes analysed: {len(d)}  ({d['enriched_in'].notna().sum()} "
          f"tissue-enriched, {d['group_enriched_in'].notna().sum()} group-enriched)")
    print(f"organ gene sets built in '{SET_MODE}' mode\n")

    ols, ftest = nested_ols(d)
    print("1. nested OLS on log10(omega)")
    for _, r in ols.iterrows():
        print(f"   {r['model']:<36} n={r['n']:>6}  R2={r['r2']:.4f}")
    print(f"   organ identity adds dR2 = {ftest['delta_r2']:.4f} "
          f"(F={ftest['f']:.1f}, df={ftest['df_diff']:.0f}, p={ftest['p']:.3g})\n")

    pc = partial_correlations(d)
    print("2. partial Spearman correlations")
    for _, r in pc.iterrows():
        print(f"   {r['pair']:<22} control: {r['control']:<22} "
              f"rho={r['rho']:+.3f}  p={r['p']:.3g}")
    print()

    matched = matched_test(d)
    print("3. matched sampling (controls matched on tau + expression)")
    print(f"   {'organ':<17}{'n':>5}{'drop':>6}{'raw':>8}{'intrin':>8}"
          f"{'compos':>8}{'rb':>7}{'q':>10}")
    for organ, r in matched.iterrows():
        star = "*" if r["q_BH"] < 0.05 else " "
        print(f"   {organ:<17}{int(r['n_matched']):>5}{int(r['n_dropped']):>6}"
              f"{r['log2_raw_vs_genome']:>+8.3f}{r['log2_intrinsic_vs_matched']:>+8.3f}"
              f"{r['log2_composition']:>+8.3f}{r['rank_biserial']:>+7.2f}"
              f"{r['q_BH']:>9.2e}{star}")

    ols.to_csv(os.path.join(DERIVED, f"t5_nested_ols_{sp_key}.tsv"), sep="\t", index=False)
    pc.to_csv(os.path.join(DERIVED, f"t5_partial_corr_{sp_key}.tsv"), sep="\t", index=False)
    matched.to_csv(os.path.join(DERIVED, f"t5_matched_{sp_key}.tsv"), sep="\t")
    p = make_figure(matched, sp_key)
    print(f"\nfigure written to: {p}")
    return d, ols, pc, matched


if __name__ == "__main__":
    report()
