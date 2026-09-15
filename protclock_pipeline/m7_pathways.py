"""
Module 7 - Pathway-level analysis.

The paper reports that the proteins rentosertib moves implicate senescence
and metabolic processes alongside the expected anti-fibrotic programme. This
module runs two complementary tests over the pathway sets M1 defined:

  OVER-REPRESENTATION  Fisher exact test of (in set) x (drug-altered).
                       Depends entirely on the FDR cut used to call
                       "altered", so it inherits M6's power problem.

  MEAN SHIFT           Welch t-test of the per-protein effect size inside
                       the set against outside it. Uses every protein and
                       never thresholds, so it keeps the power that
                       over-representation throws away. This is the
                       competitive form of gene-set testing and it is the
                       more trustworthy of the two here.

Both are reported because they can disagree, and when they do it is almost
always over-representation being underpowered rather than the mean shift
being wrong.

A SELECTION BIAS worth stating plainly: proteins with larger planted effects
are easier to detect, and in this simulation the aging-associated responders
carry the larger effects. Any enrichment computed on a THRESHOLDED set is
therefore biased upward relative to the truth. M6 shows this directly - the
paired contrast recovers an odds ratio near 2.2 where the planted value is
about 1.73. The same bias applies to the paper's 1.74 to an unknown degree.
"""
import numpy as np
import pandas as pd
from scipy import stats

from protclock_pipeline.m6_enrichment import benjamini_hochberg


def over_representation(altered, set_idx, n_total):
    """Fisher exact test of set membership against the altered set."""
    in_set = np.zeros(n_total, dtype=bool)
    in_set[set_idx] = True

    a = int(np.sum(in_set & altered))
    b = int(np.sum(in_set & ~altered))
    c = int(np.sum(~in_set & altered))
    d = int(np.sum(~in_set & ~altered))
    if min(a, b, c, d) == 0:
        return {"odds_ratio": float("nan"), "p": float("nan"),
                "n_set": len(set_idx), "n_hit": a}
    odds, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
    return {"odds_ratio": float(odds), "p": float(p),
            "n_set": len(set_idx), "n_hit": a}


def aging_aligned_shift(effect, age_slope, set_idx, n_total):
    """
    Does the drug move a set's proteins AGAINST their own aging direction?

    A plain mean shift is the wrong test here. Within any aging-related set
    some proteins rise with age and others fall, so a drug that rejuvenates
    every one of them produces shifts of both signs whose mean is close to
    zero. M7's first run showed exactly that: senescence came out at Cohen's
    d = -0.07, p = 0.91, despite being built from aging proteins.

    Projecting each effect onto the protein's own aging direction fixes it:

        aligned = effect_npx * sign(age_slope)

    Negative aligned values mean the protein moved the way it moves in a
    YOUNGER person. The age slope is taken from the reference cohort, so this
    stays a discovered quantity and never touches the planted truth.
    """
    in_set = np.zeros(n_total, dtype=bool)
    in_set[set_idx] = True
    aligned = effect * np.sign(age_slope)
    a, b = aligned[in_set], aligned[~in_set]
    t, p = stats.ttest_ind(a, b, equal_var=False)
    pooled = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    return {"aligned_in": float(a.mean()), "aligned_out": float(b.mean()),
            "diff": float(a.mean() - b.mean()),
            "cohens_d": float((a.mean() - b.mean()) / pooled) if pooled else
            float("nan"),
            "t": float(t), "p": float(p)}


def mean_shift(effect, set_idx, n_total):
    """
    Welch test of the effect size inside a set against outside it.

    Returns the difference in mean NPX effect and a standardised effect size,
    so a set can be large-and-weak or small-and-strong and still be read off.
    """
    in_set = np.zeros(n_total, dtype=bool)
    in_set[set_idx] = True
    a, b = effect[in_set], effect[~in_set]
    t, p = stats.ttest_ind(a, b, equal_var=False)
    pooled = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    return {"mean_in": float(a.mean()), "mean_out": float(b.mean()),
            "diff": float(a.mean() - b.mean()),
            "cohens_d": float((a.mean() - b.mean()) / pooled) if pooled else
            float("nan"),
            "t": float(t), "p": float(p), "n_set": int(in_set.sum())}


def run(pre, differential, aging):
    """All three tests for every pathway set, BH corrected across sets."""
    sets = pre.gene_sets_reindexed()
    altered = differential["altered"].to_numpy()
    effect = differential["effect_npx"].to_numpy()
    age_slope = aging["age_slope"].to_numpy()
    n_total = len(differential)

    rows = []
    for name, idx in sets.items():
        if len(idx) < 5:
            continue
        ora = over_representation(altered, idx, n_total)
        ms = mean_shift(effect, idx, n_total)
        al = aging_aligned_shift(effect, age_slope, idx, n_total)
        rows.append({"pathway": name, "n_set": ora["n_set"],
                     "n_altered_in_set": ora["n_hit"],
                     "ora_odds_ratio": ora["odds_ratio"], "ora_p": ora["p"],
                     "shift_npx": ms["diff"], "shift_p": ms["p"],
                     "aligned_npx": al["diff"], "aligned_d": al["cohens_d"],
                     "aligned_p": al["p"]})

    tab = pd.DataFrame(rows)
    if len(tab):
        tab["ora_q"] = benjamini_hochberg(tab["ora_p"].fillna(1.0))
        tab["shift_q"] = benjamini_hochberg(tab["shift_p"])
        tab["aligned_q"] = benjamini_hochberg(tab["aligned_p"])
    return tab.sort_values("aligned_p").reset_index(drop=True)


def report(pre=None, differential=None, aging=None):
    from protclock_pipeline import m2_preprocess, m6_enrichment
    if pre is None:
        pre = m2_preprocess.run(n_reference=2500)
    if differential is None:
        differential = m6_enrichment.differential_abundance(
            pre, contrast="paired_treated")
    if aging is None:
        aging = m6_enrichment.aging_associated(pre)

    tab = run(pre, differential, aging)

    print("M7  PATHWAY ANALYSIS")
    print("-" * 68)
    print("  aligned = effect projected onto each protein's aging direction;")
    print("  NEGATIVE means the set moved the way it moves in a younger"
          " person.\n")
    cols = ["pathway", "n_set", "n_altered_in_set", "ora_odds_ratio", "ora_q",
            "shift_npx", "shift_q", "aligned_npx", "aligned_d", "aligned_q"]
    print(tab[cols].round(4).to_string(index=False))

    sig_al = tab[(tab["aligned_q"] < 0.05) & (tab["aligned_npx"] < 0)]
    sig_shift = tab[tab["shift_q"] < 0.05]["pathway"].tolist()
    sig_ora = tab[tab["ora_q"] < 0.05]["pathway"].tolist()
    print(f"\n  rejuvenated by aligned shift     : "
          f"{', '.join(sig_al['pathway']) if len(sig_al) else 'none'}")
    print(f"  significant by plain mean shift  : "
          f"{', '.join(sig_shift) if sig_shift else 'none'}"
          f"   <- cancels, see aging_aligned_shift docstring")
    print(f"  significant by over-representation: "
          f"{', '.join(sig_ora) if sig_ora else 'none'}")
    print("\n  PLANTED truth: senescence_SASP and metabolic_process were")
    print("  drawn from aging proteins, fibrosis_ECM overlaps the drug")
    print("  target, immune_activation is a negative control and should")
    print("  NOT come out significant.")
    return tab


if __name__ == "__main__":
    report()
