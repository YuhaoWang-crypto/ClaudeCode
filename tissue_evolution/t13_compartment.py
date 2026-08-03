"""
T13 - testing the "interface vs core" hypothesis.

T11 found the same pattern inside every organ: the slowest proteins are the
structural, synaptic, contractile and catalytic core, and the fastest are
secreted or cell-surface. That suggests an organ's apparent evolutionary rate
might not be a property of the organ at all, but of HOW MUCH OF ITS PROTEOME
FACES OUTWARD - liver and spleen deploy plasma proteins and immune receptors,
brain and muscle deploy ion channels and sarcomeres.

Hypothesis (H1): organ rate differences are mostly explained by the
secreted/membrane vs intracellular composition of each organ's gene set, not by
anything tissue-specific.

Prediction that makes it falsifiable: adding subcellular compartment to the T5
controls should COLLAPSE the intrinsic organ effects. If an organ's intrinsic
effect survives exact compartment matching, H1 does not explain that organ.

RESULT: H1 IS FALSIFIED, and the premise it was built on is fine.
  * The premise holds. Secreted 0.170 > Secretory 0.137 > Membrane 0.128 >
    Intracellular 0.108; interface proteins evolve 1.26x faster than core ones
    (+0.337 log2, p=4.6e-60).
  * Organs really do differ enormously in composition - pancreas 0.80 and
    stomach 0.78 interface, against 0.22 for heart and 0.23 for skeletal muscle,
    on a genomic background of 0.355.
  * And it explains none of the organ effect. Compartment absorbs -0.2% of the
    organ term in the nested regression (organ dR2 0.2982 alone, 0.2988 given
    compartment), and every organ with a robust intrinsic effect keeps it after
    exact compartment matching: brain -1.271 -> -1.324, cerebellum -1.198 ->
    -1.256, testis +0.648 -> +0.732, spleen +0.510 -> +0.599.

Why it fails is visible in two rows. TESTIS has a BELOW-average interface
fraction (0.277 vs 0.355) and is among the fastest organs; BRAIN has an
ABOVE-average one (0.556) and is the slowest. The two extremes of the rate
ranking both point the wrong way, so no monotone function of compartment
composition can order the organs. Across all 23 organs the correlation between
interface fraction and median omega is only rho=+0.40 (p=0.056).

The residue worth keeping: compartment is a real gene-level rate factor that is
simply ORTHOGONAL to the organ effect - which means the tissue effect found in
T5 is not a repackaging of protein topology either.

Compartment classes, from UniProt (reviewed human), in this order:
  1. any transmembrane feature              -> Membrane
  2. else "Secreted" in subcellular location -> Secreted
  3. else a signal peptide                   -> Secretory (enters the pathway,
                                                resident in ER/Golgi/lysosome)
  4. otherwise                               -> Intracellular
Order matters: a protein with both a signal peptide and a TM helix is a
membrane protein, not a secreted one.

The binary collapse used for the headline test is
    interface = Membrane + Secreted + Secretory,  core = Intracellular
because the hypothesis is about facing the outside world, not about topology.
"""
import os
import argparse

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree

from .common import (RAW, DERIVED, FIGDIR, apply_figure_style, PALETTE,
                     download, log)
from . import t4_tissue_rates, t5_confounders

UNIPROT_URL = (
    "https://rest.uniprot.org/uniprotkb/stream?"
    "query=organism_id:9606+AND+reviewed:true"
    "&fields=accession,gene_primary,cc_subcellular_location,ft_signal,ft_transmem"
    "&format=tsv"
)
CLASSES = ["Secreted", "Secretory", "Membrane", "Intracellular"]
INTERFACE = {"Secreted", "Secretory", "Membrane"}
MIN_ENRICHED = 40


def fetch_annotation():
    path = os.path.join(RAW, "uniprot_subcellular_human.tsv")
    if not os.path.exists(path) or os.path.getsize(path) < 10000:
        log("downloading UniProt subcellular annotation (reviewed human)")
        download(UNIPROT_URL, path, timeout=1800)
    return pd.read_csv(path, sep="\t")


def classify(df):
    loc = df["Subcellular location [CC]"].fillna("")
    sig = df["Signal peptide"].fillna("")
    tm = df["Transmembrane"].fillna("")
    out = np.full(len(df), "Intracellular", dtype=object)
    out[sig.str.contains("SIGNAL", regex=False).to_numpy()] = "Secretory"
    out[loc.str.contains("Secreted", regex=False).to_numpy()] = "Secreted"
    out[tm.str.contains("TRANSMEM", regex=False).to_numpy()] = "Membrane"
    return pd.Series(out, index=df.index)


def gene_compartments():
    """ENSG -> compartment class, via the BioMart accession map cached by T8."""
    ann = fetch_annotation()
    ann["compartment"] = classify(ann)
    acc2class = dict(zip(ann["Entry"], ann["compartment"]))

    from .t8_structure import uniprot_map
    ensg2acc = uniprot_map("hsapiens_gene_ensembl", "uniprot_human.tsv")
    out = {g: acc2class[a] for g, a in ensg2acc.items() if a in acc2class}
    return pd.Series(out, name="compartment")


def analysis_frame(sp_key=t4_tissue_rates.PRIMARY):
    d = t5_confounders.analysis_frame(sp_key)
    comp = gene_compartments()
    d = d.join(comp, how="inner")
    d["interface"] = d["compartment"].isin(INTERFACE)
    return d


# --------------------------------------------------------------------------
# A. does compartment predict rate at all?
# --------------------------------------------------------------------------
def compartment_rates(d):
    rows = []
    for c in CLASSES:
        s = d[d["compartment"] == c]
        rows.append(dict(compartment=c, n=len(s),
                         median_omega=s["omega"].median(),
                         median_tau=s["tau"].median(),
                         median_log_expr=s["log_expr"].median()))
    res = pd.DataFrame(rows).set_index("compartment")
    a = d.loc[d["interface"], "omega"]
    b = d.loc[~d["interface"], "omega"]
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    eff = dict(n_interface=len(a), n_core=len(b),
               median_interface=float(a.median()), median_core=float(b.median()),
               log2_ratio=float(np.log2(a.median() / b.median())),
               rank_biserial=float(2 * u / (len(a) * len(b)) - 1), p=float(p))
    return res, eff


def organ_composition(d):
    rows = []
    genome = d["interface"].mean()
    for organ in t5_confounders.organ_list(d):
        sub = d[t5_confounders.organ_membership(d, organ)]
        if len(sub) < MIN_ENRICHED:
            continue
        rec = dict(organ=organ, n=len(sub),
                   frac_interface=float(sub["interface"].mean()),
                   median_omega=float(sub["omega"].median()))
        for c in CLASSES:
            rec[f"frac_{c}"] = float((sub["compartment"] == c).mean())
        rows.append(rec)
    res = pd.DataFrame(rows).set_index("organ")
    res["frac_interface_vs_genome"] = res["frac_interface"] - genome
    return res.sort_values("frac_interface", ascending=False), genome


# --------------------------------------------------------------------------
# B. nested OLS: how much of the organ term does compartment absorb?
# --------------------------------------------------------------------------
def nested_ols(d):
    import statsmodels.formula.api as smf
    sub = d[d["enriched_in"].notna() | d["group_enriched_in"].notna()].copy()
    sub["organ_set"] = sub["enriched_in"].fillna(sub["group_enriched_in"])
    base = "log_omega ~ log_expr + tau + log_codons"
    m_base = smf.ols(base, data=sub).fit()
    m_comp = smf.ols(base + " + C(compartment)", data=sub).fit()
    m_org = smf.ols(base + " + C(organ_set)", data=sub).fit()
    m_both = smf.ols(base + " + C(compartment) + C(organ_set)", data=sub).fit()
    return dict(
        n=int(m_base.nobs),
        r2_base=float(m_base.rsquared),
        r2_compartment=float(m_comp.rsquared),
        r2_organ=float(m_org.rsquared),
        r2_both=float(m_both.rsquared),
        dr2_organ_alone=float(m_org.rsquared - m_base.rsquared),
        dr2_organ_given_compartment=float(m_both.rsquared - m_comp.rsquared),
        dr2_compartment_alone=float(m_comp.rsquared - m_base.rsquared),
    )


# --------------------------------------------------------------------------
# C. the decisive test: matched sampling with compartment held fixed
# --------------------------------------------------------------------------
def matched_control_strat(d, organ, strata="compartment", seed=0, caliper=0.25):
    """T5's nearest-neighbour matching, but controls must ALSO come from the
    same compartment class - an exact match on the mediator."""
    member = t5_confounders.organ_membership(d, organ)
    case_all, pool_all = d[member], d[~member]
    if len(case_all) < MIN_ENRICHED:
        return None
    feats = ["tau", "log_expr"]
    mu, sd = d[feats].mean(), d[feats].std(ddof=0)
    rng = np.random.default_rng(seed)
    cases, ctrls, dropped = [], [], 0
    for level in d[strata].dropna().unique():
        case = case_all[case_all[strata] == level]
        pool = pool_all[pool_all[strata] == level]
        if len(case) == 0:
            continue
        if len(pool) == 0:
            dropped += len(case)
            continue
        C = ((case[feats] - mu) / sd).to_numpy()
        P = ((pool[feats] - mu) / sd).to_numpy()
        tree = cKDTree(P)
        k = min(len(P), 200)
        dists, nbrs = tree.query(C, k=k)
        used, ci, pi = set(), [], []
        for i in rng.permutation(len(C)):
            for dist, j in zip(np.atleast_1d(dists[i]), np.atleast_1d(nbrs[i])):
                if j in used or dist > caliper:
                    continue
                used.add(int(j))
                ci.append(i)
                pi.append(int(j))
                break
        dropped += len(case) - len(ci)
        cases.append(case.iloc[ci])
        ctrls.append(pool.iloc[pi])
    if not cases:
        return None
    case = pd.concat(cases)
    ctrl = pd.concat(ctrls)
    if len(case) < MIN_ENRICHED:
        return None
    return case, ctrl, dropped


def mediation_test(d, seed=0):
    """Intrinsic organ effect, before and after holding compartment fixed."""
    genome_med = d["omega"].median()
    rows = []
    for organ in t5_confounders.organ_list(d):
        plain = t5_confounders.matched_control(d, organ, seed=seed)
        strat = matched_control_strat(d, organ, seed=seed)
        if plain is None or strat is None:
            continue
        c0, k0, _ = plain
        c1, k1, drop1 = strat
        i0 = float(np.log2(c0["omega"].median() / k0["omega"].median()))
        i1 = float(np.log2(c1["omega"].median() / k1["omega"].median()))
        u, p = stats.mannwhitneyu(c1["omega"], k1["omega"], alternative="two-sided")
        rows.append(dict(
            organ=organ,
            raw=float(np.log2(c0["omega"].median() / genome_med)),
            intrinsic_tau_expr=i0,
            intrinsic_plus_compartment=i1,
            absorbed_by_compartment=i0 - i1,
            # a ratio against a near-zero baseline is noise, not a percentage:
            # an organ whose tau/expression-matched effect is already ~0 has
            # nothing left for compartment to absorb
            frac_absorbed=float(1 - i1 / i0) if abs(i0) >= 0.15 else np.nan,
            n_matched=len(c1), n_dropped=drop1, p=float(p)))
    res = pd.DataFrame(rows).set_index("organ")
    res["q"] = t5_confounders._bh(res["p"].to_numpy())
    return res.sort_values("intrinsic_tau_expr")


def make_figure(comp_rates, orgcomp, med, genome_frac, d):
    plt = apply_figure_style()
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.6))

    vals = [d.loc[d["compartment"] == c, "omega"].dropna() for c in CLASSES]
    bp = axes[0].boxplot(vals, tick_labels=CLASSES, patch_artist=True,
                         showfliers=False)
    for patch, col in zip(bp["boxes"], [PALETTE[1], PALETTE[4], PALETTE[0], PALETTE[7]]):
        patch.set_facecolor(col)
        patch.set_alpha(0.8)
    for m in bp["medians"]:
        m.set_color("k")
    axes[0].set_ylabel(r"$\omega$")
    axes[0].set_title("compartment predicts rate")
    axes[0].tick_params(axis="x", rotation=25)
    for lab in axes[0].get_xticklabels():
        lab.set_ha("right")

    o = orgcomp.sort_values("frac_interface")
    axes[1].scatter(o["frac_interface"], o["median_omega"], s=34, color=PALETTE[0])
    rho, p = stats.spearmanr(o["frac_interface"], o["median_omega"])
    for organ in o.index:
        axes[1].annotate(organ, (o.loc[organ, "frac_interface"],
                                 o.loc[organ, "median_omega"]),
                         fontsize=6, xytext=(3, 2), textcoords="offset points")
    axes[1].axvline(genome_frac, color="k", ls="--", lw=0.9, label="genome")
    axes[1].set_xlabel("fraction of the organ's genes that are interface")
    axes[1].set_ylabel(r"organ median $\omega$")
    axes[1].set_title(f"across organs: " r"$\rho$" f"={rho:+.2f} (p={p:.1g})")
    axes[1].legend(fontsize=7)

    m = med.sort_values("intrinsic_tau_expr")
    y = np.arange(len(m))
    axes[2].barh(y - 0.2, m["intrinsic_tau_expr"], height=0.4,
                 color=PALETTE[7], label=r"control $\tau$ + expression (T5)")
    axes[2].barh(y + 0.2, m["intrinsic_plus_compartment"], height=0.4,
                 color=PALETTE[1], label="+ compartment held fixed")
    axes[2].set_yticks(y)
    axes[2].set_yticklabels(m.index, fontsize=7)
    axes[2].axvline(0, color="k", lw=0.9)
    axes[2].set_xlabel(r"intrinsic $\log_2$ fold-change vs matched controls")
    axes[2].set_title("does compartment absorb\nthe organ effect?")
    axes[2].grid(axis="y", alpha=0)
    axes[2].legend(fontsize=6.5, loc="lower right")

    p = os.path.join(FIGDIR, "t13_compartment.png")
    fig.savefig(p)
    plt.close(fig)
    return p


def report(sp_key=t4_tissue_rates.PRIMARY):
    print("=" * 72)
    print("T13 - H1: is organ rate just interface-vs-core protein composition?")
    print("=" * 72)
    d = analysis_frame(sp_key)
    print(f"genes with rate + compartment annotation: {len(d)}")
    print(f"  {d['compartment'].value_counts().to_dict()}\n")

    comp_rates, eff = compartment_rates(d)
    print("A. does compartment predict rate?")
    print(f"   {'compartment':<16}{'n':>7}{'median w':>11}{'median tau':>12}")
    for c, r in comp_rates.iterrows():
        print(f"   {c:<16}{int(r['n']):>7}{r['median_omega']:>11.4f}"
              f"{r['median_tau']:>12.3f}")
    print(f"   interface vs core: {eff['median_interface']:.4f} vs "
          f"{eff['median_core']:.4f}  "
          f"({eff['log2_ratio']:+.3f} log2, rb={eff['rank_biserial']:+.2f}, "
          f"p={eff['p']:.2g})\n")

    orgcomp, genome_frac = organ_composition(d)
    print(f"B. organ composition (genome-wide interface fraction = {genome_frac:.3f})")
    print(f"   {'organ':<18}{'n':>5}{'interface':>11}{'vs genome':>11}{'median w':>10}")
    for organ, r in orgcomp.iterrows():
        print(f"   {organ:<18}{int(r['n']):>5}{r['frac_interface']:>11.3f}"
              f"{r['frac_interface_vs_genome']:>+11.3f}{r['median_omega']:>10.4f}")
    rho, p = stats.spearmanr(orgcomp["frac_interface"], orgcomp["median_omega"])
    print(f"   Spearman(interface fraction, organ median omega) = {rho:+.3f} "
          f"(p={p:.2g}, n={len(orgcomp)} organs)\n")

    ols = nested_ols(d)
    print("C. nested OLS on log10(omega), genes assigned to an organ "
          f"(n={ols['n']})")
    print(f"   base (expr + tau + length)          R2={ols['r2_base']:.4f}")
    print(f"   + compartment                       R2={ols['r2_compartment']:.4f}"
          f"   (dR2={ols['dr2_compartment_alone']:.4f})")
    print(f"   + organ identity                    R2={ols['r2_organ']:.4f}"
          f"   (dR2={ols['dr2_organ_alone']:.4f})")
    print(f"   + both                              R2={ols['r2_both']:.4f}")
    print(f"   organ dR2 alone                     {ols['dr2_organ_alone']:.4f}")
    print(f"   organ dR2 GIVEN compartment         {ols['dr2_organ_given_compartment']:.4f}")
    absorbed = 1 - ols["dr2_organ_given_compartment"] / ols["dr2_organ_alone"]
    print(f"   -> compartment absorbs {absorbed:.1%} of the organ term\n")

    med = mediation_test(d)
    print("D. decisive test: matched sampling with compartment held fixed")
    print("   ('absorbed' shown only where the tau+expr effect exceeds 0.15 log2;")
    print("    below that the ratio is dividing by noise)")
    print(f"   {'organ':<18}{'raw':>8}{'tau+expr':>10}{'+compart':>10}"
          f"{'absorbed':>10}{'n':>6}{'q':>10}")
    for organ, r in med.iterrows():
        star = "*" if r["q"] < 0.05 else " "
        fa = r["frac_absorbed"]
        fa_s = f"{fa:>9.0%}" if np.isfinite(fa) else f"{'-':>9}"
        print(f"   {organ:<18}{r['raw']:>+8.3f}{r['intrinsic_tau_expr']:>+10.3f}"
              f"{r['intrinsic_plus_compartment']:>+10.3f}{fa_s}"
              f"{int(r['n_matched']):>6}{r['q']:>10.2e}{star}")

    comp_rates.to_csv(os.path.join(DERIVED, "t13_compartment_rates.tsv"), sep="\t")
    orgcomp.to_csv(os.path.join(DERIVED, "t13_organ_composition.tsv"), sep="\t")
    med.to_csv(os.path.join(DERIVED, "t13_mediation.tsv"), sep="\t")
    p = make_figure(comp_rates, orgcomp, med, genome_frac, d)
    print(f"\nfigure written to: {p}")
    return d, comp_rates, orgcomp, ols, med


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", default=t4_tissue_rates.PRIMARY)
    a = ap.parse_args()
    report(sp_key=a.species)
