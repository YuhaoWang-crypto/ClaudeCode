"""
T15 - pathways: how a shared pathway's components get allocated across organs,
and whether they evolve as a module.

Two questions that only make sense at pathway level, and that the gene-level
modules (T5, T13, T14) cannot ask.

1. ALLOCATION. A pathway is a fixed set of jobs, but its components are not
   equally distributed across tissues: some members are used everywhere, some
   are the tissue-restricted version. Comparing those two groups WITHIN THE SAME
   PATHWAY is the cleanest control available for "does tissue restriction itself
   raise the rate" - pathway identity, and therefore function, biochemistry and
   complex membership, is held fixed by construction. Paired across pathways
   with a Wilcoxon signed-rank test, so a handful of large pathways cannot carry
   the result.

2. EMERGENCE / MODULARITY. Do a pathway's members speed up and slow down
   TOGETHER across mammalian lineages? That is evolutionary rate covariation
   (ERC; Clark, Wolfe & Singh 2012), computed here from the T12 RER matrix
   (4195 genes x 240 branches x 121 mammals). If pathways are coevolving
   modules, within-pathway ERC exceeds the background. And the allocation
   question can then be asked dynamically: when a pathway's members sit in
   DIFFERENT organs, do they still coevolve, or has the pathway fragmented into
   organ-specific modules?

The null that matters for ERC: two genes with similar overall evolutionary rates
correlate across branches for reasons that have nothing to do with shared
function, and pathway members are enriched for similar rates. So the background
here is not random pairs, it is RATE-MATCHED random pairs, drawn from the same
joint distribution of gene rates as the within-pathway pairs. Reporting a raw
within-pathway ERC against an unmatched background would mostly measure that
confound - the same class of error the T12 null exposed.
"""
import os
import argparse

import numpy as np
import pandas as pd
from scipy import stats

from .common import DERIVED, FIGDIR, apply_figure_style, PALETTE, log
from . import t4_tissue_rates, t5_confounders, t12_rer, t14_mediators

MIN_PATHWAY = 8          # members present in the rate table
MAX_PATHWAY = 300        # above this a "pathway" is a Reactome super-category
MIN_PER_GROUP = 3        # tissue-enriched / broad members needed for a pair
N_NULL_PAIRS = 200000


# --------------------------------------------------------------------------
# pathway membership
# --------------------------------------------------------------------------
def pathway_members():
    """pathway id -> (name, set of ENSG)."""
    cache = os.path.join(DERIVED, "t15_pathway_members.tsv.gz")
    if os.path.exists(cache):
        df = pd.read_csv(cache, sep="\t")
        return {r.pathway: (r.name, set(r.genes.split(",")))
                for r in df.itertuples()}
    if not os.path.exists(t14_mediators.REACTOME):
        from .common import download
        download(t14_mediators.REACTOME_URL, t14_mediators.REACTOME, timeout=1800)
    e2g = t14_mediators.entrez_to_ensg()
    members, names = {}, {}
    with open(t14_mediators.REACTOME) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6 or f[5] != "Homo sapiens":
                continue
            pid, pname = f[1], f[3].strip()
            names[pid] = pname
            for g in e2g.get(f[0], ()):
                members.setdefault(pid, set()).add(g)
    out = {p: (names[p], g) for p, g in members.items()}
    pd.DataFrame([dict(pathway=p, name=n, genes=",".join(sorted(g)))
                  for p, (n, g) in out.items()]).to_csv(
        cache, sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
    log(f"Reactome: {len(out)} human pathways")
    return out


def usable_pathways(d, members):
    out = {}
    for pid, (name, genes) in members.items():
        present = [g for g in genes if g in d.index]
        if MIN_PATHWAY <= len(present) <= MAX_PATHWAY:
            out[pid] = (name, present)
    return out


# --------------------------------------------------------------------------
# 1. allocation: tissue-restricted vs broad members of the SAME pathway
# --------------------------------------------------------------------------
HISTONE = ("H1-", "H2A", "H2B", "H3-", "H3C", "H4C", "H2AC", "H2BC", "H1F")


def _is_histone(sym):
    return isinstance(sym, str) and sym.startswith(HISTONE)


def allocation_test(d, paths, min_per_group=MIN_PER_GROUP, drop_histones=False):
    """Within each pathway, median omega of tissue-enriched vs broad members.

    NOT independent of tau: "tissue-enriched" is a high-tau class by definition,
    so this cannot show an effect BEYOND tau. What it does show is that the
    tau-rate relation holds INSIDE pathways - it is not an artefact of
    tissue-specific genes sitting in different, faster pathways than broadly
    expressed ones. The tau gap between the two groups is reported alongside so
    that is visible rather than implied.
    """
    is_enriched = (d["enriched_in"].notna() | d["group_enriched_in"].notna())
    rows = []
    for pid, (name, genes) in paths.items():
        sub = d.loc[genes]
        if drop_histones:
            sub = sub[~sub["symbol"].map(_is_histone)]
        a = sub.loc[is_enriched.reindex(sub.index).fillna(False), "omega"].dropna()
        b = sub.loc[~is_enriched.reindex(sub.index).fillna(False), "omega"].dropna()
        if len(a) < min_per_group or len(b) < min_per_group:
            continue
        rows.append(dict(pathway=pid, name=name, n_enriched=len(a), n_broad=len(b),
                         omega_enriched=float(a.median()),
                         omega_broad=float(b.median()),
                         log2_diff=float(np.log2(a.median() / b.median())),
                         tau_enriched=float(sub.loc[a.index, "tau"].median()),
                         tau_broad=float(sub.loc[b.index, "tau"].median())))
    res = pd.DataFrame(rows).set_index("pathway")
    if len(res):
        w = stats.wilcoxon(res["log2_diff"], alternative="two-sided")
        stat = dict(n_pathways=len(res),
                    median_log2_diff=float(res["log2_diff"].median()),
                    frac_positive=float((res["log2_diff"] > 0).mean()),
                    median_tau_gap=float((res["tau_enriched"]
                                          - res["tau_broad"]).median()),
                    p=float(w.pvalue))
    else:
        stat = {}
    return res.sort_values("log2_diff", ascending=False), stat


def pathway_fixed_effects(d, paths):
    """Is the tau-rate slope a WITHIN-pathway relation or a between-pathway one?

    Compares the tau coefficient with and without pathway fixed effects. If
    adding C(pathway) leaves the slope intact, tissue restriction raises the rate
    inside pathways, so it is not composition across pathways. Genes in several
    pathways are assigned to their smallest one, so each gene appears once and
    the fixed effect is identified.
    """
    import statsmodels.formula.api as smf
    smallest = {}
    for pid, (name, genes) in paths.items():
        for g in genes:
            if g not in smallest or len(paths[smallest[g]][1]) > len(genes):
                smallest[g] = pid
    sub = d.loc[[g for g in smallest if g in d.index]].copy()
    sub["pathway"] = [smallest[g] for g in sub.index]
    sub = sub.dropna(subset=["log_omega", "tau", "log_expr", "log_codons"])
    base = smf.ols("log_omega ~ tau + log_expr + log_codons", data=sub).fit()
    fe = smf.ols("log_omega ~ tau + log_expr + log_codons + C(pathway)",
                 data=sub).fit()
    return dict(n=int(base.nobs), n_pathways=sub["pathway"].nunique(),
                tau_coef=float(base.params["tau"]),
                tau_coef_fe=float(fe.params["tau"]),
                tau_p_fe=float(fe.pvalues["tau"]),
                retained=float(fe.params["tau"] / base.params["tau"]),
                r2=float(base.rsquared), r2_fe=float(fe.rsquared))


def dispersion_test(d, paths):
    """Does a pathway spread across more organs go with a wider rate spread?"""
    rows = []
    for pid, (name, genes) in paths.items():
        sub = d.loc[genes].dropna(subset=["omega"])
        if len(sub) < MIN_PATHWAY:
            continue
        organs = set()
        for v in sub["enriched_in"].dropna():
            organs.add(v)
        for v in sub["group_enriched_in"].dropna():
            organs.update(v.split(","))
        lo = np.log2(sub["omega"].clip(lower=1e-4))
        rows.append(dict(pathway=pid, name=name, n=len(sub),
                         n_organs=len(organs),
                         frac_enriched=float((sub["enriched_in"].notna()
                                              | sub["group_enriched_in"].notna()).mean()),
                         omega_median=float(sub["omega"].median()),
                         omega_iqr_log2=float(lo.quantile(0.75) - lo.quantile(0.25))))
    res = pd.DataFrame(rows).set_index("pathway")
    rho, p = stats.spearmanr(res["n_organs"], res["omega_iqr_log2"])
    return res, dict(rho=float(rho), p=float(p), n=len(res))


# --------------------------------------------------------------------------
# 2. ERC from the T12 RER matrix
# --------------------------------------------------------------------------
def load_rer():
    B, genes, order, splits, ncod = t12_rer.build_matrix()
    RER = t12_rer.compute_rer(B)
    sym = [g.split(".", 2)[1] if g.count(".") >= 2 else g for g in genes]
    Z = (RER - RER.mean(axis=1, keepdims=True))
    Z /= (Z.std(axis=1, keepdims=True) + 1e-12)
    rate = B.mean(axis=1)
    return Z, np.array(sym), rate, RER.shape[1]


def _rate_bins(rate, n_bins=20):
    edges = np.quantile(rate, np.linspace(0, 1, n_bins + 1))
    return np.clip(np.digitize(rate, edges[1:-1]), 0, n_bins - 1)


def erc_analysis(d, paths, seed=0):
    """Within-pathway ERC against a RATE-MATCHED random background."""
    Z, sym, rate, n_branch = load_rer()
    sym2row = {}
    for i, s in enumerate(sym):
        sym2row.setdefault(s, i)
    ensg2sym = d["symbol"].dropna().to_dict()
    bins = _rate_bins(rate)

    pairs, pair_bins = [], []
    for pid, (name, genes) in paths.items():
        rows = sorted({sym2row[ensg2sym[g]] for g in genes
                       if g in ensg2sym and ensg2sym[g] in sym2row})
        if len(rows) < MIN_PER_GROUP:
            continue
        for i in range(len(rows)):
            for j in range(i + 1, len(rows)):
                pairs.append((rows[i], rows[j]))
                pair_bins.append((bins[rows[i]], bins[rows[j]]))
    if len(pairs) < 100:
        return None
    pairs = np.array(pairs)
    obs = (Z[pairs[:, 0]] * Z[pairs[:, 1]]).sum(axis=1) / n_branch

    # rate-matched null: sample random pairs with the SAME joint rate-bin
    # composition as the observed within-pathway pairs
    rng = np.random.default_rng(seed)
    by_bin = {b: np.flatnonzero(bins == b) for b in range(bins.max() + 1)}
    pb = np.array(pair_bins)
    idx = rng.choice(len(pb), size=min(N_NULL_PAIRS, len(pb) * 20))
    null_a, null_b = [], []
    for k in idx:
        ba, bb = pb[k]
        if len(by_bin[ba]) == 0 or len(by_bin[bb]) == 0:
            continue
        ia = by_bin[ba][rng.integers(len(by_bin[ba]))]
        ib = by_bin[bb][rng.integers(len(by_bin[bb]))]
        if ia != ib:
            null_a.append(ia)
            null_b.append(ib)
    null_a, null_b = np.array(null_a), np.array(null_b)
    nul = (Z[null_a] * Z[null_b]).sum(axis=1) / n_branch

    u, p = stats.mannwhitneyu(obs, nul, alternative="two-sided")
    return dict(n_pairs=len(obs), n_null=len(nul),
                median_within=float(np.median(obs)),
                median_null=float(np.median(nul)),
                delta=float(np.median(obs) - np.median(nul)),
                rank_biserial=float(2 * u / (len(obs) * len(nul)) - 1),
                p=float(p)), obs, nul


def erc_permutation_null(d, paths, n_perm=200, seed=0):
    """Shuffle which genes belong to which pathway, keeping pathway sizes.

    The rate-matched null in `erc_analysis` fixes the rate confound but not the
    dependence between pairs: 4195 genes generate 240k pairs, each gene appears
    in hundreds of them, so the Mann-Whitney p-value is meaningless (it assumes
    independent observations and will return ~0 for any effect at all). This
    null preserves the pair structure exactly and asks whether the observed
    median exceeds what random gene sets of the same sizes produce.
    """
    Z, sym, rate, n_branch = load_rer()
    sym2row = {}
    for i, s_ in enumerate(sym):
        sym2row.setdefault(s_, i)
    ensg2sym = d["symbol"].dropna().to_dict()

    sets = []
    for pid, (name, genes) in paths.items():
        rows = sorted({sym2row[ensg2sym[g]] for g in genes
                       if g in ensg2sym and ensg2sym[g] in sym2row})
        if len(rows) >= MIN_PER_GROUP:
            sets.append(rows)
    if not sets:
        return None

    def median_erc(list_of_rows):
        vals = []
        for rows in list_of_rows:
            r = np.array(rows)
            if len(r) < 2:
                continue
            sub = Z[r]
            C = (sub @ sub.T) / n_branch
            iu = np.triu_indices(len(r), k=1)
            vals.append(C[iu])
        return float(np.median(np.concatenate(vals))) if vals else np.nan

    obs = median_erc(sets)
    rng = np.random.default_rng(seed)
    n_genes = Z.shape[0]
    null = []
    for _ in range(n_perm):
        shuffled = [rng.choice(n_genes, size=len(r), replace=False) for r in sets]
        null.append(median_erc(shuffled))
    null = np.array(null)
    return dict(observed=obs, null_median=float(np.median(null)),
                null_p95=float(np.percentile(null, 95)),
                null_max=float(null.max()),
                empirical_p=float((null >= obs).mean()), n_perm=len(null),
                n_sets=len(sets))


def erc_same_vs_cross_organ(d, paths, seed=0, n_perm=200):
    """Within a pathway, do members sharing an organ coevolve more tightly than
    members allocated to different organs?

    Carries its own permutation null for the same reason as the ERC test: the
    same-organ and different-organ pair sets are built from a few hundred genes,
    so pairs are not independent and a Mann-Whitney p near 0.02 means nothing on
    its own. The null shuffles the ORGAN LABELS across the annotated genes,
    keeping the pathway structure and the number of genes per organ fixed, and
    asks how often that produces a same-vs-cross gap this large."""
    Z, sym, rate, n_branch = load_rer()
    sym2row = {}
    for i, s in enumerate(sym):
        sym2row.setdefault(s, i)
    ensg2sym = d["symbol"].dropna().to_dict()

    def organs_of(g):
        o = set()
        v = d.at[g, "enriched_in"] if g in d.index else None
        if isinstance(v, str):
            o.add(v)
        v = d.at[g, "group_enriched_in"] if g in d.index else None
        if isinstance(v, str):
            o.update(v.split(","))
        return o

    same, cross = [], []
    for pid, (name, genes) in paths.items():
        info = [(g, sym2row[ensg2sym[g]], organs_of(g)) for g in genes
                if g in ensg2sym and ensg2sym[g] in sym2row]
        info = [x for x in info if x[2]]
        if len(info) < 4:
            continue
        for i in range(len(info)):
            for j in range(i + 1, len(info)):
                ri, rj = info[i][1], info[j][1]
                if ri == rj:
                    continue
                e = float((Z[ri] * Z[rj]).sum() / n_branch)
                (same if info[i][2] & info[j][2] else cross).append(e)
    if len(same) < 50 or len(cross) < 50:
        return None
    u, p = stats.mannwhitneyu(same, cross, alternative="two-sided")
    obs_delta = float(np.median(same) - np.median(cross))

    # permutation null: reassign organ labels across the annotated genes
    per_path = []
    for pid, (name, genes) in paths.items():
        info = [(sym2row[ensg2sym[g]], organs_of(g)) for g in genes
                if g in ensg2sym and ensg2sym[g] in sym2row]
        info = [x for x in info if x[1]]
        if len(info) >= 4:
            per_path.append(info)
    labels = [o for info in per_path for _, o in info]
    rng = np.random.default_rng(seed)
    null = []
    for _ in range(n_perm):
        perm = list(labels)
        rng.shuffle(perm)
        k, s_, c_ = 0, [], []
        for info in per_path:
            rows = [r for r, _ in info]
            labs = perm[k:k + len(info)]
            k += len(info)
            for i in range(len(rows)):
                for j in range(i + 1, len(rows)):
                    if rows[i] == rows[j]:
                        continue
                    e = float((Z[rows[i]] * Z[rows[j]]).sum() / n_branch)
                    (s_ if labs[i] & labs[j] else c_).append(e)
        if len(s_) >= 20 and len(c_) >= 20:
            null.append(np.median(s_) - np.median(c_))
    null = np.array(null)
    return dict(n_same=len(same), n_cross=len(cross),
                median_same=float(np.median(same)),
                median_cross=float(np.median(cross)),
                delta=obs_delta,
                rank_biserial=float(2 * u / (len(same) * len(cross)) - 1),
                p=float(p),
                null_median=float(np.median(null)) if len(null) else np.nan,
                null_p95=float(np.percentile(null, 95)) if len(null) else np.nan,
                empirical_p=float((null >= obs_delta).mean()) if len(null) else np.nan,
                n_perm=len(null))


def make_figure(alloc, alloc_stat, disp, disp_stat, erc, obs, nul, samecross):
    plt = apply_figure_style()
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4))

    axes[0].hist(alloc["log2_diff"], bins=50, color=PALETTE[0], alpha=0.85)
    axes[0].axvline(0, color="k", lw=1.0)
    axes[0].axvline(alloc_stat["median_log2_diff"], color=PALETTE[1], lw=1.5,
                    label=f"median {alloc_stat['median_log2_diff']:+.3f}")
    axes[0].set_xlabel(r"$\log_2$($\omega$ tissue-enriched / $\omega$ broad),"
                       "\nwithin the same pathway")
    axes[0].set_ylabel("pathways")
    axes[0].set_title(f"allocation: {alloc_stat['n_pathways']} pathways,\n"
                      f"{alloc_stat['frac_positive']:.0%} positive, "
                      f"p={alloc_stat['p']:.1g}")
    axes[0].legend(fontsize=7)

    if erc is not None:
        bins = np.linspace(min(nul.min(), obs.min()), max(nul.max(), obs.max()), 70)
        axes[1].hist(nul, bins=bins, density=True, alpha=0.65, color="#999999",
                     label=f"rate-matched null ({erc['median_null']:+.4f})")
        axes[1].hist(obs, bins=bins, density=True, alpha=0.65, color=PALETTE[1],
                     label=f"within pathway ({erc['median_within']:+.4f})")
        axes[1].set_xlabel("evolutionary rate covariation (ERC)")
        axes[1].set_ylabel("density")
        axes[1].set_title(f"do pathway members coevolve?\n"
                          f"rb={erc['rank_biserial']:+.3f}, p={erc['p']:.1g}")
        axes[1].legend(fontsize=6.5)

    axes[2].scatter(disp["n_organs"], disp["omega_iqr_log2"], s=8, alpha=0.35,
                    color=PALETTE[2])
    g = disp.groupby("n_organs")["omega_iqr_log2"].median()
    g = g[g.index <= 12]
    axes[2].plot(g.index, g.values, "o-", color=PALETTE[1], lw=1.6)
    axes[2].set_xlabel("distinct organs the pathway's members are enriched in")
    axes[2].set_ylabel(r"within-pathway $\log_2\omega$ IQR")
    axes[2].set_xlim(-0.5, 12.5)
    axes[2].set_title("organ spread vs rate spread\n"
                      r"$\rho$" f"={disp_stat['rho']:+.2f} (p={disp_stat['p']:.1g})")

    p = os.path.join(FIGDIR, "t15_pathway.png")
    fig.savefig(p)
    plt.close(fig)
    return p


def report(sp_key=t4_tissue_rates.PRIMARY):
    print("=" * 72)
    print("T15 - pathway allocation across organs, and pathway-level coevolution")
    print("=" * 72)
    d = t5_confounders.analysis_frame(sp_key)
    members = pathway_members()
    paths = usable_pathways(d, members)
    print(f"Reactome pathways with {MIN_PATHWAY}-{MAX_PATHWAY} members "
          f"in the rate table: {len(paths)}\n")

    alloc, stat = allocation_test(d, paths)
    print("1. ALLOCATION - within the same pathway, tissue-enriched vs broad members")
    print(f"   pathways testable            : {stat['n_pathways']}")
    print(f"   median log2(enriched/broad)  : {stat['median_log2_diff']:+.4f}"
          f"  ({2 ** stat['median_log2_diff']:.2f}x)")
    print(f"   pathways with a positive shift: {stat['frac_positive']:.1%}")
    print(f"   median tau gap (enriched-broad): {stat['median_tau_gap']:+.3f}")
    print(f"   Wilcoxon signed-rank p       : {stat['p']:.3g}")
    print("   -> pathway identity is held fixed, so this rules out the effect")
    print("      being tissue-specific genes sitting in different, faster")
    print("      pathways. It is NOT independent of tau - see the tau gap above")
    print("      and the fixed-effects model below.")

    print("\n   robustness:")
    for label, kw in [("stricter groups (>=5 per side)", dict(min_per_group=5)),
                      ("histones excluded", dict(drop_histones=True)),
                      ("both", dict(min_per_group=5, drop_histones=True))]:
        a2, s2 = allocation_test(d, paths, **kw)
        if s2:
            print(f"      {label:<32} n={s2['n_pathways']:>4}  "
                  f"median {s2['median_log2_diff']:+.3f}  "
                  f"{s2['frac_positive']:.0%} positive  p={s2['p']:.2g}")

    fe = pathway_fixed_effects(d, paths)
    print(f"\n   tau slope on log10(omega), n={fe['n']} genes in "
          f"{fe['n_pathways']} pathways:")
    print(f"      without pathway fixed effects : {fe['tau_coef']:+.4f}")
    print(f"      with    pathway fixed effects : {fe['tau_coef_fe']:+.4f} "
          f"(p={fe['tau_p_fe']:.2g})")
    print(f"      retained within pathways      : {fe['retained']:.0%}")
    print("\n   largest shifts (tissue-restricted members far faster):")
    for pid, r in alloc.head(8).iterrows():
        print(f"      {r['name'][:52]:<54} {r['log2_diff']:+.2f}  "
              f"({r['n_enriched']}v{r['n_broad']})")
    print("   largest reverse shifts:")
    for pid, r in alloc.tail(5).iloc[::-1].iterrows():
        print(f"      {r['name'][:52]:<54} {r['log2_diff']:+.2f}  "
              f"({r['n_enriched']}v{r['n_broad']})")

    disp, disp_stat = dispersion_test(d, paths)
    print(f"\n2. DISPERSION - organ spread vs rate spread across {disp_stat['n']} pathways")
    print(f"   Spearman(n organs, within-pathway log2 omega IQR) = "
          f"{disp_stat['rho']:+.3f} (p={disp_stat['p']:.2g})")

    print("\n3. ERC - do pathway members coevolve across 121 mammals?")
    erc_out = erc_analysis(d, paths)
    if erc_out is None:
        print("   not enough overlap with the T12 RER matrix")
        erc, obs, nul = None, None, None
    else:
        erc, obs, nul = erc_out
        print(f"   within-pathway pairs         : {erc['n_pairs']}")
        print(f"   median ERC within pathway    : {erc['median_within']:+.4f}")
        print(f"   median ERC, rate-matched null: {erc['median_null']:+.4f}")
        print(f"   difference                   : {erc['delta']:+.4f}  "
              f"(rb={erc['rank_biserial']:+.3f}, p={erc['p']:.2g})")
        print("   the null is rate-matched on purpose: genes of similar overall")
        print("   rate covary for reasons unrelated to shared function.")
        print(f"   NOTE: p={erc['p']:.2g} is NOT usable - {erc['n_pairs']} pairs are")
        print("   built from 4195 genes, so they are heavily non-independent and")
        print("   Mann-Whitney will call any effect significant. The effect size")
        print(f"   (rank-biserial {erc['rank_biserial']:+.3f}) is the honest reading,")
        print("   and the permutation null below is the honest test.")

    perm = erc_permutation_null(d, paths)
    if perm:
        print(f"\n   permutation null ({perm['n_perm']} shuffles of pathway "
              f"membership over {perm['n_sets']} sets, sizes preserved):")
        print(f"      observed median ERC : {perm['observed']:+.5f}")
        print(f"      null median         : {perm['null_median']:+.5f}")
        print(f"      null 95th pct       : {perm['null_p95']:+.5f}")
        print(f"      empirical p         : {perm['empirical_p']:.3f}")
        if perm["empirical_p"] > 0.05:
            print("      -> pathway coevolution is NOT distinguishable from random")
            print("         gene sets of the same size.")

    sc = erc_same_vs_cross_organ(d, paths)
    print("\n4. FRAGMENTATION - within a pathway, same-organ vs different-organ pairs")
    if sc is None:
        print("   too few annotated pairs")
    else:
        print(f"   same-organ pairs      : {sc['n_same']}, median ERC {sc['median_same']:+.4f}")
        print(f"   different-organ pairs : {sc['n_cross']}, median ERC {sc['median_cross']:+.4f}")
        print(f"   difference            : {sc['delta']:+.4f} "
              f"(rb={sc['rank_biserial']:+.3f}, Mann-Whitney p={sc['p']:.2g} "
              f"- not usable, pairs are non-independent)")
        if np.isfinite(sc.get("empirical_p", np.nan)):
            print(f"   permutation null ({sc['n_perm']} shuffles of organ labels):")
            print(f"      null median gap {sc['null_median']:+.4f}, "
                  f"95th pct {sc['null_p95']:+.4f}")
            print(f"      empirical p     {sc['empirical_p']:.3f}")
            if sc["empirical_p"] > 0.05:
                print("      -> organ-fragmented coevolution is NOT supported once")
                print("         the pair dependence is accounted for.")

    alloc.to_csv(os.path.join(DERIVED, "t15_allocation.tsv"), sep="\t")
    disp.to_csv(os.path.join(DERIVED, "t15_dispersion.tsv"), sep="\t")
    pd.DataFrame([{**(erc or {}), **{f"sc_{k}": v for k, v in (sc or {}).items()}}]
                 ).to_csv(os.path.join(DERIVED, "t15_erc.tsv"), sep="\t", index=False)
    if erc is not None:
        p = make_figure(alloc, stat, disp, disp_stat, erc, obs, nul, sc)
        print(f"\nfigure written to: {p}")
    return d, alloc, disp, erc, sc


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", default=t4_tissue_rates.PRIMARY)
    a = ap.parse_args()
    report(sp_key=a.species)
