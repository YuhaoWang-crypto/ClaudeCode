"""
T14 - what, if anything, explains the organ effect that survived T5 and T13?

The tissue effect has now survived three controls: tissue-specificity tau,
peak expression level, and subcellular compartment. This module throws the
three remaining standard gene-level rate determinants at it, in the same
mediator design, so each one gets a number rather than a hand-wave:

  * PPI DEGREE (STRING v12.0, combined_score >= 700). Highly connected proteins
    evolve slowly (Fraser et al. 2002) - every interaction is a surface under
    constraint. Organs plausibly differ: a synapse is a dense machine, a
    secreted plasma protein binds little.
  * PATHWAY PLEIOTROPY (Reactome, all levels). Number of distinct pathways a
    gene sits in. This is FUNCTIONAL pleiotropy, which tau does not capture -
    tau counts tissues, and a gene can be used in one tissue but twenty
    processes, or in thirty tissues for a single job.
  * GENE AGE, proxied by how many of the four ladder species retain a 1:1
    orthologue (0-4, i.e. survives to macaque / mouse / opossum / chicken).
    Young genes evolve fast; if an organ's gene set is young, its rate is a
    gene-age effect wearing a tissue costume.

Design is identical to T13 so the numbers are comparable: nested OLS for the
variance the mediator absorbs from the organ term, plus matched sampling with
the mediator added to (tau, peak expression) as a matching covariate. A mediator
that explains the organ effect must SHRINK the intrinsic effects; one that
merely correlates with rate will not.

Continuous mediators are matched on rather than stratified (unlike compartment,
which is categorical), so they enter the k-d tree alongside tau and expression.
"""
import os
import gzip
import argparse

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree

from .common import (RAW, DERIVED, FIGDIR, apply_figure_style, PALETTE,
                     biomart_query, download, log)
from . import t4_tissue_rates, t5_confounders, t7_esm

STRING_LINKS = os.path.join(RAW, "string_links.txt.gz")
STRING_INFO = os.path.join(RAW, "string_info.txt.gz")
REACTOME = os.path.join(RAW, "reactome_ncbi.txt")
STRING_LINKS_URL = ("https://stringdb-downloads.org/download/protein.links.v12.0/"
                    "9606.protein.links.v12.0.txt.gz")
STRING_INFO_URL = ("https://stringdb-downloads.org/download/protein.info.v12.0/"
                   "9606.protein.info.v12.0.txt.gz")
REACTOME_URL = "https://reactome.org/download/current/NCBI2Reactome_All_Levels.txt"

STRING_MIN_SCORE = 700
MIN_ENRICHED = 40
MEDIATORS = ["log_degree", "log_pathways", "n_ortholog_species"]


# --------------------------------------------------------------------------
# mediator 1: STRING PPI degree
# --------------------------------------------------------------------------
def string_degree():
    """gene symbol -> number of STRING partners at combined_score >= 700."""
    cache = os.path.join(DERIVED, "t14_string_degree.tsv")
    if os.path.exists(cache):
        s = pd.read_csv(cache, sep="\t", index_col=0)["degree"]
        return s
    for path, url in ((STRING_LINKS, STRING_LINKS_URL),
                      (STRING_INFO, STRING_INFO_URL)):
        if not os.path.exists(path):
            log(f"downloading {os.path.basename(path)}")
            download(url, path, timeout=1800)

    ensp2sym = {}
    with gzip.open(STRING_INFO, "rt") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                ensp2sym[f[0]] = f[1]

    deg = {}
    with gzip.open(STRING_LINKS, "rt") as fh:
        next(fh)
        for line in fh:
            a, b, sc = line.split()
            if int(sc) < STRING_MIN_SCORE:
                continue
            # the file lists each undirected edge in both orientations, so
            # counting one endpoint per line already gives the true degree
            sa = ensp2sym.get(a)
            if sa:
                deg[sa] = deg.get(sa, 0) + 1
    s = pd.Series(deg, name="degree")
    s.to_frame().to_csv(cache, sep="\t")
    log(f"STRING degree for {len(s)} symbols at score>={STRING_MIN_SCORE}")
    return s


# --------------------------------------------------------------------------
# mediator 2: Reactome pathway count
# --------------------------------------------------------------------------
def entrez_to_ensg():
    cache = os.path.join(RAW, "biomart_entrez.tsv")
    if not os.path.exists(cache):
        txt = biomart_query("hsapiens_gene_ensembl",
                            ["ensembl_gene_id", "entrezgene_id"])
        open(cache, "w").write(txt)
    out = {}
    for line in open(cache).read().rstrip("\n").split("\n")[1:]:
        f = line.split("\t")
        if len(f) >= 2 and f[1]:
            out.setdefault(f[1].split(".")[0], set()).add(f[0])
    return out


def pathway_counts():
    """ENSG -> number of distinct Reactome pathways (all levels)."""
    cache = os.path.join(DERIVED, "t14_pathway_counts.tsv")
    if os.path.exists(cache):
        return pd.read_csv(cache, sep="\t", index_col=0)["n_pathways"]
    if not os.path.exists(REACTOME):
        log("downloading Reactome NCBI2Reactome_All_Levels")
        download(REACTOME_URL, REACTOME, timeout=1800)
    e2g = entrez_to_ensg()
    counts = {}
    with open(REACTOME) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6 or f[5] != "Homo sapiens":
                continue
            for g in e2g.get(f[0], ()):
                counts.setdefault(g, set()).add(f[1])
    s = pd.Series({g: len(v) for g, v in counts.items()}, name="n_pathways")
    s.to_frame().to_csv(cache, sep="\t")
    log(f"Reactome pathway membership for {len(s)} genes")
    return s


# --------------------------------------------------------------------------
# assembly
# --------------------------------------------------------------------------
def analysis_frame(sp_key=t4_tissue_rates.PRIMARY):
    d = t5_confounders.analysis_frame(sp_key)
    deg = string_degree()
    d["degree"] = d["symbol"].map(deg)
    d["log_degree"] = np.log10(d["degree"].fillna(0.0) + 1.0)
    d["n_pathways"] = pathway_counts().reindex(d.index)
    d["log_pathways"] = np.log10(d["n_pathways"].fillna(0.0) + 1.0)
    breadth, _ = t7_esm.ortholog_breadth()
    d["n_ortholog_species"] = breadth.reindex(d.index).fillna(0.0)
    return d


def mediator_marginals(d):
    rows = []
    for m in MEDIATORS:
        sub = d.dropna(subset=[m, "log_omega"])
        rho, p = stats.spearmanr(sub[m], sub["log_omega"])
        rp, _ = stats.spearmanr(sub[m], sub["tau"])
        rows.append(dict(mediator=m, n=len(sub), rho_with_omega=rho, p=p,
                         rho_with_tau=rp))
    return pd.DataFrame(rows).set_index("mediator")


def organ_profiles(d):
    rows = []
    for organ in t5_confounders.organ_list(d):
        sub = d[t5_confounders.organ_membership(d, organ)]
        if len(sub) < MIN_ENRICHED:
            continue
        rec = dict(organ=organ, n=len(sub), median_omega=sub["omega"].median())
        for m in MEDIATORS:
            rec[m] = float(sub[m].median())
        rows.append(rec)
    return pd.DataFrame(rows).set_index("organ")


# --------------------------------------------------------------------------
# nested OLS
# --------------------------------------------------------------------------
def nested_ols(d):
    import statsmodels.formula.api as smf
    sub = d[d["enriched_in"].notna() | d["group_enriched_in"].notna()].copy()
    sub["organ_set"] = sub["enriched_in"].fillna(sub["group_enriched_in"])
    sub = sub.dropna(subset=MEDIATORS)
    base = "log_omega ~ log_expr + tau + log_codons"
    out = []
    m_base = smf.ols(base, data=sub).fit()
    m_org = smf.ols(base + " + C(organ_set)", data=sub).fit()
    dr2_alone = m_org.rsquared - m_base.rsquared
    for label, extra in ([(m, f" + {m}") for m in MEDIATORS]
                         + [("all three", " + " + " + ".join(MEDIATORS))]):
        m_med = smf.ols(base + extra, data=sub).fit()
        m_both = smf.ols(base + extra + " + C(organ_set)", data=sub).fit()
        dr2_given = m_both.rsquared - m_med.rsquared
        out.append(dict(mediator=label,
                        r2_with_mediator=float(m_med.rsquared),
                        dr2_mediator=float(m_med.rsquared - m_base.rsquared),
                        dr2_organ_given=float(dr2_given),
                        frac_absorbed=float(1 - dr2_given / dr2_alone)))
    return pd.DataFrame(out).set_index("mediator"), float(dr2_alone), int(m_base.nobs)


# --------------------------------------------------------------------------
# matched sampling with extra covariates
# --------------------------------------------------------------------------
def matched_control_extra(d, organ, extra=(), seed=0, caliper=0.25):
    """T5 matching with `extra` continuous covariates added to the k-d tree.

    The caliper is a Euclidean distance in standardised units over ALL matching
    covariates, so adding a covariate makes matching strictly harder; the number
    of unmatched cases is returned so a shrinking effect cannot be confused with
    a shrinking sample.
    """
    member = t5_confounders.organ_membership(d, organ)
    feats = ["tau", "log_expr"] + list(extra)
    dd = d.dropna(subset=feats)
    member = member.reindex(dd.index).fillna(False)
    case, pool = dd[member], dd[~member]
    if len(case) < MIN_ENRICHED or len(pool) < len(case):
        return None
    mu, sd = dd[feats].mean(), dd[feats].std(ddof=0).replace(0, 1.0)
    C = ((case[feats] - mu) / sd).to_numpy()
    P = ((pool[feats] - mu) / sd).to_numpy()
    tree = cKDTree(P)
    k = min(len(P), 300)
    dists, nbrs = tree.query(C, k=k)
    rng = np.random.default_rng(seed)
    used, ci, pi = set(), [], []
    cal = caliper * np.sqrt(len(feats) / 2.0)   # keep per-covariate tightness
    for i in rng.permutation(len(C)):
        for dist, j in zip(np.atleast_1d(dists[i]), np.atleast_1d(nbrs[i])):
            if j in used or dist > cal:
                continue
            used.add(int(j))
            ci.append(i)
            pi.append(int(j))
            break
    if len(ci) < MIN_ENRICHED:
        return None
    return case.iloc[ci], pool.iloc[pi], len(case) - len(ci)


def mediation_test(d, seed=0):
    rows = []
    for organ in t5_confounders.organ_list(d):
        base = matched_control_extra(d, organ, extra=(), seed=seed)
        if base is None:
            continue
        c0, k0, _ = base
        i0 = float(np.log2(c0["omega"].median() / k0["omega"].median()))
        rec = dict(organ=organ, n_base=len(c0), intrinsic_base=i0)
        for m in MEDIATORS + ["all"]:
            extra = MEDIATORS if m == "all" else [m]
            r = matched_control_extra(d, organ, extra=extra, seed=seed)
            if r is None:
                rec[f"intrinsic_{m}"] = np.nan
                continue
            c1, k1, drop = r
            rec[f"intrinsic_{m}"] = float(
                np.log2(c1["omega"].median() / k1["omega"].median()))
            rec[f"n_{m}"] = len(c1)
            if m == "all":
                u, p = stats.mannwhitneyu(c1["omega"], k1["omega"],
                                          alternative="two-sided")
                rec["p_all"] = float(p)
        rows.append(rec)
    res = pd.DataFrame(rows).set_index("organ")
    if "p_all" in res:
        ok = res["p_all"].notna()
        q = np.full(len(res), np.nan)
        q[ok.to_numpy()] = t5_confounders._bh(res.loc[ok, "p_all"].to_numpy())
        res["q_all"] = q
    return res.sort_values("intrinsic_base")


def make_figure(marg, prof, med, d):
    plt = apply_figure_style()
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.6))

    for i, m in enumerate(MEDIATORS):
        sub = d.dropna(subset=[m, "omega"])
        try:
            q = pd.qcut(sub[m], 6, duplicates="drop")
        except ValueError:
            continue
        g = sub.groupby(q, observed=True)["omega"].median()
        axes[0].plot(range(len(g)), g.values, "o-", color=PALETTE[i],
                     label=f"{m} (" r"$\rho$" f"={marg.loc[m, 'rho_with_omega']:+.2f})")
    axes[0].set_xlabel("mediator sextile (low to high)")
    axes[0].set_ylabel(r"median $\omega$")
    axes[0].set_title("each mediator does predict rate")
    axes[0].legend(fontsize=6.5)

    p = prof.sort_values("median_omega")
    axes[1].scatter(p["log_degree"], p["median_omega"], s=32, color=PALETTE[0])
    rho, pv = stats.spearmanr(p["log_degree"], p["median_omega"])
    for organ in p.index:
        axes[1].annotate(organ, (p.loc[organ, "log_degree"],
                                 p.loc[organ, "median_omega"]),
                         fontsize=6, xytext=(3, 2), textcoords="offset points")
    axes[1].set_xlabel(r"organ median $\log_{10}$(STRING degree + 1)")
    axes[1].set_ylabel(r"organ median $\omega$")
    axes[1].set_title(f"organ connectivity vs rate\n" r"$\rho$"
                      f"={rho:+.2f} (p={pv:.1g})")

    m = med.sort_values("intrinsic_base")
    y = np.arange(len(m))
    axes[2].barh(y - 0.2, m["intrinsic_base"], height=0.4, color=PALETTE[7],
                 label=r"$\tau$ + expression only")
    axes[2].barh(y + 0.2, m["intrinsic_all"], height=0.4, color=PALETTE[1],
                 label="+ degree, pathways, age")
    axes[2].set_yticks(y)
    axes[2].set_yticklabels(m.index, fontsize=7)
    axes[2].axvline(0, color="k", lw=0.9)
    axes[2].set_xlabel(r"intrinsic $\log_2$ FC vs matched controls")
    axes[2].set_title("do the mediators absorb\nthe organ effect?")
    axes[2].grid(axis="y", alpha=0)
    axes[2].legend(fontsize=6.5, loc="lower right")

    path = os.path.join(FIGDIR, "t14_mediators.png")
    fig.savefig(path)
    plt.close(fig)
    return path


def report(sp_key=t4_tissue_rates.PRIMARY):
    print("=" * 72)
    print("T14 - PPI degree, pathway pleiotropy and gene age as mediators")
    print("=" * 72)
    d = analysis_frame(sp_key)
    print(f"genes analysed: {len(d)}")
    print(f"  with STRING degree   : {d['degree'].notna().sum()}")
    print(f"  with Reactome pathway: {d['n_pathways'].notna().sum()}\n")

    marg = mediator_marginals(d)
    print("A. does each mediator predict rate on its own?")
    print(f"   {'mediator':<22}{'n':>7}{'rho vs omega':>14}{'rho vs tau':>13}{'p':>11}")
    for m, r in marg.iterrows():
        print(f"   {m:<22}{int(r['n']):>7}{r['rho_with_omega']:>+14.3f}"
              f"{r['rho_with_tau']:>+13.3f}{r['p']:>11.2g}")

    prof = organ_profiles(d)
    print("\nB. organ medians on each mediator")
    print(f"   {'organ':<18}{'n':>5}{'omega':>9}{'log_deg':>9}{'log_path':>10}{'age':>7}")
    for organ, r in prof.sort_values("median_omega", ascending=False).iterrows():
        print(f"   {organ:<18}{int(r['n']):>5}{r['median_omega']:>9.4f}"
              f"{r['log_degree']:>9.2f}{r['log_pathways']:>10.2f}"
              f"{r['n_ortholog_species']:>7.1f}")

    ols, dr2_alone, n = nested_ols(d)
    print(f"\nC. nested OLS (n={n}); organ dR2 with no mediator = {dr2_alone:.4f}")
    print(f"   {'mediator':<22}{'its own dR2':>13}{'organ dR2 given':>17}{'absorbed':>10}")
    for m, r in ols.iterrows():
        print(f"   {m:<22}{r['dr2_mediator']:>13.4f}{r['dr2_organ_given']:>17.4f}"
              f"{r['frac_absorbed']:>9.1%}")

    med = mediation_test(d)
    print("\nD. matched sampling, mediators added as matching covariates")
    print(f"   {'organ':<18}{'base':>8}{'+degree':>9}{'+paths':>9}{'+age':>8}"
          f"{'+all':>8}{'n':>6}{'q':>10}")
    def _f(v, w, prec=3, sign=True):
        if not np.isfinite(v):
            return f"{'-':>{w}}"
        return f"{v:>+{w}.{prec}f}" if sign else f"{v:>{w}.{prec}f}"

    for organ, r in med.iterrows():
        q = r.get("q_all", np.nan)
        n_all = r.get("n_all", np.nan)
        star = "*" if np.isfinite(q) and q < 0.05 else " "
        qs = f"{q:>10.2e}" if np.isfinite(q) else f"{'-':>10}"
        ns = f"{int(n_all):>6}" if np.isfinite(n_all) else f"{'-':>6}"
        print(f"   {organ:<18}{_f(r['intrinsic_base'], 8)}"
              f"{_f(r['intrinsic_log_degree'], 9)}"
              f"{_f(r['intrinsic_log_pathways'], 9)}"
              f"{_f(r['intrinsic_n_ortholog_species'], 8)}"
              f"{_f(r['intrinsic_all'], 8)}{ns}{qs}{star}")

    marg.to_csv(os.path.join(DERIVED, "t14_mediator_marginals.tsv"), sep="\t")
    prof.to_csv(os.path.join(DERIVED, "t14_organ_profiles.tsv"), sep="\t")
    ols.to_csv(os.path.join(DERIVED, "t14_nested_ols.tsv"), sep="\t")
    med.to_csv(os.path.join(DERIVED, "t14_mediation.tsv"), sep="\t")
    p = make_figure(marg, prof, med, d)
    print(f"\nfigure written to: {p}")
    return d, marg, prof, ols, med


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", default=t4_tissue_rates.PRIMARY)
    a = ap.parse_args()
    report(sp_key=a.species)
