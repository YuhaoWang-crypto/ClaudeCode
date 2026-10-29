"""Reproduce the paper's headline numbers from our recomputed scores.

Every number printed here is computed from public data; the paper's value is
shown next to it so agreement (or not) is explicit. Results go to
gismo/data/benchmark.json and figures/gismo_benchmark.png.

Not reproducible from public data (skipped, stated in the report):
  * UK Biobank exome CNV deletion/duplication rates (Fig. 2g,h)
  * LDSC partitioned heritability (needs GWAS sumstats; authors' outputs ship in
    the Zenodo archive and are summarised by `heritability_from_archive`).
"""
import gzip
import json

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, spearmanr

from gismo import paths
from gismo.core import load_scores

PAPER = {  # value reported in Liao et al. 2026
    "rho_gismo_loeuf": 0.38, "rho_mis_loeuf": 0.56,
    "rho_gismo_moeuf": 0.33, "rho_mis_moeuf": 0.61,
    "rho_gismo_phaplo": -0.32, "rho_mis_phaplo": -0.44,
    "rho_gismo_ptriplo": -0.32, "rho_mis_ptriplo": -0.50,
    "rho_mis_alphamissense": 0.84,
    "rho_gismo_cds": -0.12, "rho_mis_cds": -0.12, "rho_loeuf_cds": -0.50,
    "OR_d1_mouse_het_lethal_gismo": 3.43, "OR_d1_mouse_het_lethal_mis": 5.78,
    "OR_d1_mouse_lethal_gismo": 8.21, "OR_d1_mouse_lethal_mis": 13.1,
    "OR_d1_cell_essential_gismo": 4.34, "OR_d1_cell_essential_mis": 10.70,
    "median_cds_top15_gismo": 1502, "median_cds_top15_mis": 1302,
    "median_cds_loeuf_lt_0.35": 2267,
    "OR_d1_recessive_gismo": 3.37, "OR_d1_recessive_mis": 0.91,
    "OR_d1_dominant_gismo": 9.92, "OR_d1_dominant_mis": 4.38,
    "median_losses_per_gene": 18, "median_losses_per_species": 1130,
    "median_total_missense_per_species": 610603,
    "median_total_synonymous_per_species": 1084983,
    "median_mean_missense": 27.7, "median_mean_synonymous": 53.4,
    "n_zero_loss_genes": 6,
}


# Approximate bar heights read off Fig. 4c,d (⚠️ visual digitisation, ±1%).
FIG4_APPROX = {
    "recessive_gismo": [12.8, 12.8, 13.4, 12.7, 12.4, 11.0, 11.0, 9.2, 5.1, 1.7],
    "recessive_mis": [3.2, 5.7, 10.3, 12.3, 12.3, 14.8, 11.7, 12.9, 9.0, 5.6],
    "dominant_gismo": [15.4, 14.0, 10.5, 11.7, 9.8, 9.2, 8.3, 10.9, 8.5, 4.5],
    "dominant_mis": [15.4, 14.6, 12.7, 12.0, 10.8, 8.0, 7.9, 8.1, 4.4, 3.6],
}


def _gene_list(name, upper=False):
    s = pd.read_csv(paths.ref(name), sep="\t", header=None, comment="#")[0].astype(str)
    return set(s.str.upper() if upper else s)


def load_reference():
    v2 = pd.read_csv(paths.ref("gnomad_v2.1.1_lof.txt.bgz"), sep="\t", compression="gzip",
                     usecols=["gene", "oe_lof_upper", "oe_mis_upper", "cds_length", "chromosome"])
    v2 = v2.dropna(subset=["oe_lof_upper"]).sort_values("oe_lof_upper").drop_duplicates("gene")
    v2 = v2.rename(columns={"oe_lof_upper": "loeuf", "oe_mis_upper": "moeuf"})
    v4 = pd.read_csv(paths.ref("gnomad_v4.1_constraint.tsv"), sep="\t",
                     usecols=["gene", "mane_select", "lof.oe_ci.upper"], low_memory=False)
    v4 = v4[v4.mane_select.astype(str).str.lower() == "true"].drop_duplicates("gene")
    v4 = v4.rename(columns={"lof.oe_ci.upper": "loeuf_v4"})[["gene", "loeuf_v4"]]
    ph = pd.read_csv(paths.ref("phaplo_ptriplo.tsv.gz"), sep="\t").rename(columns={"#gene": "gene"})
    am = pd.read_csv(paths.ref("am_gene.tsv.gz"), sep="\t", comment="#")
    am["transcript_id"] = am.transcript_id.str.split(".").str[0]
    return v2, v4, ph, am


def merged():
    d = load_scores()
    v2, v4, ph, am = load_reference()
    d = d.merge(v2, left_on="symbol", right_on="gene", how="left").drop(columns="gene")
    d = d.merge(v4, left_on="symbol", right_on="gene", how="left").drop(columns="gene")
    d = d.merge(ph, left_on="symbol", right_on="gene", how="left").drop(columns="gene")
    d = d.merge(am, on="transcript_id", how="left")
    return d


def _rho(a, b):
    m = a.notna() & b.notna()
    return round(float(spearmanr(a[m], b[m])[0]), 3), int(m.sum())


def _or(dec, members, top=(1,), bottom=None):
    """Fisher OR for a gene set: genes in `top` deciles vs `bottom` deciles
    (default: all other deciles)."""
    ok = dec.notna()
    if bottom is not None:
        ok &= dec.isin(tuple(top) + tuple(bottom))
    x = dec[ok].isin(top).to_numpy()
    y = members[ok].to_numpy()
    t = [[(x & y).sum(), (x & ~y).sum()], [(~x & y).sum(), (~x & ~y).sum()]]
    o, p = fisher_exact(t)
    return round(float(o), 2), float(p), int(y.sum())


def _logit_or(dec, members):
    """OR per decile step from logistic regression (membership ~ decile),
    scaled to decile 1 vs decile 10 (exp(9*beta))."""
    from scipy.optimize import minimize
    ok = dec.notna()
    x = (10 - dec[ok].to_numpy(float)) / 9.0
    y = members[ok].to_numpy(float)

    def nll(b):
        z = b[0] + b[1] * x
        return np.sum(np.logaddexp(0, z) - y * z)
    b = minimize(nll, [np.log(y.mean() / (1 - y.mean())), 0.0], method="BFGS").x
    return round(float(np.exp(b[1])), 2)


def descriptive():
    """Per-gene / per-species loss and mis/syn counts (Fig. 1)."""
    from gismo.core import load_loss_matrix, load_mis_syn, load_species
    root = paths.zenodo_root()
    sp = load_species(root)
    loss, sp_used, _ = load_loss_matrix(root, sp)
    mis, syn = load_mis_syn(root, sp_used)
    return {
        "n_genes": int(loss.shape[0]), "n_species": int(loss.shape[1]),
        "median_losses_per_gene": float(np.median(loss.sum(axis=1))),
        # paper's "1,130 per species" = TOGA summary lost + uncertain_loss;
        # the one2zero count used inside GISMO is higher
        "median_losses_per_species": float((sp_used.lost + sp_used.uncertain_loss).median()),
        "median_one2zero_per_species": float(np.median(loss.sum(axis=0))),
        "median_mean_missense": round(float(mis.mean(axis=0).median()), 1),
        "median_mean_synonymous": round(float(syn.mean(axis=0).median()), 1),
        "median_total_missense_per_species": float(mis.sum(axis=1).median()),
        "median_total_synonymous_per_species": float(syn.sum(axis=1).median()),
        "n_zero_loss_genes": int((loss.sum(axis=1) == 0).sum()),
    }


def run(fig=True):
    d = merged()
    r = {"descriptive": descriptive()}

    # --- correlations (Supplementary Fig. 5, Fig. 2f, Fig. 3a)
    c = {}
    for k, col in [("loeuf", "loeuf"), ("moeuf", "moeuf"), ("phaplo", "pHaplo"),
                   ("ptriplo", "pTriplo"), ("cds", "cds_length"), ("loeuf_v4", "loeuf_v4")]:
        c[f"rho_gismo_{k}"] = _rho(d.gismo, d[col])
        c[f"rho_mis_{k}"] = _rho(d.gismo_mis, d[col])
    c["rho_mis_alphamissense"] = _rho(d.gismo_mis, d.mean_am_pathogenicity)
    c["rho_gismo_alphamissense"] = _rho(d.gismo, d.mean_am_pathogenicity)
    c["rho_loeuf_cds"] = _rho(d.loeuf, d.cds_length)
    c["rho_gismo_mis"] = _rho(d.gismo, d.gismo_mis)
    r["correlations"] = c

    # --- decile enrichment (Fig. 2a-e)
    sym = d.symbol.astype(str)
    sets = {
        "mouse_het_lethal": sym.isin(_gene_list("mouse_het_lethal.tsv", upper=True)),
        "mouse_lethal": sym.isin(_gene_list("mouse_lethal.tsv", upper=True)),
        "cell_essential": sym.isin(_gene_list("CEGv2.tsv")),
        "nonessential": sym.isin(_gene_list("NEGv1.tsv")),
        "olfactory": sym.isin(_gene_list("olfactory_receptors.tsv")),
        "recessive": sym.isin(_gene_list("all_ar.tsv")),
        "dominant": sym.isin(_gene_list("all_ad.tsv")),
    }
    enr, pct = {}, {}
    for name, mem in sets.items():
        for tag, col in [("gismo", "gismo_decile"), ("mis", "gismo_mis_decile")]:
            # The paper does not state its OR definition; report several.
            enr[f"OR_d1_{name}_{tag}"] = _or(d[col], mem)                       # d1 vs d2-10
            enr[f"OR_d1vd10_{name}_{tag}"] = _or(d[col], mem, (1,), (10,))
            enr[f"OR_d12vd910_{name}_{tag}"] = _or(d[col], mem, (1, 2), (9, 10))
            enr[f"OR_logit_{name}_{tag}"] = _logit_or(d[col], mem)
            ok = d[col].notna()
            # Fig. 2 / Fig. 4c,d y-axis: % of the gene set falling in each decile
            pct[f"{name}_{tag}"] = (mem[ok].groupby(d.loc[ok, col]).sum()
                                    / mem[ok].sum() * 100).round(2).tolist()
    r["enrichment"], r["pct_by_decile"] = enr, pct

    # --- short constrained genes (Fig. 3b): top 15% most constrained
    s = {}
    for tag, col in [("gismo", "gismo"), ("mis", "gismo_mis")]:
        x = d.dropna(subset=[col, "cds_length"])
        top = x[x[col] <= x[col].quantile(0.15)]
        s[f"median_cds_top15_{tag}"] = (float(top.cds_length.median()), len(top))
    lo = d.dropna(subset=["loeuf", "cds_length"])
    s["median_cds_loeuf_lt_0.35"] = (float(lo[lo.loeuf < 0.35].cds_length.median()),
                                     int((lo.loeuf < 0.35).sum()))
    both = d.dropna(subset=["loeuf", "gismo", "cds_length"])
    gi = both.gismo <= both.gismo.quantile(0.15)
    s["median_cds_gismo_int_loeuf_tol"] = (float(both[gi & (both.loeuf > 0.35)].cds_length.median()),
                                           int((gi & (both.loeuf > 0.35)).sum()))
    s["median_cds_gismo_int_loeuf_int"] = (float(both[gi & (both.loeuf < 0.35)].cds_length.median()),
                                           int((gi & (both.loeuf < 0.35)).sum()))
    r["short_genes"] = s

    # --- recessive vs dominant: GISMO captures homozygous selection
    # Rule from the paper: GISMO low & LOEUF high => recessive-like.
    q = d.dropna(subset=["gismo_decile", "loeuf"]).copy()
    q["class"] = np.select(
        [(q.gismo_decile <= 2) & (q.loeuf < 0.35), (q.gismo_decile <= 2) & (q.loeuf >= 0.6)],
        ["both_constrained", "gismo_only"], "other")
    rec, dom = sets["recessive"][q.index], sets["dominant"][q.index]
    r["recessive_rule"] = {k: {"n": int((q["class"] == k).sum()),
                               "pct_AR": round(float(rec[q["class"] == k].mean() * 100), 2),
                               "pct_AD": round(float(dom[q["class"] == k].mean() * 100), 2)}
                           for k in ["both_constrained", "gismo_only", "other"]}

    # --- X-chromosome discordance (LOEUF constrained, GISMO least constrained)
    x = d.dropna(subset=["loeuf", "gismo_decile"])
    disc = (x.loeuf < 0.35) & (x.gismo_decile >= 9)
    isx = x.chrom.astype(str).isin(["X", "chrX"])
    r["x_chrom"] = {"pct_X_discordant": round(float(isx[disc].mean() * 100), 2),
                    "pct_X_all": round(float(isx.mean() * 100), 2),
                    "n_discordant": int(disc.sum()),
                    "OR": _or(pd.Series(np.where(disc, 1, 2), index=x.index), isx)}

    r["gc_bias"] = gc_bias(d, sets)
    r["heritability_archive"] = heritability_from_archive()
    r["fig4_vs_recomputed"] = {
        k: {"paper_approx": v, "recomputed": pct[k],
            "max_abs_diff": round(float(np.max(np.abs(np.array(v) - np.array(pct[k])))), 2)}
        for k, v in FIG4_APPROX.items()}
    r["paper"] = PAPER

    out = paths.PKG_DATA / "benchmark.json"
    out.write_text(json.dumps(r, indent=1, default=float))
    if fig:
        plot(r)
    return r


def gc_bias(d, sets):
    """GC-content artefact in GISMO (not reported in the paper)."""
    if "gc_cds" not in d or d.gc_cds.isna().all():
        return None
    ok = d.gismo.notna() & d.gc_cds.notna()
    ml = sets["mouse_lethal"]
    out = {"rho_gc_gismo_all": _rho(d.gismo[ok], d.gc_cds[ok]),
           "rho_gc_gismo_mouse_lethal": _rho(d.gismo[ok & ml], d.gc_cds[ok & ml]),
           "rho_gc_gismo_loeuf_lt_0.35": _rho(d.gismo[ok & (d.loeuf < 0.35)],
                                              d.gc_cds[ok & (d.loeuf < 0.35)]),
           "rho_gc_gismo_hq": _rho(d.gismo_hq[ok], d.gc_cds[ok]),
           "rho_gc_gismo_mis": _rho(d.gismo_mis, d.gc_cds)}
    q = pd.qcut(d.gc_cds[ok], 5, labels=False) + 1
    out["mouse_lethal_pct_in_d9_10_by_gc_quintile"] = [
        round(float((d.gismo_decile[ok & ml & (q == k)] >= 9).mean() * 100), 1)
        for k in range(1, 6)]
    hi = ok & (d.gismo_decile >= 8)
    flag = hi & (d.gc_cds_pct >= 0.8)
    rest = hi & (d.gc_cds_pct < 0.8)
    out["flag_d8_10_gc_top20"] = {
        k: {"n": int(m.sum()), "pct_mouse_lethal": round(float(ml[m].mean() * 100), 1),
            "pct_loeuf_lt_0.35": round(float((d.loeuf[m] < 0.35).mean() * 100), 1),
            "pct_phaplo_gt_0.86": round(float((d.pHaplo[m] > 0.86).mean() * 100), 1)}
        for k, m in [("flagged", flag), ("other_d8_10", rest)]}
    return out


def heritability_from_archive():
    """Mean LDSC enrichment per decile, from the authors' shipped outputs."""
    root = paths.zenodo_root() / "partitioned-heritability"
    res = {}
    for tag in ["GISMO", "GISMO-mis"]:
        f = root / tag / "outputs" / "all-traits_merged_result.tsv"
        if not f.exists():
            continue
        t = pd.read_csv(f, sep="\t")
        t = t[t.group.astype(str).str.match(r"d\d+$")]
        m = t.groupby("group").Enrichment.mean()
        m.index = m.index.str[1:].astype(int)
        res[tag] = {"n_traits": int(t.trait.nunique()),
                    "mean_enrichment_by_decile": m.sort_index().round(3).tolist()}
    return res


def plot(r):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    names = ["mouse_het_lethal", "mouse_lethal", "cell_essential", "nonessential",
             "olfactory", "recessive", "dominant"]
    fig, ax = plt.subplots(2, 4, figsize=(17, 7.5))
    x = np.arange(1, 11)
    for a, n in zip(ax.flat, names):
        a.bar(x - 0.2, r["pct_by_decile"][f"{n}_gismo"], 0.4, color="#b5838d", label="GISMO")
        a.bar(x + 0.2, r["pct_by_decile"][f"{n}_mis"], 0.4, color="#6a9fb5", label="GISMO-mis")
        a.set_title(n.replace("_", " ")); a.set_xticks(x); a.set_xlabel("decile (1 = most constrained)")
        a.set_ylabel("% of gene set in decile")
    ax[0, 0].legend()
    c = r["correlations"]; p = r["paper"]
    keys = ["loeuf", "moeuf", "phaplo", "ptriplo", "cds"]
    a = ax[1, 3]
    ours = [c[f"rho_gismo_{k}"][0] for k in keys] + [c[f"rho_mis_{k}"][0] for k in keys]
    theirs = [p[f"rho_gismo_{k}"] for k in keys] + [p[f"rho_mis_{k}"] for k in keys]
    a.scatter(theirs, ours, c=["#b5838d"] * 5 + ["#6a9fb5"] * 5)
    for t, o, k in zip(theirs, ours, keys * 2):
        a.annotate(k, (t, o), fontsize=7)
    a.plot([-0.7, 0.7], [-0.7, 0.7], "k--", lw=0.8)
    a.set_xlabel("paper Spearman ρ"); a.set_ylabel("recomputed ρ")
    a.set_title("correlations: paper vs recomputed")
    fig.suptitle("GISMO / GISMO-mis recomputed from Zenodo inputs (Liao et al. 2026)")
    fig.tight_layout()
    paths.FIG_DIR.mkdir(exist_ok=True)
    fig.savefig(paths.FIG_DIR / "gismo_benchmark.png", dpi=110)


def table(r):
    """Paper-vs-recomputed rows for the report."""
    flat = {}
    for sec in ["descriptive", "correlations", "enrichment", "short_genes"]:
        for k, v in r[sec].items():
            flat[k] = v[0] if isinstance(v, (list, tuple)) else v
    rows = []
    for k, pv in PAPER.items():
        if k.startswith("OR_"):
            # OR_d1_<set>_<tag> or OR_<set>_<tag>: show every definition
            core = k.replace("OR_d1_", "").replace("OR_", "")
            rows.append((k, pv, " / ".join(str(flat.get(f"OR_{v}_{core}")) for v in
                                             ["d1", "d1vd10", "d12vd910", "logit"])))
        elif k in flat:
            rows.append((k, pv, flat[k]))
    return pd.DataFrame(rows, columns=["metric", "paper",
                                       "recomputed (OR: d1-vs-rest / d1-vs-d10 / d12-vs-d910 / logit)"])


if __name__ == "__main__":
    res = run()
    pd.set_option("display.width", 200)
    print(table(res).to_string(index=False))
    for k in ["recessive_rule", "x_chrom", "heritability_archive"]:
        print(k, json.dumps(res[k], default=float))
    print("gc_bias", json.dumps(res["gc_bias"], default=float))
    for k, v in res["fig4_vs_recomputed"].items():
        print(k, "max|Δ| vs Fig.4 =", v["max_abs_diff"])
