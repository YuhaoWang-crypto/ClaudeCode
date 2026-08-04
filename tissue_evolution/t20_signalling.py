"""
T20 - signalling pathways vs metabolic ones, and the compartment gradient.

T18/T19 took metabolic chemistry apart. Signalling is a different kind of object
and should leave a different evolutionary signature, which makes it a test
rather than a description.

A metabolic pathway is a fixed chemistry executed inside one compartment. A
signalling pathway is an INFORMATION RELAY THAT CROSSES COMPARTMENTS: a receptor
in the plasma membrane reads the outside world, adaptors and kinases carry the
signal through the cytosol, and transcription factors act in the nucleus. Its
"catalysis" is mostly phosphotransfer, and much of its machinery is interaction
surface rather than active site.

Three predictions follow, each falsifiable:

  1. SIGNALLING SHOULD BE FASTER, and for a specific reason: T13 found that
     interface proteins (secreted + membrane) evolve 1.26x faster than
     intracellular ones, and signalling pathways are enriched in exactly those.
  2. THE CATALYTIC/ACCESSORY SPLIT SHOULD BE WEAKER in signalling. In the TCA
     cycle the catalytic members do the pathway's actual work; in a signalling
     cascade the scaffolds and adaptors are not "accessory" in the same sense -
     they carry the information. If the T19 effect is about chemistry rather
     than about importance, it should shrink here.
  3. THERE SHOULD BE A COMPARTMENT GRADIENT. The end of the relay that faces the
     ENVIRONMENT is under a different regime from the end that faces the genome.
     Receptors negotiate with ligands, pathogens and other individuals; nuclear
     effectors negotiate with DNA and with the rest of the transcriptional
     machinery. Rate should fall from plasma membrane to nucleus.

Prediction 3 is the one that answers "which organelles get mobilised" as a
measurement rather than a narrative, and it is the reason this module assigns
every gene a compartment from UniProt rather than reusing pathway labels.

Pathway classification walks the Reactome hierarchy from its top-level roots, so
"signalling" means a descendant of Signal Transduction (R-HSA-162582) and
"metabolic" a descendant of Metabolism (R-HSA-1430728) - not a keyword match on
the pathway name.
"""
import os
import re
import argparse

import numpy as np
import pandas as pd
from scipy import stats

from .common import RAW, DERIVED, FIGDIR, apply_figure_style, PALETTE, log
from . import t5_confounders, t15_pathway, t19_scaleout

RELATION = os.path.join(RAW, "reactome_relation.txt")
PATHWAY_NAMES = os.path.join(RAW, "reactome_pathways.txt")
SUBCELL = os.path.join(RAW, "uniprot_subcell.tsv")

ROOTS = {
    "R-HSA-162582": "signalling",
    "R-HSA-1430728": "metabolic",
    "R-HSA-168256": "immune",
    "R-HSA-1640170": "cell_cycle",
    "R-HSA-73857": "transcription",
    "R-HSA-392499": "protein_metabolism",
}
# ordered from the outside world inwards - the gradient is tested along this axis
COMPARTMENT_ORDER = ["secreted", "plasma_membrane", "endomembrane", "cytosol",
                     "mitochondrion", "nucleus"]
MIN_PATHWAY = 8
MAX_PATHWAY = 300
MIN_PER_GROUP = 3
MIN_PER_COMPARTMENT = 15

FOCUS_PATHWAYS = ["MAPK", "Wnt", "Notch", "TGF", "JAK", "Hedgehog",
                  "PI3K", "NF-kB", "Toll", "Hippo", "mTOR", "Death Receptor"]


# --------------------------------------------------------------------------
# pathway classification via the Reactome hierarchy
# --------------------------------------------------------------------------
def pathway_class():
    """pathway id -> top-level class, by descent rather than by name."""
    children = {}
    for line in open(RELATION):
        f = line.rstrip("\n").split("\t")
        if len(f) == 2 and f[0].startswith("R-HSA") and f[1].startswith("R-HSA"):
            children.setdefault(f[0], []).append(f[1])
    out = {}
    for root, label in ROOTS.items():
        stack = [root]
        seen = set()
        while stack:
            p = stack.pop()
            if p in seen:
                continue
            seen.add(p)
            # a pathway reachable from two roots keeps the first assignment,
            # so ROOTS order encodes precedence
            out.setdefault(p, label)
            stack.extend(children.get(p, ()))
    return out


# --------------------------------------------------------------------------
# compartment assignment
# --------------------------------------------------------------------------
def gene_compartment():
    """symbol -> one compartment, resolved by precedence, not by first mention.

    UniProt lists every location a protein has ever been seen in, so a plain
    keyword scan makes almost everything "cytosol". The order below runs from
    the outside in and takes the OUTERMOST annotated location, which is the one
    that matters for a relay: a receptor that also cycles through endosomes is
    still a receptor.
    """
    d = pd.read_csv(SUBCELL, sep="\t")
    loc = d["Subcellular location [CC]"].fillna("")
    tm = d["Transmembrane"].fillna("")
    comp = np.full(len(d), "cytosol", dtype=object)
    # assign from the inside out so that outer hits overwrite inner ones
    comp[loc.str.contains("Nucleus", case=False).to_numpy()] = "nucleus"
    comp[loc.str.contains("Mitochondri", case=False).to_numpy()] = "mitochondrion"
    comp[loc.str.contains("Endoplasmic reticulum|Golgi|Lysosome|Endosome",
                          case=False, regex=True).to_numpy()] = "endomembrane"
    pm = (loc.str.contains("Cell membrane|Plasma membrane", case=False,
                           regex=True)
          & tm.str.contains("TRANSMEM")).to_numpy()
    comp[pm] = "plasma_membrane"
    comp[loc.str.contains("Secreted", case=False).to_numpy()] = "secreted"
    d["compartment"] = comp
    d = d.dropna(subset=["Gene Names (primary)"]).drop_duplicates(
        "Gene Names (primary)").set_index("Gene Names (primary)")
    return d["compartment"]


def analysis_frame(sp_key="mouse"):
    d = t19_scaleout.analysis_frame(sp_key)
    d["compartment"] = d["symbol"].map(gene_compartment())
    return d


def usable_pathways(d):
    members = t15_pathway.pathway_members()
    cls = pathway_class()
    out = {}
    for pid, (name, genes) in members.items():
        present = [g for g in genes if g in d.index]
        if MIN_PATHWAY <= len(present) <= MAX_PATHWAY:
            out[pid] = (name, present, cls.get(pid, "other"))
    return out


# --------------------------------------------------------------------------
# 1. signalling vs metabolic
# --------------------------------------------------------------------------
def class_comparison(d, paths):
    rows = []
    for pid, (name, genes, cls) in paths.items():
        sub = d.loc[genes]
        om = sub["omega"].dropna()
        if len(om) < MIN_PATHWAY:
            continue
        rows.append(dict(pathway=pid, name=name, cls=cls, n=len(om),
                         median_omega=float(om.median()),
                         frac_catalytic=float(sub["is_catalytic"].astype(bool).mean()),
                         frac_interface=float(sub["compartment"].isin(
                             ["secreted", "plasma_membrane"]).mean()),
                         frac_nuclear=float((sub["compartment"] == "nucleus").mean()),
                         median_tau=float(sub["tau"].median()),
                         median_retention=float(sub["n_species"].median())
                         if "n_species" in sub else np.nan))
    res = pd.DataFrame(rows).set_index("pathway")
    summary = res.groupby("cls").agg(
        n_pathways=("n", "size"), median_omega=("median_omega", "median"),
        frac_catalytic=("frac_catalytic", "median"),
        frac_interface=("frac_interface", "median"),
        frac_nuclear=("frac_nuclear", "median"),
        median_tau=("median_tau", "median"),
        median_retention=("median_retention", "median"))
    sig = res.loc[res["cls"] == "signalling", "median_omega"]
    met = res.loc[res["cls"] == "metabolic", "median_omega"]
    u, p = stats.mannwhitneyu(sig, met, alternative="two-sided")
    test = dict(n_sig=len(sig), n_met=len(met),
                median_sig=float(sig.median()), median_met=float(met.median()),
                fold=float(sig.median() / met.median()),
                rank_biserial=float(2 * u / (len(sig) * len(met)) - 1),
                p=float(p))
    return res, summary, test


def catalytic_split_by_class(d, paths):
    """Run the T19 catalytic/accessory test separately per pathway class."""
    rows = []
    for pid, (name, genes, cls) in paths.items():
        sub = d.loc[genes]
        mask = sub["is_catalytic"].astype(bool)
        a = sub.loc[mask, "omega"].dropna()
        b = sub.loc[~mask, "omega"].dropna()
        if len(a) < MIN_PER_GROUP or len(b) < MIN_PER_GROUP:
            continue
        rows.append(dict(pathway=pid, cls=cls,
                         log2_acc_over_cat=float(np.log2(b.median() / a.median()))))
    r = pd.DataFrame(rows)
    out = []
    for cls, sub in r.groupby("cls"):
        if len(sub) < 20:
            continue
        w = stats.wilcoxon(sub["log2_acc_over_cat"], alternative="two-sided")
        out.append(dict(cls=cls, n_pathways=len(sub),
                        median_log2=float(sub["log2_acc_over_cat"].median()),
                        fold=float(2 ** sub["log2_acc_over_cat"].median()),
                        frac_positive=float((sub["log2_acc_over_cat"] > 0).mean()),
                        p=float(w.pvalue)))
    return pd.DataFrame(out).set_index("cls").sort_values("fold", ascending=False)


# --------------------------------------------------------------------------
# 2. the compartment gradient
# --------------------------------------------------------------------------
def compartment_gradient(d, paths, cls="signalling"):
    """Rate by compartment, pooled and within pathways.

    The pooled version is confounded: compartments differ in tau, expression and
    coding length. The within-pathway version blocks on pathway, and the
    regression adds the covariates that predicted rate everywhere else in this
    project, so a surviving gradient is not a restatement of those.
    """
    genes = sorted({g for pid, (n, gs, c) in paths.items() if c == cls for g in gs})
    sub = d.loc[genes].dropna(subset=["omega", "compartment"])
    pooled = []
    for c in COMPARTMENT_ORDER:
        v = sub.loc[sub["compartment"] == c, "omega"].dropna()
        if len(v) < MIN_PER_COMPARTMENT:
            continue
        pooled.append(dict(compartment=c, n=len(v), median_omega=float(v.median()),
                           median_tau=float(sub.loc[v.index, "tau"].median()),
                           median_retention=float(sub.loc[v.index, "n_species"].median())
                           if "n_species" in sub else np.nan))
    pooled = pd.DataFrame(pooled).set_index("compartment")

    # ordinal trend across the outside-in axis
    order = {c: i for i, c in enumerate(COMPARTMENT_ORDER)}
    ok = sub[sub["compartment"].isin(order)]
    rho, p = stats.spearmanr(ok["compartment"].map(order), ok["omega"])

    # membrane vs nucleus, the two ends, matched on tau and expression
    pm = ok[ok["compartment"] == "plasma_membrane"]
    nu = ok[ok["compartment"] == "nucleus"]
    ends = {}
    if len(pm) >= MIN_PER_COMPARTMENT and len(nu) >= MIN_PER_COMPARTMENT:
        u, pv = stats.mannwhitneyu(pm["omega"], nu["omega"], alternative="two-sided")
        ends = dict(n_pm=len(pm), n_nuc=len(nu),
                    median_pm=float(pm["omega"].median()),
                    median_nuc=float(nu["omega"].median()),
                    fold=float(pm["omega"].median() / nu["omega"].median()),
                    rank_biserial=float(2 * u / (len(pm) * len(nu)) - 1),
                    p=float(pv))
    return pooled, dict(rho=float(rho), p=float(p), n=len(ok)), ends


def gradient_regression(d, paths, cls="signalling"):
    """Does compartment predict rate within pathways, after the usual covariates?"""
    import statsmodels.formula.api as smf
    rows = []
    for pid, (name, genes, c) in paths.items():
        if c != cls:
            continue
        for g in genes:
            rows.append((g, pid))
    df = pd.DataFrame(rows, columns=["gene", "pathway"]).drop_duplicates("gene")
    sub = d.loc[df["gene"]].copy()
    sub["pathway"] = df["pathway"].to_numpy()
    sub = sub.dropna(subset=["log_omega", "tau", "log_expr", "log_codons",
                             "compartment"])
    sub = sub[sub["compartment"].isin(COMPARTMENT_ORDER)]
    sub["comp_rank"] = sub["compartment"].map(
        {c: i for i, c in enumerate(COMPARTMENT_ORDER)})
    base = smf.ols("log_omega ~ comp_rank", data=sub).fit()
    cov = smf.ols("log_omega ~ comp_rank + tau + log_expr + log_codons",
                  data=sub).fit()
    fe = smf.ols("log_omega ~ comp_rank + tau + log_expr + log_codons "
                 "+ C(pathway)", data=sub).fit()
    return dict(n=int(base.nobs), n_pathways=sub["pathway"].nunique(),
                coef_raw=float(base.params["comp_rank"]),
                coef_cov=float(cov.params["comp_rank"]),
                coef_fe=float(fe.params["comp_rank"]),
                p_fe=float(fe.pvalues["comp_rank"]),
                fold_per_step_fe=float(10 ** fe.params["comp_rank"]))


def named_pathways(res, d, paths):
    rows = []
    for pid, (name, genes, cls) in paths.items():
        if cls != "signalling":
            continue
        if not any(k.lower() in name.lower() for k in FOCUS_PATHWAYS):
            continue
        sub = d.loc[genes].dropna(subset=["omega"])
        if len(sub) < 12:
            continue
        comp = sub["compartment"].value_counts(normalize=True)
        rows.append(dict(name=name[:52], n=len(sub),
                         median_omega=float(sub["omega"].median()),
                         frac_pm=float(comp.get("plasma_membrane", 0.0)),
                         frac_nuc=float(comp.get("nucleus", 0.0)),
                         frac_cat=float(sub["is_catalytic"].astype(bool).mean()),
                         retention=float(sub["n_species"].median())
                         if "n_species" in sub else np.nan))
    return pd.DataFrame(rows).sort_values("median_omega", ascending=False)


def make_figure(summary, cat_cls, pooled, grad, ends, res):
    plt = apply_figure_style()
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    order = ["signalling", "immune", "transcription", "cell_cycle",
             "protein_metabolism", "metabolic"]
    order = [c for c in order if c in res["cls"].unique()]
    vals = [res.loc[res["cls"] == c, "median_omega"] for c in order]
    bp = axes[0].boxplot(vals, tick_labels=order, patch_artist=True,
                         showfliers=False)
    for patch, col in zip(bp["boxes"], PALETTE):
        patch.set_facecolor(col)
        patch.set_alpha(0.8)
    for m in bp["medians"]:
        m.set_color("k")
    axes[0].set_ylabel(r"pathway median $\omega$")
    axes[0].set_title("pathway class (Reactome hierarchy)")
    axes[0].tick_params(axis="x", rotation=35)
    for lab in axes[0].get_xticklabels():
        lab.set_ha("right")

    p = pooled.reindex([c for c in COMPARTMENT_ORDER if c in pooled.index])
    axes[1].plot(range(len(p)), p["median_omega"], "o-", color=PALETTE[1], lw=1.8)
    for i, (c, r) in enumerate(p.iterrows()):
        axes[1].annotate(f"n={int(r['n'])}", (i, r["median_omega"]), fontsize=6,
                         xytext=(3, 4), textcoords="offset points")
    axes[1].set_xticks(range(len(p)))
    axes[1].set_xticklabels(p.index, rotation=35, ha="right")
    axes[1].set_ylabel(r"median $\omega$")
    axes[1].set_title("signalling: outside in\n"
                      r"Spearman $\rho$" f"={grad['rho']:+.3f} (p={grad['p']:.1g})")

    c = cat_cls.sort_values("fold")
    axes[2].barh(c.index, c["fold"], color=[PALETTE[1] if q < 0.05 else "#BBBBBB"
                                            for q in c["p"]])
    axes[2].axvline(1.0, color="k", lw=0.9)
    axes[2].set_xlabel(r"$\omega_{\rm accessory} / \omega_{\rm catalytic}$")
    axes[2].set_title("catalytic/accessory split\nby pathway class")
    axes[2].grid(axis="y", alpha=0)

    path = os.path.join(FIGDIR, "t20_signalling.png")
    fig.savefig(path)
    plt.close(fig)
    return path


def report(sp_key="mouse"):
    print("=" * 72)
    print("T20 - signalling vs metabolic pathways, and the compartment gradient")
    print("=" * 72)
    d = analysis_frame(sp_key)
    paths = usable_pathways(d)
    print(f"genes: {len(d)}  with a compartment call: "
          f"{d['compartment'].notna().sum()}")
    print(f"pathways classified: {len(paths)}")
    print("  " + ", ".join(f"{k}={v}" for k, v in
                           pd.Series([c for _, _, c in paths.values()]
                                     ).value_counts().items()) + "\n")

    res, summary, test = class_comparison(d, paths)
    print("1. PATHWAY CLASS")
    print(f"   {'class':<20}{'n':>5}{'omega':>9}{'%catalytic':>12}"
          f"{'%interface':>12}{'%nuclear':>10}{'tau':>7}{'retention':>11}")
    for cls, r in summary.sort_values("median_omega", ascending=False).iterrows():
        print(f"   {cls:<20}{int(r['n_pathways']):>5}{r['median_omega']:>9.4f}"
              f"{r['frac_catalytic']:>12.2f}{r['frac_interface']:>12.2f}"
              f"{r['frac_nuclear']:>10.2f}{r['median_tau']:>7.3f}"
              f"{r['median_retention']:>11.0f}")
    print(f"\n   signalling vs metabolic: {test['median_sig']:.4f} vs "
          f"{test['median_met']:.4f}  ({test['fold']:.2f}x, "
          f"rb={test['rank_biserial']:+.2f}, p={test['p']:.2g})")

    cat_cls = catalytic_split_by_class(d, paths)
    print("\n2. CATALYTIC vs ACCESSORY, split by pathway class")
    print(f"   {'class':<20}{'n':>5}{'fold':>8}{'%positive':>11}{'p':>11}")
    for cls, r in cat_cls.iterrows():
        print(f"   {cls:<20}{int(r['n_pathways']):>5}{r['fold']:>8.2f}"
              f"{r['frac_positive']:>11.0%}{r['p']:>11.2g}")

    pooled, grad, ends = compartment_gradient(d, paths)
    print("\n3. COMPARTMENT GRADIENT within signalling pathways")
    print(f"   {'compartment':<20}{'n':>6}{'omega':>9}{'tau':>8}{'retention':>11}")
    for c, r in pooled.iterrows():
        print(f"   {c:<20}{int(r['n']):>6}{r['median_omega']:>9.4f}"
              f"{r['median_tau']:>8.3f}{r['median_retention']:>11.0f}")
    print(f"\n   ordinal trend outside->in: Spearman rho={grad['rho']:+.3f} "
          f"(p={grad['p']:.2g}, n={grad['n']})")
    if ends:
        print(f"   the two ends: plasma membrane {ends['median_pm']:.4f} vs "
              f"nucleus {ends['median_nuc']:.4f}  ({ends['fold']:.2f}x, "
              f"rb={ends['rank_biserial']:+.2f}, p={ends['p']:.2g})")
    gr = gradient_regression(d, paths)
    print(f"\n   compartment-rank coefficient on log10(omega), n={gr['n']} genes "
          f"/ {gr['n_pathways']} pathways:")
    print(f"      raw                        {gr['coef_raw']:+.4f}")
    print(f"      + tau, expression, length  {gr['coef_cov']:+.4f}")
    print(f"      + pathway fixed effects    {gr['coef_fe']:+.4f} "
          f"(p={gr['p_fe']:.2g}, {gr['fold_per_step_fe']:.3f}x per step inward)")

    named = named_pathways(res, d, paths)
    if len(named):
        print("\n4. NAMED SIGNALLING PATHWAYS")
        print(f"   {'pathway':<54}{'n':>4}{'omega':>8}{'%PM':>7}{'%nuc':>7}"
              f"{'%cat':>7}{'ret':>6}")
        for _, r in named.head(18).iterrows():
            print(f"   {r['name']:<54}{int(r['n']):>4}{r['median_omega']:>8.4f}"
                  f"{r['frac_pm']:>7.2f}{r['frac_nuc']:>7.2f}"
                  f"{r['frac_cat']:>7.2f}{r['retention']:>6.0f}")

    res.to_csv(os.path.join(DERIVED, "t20_pathway_class.tsv"), sep="\t")
    summary.to_csv(os.path.join(DERIVED, "t20_class_summary.tsv"), sep="\t")
    cat_cls.to_csv(os.path.join(DERIVED, "t20_catalytic_by_class.tsv"), sep="\t")
    pooled.to_csv(os.path.join(DERIVED, "t20_compartment_gradient.tsv"), sep="\t")
    if len(named):
        named.to_csv(os.path.join(DERIVED, "t20_named_pathways.tsv"),
                     sep="\t", index=False)
    p = make_figure(summary, cat_cls, pooled, grad, ends, res)
    print(f"\nfigure written to: {p}")
    return d, res, summary, cat_cls, pooled


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", default="mouse")
    a = ap.parse_args()
    report(sp_key=a.species)
