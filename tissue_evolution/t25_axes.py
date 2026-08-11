"""
T25 - the three-axis map, with decisiveness as the fourth channel.

The proposal is to place every candidate protein in a space of three weights:

    axis 1   tissue / cell division of labour
    axis 2   species / ecology
    axis 3   coevolution and compensation

and to stop treating expression difference as evidence of importance -
importance is defined by perturbation, C_i^J = dlnJ/dln a_i. So C^J is not a
fourth axis here, it is the size and colour of the point: position says what
kind of variation a protein carries, the marker says whether any of it moves
the flux.

HOW EACH AXIS IS MEASURED, and what each one is NOT.

axis 1, division of labour = tau (Yanai tissue-specificity index) over 32 GTEx
    organs, plus, for paralogue families, how strongly a gene partitions
    organs against its closest paralogue. tau alone would call a gene
    "generalist" even when it and its paralogue split the body cleanly between
    them, because each one is broadly expressed within its own half. The
    partition term is a Spearman anticorrelation of the two organ profiles, so
    SUCLA2 (muscle/brain, ATP-forming) against SUCLG2 (liver/kidney,
    GTP-forming) scores high even though neither is tau-extreme.
    NOT a measure of how much the tissue needs the gene.

axis 2, species / ecology = omega (dN/dS) against mouse, z-scored, combined
    with how much that rate shifts across the four divergence depths. A
    protein that is uniformly fast everywhere is not responding to ecology; a
    protein whose rate changes between lineages might be.
    NOT a test of positive selection - omega here is genome-wide-scaled, and
    T17 already showed that separating adaptation from relaxation needs
    polymorphism data, not divergence alone.

axis 3, coevolution and compensation = the T24 interface excess, i.e. how much
    more a protein's contacting residues prefer their own species' partner
    than equally-variable NON-contacting residues in the same two proteins do.
    Using the raw native-minus-mismatched score here would put every protein
    high, because relatives resemble relatives; the excess over the matched
    control is the only part of that which is about the interface.
    NOT available for proteins with no solved interface in the panel - those
    are plotted at the origin of this axis and flagged, not imputed.

marker, decisiveness = C_i^J from T23, the median over 150 sampled parameter
    sets. Only the ten reactions in the kinetic model have one; everything
    else is drawn hollow. An enzyme drawn hollow is not "zero control", it is
    "not modelled", and the two must not be confused.

The point of drawing it this way is that the axes are nearly independent of
the marker. If tissue-specificity or evolutionary rate predicted control, the
big markers would cluster; the map is worth making precisely because they do
not.
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

from .common import DERIVED, FIGDIR, apply_figure_style, PALETTE, log

# the four nominated groups, plus the TCA catalytic core they sit in
GROUPS = {
    "SUCL (ATP/GTP split)": ["SUCLA2", "SUCLG2", "SUCLG1"],
    "redox tolerance": ["SIRT3", "FXN", "NNT", "SOD2"],
    "complex II + quinone": ["SDHA", "SDHB", "SDHC", "SDHD",
                             "SDHAF1", "SDHAF2", "SDHAF3", "SDHAF4"],
    "mito/nuclear pairing": ["MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB",
                             "MT-ATP6", "MT-ND1", "MT-ND2",
                             "COX4I1", "COX5A", "COX5B", "COX6B1", "COX6C",
                             "COX7A2", "COX8A", "UQCRC1", "UQCRC2",
                             "UQCRFS1", "UQCRB", "UQCRQ", "CYC1",
                             "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS7",
                             "NDUFS8", "NDUFV1", "NDUFA9"],
    "TCA core": ["CS", "ACO2", "IDH2", "IDH3A", "IDH3B", "IDH3G", "OGDH",
                 "DLST", "DLD", "FH", "MDH2", "ACO1", "IDH1", "MDH1",
                 "GOT1", "GOT2", "OGDHL"],
}

# reaction in the T23 model -> the gene whose activity it stands for
MODEL_GENE = {"CS": "CS", "ACO2": "ACO2", "IDH3": "IDH3A", "OGDH": "OGDH",
              "SUCL": "SUCLA2", "SDH": "SDHA", "FH": "FH", "MDH2": "MDH2",
              "CI": "NDUFS1", "OX": "MT-CO1"}


def expression():
    f = pd.read_csv(os.path.join(DERIVED, "t1_gene_expression_features.tsv.gz"),
                    sep="\t")
    return f.dropna(subset=["symbol"]).drop_duplicates("symbol").set_index("symbol")


def organ_profiles():
    m = pd.read_csv(os.path.join(DERIVED, "t1_organ_median_tpm.tsv.gz"),
                    sep="\t", index_col=0)
    f = pd.read_csv(os.path.join(DERIVED, "t1_gene_expression_features.tsv.gz"),
                    sep="\t").dropna(subset=["symbol"])
    sym = f.drop_duplicates("Name").set_index("Name")["symbol"]
    m = m.join(sym, how="inner").dropna(subset=["symbol"])
    return m.groupby("symbol").mean(numeric_only=True)


def paralogue_partition(prof, families):
    """How cleanly do members of a family split the organs between them?

    Spearman anticorrelation of log organ profiles: -1 means a perfect split
    (wherever one is high the other is low), +1 means redundancy.
    """
    out = {}
    for fam in families:
        present = [g for g in fam if g in prof.index]
        if len(present) < 2:
            continue
        L = np.log1p(prof.loc[present])
        for g in present:
            others = [h for h in present if h != g]
            rhos = [stats.spearmanr(L.loc[g], L.loc[h]).statistic
                    for h in others]
            out[g] = float(-np.nanmin(rhos))   # strongest anticorrelation
    return out


def rates():
    """omega against mouse, and how much the rate moves across depths."""
    import glob
    frames = {}
    for f in glob.glob(os.path.join(DERIVED, "t3_dnds_human_*.tsv.gz")):
        sp = os.path.basename(f).replace("t3_dnds_human_", "").replace(".tsv.gz", "")
        d = pd.read_csv(f, sep="\t")
        d = d[d["n_codons"] >= 100]
        frames[sp] = d.set_index("human_gene")[["omega", "dN"]]
    feat = pd.read_csv(os.path.join(DERIVED, "t1_gene_expression_features.tsv.gz"),
                       sep="\t").dropna(subset=["symbol"])
    g2s = feat.drop_duplicates("Name").set_index("Name")["symbol"]

    out = {}
    for sp, d in frames.items():
        s = d.join(g2s, how="inner").dropna(subset=["symbol"])
        s = s.drop_duplicates("symbol").set_index("symbol")
        # rate relative to the genome median at that depth, so depths are
        # comparable despite wildly different absolute divergence.
        # dN == 0 (protein identical to the orthologue) would give -inf and
        # poison the mean and standard deviation of the whole axis, so those
        # are dropped rather than clipped - "zero substitutions" is a censored
        # observation, not a very small rate.
        rel = np.log2(s["dN"].where(s["dN"] > 0) / s["dN"][s["dN"] > 0].median())
        out[sp] = rel.replace([np.inf, -np.inf], np.nan)
    R = pd.DataFrame(out)
    om = frames["mouse"].join(g2s, how="inner").dropna(subset=["symbol"])
    om = om.drop_duplicates("symbol").set_index("symbol")["omega"]
    return pd.DataFrame(dict(omega=om, rel_rate=R.mean(axis=1),
                             rate_shift=R.std(axis=1), n_depths=R.notna().sum(axis=1)))


def coevolution():
    p = os.path.join(DERIVED, "t24_swap.tsv")
    if not os.path.exists(p):
        return pd.Series(dtype=float), pd.Series(dtype=float)
    d = pd.read_csv(p, sep="\t")
    d = d[d["excess"].notna()]
    if d.empty:
        return pd.Series(dtype=float), pd.Series(dtype=float)
    long = pd.concat([
        d[["gene_a", "excess", "n_informative"]].rename(columns={"gene_a": "g"}),
        d[["gene_b", "excess", "n_informative"]].rename(columns={"gene_b": "g"})])
    return (long.groupby("g")["excess"].mean(),
            long.groupby("g")["n_informative"].sum())


def control():
    p = os.path.join(DERIVED, "t23_control_summary.tsv")
    if not os.path.exists(p):
        return pd.Series(dtype=float)
    d = pd.read_csv(p, sep="\t", index_col=0)
    return pd.Series({MODEL_GENE[k]: v for k, v in d["median"].items()
                      if k in MODEL_GENE})


def build():
    genes = sorted({g for v in GROUPS.values() for g in v})
    grp = {g: k for k, v in GROUPS.items() for g in v}

    ex = expression()
    prof = organ_profiles()
    part = paralogue_partition(prof, [
        ["SUCLA2", "SUCLG2"], ["ACO1", "ACO2"], ["IDH1", "IDH2", "IDH3A"],
        ["MDH1", "MDH2"], ["GOT1", "GOT2"], ["OGDH", "OGDHL"],
    ])
    rt = rates()
    coev, coev_n = coevolution()
    cj = control()

    d = pd.DataFrame(index=genes)
    d["group"] = [grp[g] for g in genes]
    d["tau"] = ex["tau"].reindex(genes)
    d["max_tpm"] = ex["max_tpm"].reindex(genes)
    d["top_organ"] = ex["top_organ"].reindex(genes)
    d["partition"] = pd.Series(part).reindex(genes)
    d["omega"] = rt["omega"].reindex(genes)
    d["rel_rate"] = rt["rel_rate"].reindex(genes)
    d["rate_shift"] = rt["rate_shift"].reindex(genes)
    d["coev_excess"] = coev.reindex(genes)
    d["coev_n"] = coev_n.reindex(genes)
    d["C_J"] = cj.reindex(genes)

    def z(s):
        return (s - s.mean()) / s.std(ddof=0)

    # axis 1: specificity, raised by a clean paralogue split
    d["axis_tissue"] = z(d["tau"].fillna(d["tau"].median())) \
        + 0.5 * z(d["partition"].fillna(0.0))
    # axis 2: fast, and unevenly fast across lineages
    d["axis_species"] = z(d["rel_rate"].fillna(d["rel_rate"].median())) \
        + 0.5 * z(d["rate_shift"].fillna(d["rate_shift"].median()))
    # axis 3: interface excess only; no imputation
    d["axis_coev"] = d["coev_excess"]
    return d


def make_figure(d):
    plt = apply_figure_style()
    fig = plt.figure(figsize=(15, 4.6))
    groups = list(GROUPS)
    colour = {g: PALETTE[i % len(PALETTE)] for i, g in enumerate(groups)}

    def scatter(ax, xcol, ycol, xlab, ylab):
        for g in groups:
            s = d[d["group"] == g]
            if s.empty:
                continue
            has = s["C_J"].notna()
            ax.scatter(s.loc[~has, xcol], s.loc[~has, ycol], s=22,
                       facecolors="none", edgecolors=colour[g], lw=0.9,
                       label=g)
            if has.any():
                ax.scatter(s.loc[has, xcol], s.loc[has, ycol],
                           s=40 + 900 * s.loc[has, "C_J"].clip(0),
                           color=colour[g], alpha=0.85, edgecolors="k", lw=0.5)
        for k, r in d.iterrows():
            if pd.notna(r[xcol]) and pd.notna(r[ycol]) and \
                    (pd.notna(r["C_J"]) or abs(r[ycol]) > 0.4):
                ax.annotate(k, (r[xcol], r[ycol]), fontsize=5.5,
                            xytext=(3, 2), textcoords="offset points")
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)

    ax = fig.add_subplot(1, 3, 1)
    scatter(ax, "axis_tissue", "axis_species",
            "axis 1: tissue / cell division of labour",
            "axis 2: species / ecology")
    ax.set_title("filled = has a control coefficient\narea proportional to $C^J$")
    ax.legend(fontsize=5.5, loc="upper left")

    ax = fig.add_subplot(1, 3, 2)
    scatter(ax, "axis_tissue", "axis_coev",
            "axis 1: tissue / cell division of labour",
            "axis 3: coevolution excess over matched control")
    ax.axhline(0, color="k", lw=0.8, ls="--")
    ax.set_title("above the dashed line = interface carries\nmore than phylogeny alone")

    ax = fig.add_subplot(1, 3, 3)
    s = d.dropna(subset=["C_J"])
    x = np.arange(len(s))
    s = s.sort_values("C_J", ascending=False)
    ax.bar(x, s["C_J"], color=[colour[g] for g in s["group"]])
    ax.set_xticks(x)
    ax.set_xticklabels(s.index, rotation=60, ha="right", fontsize=6)
    ax.set_ylabel(r"$C_i^{J}$ (median over 150 parameter sets)")
    ax.set_title("decisiveness is concentrated;\nexpression and rate do not predict it")

    p = os.path.join(FIGDIR, "t25_three_axes.png")
    fig.savefig(p)
    plt.close(fig)
    return p


def report():
    print("=" * 72)
    print("T25 - three axes, with decisiveness defined by perturbation")
    print("=" * 72)
    d = build()

    print(f"\n{'gene':<9}{'group':<22}{'tau':>6}{'part':>6}{'relrate':>8}"
          f"{'shift':>7}{'coev':>8}{'C^J':>8}")
    for g, r in d.sort_values(["group", "gene" if "gene" in d else "axis_tissue"]).iterrows():
        cj = "  -  " if pd.isna(r["C_J"]) else f"{r['C_J']:.3f}"
        cv = "   -  " if pd.isna(r["coev_excess"]) else f"{r['coev_excess']:+.3f}"
        print(f"{g:<9}{r['group']:<22}{r['tau']:>6.2f}"
              f"{(r['partition'] if pd.notna(r['partition']) else np.nan):>6.2f}"
              f"{r['rel_rate']:>8.2f}{r['rate_shift']:>7.2f}{cv:>8}{cj:>8}")

    print("\ndoes any axis predict decisiveness?")
    s = d.dropna(subset=["C_J"])
    for col, lab in [("axis_tissue", "axis 1 tissue"),
                     ("axis_species", "axis 2 species"),
                     ("tau", "tau alone"), ("rel_rate", "rate alone"),
                     ("max_tpm", "peak expression")]:
        ss = s.dropna(subset=[col])
        if len(ss) < 5:
            continue
        rho = stats.spearmanr(ss[col], ss["C_J"])
        print(f"  Spearman({lab:<16}, C^J) = {rho.statistic:+.3f}  "
              f"p={rho.pvalue:.3f}   n={len(ss)}")
    print("  n is 10 - the number of reactions in the kinetic model - so these")
    print("  correlations are descriptive. They cannot establish independence;")
    print("  they can only fail to show the expected dependence.")

    d.to_csv(os.path.join(DERIVED, "t25_axes.tsv"), sep="\t")
    p = make_figure(d)
    print(f"\nfigure written to: {p}")
    return d


if __name__ == "__main__":
    report()
