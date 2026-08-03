"""
T4 - per-organ evolutionary rate, three different ways, plus orthologue coverage.

"How fast does a tissue evolve" is not one number, and the three reasonable
definitions do not agree. This module computes all three and keeps them apart:

  (a) ENRICHED-SET rate - median omega over the genes that are tissue-enriched
      in that organ. This is what most papers report. It measures the organ's
      *identity genes*, and it is heavily confounded by tau (an enriched gene is
      a high-tau gene by construction).

  (b) EXPRESSION-WEIGHTED rate - sum_g TPM_g * omega_g / sum_g TPM_g over every
      expressed gene. This measures the average constraint on the proteome the
      organ actually deploys. It is dominated by a handful of very highly
      expressed genes (ALB in liver, MYH7 in heart), so a log-weighted variant
      is reported alongside it.

  (c) BREADTH rate - unweighted median omega over all genes above 1 TPM in the
      organ. Every expressed gene counts once.

ORTHOLOGUE COVERAGE (the fourth output) is not a rate at all but answers a
question the rates cannot: what fraction of an organ's NUCLEAR PROTEIN-CODING
expression mass sits on genes with NO 1:1 orthologue in the compared species?
Those genes are exactly the lineage-specific / duplicated / fast-diverging ones
that a dN/dS analysis must discard. An organ with low coverage has a
transcriptome that is partly evolutionarily novel, and its rate estimates are
correspondingly biased DOWN, because the fastest-changing part of its proteome
is unmeasurable by this method. Reporting rates without coverage would hide that
bias.

The mitochondrial genome is excluded from that denominator and reported on its
own line - see `t1_expression.coding_universe` for why leaving it in silently
turns the coverage column into a ranking of oxidative capacity.
"""
import os

import numpy as np
import pandas as pd

from .common import DERIVED, FIGDIR, SPECIES, apply_figure_style, PALETTE, log
from . import t1_expression

MIN_TPM = 1.0
PRIMARY = "mouse"          # best signal-to-saturation trade-off (see T6)


def load_rates(sp_key):
    path = os.path.join(DERIVED, f"t3_dnds_human_{sp_key}.tsv.gz")
    if not os.path.exists(path):
        raise FileNotFoundError(f"{path} - run tissue_evolution.t3_dnds first")
    d = pd.read_csv(path, sep="\t")
    d = d[d["n_codons"] >= 100]
    return d.set_index("human_gene")


def merged_table(sp_key):
    """Gene-level table: expression features + evolutionary rates, inner join."""
    feat, mat = t1_expression.build()
    rates = load_rates(sp_key)
    df = feat.join(rates[["dN", "dS", "omega", "identity", "n_codons"]], how="inner")
    # Two expression covariates, because they are NOT interchangeable here.
    # mean TPM across organs is close to a restatement of breadth, so it is
    # ~collinear with tau (Spearman about -0.74) and "controlling for expression"
    # with it silently controls for tau twice over. PEAK TPM is the level the
    # gene actually reaches where it is used, and is far less entangled with
    # tau, so it is the default covariate for matching and regression.
    df["log_expr_mean"] = np.log10(df["mean_tpm"] + 1.0)
    df["log_expr"] = np.log10(df["max_tpm"] + 1.0)
    df["log_codons"] = np.log10(df["n_codons"])
    return df, mat


def orthologue_coverage(mat, rates_index):
    """Per-organ fraction of expression mass (and of expressed genes) that sits
    on genes with a usable 1:1 orthologue.

    The denominator is NUCLEAR PROTEIN-CODING genes only. Computed over all
    expressed genes instead, this statistic is dominated by the 13 mitochondrial
    protein genes plus MT-RNR1/2, which alone carry 28-65% of a tissue's TPM and
    can never have a nuclear-code orthologue score - it would rank organs by
    mitochondrial content while appearing to rank them by evolutionary novelty.
    `frac_expression_mitochondrial` is reported separately so that signal is
    visible rather than smuggled into the coverage number.
    """
    nuclear, mito = t1_expression.coding_universe()
    is_nuc = mat.index.isin(nuclear)
    is_mt = mat.index.isin(mito)
    has = mat.index.isin(rates_index)
    rows = []
    for organ in mat.columns:
        tpm = mat[organ].to_numpy(dtype=float)
        expressed = tpm >= MIN_TPM
        base = expressed & is_nuc
        mass_total = tpm[base].sum()
        mass_cov = tpm[base & has].sum()
        all_mass = tpm[expressed].sum()
        rows.append(dict(
            organ=organ,
            n_expressed_coding=int(base.sum()),
            frac_genes_with_ortholog=float((base & has).sum() / max(base.sum(), 1)),
            frac_expression_with_ortholog=float(mass_cov / mass_total) if mass_total else np.nan,
            frac_expression_mitochondrial=float(tpm[expressed & is_mt].sum() / all_mass)
            if all_mass else np.nan,
        ))
    return pd.DataFrame(rows).set_index("organ")


def per_organ_rates(df, mat, min_genes=25):
    """Compute estimators (a), (b), (c) for every organ."""
    common = mat.index.intersection(df.index)
    sub_mat = mat.loc[common]
    om = df.loc[common, "omega"]
    dn = df.loc[common, "dN"]
    ok = om.notna().to_numpy()

    rows = []
    for organ in mat.columns:
        tpm = sub_mat[organ].to_numpy(dtype=float)

        # (a) enriched-set
        enr = df[(df["enriched_in"] == organ) & df["omega"].notna()]
        a_med = enr["omega"].median() if len(enr) >= min_genes else np.nan
        a_dn = enr["dN"].median() if len(enr) >= min_genes else np.nan

        # (b) expression-weighted over genes with a usable omega
        w = tpm.copy()
        w[~ok] = 0.0
        b_lin = float(np.average(om.to_numpy()[ok], weights=w[ok])) if w[ok].sum() > 0 else np.nan
        wl = np.log2(tpm + 1.0)
        wl[~ok] = 0.0
        b_log = float(np.average(om.to_numpy()[ok], weights=wl[ok])) if wl[ok].sum() > 0 else np.nan

        # (c) breadth
        sel = ok & (tpm >= MIN_TPM)
        c_med = float(np.median(om.to_numpy()[sel])) if sel.sum() >= min_genes else np.nan

        rows.append(dict(organ=organ,
                         n_enriched=int(len(enr)),
                         omega_enriched=a_med, dN_enriched=a_dn,
                         omega_expr_weighted=b_lin, omega_logexpr_weighted=b_log,
                         omega_breadth=c_med, n_breadth=int(sel.sum()),
                         median_tau_enriched=enr["tau"].median() if len(enr) else np.nan))
    return pd.DataFrame(rows).set_index("organ")


def make_figure(res, cov, sp_key):
    plt = apply_figure_style()
    fig, axes = plt.subplots(1, 3, figsize=(14, 5.0))

    d = res.dropna(subset=["omega_enriched"]).sort_values("omega_enriched")
    axes[0].barh(d.index, d["omega_enriched"], color=PALETTE[0])
    axes[0].set_xlabel(r"median $\omega$ (tissue-enriched genes)")
    axes[0].set_title(f"(a) identity genes, human-{sp_key}")
    axes[0].grid(axis="y", alpha=0)

    d2 = res.dropna(subset=["omega_breadth"]).sort_values("omega_breadth")
    axes[1].barh(d2.index, d2["omega_breadth"], color=PALETTE[1])
    axes[1].set_xlabel(r"median $\omega$ (all genes $\geq$1 TPM)")
    axes[1].set_title("(c) whole expressed transcriptome")
    axes[1].grid(axis="y", alpha=0)

    d3 = cov.sort_values("frac_expression_with_ortholog")
    axes[2].barh(d3.index, d3["frac_expression_with_ortholog"], color=PALETTE[2])
    axes[2].set_xlabel("fraction of expression mass\nwith a 1:1 orthologue")
    axes[2].set_title("orthologue coverage\n(low = novel transcriptome)")
    axes[2].grid(axis="y", alpha=0)

    fig.suptitle(f"Per-organ evolutionary rate, three estimators (human vs {sp_key})",
                 y=1.01)
    p = os.path.join(FIGDIR, f"t4_organ_rates_{sp_key}.png")
    fig.savefig(p)
    plt.close(fig)
    return p


def build(sp_key=PRIMARY):
    df, mat = merged_table(sp_key)
    res = per_organ_rates(df, mat)
    cov = orthologue_coverage(mat, set(df.index))
    res = res.join(cov)
    res.to_csv(os.path.join(DERIVED, f"t4_organ_rates_{sp_key}.tsv"), sep="\t")
    df.to_csv(os.path.join(DERIVED, f"t4_gene_table_{sp_key}.tsv.gz"), sep="\t")
    return df, mat, res


def report(sp_key=PRIMARY):
    print("=" * 72)
    print(f"T4 - per-organ evolutionary rate, three estimators (human vs {sp_key})")
    print("=" * 72)
    df, mat, res = build(sp_key)
    print(f"genes with expression AND a usable 1:1 orthologue: {len(df)}")
    print(f"of those, usable omega (dS unsaturated)          : {df['omega'].notna().sum()}")
    print(f"genome-wide median omega                          : {df['omega'].median():.4f}\n")

    show = res.dropna(subset=["omega_enriched"]).sort_values("omega_enriched",
                                                             ascending=False)
    print(f"{'organ':<17}{'n_enr':>6}{'w_enr':>8}{'tau_enr':>9}"
          f"{'w_breadth':>11}{'w_exprwt':>10}{'cov_expr':>10}{'mt_frac':>9}")
    for organ, r in show.iterrows():
        print(f"{organ:<17}{int(r['n_enriched']):>6}{r['omega_enriched']:>8.4f}"
              f"{r['median_tau_enriched']:>9.3f}{r['omega_breadth']:>11.4f}"
              f"{r['omega_expr_weighted']:>10.4f}"
              f"{r['frac_expression_with_ortholog']:>10.3f}"
              f"{r['frac_expression_mitochondrial']:>9.3f}")

    print("\nspread across organs:")
    for col, label in [("omega_enriched", "(a) enriched-set"),
                       ("omega_breadth", "(c) breadth"),
                       ("omega_expr_weighted", "(b) expression-weighted")]:
        v = res[col].dropna()
        print(f"  {label:<26} min {v.min():.4f}  max {v.max():.4f}  "
              f"max/min {v.max() / v.min():.2f}x")

    p = make_figure(res, res[["frac_expression_with_ortholog"]].dropna(), sp_key)
    print(f"\nfigure written to: {p}")
    return df, mat, res


if __name__ == "__main__":
    report()
