"""
T11 - is an organ's rate carried by a few dominant proteins, or by the bulk?

T4/T5 report medians, which deliberately hide this. "Testis is fast" could mean
either of two completely different biological situations:

  * a CONCENTRATED effect - a handful of extremely fast proteins (say, sperm
    surface and seminal proteins under sperm competition) dragging a summary
    statistic, while the other 95% of the organ's genes are ordinary; or
  * a DIFFUSE effect - the entire tissue-specific gene set shifted, which would
    mean something about the tissue as a whole rather than about a pathway.

They imply different mechanisms and different follow-up experiments, so this
module measures which one it is, three ways:

  1. LORENZ / GINI on each gene's contribution to the organ's excess rate. How
     unequal is the contribution?
  2. TRIMMED RE-RANKING. Recompute every organ's rate after deleting its top
     5 / 10 / 25% fastest genes. If the ranking survives, the effect is diffuse.
  3. NAMED GENES. The actual fastest and slowest proteins per organ, so the
     result can be checked against known biology instead of taken on trust.

Contribution of gene g in organ O is defined as log2(omega_g / median omega of
the genome), i.e. how much that gene pulls the organ away from the genomic
background. Genes below the background contribute negatively, which is what
makes the Lorenz curve informative rather than trivially monotone.
"""
import os

import numpy as np
import pandas as pd

from .common import DERIVED, FIGDIR, apply_figure_style, PALETTE, log
from . import t4_tissue_rates, t5_confounders

TRIM_LEVELS = [0.05, 0.10, 0.25]
MIN_GENES = 40


def frame(sp_key=t4_tissue_rates.PRIMARY):
    d = t5_confounders.analysis_frame(sp_key)
    d = d[d["omega"].notna()].copy()
    d["log2_vs_genome"] = np.log2(d["omega"] / d["omega"].median())
    return d


def gini(x):
    """Gini coefficient of |contribution| - a pure concentration measure."""
    v = np.sort(np.abs(np.asarray(x, dtype=float)))
    n = len(v)
    if n == 0 or v.sum() == 0:
        return np.nan
    return float((2.0 * np.arange(1, n + 1) - n - 1).dot(v) / (n * v.sum()))


def concentration(d):
    """Per organ: Gini, and the share of the total excess carried by the top decile."""
    rows = []
    for organ in t5_confounders.organ_list(d):
        sub = d[t5_confounders.organ_membership(d, organ)]
        if len(sub) < MIN_GENES:
            continue
        c = sub["log2_vs_genome"].to_numpy()
        order = np.argsort(-np.abs(c))
        k = max(1, int(round(0.10 * len(c))))
        top_share = np.abs(c[order[:k]]).sum() / np.abs(c).sum()
        rows.append(dict(organ=organ, n=len(sub),
                         median_log2=float(np.median(c)),
                         gini=gini(c),
                         top10pct_share=float(top_share),
                         frac_above_genome=float((c > 0).mean())))
    return pd.DataFrame(rows).set_index("organ")


def trimmed_ranking(d):
    """Recompute each organ's median omega after deleting its fastest genes."""
    rows = []
    for organ in t5_confounders.organ_list(d):
        sub = d[t5_confounders.organ_membership(d, organ)]
        if len(sub) < MIN_GENES:
            continue
        w = np.sort(sub["omega"].to_numpy())
        rec = dict(organ=organ, n=len(sub), omega_full=float(np.median(w)))
        for t in TRIM_LEVELS:
            keep = w[:max(1, int(round(len(w) * (1 - t))))]
            rec[f"omega_trim{int(t * 100)}"] = float(np.median(keep))
        rows.append(rec)
    res = pd.DataFrame(rows).set_index("organ")
    for t in TRIM_LEVELS:
        res[f"rank_trim{int(t * 100)}"] = res[f"omega_trim{int(t * 100)}"].rank(ascending=False)
    res["rank_full"] = res["omega_full"].rank(ascending=False)
    return res.sort_values("omega_full", ascending=False)


def named_genes(d, organs, n=12):
    """Fastest and slowest named proteins per organ."""
    out = {}
    for organ in organs:
        sub = d[t5_confounders.organ_membership(d, organ)].dropna(subset=["omega"])
        if len(sub) < MIN_GENES:
            continue
        s = sub.sort_values("omega", ascending=False)
        out[organ] = dict(
            fast=s[["symbol", "omega", "dN", "tau"]].head(n),
            slow=s[["symbol", "omega", "dN", "tau"]].tail(n).iloc[::-1])
    return out


def make_figure(conc, trim, d):
    plt = apply_figure_style()
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4))

    # Lorenz curves for a few contrasting organs
    for i, organ in enumerate(["Testis", "Spleen", "Brain", "Cerebellum",
                               "Skeletal_muscle"]):
        sub = d[t5_confounders.organ_membership(d, organ)]
        if len(sub) < MIN_GENES:
            continue
        v = np.sort(np.abs(sub["log2_vs_genome"].to_numpy()))
        cum = np.concatenate([[0], np.cumsum(v) / v.sum()])
        x = np.linspace(0, 1, len(cum))
        axes[0].plot(x, cum, lw=1.5, color=PALETTE[i % len(PALETTE)],
                     label=f"{organ} (G={gini(sub['log2_vs_genome']):.2f})")
    axes[0].plot([0, 1], [0, 1], "k--", lw=0.8, label="perfectly diffuse")
    axes[0].set_xlabel("fraction of genes (ranked by |contribution|)")
    axes[0].set_ylabel("cumulative share of the organ's excess")
    axes[0].set_title("is the effect concentrated?")
    axes[0].legend(fontsize=6.5, loc="upper left")

    t = trim.dropna()
    axes[1].scatter(t["omega_full"], t["omega_trim10"], s=30, color=PALETTE[0])
    lim = [min(t["omega_trim10"].min(), t["omega_full"].min()) * 0.9,
           t["omega_full"].max() * 1.05]
    axes[1].plot(lim, lim, "k--", lw=0.8)
    for organ in t.index:
        axes[1].annotate(organ, (t.loc[organ, "omega_full"],
                                 t.loc[organ, "omega_trim10"]),
                         fontsize=6, xytext=(3, 2), textcoords="offset points")
    axes[1].set_xlabel(r"median $\omega$, all genes")
    axes[1].set_ylabel(r"median $\omega$ after deleting the fastest 10%")
    axes[1].set_title("ranking survives trimming?")

    c = conc.sort_values("top10pct_share")
    axes[2].barh(c.index, c["top10pct_share"], color=PALETTE[1])
    axes[2].axvline(0.10, color="k", ls="--", lw=0.9,
                    label="0.10 = perfectly diffuse")
    axes[2].set_xlabel("share of the organ's total excess\ncarried by its top 10% of genes")
    axes[2].set_title("concentration by organ")
    axes[2].grid(axis="y", alpha=0)
    axes[2].legend(fontsize=7)

    p = os.path.join(FIGDIR, "t11_drivers.png")
    fig.savefig(p)
    plt.close(fig)
    return p


def report(sp_key=t4_tissue_rates.PRIMARY, n_named=12):
    print("=" * 72)
    print(f"T11 - concentrated or diffuse? which proteins drive it (human-{sp_key})")
    print("=" * 72)
    d = frame(sp_key)
    conc = concentration(d)
    trim = trimmed_ranking(d)

    print("1. concentration of each organ's excess rate")
    print("   (top10% share = 0.10 would mean perfectly diffuse; 1.0 = one gene does it all)")
    print(f"   {'organ':<18}{'n':>5}{'gini':>8}{'top10%':>9}{'>genome':>9}")
    for organ, r in conc.sort_values("top10pct_share", ascending=False).iterrows():
        print(f"   {organ:<18}{int(r['n']):>5}{r['gini']:>8.3f}"
              f"{r['top10pct_share']:>9.3f}{r['frac_above_genome']:>9.1%}")

    print("\n2. does the ranking survive deleting the fastest genes?")
    print(f"   {'organ':<18}{'full':>9}{'-5%':>9}{'-10%':>9}{'-25%':>9}"
          f"{'rank full':>11}{'rank -25%':>11}")
    for organ, r in trim.iterrows():
        print(f"   {organ:<18}{r['omega_full']:>9.4f}{r['omega_trim5']:>9.4f}"
              f"{r['omega_trim10']:>9.4f}{r['omega_trim25']:>9.4f}"
              f"{int(r['rank_full']):>11}{int(r['rank_trim25']):>11}")
    rho = trim["rank_full"].corr(trim["rank_trim25"], method="spearman")
    print(f"   Spearman(rank full, rank after deleting fastest 25%) = {rho:+.3f}")

    focus = ["Testis", "Spleen", "Blood", "Lung", "Brain", "Cerebellum",
             "Skeletal_muscle", "Liver"]
    named = named_genes(d, focus, n=n_named)
    print("\n3. the actual proteins")
    for organ, tables in named.items():
        print(f"\n   --- {organ} ---")
        fast = ", ".join(f"{r.symbol}({r.omega:.2f})"
                         for r in tables["fast"].itertuples() if isinstance(r.symbol, str))
        slow = ", ".join(f"{r.symbol}({r.omega:.3f})"
                         for r in tables["slow"].itertuples() if isinstance(r.symbol, str))
        print(f"   fastest: {fast}")
        print(f"   slowest: {slow}")

    conc.to_csv(os.path.join(DERIVED, f"t11_concentration_{sp_key}.tsv"), sep="\t")
    trim.to_csv(os.path.join(DERIVED, f"t11_trimmed_{sp_key}.tsv"), sep="\t")
    allnamed = []
    for organ, tables in named.items():
        for kind in ("fast", "slow"):
            t = tables[kind].copy()
            t["organ"], t["end"] = organ, kind
            allnamed.append(t)
    if allnamed:
        pd.concat(allnamed).to_csv(
            os.path.join(DERIVED, f"t11_named_genes_{sp_key}.tsv"), sep="\t")
    p = make_figure(conc, trim, d)
    print(f"\nfigure written to: {p}")
    return d, conc, trim, named


if __name__ == "__main__":
    report()
