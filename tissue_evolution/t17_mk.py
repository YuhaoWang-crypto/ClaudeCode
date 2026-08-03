"""
T17 - McDonald-Kreitman: is testis under POSITIVE SELECTION or RELAXED CONSTRAINT?

Every rate result so far has had the same irreducible ambiguity. A high dN/dS is
equally consistent with adaptation (new variants driven to fixation) and with
relaxed constraint (deleterious variants no longer removed). Pairwise counting
cannot tell them apart, and I have said so at every step. The MK test can,
because it adds a second, independent axis: POLYMORPHISM within humans.

    Dn, Ds  fixed differences between species  (T3 counts, human vs macaque)
    Pn, Ps  polymorphisms segregating in humans (gnomAD v2.1.1 constraint)

Under strict neutrality Dn/Ds = Pn/Ps. The two mechanisms separate cleanly:

  * POSITIVE SELECTION drives amino-acid changes to fixation fast, so they show
    up in D but barely in P:   Dn/Ds > Pn/Ps,  alpha = 1 - (Ps*Dn)/(Pn*Ds) > 0.
  * RELAXED CONSTRAINT lets mildly deleterious amino-acid changes both segregate
    AND fix, so BOTH ratios rise together and alpha stays near 0.

That is exactly the discrimination the testis question needs, and it is why this
module exists rather than another dN/dS variant.

Divergence uses human-MACAQUE, not human-mouse: MK assumes the fixed differences
are single hits, and at dS=0.085 (T6) that is nearly true, while human-mouse
dS=0.625 is not.

THE BIAS THAT MATTERS, stated before the results
------------------------------------------------
gnomAD's obs_mis / obs_syn counts include every variant down to singletons.
Slightly deleterious amino-acid variants segregate at low frequency and are
over-represented there, which inflates Pn and biases alpha DOWNWARD - the
standard fix is a minor-allele-frequency cutoff (Fay-Wyckoff-Wu), and the
per-gene constraint file does not carry frequency-stratified counts. So absolute
alpha here is a LOWER BOUND and negative values are expected genome-wide. The
comparison BETWEEN organs is the interpretable part, and even that assumes the
bias is similar across gene classes, which is not guaranteed because the
distribution of fitness effects differs between them.

DoS = Dn/(Dn+Ds) - Pn/(Pn+Ps) (Stoletzki & Eyre-Walker 2011) is reported
alongside alpha because it is far better behaved when counts are small; alpha is
a ratio of ratios and explodes on small Ps or Pn.
"""
import os
import gzip
import argparse

import numpy as np
import pandas as pd
from scipy import stats

from .common import (RAW, DERIVED, FIGDIR, apply_figure_style, PALETTE,
                     download, log)
from . import t4_tissue_rates, t5_confounders

GNOMAD = os.path.join(RAW, "gnomad_constraint.txt.bgz")
GNOMAD_URL = ("https://storage.googleapis.com/gcp-public-data--gnomad/release/"
              "2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz")
DIVERGENCE = "macaque"          # unsaturated; see module docstring
MIN_ENRICHED = 40
MIN_COUNTS = 5                  # per-gene floor on each cell of the 2x2 table


def load_gnomad():
    if not os.path.exists(GNOMAD):
        log("downloading gnomAD v2.1.1 constraint metrics")
        download(GNOMAD_URL, GNOMAD, timeout=1800)
    with gzip.open(GNOMAD, "rt") as fh:
        d = pd.read_csv(fh, sep="\t", low_memory=False)
    d = d[d["gene_id"].notna()]
    d["gene_id"] = d["gene_id"].str.split(".").str[0]
    # one canonical transcript per gene: keep the row with the most synonymous
    # sites, i.e. the longest usable CDS
    d = d.sort_values("possible_syn", ascending=False).drop_duplicates("gene_id")
    return d.set_index("gene_id")[["obs_mis", "obs_syn", "possible_mis",
                                   "possible_syn", "exp_mis", "exp_syn",
                                   "oe_mis", "oe_lof", "pLI"]]


def analysis_frame(sp_key=DIVERGENCE):
    d = t5_confounders.analysis_frame(sp_key)
    div = pd.read_csv(os.path.join(DERIVED, f"t3_dnds_human_{sp_key}.tsv.gz"),
                      sep="\t").set_index("human_gene")
    if "Nd" not in div.columns:
        raise RuntimeError("t3 output lacks Nd/Sd counts - rerun t3_dnds")
    d = d.join(div[["Nd", "Sd", "N_sites", "S_sites"]], how="inner")
    d = d.join(load_gnomad(), how="inner")
    d = d.rename(columns={"obs_mis": "Pn", "obs_syn": "Ps"})
    d = d.dropna(subset=["Nd", "Sd", "Pn", "Ps"])
    return d


# --------------------------------------------------------------------------
# statistics
# --------------------------------------------------------------------------
def mk_stats(Dn, Ds, Pn, Ps):
    """alpha, neutrality index and DoS from 2x2 counts (arrays or scalars)."""
    Dn, Ds, Pn, Ps = (np.asarray(x, dtype=float) for x in (Dn, Ds, Pn, Ps))
    with np.errstate(divide="ignore", invalid="ignore"):
        # NI = (Pn/Ps) / (Dn/Ds); alpha = 1 - NI. Writing alpha with the
        # numerator and denominator swapped gives 1 - 1/NI, which has the right
        # sign only when NI < 1 and silently reverses the organ ordering.
        ni = (Pn * Ds) / (Ps * Dn)
        alpha = 1.0 - ni
        dos = Dn / (Dn + Ds) - Pn / (Pn + Ps)
    return alpha, ni, dos


def pooled_mk(sub):
    """Organ-level MK on POOLED counts.

    Pooling rather than averaging per-gene alpha is deliberate: alpha is a ratio
    of ratios and its per-gene distribution has no finite mean, so a mean of
    per-gene alphas is not an estimator of anything. Pooling counts across the
    gene set is the standard organ-level estimator.
    """
    Dn, Ds = sub["Nd"].sum(), sub["Sd"].sum()
    Pn, Ps = sub["Pn"].sum(), sub["Ps"].sum()
    alpha, ni, dos = mk_stats(Dn, Ds, Pn, Ps)
    table = [[int(round(Dn)), int(round(Ds))], [int(Pn), int(Ps)]]
    odds, p = stats.fisher_exact(table)
    return dict(n_genes=len(sub), Dn=float(Dn), Ds=float(Ds),
                Pn=float(Pn), Ps=float(Ps),
                DnDs=float(Dn / Ds), PnPs=float(Pn / Ps),
                alpha=float(alpha), NI=float(ni), DoS=float(dos),
                fisher_p=float(p))


def bootstrap_alpha(sub, n_boot=2000, seed=0):
    """Gene-level bootstrap CI for the pooled alpha."""
    rng = np.random.default_rng(seed)
    idx = np.arange(len(sub))
    Dn, Ds = sub["Nd"].to_numpy(), sub["Sd"].to_numpy()
    Pn, Ps = sub["Pn"].to_numpy(), sub["Ps"].to_numpy()
    out = np.empty(n_boot)
    for b in range(n_boot):
        s = rng.choice(idx, size=len(idx), replace=True)
        a, _, _ = mk_stats(Dn[s].sum(), Ds[s].sum(), Pn[s].sum(), Ps[s].sum())
        out[b] = a
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


def organ_mk(d, n_boot=2000):
    genome = pooled_mk(d)
    rows = []
    for organ in t5_confounders.organ_list(d):
        sub = d[t5_confounders.organ_membership(d, organ)]
        if len(sub) < MIN_ENRICHED:
            continue
        rec = pooled_mk(sub)
        rec["organ"] = organ
        lo, hi = bootstrap_alpha(sub, n_boot=n_boot)
        rec["alpha_lo"], rec["alpha_hi"] = lo, hi
        # DoS per gene is well behaved enough for a rank test against the rest
        rest = d[~t5_confounders.organ_membership(d, organ)]
        ok_s = sub[(sub[["Nd", "Sd", "Pn", "Ps"]] >= MIN_COUNTS).all(axis=1)]
        ok_r = rest[(rest[["Nd", "Sd", "Pn", "Ps"]] >= MIN_COUNTS).all(axis=1)]
        if len(ok_s) >= 20 and len(ok_r) >= 20:
            _, _, dos_s = mk_stats(ok_s["Nd"], ok_s["Sd"], ok_s["Pn"], ok_s["Ps"])
            _, _, dos_r = mk_stats(ok_r["Nd"], ok_r["Sd"], ok_r["Pn"], ok_r["Ps"])
            u, p = stats.mannwhitneyu(dos_s[np.isfinite(dos_s)],
                                      dos_r[np.isfinite(dos_r)],
                                      alternative="two-sided")
            rec["n_gene_dos"] = int(np.isfinite(dos_s).sum())
            rec["median_DoS_gene"] = float(np.nanmedian(dos_s))
            rec["median_DoS_rest"] = float(np.nanmedian(dos_r))
            rec["dos_p"] = float(p)
        rows.append(rec)
    res = pd.DataFrame(rows).set_index("organ")
    if "dos_p" in res:
        ok = res["dos_p"].notna()
        q = np.full(len(res), np.nan)
        q[ok.to_numpy()] = t5_confounders._bh(res.loc[ok, "dos_p"].to_numpy())
        res["dos_q"] = q
    return res.sort_values("alpha", ascending=False), genome


def make_figure(res, genome, d):
    plt = apply_figure_style()
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8))

    s = res.sort_values("alpha")
    y = np.arange(len(s))
    axes[0].barh(y, s["alpha"], color=[PALETTE[1] if a > genome["alpha"]
                                       else PALETTE[0] for a in s["alpha"]])
    axes[0].errorbar(s["alpha"], y,
                     xerr=[s["alpha"] - s["alpha_lo"], s["alpha_hi"] - s["alpha"]],
                     fmt="none", ecolor="k", elinewidth=0.8, capsize=2)
    axes[0].axvline(genome["alpha"], color="k", ls="--", lw=1.0,
                    label=f"genome {genome['alpha']:+.2f}")
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(s.index, fontsize=7)
    axes[0].set_xlabel(r"pooled $\alpha$  (95% gene bootstrap)")
    axes[0].set_title(r"$\alpha>$ genome = more adaptive fixation")
    axes[0].grid(axis="y", alpha=0)
    axes[0].legend(fontsize=7)

    axes[1].scatter(res["PnPs"], res["DnDs"], s=36, color=PALETTE[0])
    lim = [0, max(res["PnPs"].max(), res["DnDs"].max()) * 1.1]
    axes[1].plot(lim, lim, "k--", lw=0.9, label="neutral (Dn/Ds = Pn/Ps)")
    for organ in res.index:
        axes[1].annotate(organ, (res.loc[organ, "PnPs"], res.loc[organ, "DnDs"]),
                         fontsize=6, xytext=(3, 2), textcoords="offset points")
    axes[1].set_xlabel(r"$P_N/P_S$  (polymorphism, gnomAD)")
    axes[1].set_ylabel(r"$D_N/D_S$  (divergence, human-macaque)")
    axes[1].set_title("relaxed constraint moves along the line;\n"
                      "positive selection moves above it")
    axes[1].legend(fontsize=7)

    s2 = res.dropna(subset=["median_DoS_gene"]).sort_values("median_DoS_gene")
    colours = [PALETTE[1] if q < 0.05 else "#BBBBBB" for q in s2["dos_q"]]
    axes[2].barh(s2.index, s2["median_DoS_gene"], color=colours)
    axes[2].axvline(float(np.nanmedian(s2["median_DoS_rest"])), color="k",
                    ls="--", lw=0.9, label="genome")
    axes[2].set_xlabel("median per-gene DoS")
    axes[2].set_title("DoS: better behaved than "
                      r"$\alpha$" "\nat small counts (orange q<0.05)")
    axes[2].grid(axis="y", alpha=0)
    axes[2].legend(fontsize=7)

    p = os.path.join(FIGDIR, "t17_mk.png")
    fig.savefig(p)
    plt.close(fig)
    return p


def report(n_boot=2000):
    print("=" * 72)
    print("T17 - McDonald-Kreitman: positive selection vs relaxed constraint")
    print("=" * 72)
    d = analysis_frame()
    print(f"genes with divergence (human-{DIVERGENCE}) AND gnomAD polymorphism: "
          f"{len(d)}")

    res, genome = organ_mk(d, n_boot=n_boot)
    print(f"\ngenome-wide pooled: Dn={genome['Dn']:.0f} Ds={genome['Ds']:.0f} "
          f"Pn={genome['Pn']:.0f} Ps={genome['Ps']:.0f}")
    print(f"   Dn/Ds={genome['DnDs']:.4f}  Pn/Ps={genome['PnPs']:.4f}  "
          f"alpha={genome['alpha']:+.3f}  NI={genome['NI']:.3f}  "
          f"DoS={genome['DoS']:+.4f}")
    print("   a negative genome-wide alpha is the expected consequence of")
    print("   including singletons - see the module docstring. Read the organ")
    print("   ordering, not the absolute value.\n")

    print(f"{'organ':<18}{'n':>5}{'Dn/Ds':>8}{'Pn/Ps':>8}{'alpha':>8}"
          f"{'95% CI':>18}{'DoS':>9}{'q':>10}")
    for organ, r in res.iterrows():
        ci = f"[{r['alpha_lo']:+.2f},{r['alpha_hi']:+.2f}]"
        q = r.get("dos_q", np.nan)
        qs = f"{q:>10.2e}" if np.isfinite(q) else f"{'-':>10}"
        star = "*" if np.isfinite(q) and q < 0.05 else " "
        print(f"{organ:<18}{int(r['n_genes']):>5}{r['DnDs']:>8.4f}{r['PnPs']:>8.4f}"
              f"{r['alpha']:>+8.3f}{ci:>18}{r['DoS']:>+9.4f}{qs}{star}")

    print("\nreading the two axes together:")
    g = res.assign(dn_rel=res["DnDs"] / genome["DnDs"],
                   pn_rel=res["PnPs"] / genome["PnPs"])
    print(f"{'organ':<18}{'Dn/Ds vs genome':>17}{'Pn/Ps vs genome':>17}  interpretation")
    for organ, r in g.sort_values("dn_rel", ascending=False).iterrows():
        if r["dn_rel"] > 1.15 and r["pn_rel"] > 1.15:
            tag = "both up -> relaxed constraint"
        elif r["dn_rel"] > 1.15 and r["pn_rel"] <= 1.15:
            tag = "divergence only -> adaptive"
        elif r["dn_rel"] < 0.85 and r["pn_rel"] < 0.85:
            tag = "both down -> stronger purifying"
        else:
            tag = ""
        if tag:
            print(f"{organ:<18}{r['dn_rel']:>17.2f}{r['pn_rel']:>17.2f}  {tag}")

    res.to_csv(os.path.join(DERIVED, "t17_organ_mk.tsv"), sep="\t")
    pd.DataFrame([genome]).to_csv(os.path.join(DERIVED, "t17_genome_mk.tsv"),
                                  sep="\t", index=False)
    p = make_figure(res, genome, d)
    print(f"\nfigure written to: {p}")
    return d, res, genome


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-boot", type=int, default=2000)
    a = ap.parse_args()
    report(n_boot=a.n_boot)
