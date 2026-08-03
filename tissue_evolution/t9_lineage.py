"""
T9 - lineage-specific rate asymmetry: did an organ speed up on ONE lineage?

T4-T6 compare human against other species, so every number they produce is a
property of a PAIR. They cannot say whether a difference happened on the human
side or the other side, and "does the primate brain evolve faster than the
rodent brain" is a question about sides.

The relative-rate test fixes that with an outgroup. For a gene with 1:1
orthologues in human, mouse and chicken:

    dN(human, chicken)  = human branch + shared path to the amniote ancestor
    dN(mouse,  chicken) = mouse branch + the SAME shared path

so their difference cancels the shared path exactly:

    delta = dN(human, chicken) - dN(mouse, chicken)
          = (rate on the human lineage) - (rate on the mouse lineage)
            since the human/mouse split

Two things must be subtracted before any of this means anything about organs:

  1. A GENOME-WIDE OFFSET. Rodents have shorter generation times and a higher
     per-year substitution rate, so delta is negative for most genes no matter
     what tissue they serve. The organ-level statistic is therefore always
     delta MINUS the genome-wide median, never delta itself.
  2. THE SAME COMPOSITION CONFOUND AS T5. Organ-enriched sets differ in tau and
     expression, both of which affect dN. The matched-control version of the
     test is reported next to the raw one.

What a positive result does and does not mean: dN going up on one lineage is
consistent with positive selection AND with relaxed constraint (a trait being
lost accelerates its genes just as adaptation does). Distinguishing them needs
a codon model that separates the two - RELAX, or a branch-site test - which is
outside what pairwise counting can do. This module measures the asymmetry and
says which genes to take to that test; it does not claim adaptation.
"""
import os
import csv
import argparse
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy import stats

from .common import (RAW, DERIVED, FIGDIR, apply_figure_style, PALETTE,
                     biomart_query, log)
from . import t1_expression, t3_dnds, t2_orthologs, t5_confounders

MOUSE_CHROMS = [str(i) for i in range(1, 20)] + ["X", "Y"]
MIN_GENES = 30


def fetch_mouse_chicken_orthologs():
    out = os.path.join(RAW, "orthologs_mouse_chicken.tsv")
    if os.path.exists(out):
        return out
    attrs = ["ensembl_gene_id", "ggallus_homolog_ensembl_gene",
             "ggallus_homolog_orthology_type", "ggallus_homolog_perc_id",
             "ggallus_homolog_orthology_confidence"]
    part_dir = os.path.join(RAW, "biomart_parts")
    os.makedirs(part_dir, exist_ok=True)
    rows = []
    for chrom in MOUSE_CHROMS:
        part = os.path.join(part_dir, f"mouse_chicken_chr{chrom}.tsv")
        if os.path.exists(part) and os.path.getsize(part) > 0:
            txt = open(part).read()
        else:
            txt = biomart_query("mmusculus_gene_ensembl", attrs,
                                filters=[("chromosome_name", chrom)])
            open(part, "w").write(txt)
        for line in txt.rstrip("\n").split("\n")[1:]:
            f = line.split("\t")
            if len(f) < 5 or not f[1] or f[2] != "ortholog_one2one":
                continue
            rows.append(f[:5])
        log(f"  mouse-chicken: chr{chrom} -> {len(rows)} cumulative 1:1 rows")
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["mouse_gene", "chicken_gene", "orthology_type",
                    "perc_id", "confidence"])
        w.writerows(rows)
    return out


def compute_mouse_chicken_dnds(processes=4):
    out = os.path.join(DERIVED, "t9_dnds_mouse_chicken.tsv.gz")
    if os.path.exists(out):
        return out
    ortho = fetch_mouse_chicken_orthologs()
    mm = t2_orthologs.index_longest_cds(t2_orthologs.fetch_cds("mus_musculus"))
    gg = t2_orthologs.index_longest_cds(t2_orthologs.fetch_cds("gallus_gallus"))
    log(f"  CDS indexed mouse={len(mm)} chicken={len(gg)}")

    pairs = []
    with open(ortho) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            a, b = mm.get(r["mouse_gene"]), gg.get(r["chicken_gene"])
            if a and b:
                pairs.append((r["mouse_gene"], r["chicken_gene"],
                              r["perc_id"], r["confidence"], a, b))
    log(f"  mouse-chicken pairs with both CDS: {len(pairs)}")
    with Pool(processes) as pool:
        rows = [r for r in pool.imap_unordered(t3_dnds.pair_stats, pairs,
                                               chunksize=64) if r is not None]
    df = pd.DataFrame(rows).rename(columns={"human_gene": "mouse_gene",
                                            "ortholog_gene": "chicken_gene"})
    df.to_csv(out, sep="\t", index=False)
    log(f"  mouse-chicken: {len(df)} aligned, median dN={df['dN'].median():.4f}")
    return out


def build_triplets():
    """Genes with human-chicken, mouse-chicken AND human-mouse 1:1 orthology."""
    hc = pd.read_csv(os.path.join(DERIVED, "t3_dnds_human_chicken.tsv.gz"),
                     sep="\t")
    hm = pd.read_csv(os.path.join(DERIVED, "t3_dnds_human_mouse.tsv.gz"),
                     sep="\t")
    mc = pd.read_csv(compute_mouse_chicken_dnds(), sep="\t")
    for d in (hc, hm, mc):
        d.drop(d.index[d["n_codons"] < 100], inplace=True)

    # human gene -> mouse gene, from the human-mouse 1:1 table
    link = hm.set_index("human_gene")["ortholog_gene"]
    hc = hc.set_index("human_gene")
    mc = mc.set_index("mouse_gene")

    common_h = hc.index.intersection(link.index)
    mouse_of = link.loc[common_h]
    keep = mouse_of.isin(mc.index)
    common_h = common_h[keep.to_numpy()]
    mouse_of = mouse_of.loc[common_h]

    # require the two outgroup comparisons to hit the SAME chicken gene,
    # otherwise the shared path does not cancel
    same = hc.loc[common_h, "ortholog_gene"].to_numpy() == \
        mc.loc[mouse_of, "chicken_gene"].to_numpy()
    common_h = common_h[same]
    mouse_of = mouse_of.loc[common_h]

    df = pd.DataFrame({
        "mouse_gene": mouse_of.to_numpy(),
        "dN_human_chicken": hc.loc[common_h, "dN"].to_numpy(),
        "dN_mouse_chicken": mc.loc[mouse_of, "dN"].to_numpy(),
        "dS_human_chicken": hc.loc[common_h, "dS"].to_numpy(),
        "dS_mouse_chicken": mc.loc[mouse_of, "dS"].to_numpy(),
    }, index=common_h)
    df.index.name = "human_gene"
    df["delta_dN"] = df["dN_human_chicken"] - df["dN_mouse_chicken"]
    return df.dropna(subset=["delta_dN"])


def organ_asymmetry(matched=True, seed=0):
    feat, _ = t1_expression.build()
    tri = build_triplets()
    d = feat.join(tri, how="inner")
    baseline = d["delta_dN"].median()
    d["delta_centred"] = d["delta_dN"] - baseline

    # same covariates and same organ gene sets as T5, so the two modules'
    # organ rows mean the same thing: PEAK expression (not mean, which is
    # near-collinear with tau) and tissue-enriched + group-enriched membership
    d["log_expr"] = np.log10(d["max_tpm"] + 1.0)
    d["log_expr_mean"] = np.log10(d["mean_tpm"] + 1.0)
    if "group_enriched_in" not in d.columns:
        d["group_enriched_in"] = np.nan
    d = d[d["tau"].notna() & d["log_expr"].notna()]

    rows = []
    for organ in t5_confounders.organ_list(d):
        sub = d[t5_confounders.organ_membership(d, organ)]
        if len(sub) < MIN_GENES:
            continue
        w = stats.wilcoxon(sub["delta_centred"], alternative="two-sided")
        rec = dict(organ=organ, n=len(sub),
                   median_delta=sub["delta_dN"].median(),
                   median_delta_centred=sub["delta_centred"].median(),
                   frac_human_faster=float((sub["delta_dN"] > 0).mean()),
                   p_raw=float(w.pvalue))
        if matched:
            m = t5_confounders.matched_control(
                d.assign(omega=d["delta_dN"]), organ, seed=seed)
            if m is not None:
                case, ctrl, n_drop = m
                u, p = stats.mannwhitneyu(case["delta_dN"], ctrl["delta_dN"],
                                          alternative="two-sided")
                rec.update(n_matched=len(case), n_dropped=n_drop,
                           median_delta_matched_ctrl=ctrl["delta_dN"].median(),
                           matched_effect=case["delta_dN"].median()
                           - ctrl["delta_dN"].median(),
                           p_matched=float(p))
        rows.append(rec)
    res = pd.DataFrame(rows).set_index("organ")
    res["q_raw"] = t5_confounders._bh(res["p_raw"].to_numpy())
    if "p_matched" in res:
        ok = res["p_matched"].notna()
        q = np.full(len(res), np.nan)
        q[ok.to_numpy()] = t5_confounders._bh(res.loc[ok, "p_matched"].to_numpy())
        res["q_matched"] = q
    return res.sort_values("median_delta_centred"), baseline, d


def make_figure(res, baseline, d):
    plt = apply_figure_style()
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.0))

    axes[0].hist(d["delta_dN"], bins=80, color=PALETTE[7], alpha=0.85)
    axes[0].axvline(0, color="k", lw=1.0, label="equal rates")
    axes[0].axvline(baseline, color=PALETTE[1], lw=1.4,
                    label=f"genome-wide median = {baseline:+.4f}")
    axes[0].set_xlim(np.nanpercentile(d["delta_dN"], [1, 99]))
    axes[0].set_xlabel(r"$\Delta d_N = d_N(\mathrm{human,chicken}) - d_N(\mathrm{mouse,chicken})$")
    axes[0].set_ylabel("genes")
    axes[0].set_title("mouse lineage is faster genome-wide\n(generation-time effect)")
    axes[0].legend(fontsize=7)

    s = res.sort_values("median_delta_centred")
    colours = [PALETTE[1] if q < 0.05 else "#BBBBBB" for q in s["q_raw"]]
    axes[1].barh(s.index, s["median_delta_centred"], color=colours)
    axes[1].axvline(0, color="k", lw=0.9)
    axes[1].set_xlabel(r"organ median $\Delta d_N$ - genome-wide median")
    axes[1].set_title("organ-specific lineage asymmetry\n"
                      "(right = faster on the primate lineage)")
    axes[1].grid(axis="y", alpha=0)

    p = os.path.join(FIGDIR, "t9_lineage_asymmetry.png")
    fig.savefig(p)
    plt.close(fig)
    return p


def report(processes=4):
    print("=" * 72)
    print("T9 - lineage-specific rate asymmetry (human vs mouse, chicken outgroup)")
    print("=" * 72)
    res, baseline, d = organ_asymmetry()
    print(f"three-way orthologue triplets used : {len(d)}")
    print(f"genome-wide median delta_dN        : {baseline:+.5f}")
    print(f"genes faster on the human lineage  : {(d['delta_dN'] > 0).mean():.1%}")
    print("  (a negative genome-wide median is the expected generation-time")
    print("   effect: the mouse lineage fixes more substitutions per year)\n")

    has_matched = "matched_effect" in res.columns
    hdr = f"{'organ':<18}{'n':>5}{'d_centred':>11}{'%h_fast':>9}{'q_raw':>10}"
    if has_matched:
        hdr += f"{'matched':>10}{'q_match':>10}"
    print(hdr)
    for organ, r in res.iterrows():
        line = (f"{organ:<18}{int(r['n']):>5}{r['median_delta_centred']:>+11.5f}"
                f"{r['frac_human_faster']:>8.0%}{r['q_raw']:>10.2e}")
        if has_matched:
            me = r.get("matched_effect", np.nan)
            qm = r.get("q_matched", np.nan)
            line += (f"{me:>+10.5f}" if np.isfinite(me) else f"{'-':>10}")
            line += (f"{qm:>10.2e}" if np.isfinite(qm) else f"{'-':>10}")
        print(line)

    res.to_csv(os.path.join(DERIVED, "t9_lineage_asymmetry.tsv"), sep="\t")
    p = make_figure(res, baseline, d)
    print(f"\nfigure written to: {p}")
    return res, baseline, d


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--processes", type=int, default=4)
    a = ap.parse_args()
    report(processes=a.processes)
