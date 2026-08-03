"""
T6 - does the organ rate ranking survive going deeper in the tree, and where
does the method break?

Running the same analysis at four divergence depths (human vs macaque ~29 My,
mouse ~87 My, opossum ~159 My, chicken ~319 My) tests two separate things:

  1. STABILITY. If "testis is fast, brain is slow" is a real property of the
     tissues it should reproduce at every depth. If it only appears at one
     depth it is a property of that particular comparison.

  2. WHERE THE METHOD DIES. dS accumulates ~5-10x faster than dN, so synonymous
     sites saturate first: by ~150 My a large fraction of genes have pS >= 3/4
     and the Jukes-Cantor correction is undefined. Past that point omega is not
     merely noisy, it is *selectively* missing - the genes that drop out are the
     fast ones, so the surviving omega distribution is biased downward. This
     module measures the drop-out rate instead of ignoring it, and falls back to
     dN, which saturates far later, for the deep comparisons.

dN is additionally divided by the TimeTree divergence estimate to give a rate
in substitutions per non-synonymous site per 100 My. That conversion assumes a
constant rate along each lineage, which is false in detail (rodents are fast,
primates slow); it is used to compare ORGANS WITHIN a comparison, never to
compare absolute rates across comparisons.
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

from .common import DERIVED, FIGDIR, SPECIES, apply_figure_style, PALETTE, log
from . import t1_expression, t4_tissue_rates, t5_confounders  # noqa: F401

FOCUS = ["Testis", "Liver", "Skin", "Blood", "Skeletal_muscle",
         "Brain", "Cerebellum", "Pituitary", "Kidney", "Spleen"]


def saturation_table():
    rows = []
    for sp in SPECIES:
        f = os.path.join(DERIVED, f"t3_dnds_human_{sp['key']}.tsv.gz")
        if not os.path.exists(f):
            continue
        d = pd.read_csv(f, sep="\t")
        d = d[d["n_codons"] >= 100]
        rows.append(dict(
            species=sp["key"], mya=sp["divergence_mya"], n_pairs=len(d),
            median_identity=d["identity"].median(),
            median_dN=d["dN"].median(), median_dS=d["dS"].median(),
            frac_dS_over_1=float((d["dS"] > 1).mean()),
            frac_dS_undefined=float(d["dS"].isna().mean()),
            frac_omega_usable=float(d["omega"].notna().mean()),
            median_omega=d["omega"].median()))
    return pd.DataFrame(rows).set_index("species")


def organ_rates_all_depths():
    """Per-organ median omega and dN over tissue-enriched genes, every depth."""
    feat, mat = t1_expression.build()
    frames = {}
    for sp in SPECIES:
        f = os.path.join(DERIVED, f"t3_dnds_human_{sp['key']}.tsv.gz")
        if not os.path.exists(f):
            continue
        r = pd.read_csv(f, sep="\t")
        r = r[r["n_codons"] >= 100].set_index("human_gene")
        d = feat.join(r[["dN", "dS", "omega"]], how="inner")
        if "group_enriched_in" not in d.columns:
            d["group_enriched_in"] = np.nan
        rows = []
        # same organ gene sets as T5: the strict single-organ rule leaves only
        # a handful of organs above the size floor once deep comparisons have
        # thinned the orthologue set, and a 4-point rank correlation is not a
        # measurement of anything
        for organ in t5_confounders.organ_list(d):
            sub = d[t5_confounders.organ_membership(d, organ)]
            if len(sub) < 40:
                continue
            rows.append(dict(organ=organ, n=len(sub),
                             omega=sub["omega"].median(),
                             dN=sub["dN"].median(),
                             n_omega=int(sub["omega"].notna().sum()),
                             frac_omega=float(sub["omega"].notna().mean())))
        frames[sp["key"]] = pd.DataFrame(rows).set_index("organ")
    return frames


def rank_stability(frames, col="dN"):
    keys = list(frames)
    common = None
    for k in keys:
        idx = frames[k][col].dropna().index
        common = idx if common is None else common.intersection(idx)
    M = pd.DataFrame({k: frames[k].loc[common, col] for k in keys})
    rho = pd.DataFrame(index=keys, columns=keys, dtype=float)
    for a in keys:
        for b in keys:
            rho.loc[a, b] = stats.spearmanr(M[a], M[b]).statistic
    return M, rho


def intrinsic_across_depths():
    """Run the T5 matched test at every depth; keep the intrinsic component."""
    out = {}
    for sp in SPECIES:
        f = os.path.join(DERIVED, f"t3_dnds_human_{sp['key']}.tsv.gz")
        if not os.path.exists(f):
            continue
        try:
            d = t5_confounders.analysis_frame(sp["key"])
            out[sp["key"]] = t5_confounders.matched_test(d)
        except Exception as exc:                       # noqa: BLE001
            log(f"  matched test failed for {sp['key']}: {exc}")
    return out


def make_figure(sat, M_dn, frames):
    plt = apply_figure_style()
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4))

    x = sat["mya"].to_numpy()
    axes[0].plot(x, sat["median_dS"], "o-", color=PALETTE[0], label=r"median $d_S$")
    axes[0].plot(x, sat["median_dN"], "s-", color=PALETTE[1], label=r"median $d_N$")
    axes[0].plot(x, sat["frac_omega_usable"], "^--", color=PALETTE[2],
                 label=r"fraction with usable $\omega$")
    axes[0].axhline(1.0, color="k", lw=0.7, ls=":")
    axes[0].set_xlabel("divergence from human (Mya)")
    axes[0].set_title("saturation: where the method dies")
    axes[0].legend(fontsize=7)
    for xi, name in zip(x, sat.index):
        axes[0].annotate(name, (xi, sat.loc[name, "median_dS"]), fontsize=6.5,
                         xytext=(2, 4), textcoords="offset points")

    mya = {s["key"]: s["divergence_mya"] for s in SPECIES}
    for i, organ in enumerate([o for o in FOCUS if o in M_dn.index]):
        xs = [mya[k] for k in M_dn.columns]
        axes[1].plot(xs, M_dn.loc[organ], "o-", lw=1.2, ms=3.5,
                     color=PALETTE[i % len(PALETTE)], label=organ)
    axes[1].set_xlabel("divergence from human (Mya)")
    axes[1].set_ylabel(r"median $d_N$ of tissue-enriched genes")
    axes[1].set_title(r"organ $d_N$ vs depth")
    axes[1].legend(fontsize=6, ncol=2)

    keys = list(M_dn.columns)
    ref = keys[0]
    for i, organ in enumerate([o for o in FOCUS if o in M_dn.index]):
        rel = M_dn.loc[organ] / M_dn.median(axis=0)
        axes[2].plot([mya[k] for k in keys], rel, "o-", lw=1.2, ms=3.5,
                     color=PALETTE[i % len(PALETTE)], label=organ)
    axes[2].axhline(1.0, color="k", lw=0.8)
    axes[2].set_xlabel("divergence from human (Mya)")
    axes[2].set_ylabel(r"organ $d_N$ / median organ $d_N$")
    axes[2].set_title("relative ranking is what is stable")
    axes[2].legend(fontsize=6, ncol=2)

    p = os.path.join(FIGDIR, "t6_divergence_depth.png")
    fig.savefig(p)
    plt.close(fig)
    return p


def report():
    print("=" * 72)
    print("T6 - rank stability across divergence depth, and where dN/dS breaks")
    print("=" * 72)
    sat = saturation_table()
    print(f"{'species':<10}{'Mya':>7}{'pairs':>8}{'ident':>8}{'med dN':>9}"
          f"{'med dS':>9}{'dS>1':>8}{'omega ok':>10}")
    for k, r in sat.iterrows():
        print(f"{k:<10}{r['mya']:>7.0f}{int(r['n_pairs']):>8}{r['median_identity']:>8.3f}"
              f"{r['median_dN']:>9.4f}{r['median_dS']:>9.3f}"
              f"{r['frac_dS_over_1']:>7.0%}{r['frac_omega_usable']:>10.0%}")

    frames = organ_rates_all_depths()
    if len(frames) < 2:
        print("\n(only one depth computed - rank stability needs >= 2)")
        return sat, frames, None, None

    M_dn, rho_dn = rank_stability(frames, "dN")
    M_om, rho_om = rank_stability(frames, "omega")
    print(f"\norgan-ranking Spearman rho between depths ({M_dn.shape[0]} organs)")
    print("  using dN:")
    print(rho_dn.round(3).to_string().replace("\n", "\n    "))
    print("  using omega:")
    print(rho_om.round(3).to_string().replace("\n", "\n    "))

    n_show = max(1, min(3, len(M_dn) // 2))
    print(f"\nfastest / slowest organs by median dN at each depth "
          f"({len(M_dn)} organs shared by all depths):")
    for k in M_dn.columns:
        s = M_dn[k].sort_values(ascending=False)
        print(f"  {k:<9} fastest: {', '.join(s.index[:n_show])}   "
              f"slowest: {', '.join(s.index[-n_show:][::-1])}")

    sat.to_csv(os.path.join(DERIVED, "t6_saturation.tsv"), sep="\t")
    M_dn.to_csv(os.path.join(DERIVED, "t6_organ_dN_by_depth.tsv"), sep="\t")
    M_om.to_csv(os.path.join(DERIVED, "t6_organ_omega_by_depth.tsv"), sep="\t")
    p = make_figure(sat, M_dn, frames)
    print(f"\nfigure written to: {p}")
    return sat, frames, M_dn, rho_dn


if __name__ == "__main__":
    report()
