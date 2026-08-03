"""
T10 - is a protein's TISSUE DISTRIBUTION conserved, and is that coupled to how
fast its sequence evolves?

T1-T9 treat expression as a fixed human-derived covariate. But an orthologue's
organ profile is itself an evolving trait: a gene can keep its sequence and move
organs, or hold its organ and change its sequence. Those are different
evolutionary events and only the second one shows up in dN/dS.

Data: Cardoso-Moreira et al. 2019 (Nature) developmental organ transcriptomes,
via Expression Atlas - human E-MTAB-6814 and mouse E-MTAB-6798, the same seven
organs assayed the same way in both species. Each organ is pooled over its FULL
developmental series in both species; see the note by HUMAN_ADULT for why an
adult-only comparison is not the primary analysis.

Statistic: for each 1:1 orthologue, the correlation between its human 7-organ
profile and its mouse 7-organ profile. Correlating a gene's profile WITH ITSELF
ACROSS SPECIES is deliberately chosen over comparing absolute levels: it is
invariant to per-gene scaling, to library-size normalisation and to the
between-species batch effect that makes raw cross-species TPM comparisons
unreliable. Only the SHAPE of the profile is compared.

Limits, in order of how much they matter:

  * The coupling between profile conservation and omega DOES NOT REPLICATE in
    the adult-only subset (rho = -0.25 pooled vs -0.01 adult-only). Isolating
    the variable shows this is the stage pooling, not the organ count: at the
    same four organs, pooled stages still give rho = -0.19. The most likely
    reason is noise attenuation - adult-only rests on 2-3 human and 1 mouse
    assay group per organ, against 11-21 and 11-14 when stages are pooled - but
    it cannot be separated here from the alternative that what is conserved is
    the DEVELOPMENTAL trajectory rather than adult tissue identity. Treat the
    pooled correlation as real but not yet attributable to tissue identity.
  * "Peak organ moved between species" is fragile whenever two of the seven
    organs are near-duplicates: forebrain/hindbrain swaps dominate that table
    and are annotation noise, not relocation events.
  * Seven organs is a short vector, so Spearman per gene is coarse; Pearson on
    log2(TPM+1) is reported alongside it.
  * A gene expressed in one organ in both species scores high almost
    automatically, so results are stratified by tissue class rather than pooled.
"""
import os
import re
import argparse

import numpy as np
import pandas as pd
from scipy import stats

from .common import (RAW, DERIVED, FIGDIR, apply_figure_style, PALETTE,
                     download, log, to_tsv_gz)
from . import t1_expression

ATLAS_FTP = "https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments"
HUMAN_ACC, MOUSE_ACC = "E-MTAB-6814", "E-MTAB-6798"

# the seven organs assayed in both species; Atlas labels differ slightly
ORGANS = ["forebrain", "hindbrain", "heart", "kidney", "liver", "ovary", "testis"]
HUMAN_ADULT = ["young adult", "middle adult", "elderly"]
MOUSE_ADULT = ["postnatal day 63"]
# The two species are NOT sampled at comparable post-natal depth: mouse has all
# seven organs at P63, human adult samples cover only forebrain/heart/liver/
# testis (no hindbrain, kidney or ovary). Stage-matching human to mouse would
# need the published developmental correspondence, so the primary analysis
# instead pools each organ over its FULL developmental series in both species -
# a symmetric choice that keeps all seven organs. The adult-only 4-organ
# comparison is run as a sensitivity check, not discarded.


def _fetch(acc, suffix):
    dest = os.path.join(RAW, f"{acc}-{suffix}")
    if not os.path.exists(dest) or os.path.getsize(dest) == 0:
        download(f"{ATLAS_FTP}/{acc}/{acc}-{suffix}", dest, timeout=1800)
    return dest


def group_labels(acc):
    """assay-group id -> (stage, organ) from the Atlas configuration XML."""
    path = _fetch(acc, "configuration.xml")
    txt = open(path, encoding="latin-1").read()
    out = {}
    for gid, label in re.findall(r'<assay_group id="(g\d+)" label="([^"]+)"', txt):
        if ";" not in label:
            continue
        stage, organ = [s.strip() for s in label.split(";", 1)]
        out[gid] = (stage, organ)
    return out


def _parse_quartile_cell(series):
    """Expression Atlas stores each assay group as a comma-separated five-number
    summary (min,Q1,median,Q3,max), not a single value. Take the median element.
    A group with one replicate is written as a bare number."""
    parts = series.fillna("").astype(str).str.split(",")
    med = parts.apply(lambda p: p[len(p) // 2] if p and p[0] != "" else "")
    return pd.to_numeric(med, errors="coerce").fillna(0.0)


def organ_matrix(acc, adult_stages=None):
    """gene -> median TPM per organ. `adult_stages=None` pools all stages."""
    path = _fetch(acc, "tpms.tsv")
    labels = group_labels(acc)
    header = pd.read_csv(path, sep="\t", comment="#", nrows=0)
    gene_col = header.columns[0]
    keep = {}
    for gid in header.columns:
        if gid not in labels:
            continue
        stage, organ = labels[gid]
        if (adult_stages is None or stage in adult_stages) and organ in ORGANS:
            keep.setdefault(organ, []).append(gid)
    wanted = [gene_col] + [g for cols in keep.values() for g in cols]
    df = pd.read_csv(path, sep="\t", comment="#", usecols=wanted,
                     dtype=str, low_memory=False).set_index(gene_col)
    df = df[~df.index.duplicated(keep="first")]
    parsed = pd.DataFrame({c: _parse_quartile_cell(df[c]) for c in df.columns},
                          index=df.index)
    mat = pd.DataFrame({organ: parsed[cols].median(axis=1)
                        for organ, cols in keep.items()})
    mat = mat.reindex(columns=[o for o in ORGANS if o in mat.columns])
    mat = mat[mat.notna().all(axis=1)]
    log(f"  {acc}: {mat.shape[0]} genes x {mat.shape[1]} organs "
        f"({', '.join(mat.columns)})")
    return mat


def profile_conservation(sp_key="mouse", adult_only=False):
    hs = organ_matrix(HUMAN_ACC, HUMAN_ADULT if adult_only else None)
    mm = organ_matrix(MOUSE_ACC, MOUSE_ADULT if adult_only else None)
    organs = [o for o in ORGANS if o in hs.columns and o in mm.columns]
    hs, mm = hs[organs], mm[organs]

    rates = pd.read_csv(os.path.join(DERIVED, f"t3_dnds_human_{sp_key}.tsv.gz"),
                        sep="\t")
    rates = rates[rates["n_codons"] >= 100].set_index("human_gene")

    idx = rates.index.intersection(hs.index)
    pairs = rates.loc[idx, "ortholog_gene"]
    ok = pairs.isin(mm.index)
    idx, pairs = idx[ok.to_numpy()], pairs[ok.to_numpy()]

    H = np.log2(hs.loc[idx].to_numpy(dtype=float) + 1.0)
    M = np.log2(mm.loc[pairs].to_numpy(dtype=float) + 1.0)

    # drop genes that are flat in either species: a constant vector has no
    # profile to conserve and its correlation is undefined
    live = (H.std(axis=1) > 1e-6) & (M.std(axis=1) > 1e-6)
    idx, H, M = idx[live], H[live], M[live]

    Hc = H - H.mean(axis=1, keepdims=True)
    Mc = M - M.mean(axis=1, keepdims=True)
    pear = (Hc * Mc).sum(axis=1) / (np.linalg.norm(Hc, axis=1)
                                    * np.linalg.norm(Mc, axis=1))
    spear = np.array([stats.spearmanr(h, m).statistic for h, m in zip(H, M)])

    out = pd.DataFrame({"profile_pearson": pear, "profile_spearman": spear},
                       index=idx)
    out["mouse_gene"] = pairs.loc[idx].to_numpy()
    out = out.join(rates[["dN", "dS", "omega", "identity"]])
    feat, _ = t1_expression.build()
    out = out.join(feat[["symbol", "tau", "mean_tpm", "enriched_in"]])
    out["log_expr"] = np.log10(out["mean_tpm"] + 1.0)
    # the organ in which the gene peaks, in each species independently
    out["top_organ_human"] = np.array(organs)[H.argmax(axis=1)]
    out["top_organ_mouse"] = np.array(organs)[M.argmax(axis=1)]
    out["top_organ_conserved"] = out["top_organ_human"] == out["top_organ_mouse"]
    return out, organs


def make_figure(d, organs):
    plt = apply_figure_style()
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.3))

    axes[0].hist(d["profile_pearson"], bins=60, color=PALETTE[0], alpha=0.85)
    axes[0].axvline(0, color="k", lw=0.9)
    axes[0].axvline(d["profile_pearson"].median(), color=PALETTE[1], lw=1.4,
                    label=f"median = {d['profile_pearson'].median():.3f}")
    axes[0].set_xlabel("human-mouse organ-profile correlation")
    axes[0].set_ylabel("orthologues")
    axes[0].set_title(f"is the tissue profile conserved?\n({len(organs)} organs, adults)")
    axes[0].legend(fontsize=7)

    ok = d.dropna(subset=["omega", "profile_pearson"])
    q = pd.qcut(ok["profile_pearson"], 6, duplicates="drop")
    grp = ok.groupby(q, observed=True)["omega"].median()
    axes[1].plot(range(len(grp)), grp.values, "o-", color=PALETTE[1])
    axes[1].set_xticks(range(len(grp)))
    axes[1].set_xticklabels([f"{iv.mid:.2f}" for iv in grp.index], rotation=30)
    axes[1].set_xlabel("organ-profile conservation (bin midpoint)")
    axes[1].set_ylabel(r"median $\omega$")
    rho, p = stats.spearmanr(ok["profile_pearson"], ok["omega"])
    axes[1].set_title(f"expression conservation vs sequence rate\n"
                      r"Spearman $\rho$" f"={rho:+.3f} (p={p:.1g})")

    org = (d.dropna(subset=["enriched_in"])
             .groupby("enriched_in")
             .agg(n=("profile_pearson", "size"),
                  med=("profile_pearson", "median")))
    org = org[org["n"] >= 25].sort_values("med")
    axes[2].barh(org.index, org["med"], color=PALETTE[2])
    axes[2].set_xlabel("median organ-profile conservation")
    axes[2].set_title("by GTEx tissue-enrichment class")
    axes[2].grid(axis="y", alpha=0)

    p = os.path.join(FIGDIR, "t10_expression_divergence.png")
    fig.savefig(p)
    plt.close(fig)
    return p


def report(force=False):
    print("=" * 72)
    print("T10 - organ-profile conservation vs sequence evolution (human/mouse)")
    print("=" * 72)
    out = os.path.join(DERIVED, "t10_profile_conservation.tsv.gz")
    if os.path.exists(out) and not force:
        d = pd.read_csv(out, sep="\t", index_col=0)
        organs = ORGANS
        log(f"loaded cached table for {len(d)} orthologues")
    else:
        d, organs = profile_conservation()
        to_tsv_gz(d, out)
    print(f"organs compared (all stages pooled): {len(ORGANS)}")

    print(f"1:1 orthologues with adult profiles in both species: {len(d)}")
    print(f"median profile Pearson  : {d['profile_pearson'].median():+.3f}")
    print(f"median profile Spearman : {d['profile_spearman'].median():+.3f}")
    print(f"peak organ conserved    : {d['top_organ_conserved'].mean():.1%} "
          f"of orthologues")

    ok = d.dropna(subset=["omega", "profile_pearson"])
    rho, p = stats.spearmanr(ok["profile_pearson"], ok["omega"])
    t_ok = d.dropna(subset=["profile_pearson", "tau"])
    rho_t, p_t = stats.spearmanr(t_ok["profile_pearson"], t_ok["tau"])
    print(f"\nSpearman(profile conservation, omega) = {rho:+.3f} "
          f"(p={p:.2g}, n={len(ok)})")
    print(f"Spearman(profile conservation, tau)   = {rho_t:+.3f} "
          f"(p={p_t:.2g}, n={len(t_ok)})")
    print("  -> the tau correlation is why the per-organ table below is the")
    print("     interesting one: pooled, this is largely a tau effect again.")

    print("\nby GTEx tissue-enrichment class:")
    print(f"   {'organ':<18}{'n':>5}{'prof_r':>9}{'peak_kept':>11}{'omega':>9}")
    g = d.dropna(subset=["enriched_in"]).groupby("enriched_in")
    rows = []
    for organ, sub in sorted(g, key=lambda kv: -kv[1]["profile_pearson"].median()):
        if len(sub) < 25:
            continue
        print(f"   {organ:<18}{len(sub):>5}{sub['profile_pearson'].median():>+9.3f}"
              f"{sub['top_organ_conserved'].mean():>10.0%}"
              f"{sub['omega'].median():>9.4f}")
        rows.append(dict(organ=organ, n=len(sub),
                         profile_r=sub["profile_pearson"].median(),
                         peak_conserved=sub["top_organ_conserved"].mean(),
                         omega=sub["omega"].median()))
    pd.DataFrame(rows).to_csv(os.path.join(DERIVED, "t10_organ_profile.tsv"),
                              sep="\t", index=False)

    print("\nwhere peak organ MOVED between human and mouse (top shifts):")
    moved = d[~d["top_organ_conserved"]]
    shifts = (moved.groupby(["top_organ_human", "top_organ_mouse"])
                   .size().sort_values(ascending=False).head(8))
    for (a, b), n in shifts.items():
        print(f"   human {a:<11} -> mouse {b:<11} {n:>5} orthologues")

    try:
        da, organs_a = profile_conservation(adult_only=True)
        oka = da.dropna(subset=["omega", "profile_pearson"])
        rho_a, p_a = stats.spearmanr(oka["profile_pearson"], oka["omega"])
        print(f"\nsensitivity check - adult samples only "
              f"({len(organs_a)} shared organs, n={len(da)}):")
        print(f"   median profile Pearson {da['profile_pearson'].median():+.3f}, "
              f"peak organ conserved {da['top_organ_conserved'].mean():.1%}, "
              f"Spearman(prof, omega) {rho_a:+.3f} (p={p_a:.2g})")
        to_tsv_gz(da, os.path.join(DERIVED,
                                   "t10_profile_conservation_adult.tsv.gz"))
    except Exception as exc:                              # noqa: BLE001
        print(f"\nadult-only sensitivity check failed: {exc}")

    p = make_figure(d, organs)
    print(f"\nfigure written to: {p}")
    return d


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    report(force=a.force)
