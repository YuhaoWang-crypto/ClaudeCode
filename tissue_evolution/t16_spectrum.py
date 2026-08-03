"""
T16 - is the MUTATION INPUT different for testis genes, or only the selection?

T5/T13/T14 keep finding a tissue effect that no selective covariate explains.
There is a non-selective explanation that has not been tested: the mutational
process itself. Testis-enriched genes are transcribed in the germline, where
mutations actually arise, so transcription-coupled repair, transcription-
associated mutagenesis and the post-meiotic collapse of repair capacity all act
on them specifically. If that is happening, it should show up as a shifted
mutation SPECTRUM at sites where selection is weakest - not just a shifted rate.

Design
------
Sites: FOUR-FOLD DEGENERATE third positions only. At those, every one of the
three possible changes is synonymous, so amino-acid selection is removed by
construction rather than by assumption. (Codon-usage selection remains, but it
is weak in mammals.)

Polarisation: human/mouse/chicken triplets. Where human and mouse differ at a
4-fold site and chicken matches one of them, the matching base is ancestral by
parsimony. Without this the spectrum is only symmetric classes and the most
informative signals - GC->AT versus AT->GC, and CpG deamination - are invisible.

Normalisation: every class is a rate = observed changes / ANCESTRAL OPPORTUNITY
(the count of ancestral bases in the relevant context at 4-fold sites). Raw
proportions would mostly report each gene's base composition: a GC-rich gene has
more C>T opportunity, so more C>T changes, with no difference in process.

What a positive result would mean: a spectrum shift is evidence about mutation
INPUT, which selection cannot produce. A null result means the tissue effect
stays on the selection side of the ledger.
"""
import os
import csv
import argparse
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy import stats

from .common import DERIVED, FIGDIR, apply_figure_style, PALETTE, log
from . import t3_dnds, t4_tissue_rates, t5_confounders

# 4-fold degenerate families: 3rd position free
FOURFOLD_PREFIX = {"GC", "GG", "CC", "AC", "GT", "TC", "CG", "CT"}
BASES = "ACGT"
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}
# strand-symmetric classes, keyed by the pyrimidine-ancestral representative
CLASSES = ["A>C", "A>G", "A>T", "C>A", "C>G", "C>T"]
MIN_SITES = 60
MIN_ENRICHED = 40


def _canon(anc, der):
    """Collapse a substitution and its reverse-complement into one class."""
    if anc in "AC":
        return f"{anc}>{der}"
    return f"{COMPLEMENT[anc]}>{COMPLEMENT[der]}"


def codon_align_indexed(cds_h, cds_x):
    """Codon alignment that also returns the HUMAN codon index of each column."""
    aligned = t3_dnds.codon_align(cds_h, cds_x)
    if aligned is None:
        return None
    global _A
    al = t3_dnds._ALIGNER
    p1, p2 = t3_dnds.translate(cds_h), t3_dnds.translate(cds_x)
    if p1 is None or p2 is None:
        return None
    a = al.align(p1.replace("X", "A"), p2.replace("X", "A"))[0]
    idx, cod_h, cod_x = [], [], []
    for (s1, e1), (s2, e2) in zip(a.aligned[0], a.aligned[1]):
        for k in range(e1 - s1):
            idx.append(s1 + k)
            cod_h.append(cds_h[3 * (s1 + k):3 * (s1 + k) + 3].upper())
            cod_x.append(cds_x[3 * (s2 + k):3 * (s2 + k) + 3].upper())
    return np.array(idx), cod_h, cod_x


def _load_pairs(sp_key):
    path = os.path.join(DERIVED, f"pairs_human_{sp_key}.tsv")
    out = {}
    with open(path) as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd)
        for r in rd:
            out[r[0]] = (r[4], r[5])          # human cds, ortholog cds
    return out


_G = {}


def _init(mouse, chicken):
    _G["mouse"], _G["chicken"] = mouse, chicken


def _gene_spectrum(gene):
    """Polarised 4-fold spectrum for one gene; None if unusable."""
    mm, gg = _G["mouse"].get(gene), _G["chicken"].get(gene)
    if mm is None or gg is None:
        return None
    hm = codon_align_indexed(mm[0], mm[1])
    hc = codon_align_indexed(gg[0], gg[1])
    if hm is None or hc is None:
        return None
    idx_m, cod_h_m, cod_m = hm
    idx_c, _, cod_c = hc
    chick = {i: c for i, c in zip(idx_c, cod_c)}

    obs = {c: 0 for c in CLASSES}
    obs_cpg = {c: 0 for c in CLASSES}
    opp = {b: 0 for b in BASES}
    opp_cpg = {b: 0 for b in BASES}
    n4 = 0
    for k, i in enumerate(idx_m):
        ch, cm = cod_h_m[k], cod_m[k]
        cc = chick.get(i)
        if cc is None or len(ch) < 3 or len(cm) < 3 or len(cc) < 3:
            continue
        if ch[:2] not in FOURFOLD_PREFIX or ch[:2] != cm[:2] or ch[:2] != cc[:2]:
            continue
        bh, bm, bc = ch[2], cm[2], cc[2]
        if bh not in BASES or bm not in BASES or bc not in BASES:
            continue
        n4 += 1
        # CpG context at a 4-fold third position: the site is a C followed by
        # the next codon's first base G, or a G preceded by the 2nd-position C
        nxt = cod_h_m[k + 1][0] if k + 1 < len(cod_h_m) and cod_h_m[k + 1] else ""
        is_cpg = (bh == "C" and nxt == "G") or (bh == "G" and ch[1] == "C")
        (opp_cpg if is_cpg else opp)[bh] += 1
        if bh == bm:
            continue                                   # no substitution
        # parsimony polarisation: chicken agrees with exactly one of the two
        if bc == bh:
            anc, der = bh, bm
        elif bc == bm:
            anc, der = bm, bh
        else:
            continue
        cl = _canon(anc, der)
        if cl in obs:
            (obs_cpg if is_cpg else obs)[cl] += 1
    if n4 < MIN_SITES:
        return None
    rec = dict(human_gene=gene, n_fourfold=n4)
    for c in CLASSES:
        rec[f"obs_{c}"] = obs[c]
        rec[f"obsCpG_{c}"] = obs_cpg[c]
    for b in BASES:
        rec[f"opp_{b}"] = opp[b]
        rec[f"oppCpG_{b}"] = opp_cpg[b]
    return rec


def build(processes=4, force=False):
    out = os.path.join(DERIVED, "t16_spectrum_counts.tsv.gz")
    if os.path.exists(out) and not force:
        return pd.read_csv(out, sep="\t")
    mouse = _load_pairs("mouse")
    chicken = _load_pairs("chicken")
    genes = sorted(set(mouse) & set(chicken))
    log(f"human/mouse/chicken triplets available: {len(genes)}")
    rows = []
    with Pool(processes, initializer=_init, initargs=(mouse, chicken)) as pool:
        for i, r in enumerate(pool.imap_unordered(_gene_spectrum, genes,
                                                  chunksize=32), 1):
            if r is not None:
                rows.append(r)
            if i % 2000 == 0:
                log(f"  {i}/{len(genes)} processed, {len(rows)} usable")
    df = pd.DataFrame(rows)
    df.to_csv(out, sep="\t", index=False,
              compression={"method": "gzip", "mtime": 0})
    log(f"spectrum for {len(df)} genes")
    return df


# --------------------------------------------------------------------------
# rates and organ comparison
# --------------------------------------------------------------------------
def spectrum_rates(counts):
    """Per gene: rate of each class = observed / ancestral opportunity."""
    anc_of = {c: c.split(">")[0] for c in CLASSES}
    out = pd.DataFrame(index=counts["human_gene"])
    for c in CLASSES:
        anc = anc_of[c]
        comp = COMPLEMENT[anc]
        opp = counts[f"opp_{anc}"] + counts[f"opp_{comp}"]
        out[c] = np.where(opp > 0, counts[f"obs_{c}"].to_numpy()
                          / np.maximum(opp.to_numpy(), 1), np.nan)
    tot_obs = counts[[f"obs_{c}" for c in CLASSES]].sum(axis=1).to_numpy()
    cpg_obs = counts[[f"obsCpG_{c}" for c in CLASSES]].sum(axis=1).to_numpy()
    cpg_opp = counts[[f"oppCpG_{b}" for b in BASES]].sum(axis=1).to_numpy()
    all_opp = (counts[[f"opp_{b}" for b in BASES]].sum(axis=1).to_numpy()
               + cpg_opp)
    ti = counts["obs_A>G"].to_numpy() + counts["obs_C>T"].to_numpy()
    tv = (counts["obs_A>C"] + counts["obs_A>T"]
          + counts["obs_C>A"] + counts["obs_C>G"]).to_numpy()
    gc_to_at = counts["obs_C>A"].to_numpy() + counts["obs_C>T"].to_numpy()
    at_to_gc = counts["obs_A>C"].to_numpy() + counts["obs_A>G"].to_numpy()
    with np.errstate(divide="ignore", invalid="ignore"):
        out["TiTv"] = ti / tv
        out["CpG_frac_obs"] = cpg_obs / np.maximum(tot_obs, 1)
        out["CpG_frac_opp"] = cpg_opp / np.maximum(all_opp, 1)
        out["GC_to_AT_over_AT_to_GC"] = gc_to_at / at_to_gc
        out["GC3"] = ((counts["opp_G"] + counts["opp_C"]
                       + counts["oppCpG_G"] + counts["oppCpG_C"]).to_numpy()
                      / np.maximum(all_opp, 1))
    out["n_subs"] = tot_obs
    out["n_fourfold"] = counts["n_fourfold"].to_numpy()
    return out


def organ_spectrum(d, rates, min_subs=20):
    """Per organ: pooled class rates, plus a GC3-matched comparison."""
    r = rates[rates["n_subs"] >= min_subs]
    d = d.loc[d.index.intersection(r.index)]
    r = r.loc[d.index]
    metrics = CLASSES + ["TiTv", "CpG_frac_obs", "GC_to_AT_over_AT_to_GC", "GC3"]
    rows = []
    for organ in t5_confounders.organ_list(d):
        m = t5_confounders.organ_membership(d, organ)
        sub, rest = r[m.to_numpy()], r[~m.to_numpy()]
        if len(sub) < MIN_ENRICHED:
            continue
        rec = dict(organ=organ, n=len(sub))
        for k in metrics:
            a = sub[k].replace([np.inf, -np.inf], np.nan).dropna()
            b = rest[k].replace([np.inf, -np.inf], np.nan).dropna()
            if len(a) < 20:
                continue
            u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            rec[k] = float(a.median())
            rec[f"{k}_rest"] = float(b.median())
            rec[f"{k}_rb"] = float(2 * u / (len(a) * len(b)) - 1)
            rec[f"{k}_p"] = float(p)
        rows.append(rec)
    res = pd.DataFrame(rows).set_index("organ")
    for k in metrics:
        col = f"{k}_p"
        if col in res:
            res[f"{k}_q"] = t5_confounders._bh(res[col].to_numpy())
    return res


def gc3_matched(d, rates, organ, seed=0, n_bins=10):
    """Compare one organ's spectrum to GC3-matched controls.

    GC3 is the single strongest determinant of the observed spectrum, and organ
    gene sets differ in it, so an unmatched comparison mostly measures base
    composition.
    """
    r = rates[rates["n_subs"] >= 20]
    dd = d.loc[d.index.intersection(r.index)]
    r = r.loc[dd.index]
    m = t5_confounders.organ_membership(dd, organ).to_numpy()
    if m.sum() < MIN_ENRICHED:
        return None
    gc = r["GC3"].to_numpy()
    edges = np.quantile(gc, np.linspace(0, 1, n_bins + 1))
    b = np.clip(np.digitize(gc, edges[1:-1]), 0, n_bins - 1)
    rng = np.random.default_rng(seed)
    case_idx, ctrl_idx = [], []
    for bin_id in range(n_bins):
        ci = np.flatnonzero(m & (b == bin_id))
        pi = np.flatnonzero((~m) & (b == bin_id))
        if len(ci) == 0 or len(pi) == 0:
            continue
        take = rng.choice(pi, size=min(len(ci), len(pi)), replace=False)
        case_idx.extend(ci[:len(take)])
        ctrl_idx.extend(take)
    if len(case_idx) < MIN_ENRICHED:
        return None
    out = {}
    for k in CLASSES + ["TiTv", "CpG_frac_obs", "GC_to_AT_over_AT_to_GC"]:
        a = r.iloc[case_idx][k].replace([np.inf, -np.inf], np.nan).dropna()
        c = r.iloc[ctrl_idx][k].replace([np.inf, -np.inf], np.nan).dropna()
        if len(a) < 20 or len(c) < 20:
            continue
        u, p = stats.mannwhitneyu(a, c, alternative="two-sided")
        out[k] = dict(case=float(a.median()), ctrl=float(c.median()),
                      rb=float(2 * u / (len(a) * len(c)) - 1), p=float(p))
    out["_n"] = len(case_idx)
    return out


def make_figure(res, rates, d):
    plt = apply_figure_style()
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4))

    focus = [o for o in ["Testis", "Brain", "Cerebellum", "Liver", "Blood",
                         "Skeletal_muscle"] if o in res.index]
    x = np.arange(len(CLASSES))
    for i, organ in enumerate(focus):
        vals = [res.loc[organ, c] / res.loc[organ, f"{c}_rest"] for c in CLASSES]
        axes[0].plot(x, vals, "o-", ms=4, lw=1.2, color=PALETTE[i % len(PALETTE)],
                     label=organ)
    axes[0].axhline(1.0, color="k", lw=0.9)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(CLASSES, rotation=30)
    axes[0].set_ylabel("class rate / rest of genome")
    axes[0].set_title("mutation spectrum at 4-fold sites\n(polarised by chicken)")
    axes[0].legend(fontsize=6.5, ncol=2)

    s = res.sort_values("GC3")
    axes[1].scatter(s["GC3"], s["TiTv"], s=32, color=PALETTE[0])
    for organ in s.index:
        axes[1].annotate(organ, (s.loc[organ, "GC3"], s.loc[organ, "TiTv"]),
                         fontsize=6, xytext=(3, 2), textcoords="offset points")
    rho, p = stats.spearmanr(s["GC3"], s["TiTv"])
    axes[1].set_xlabel("organ median GC3 at 4-fold sites")
    axes[1].set_ylabel("Ti/Tv")
    axes[1].set_title(f"spectrum tracks base composition\n" r"$\rho$"
                      f"={rho:+.2f} (p={p:.1g})")

    s2 = res.sort_values("CpG_frac_obs")
    colours = [PALETTE[1] if q < 0.05 else "#BBBBBB"
               for q in s2.get("CpG_frac_obs_q", pd.Series(1.0, index=s2.index))]
    axes[2].barh(s2.index, s2["CpG_frac_obs"], color=colours)
    axes[2].axvline(float(np.median(s2["CpG_frac_obs_rest"])), color="k",
                    ls="--", lw=0.9, label="rest of genome")
    axes[2].set_xlabel("fraction of 4-fold substitutions in CpG context")
    axes[2].set_title("CpG deamination share")
    axes[2].grid(axis="y", alpha=0)
    axes[2].legend(fontsize=7)

    p = os.path.join(FIGDIR, "t16_spectrum.png")
    fig.savefig(p)
    plt.close(fig)
    return p


def report(processes=4, force=False):
    print("=" * 72)
    print("T16 - mutation spectrum at 4-fold degenerate sites, by organ")
    print("=" * 72)
    counts = build(processes=processes, force=force)
    rates = spectrum_rates(counts)
    print(f"genes with a polarised 4-fold spectrum : {len(rates)}")
    print(f"median 4-fold sites per gene           : "
          f"{int(rates['n_fourfold'].median())}")
    print(f"median polarised substitutions per gene: "
          f"{int(rates['n_subs'].median())}")
    print(f"genome-wide Ti/Tv                      : "
          f"{rates['TiTv'].replace([np.inf,-np.inf],np.nan).median():.3f}")
    print(f"genome-wide CpG share of substitutions : "
          f"{rates['CpG_frac_obs'].median():.3f} "
          f"(CpG share of opportunity {rates['CpG_frac_opp'].median():.3f})")

    d = t5_confounders.analysis_frame(t4_tissue_rates.PRIMARY)
    res = organ_spectrum(d, rates)
    print(f"\nper-organ spectrum (Mann-Whitney vs rest of genome, BH across organs)")
    print(f"{'organ':<18}{'n':>5}{'Ti/Tv':>8}{'CpG':>8}{'GC>AT/AT>GC':>13}"
          f"{'GC3':>7}{'CpG q':>10}{'TiTv q':>10}")
    for organ, r in res.sort_values("CpG_frac_obs", ascending=False).iterrows():
        print(f"{organ:<18}{int(r['n']):>5}{r['TiTv']:>8.3f}"
              f"{r['CpG_frac_obs']:>8.3f}{r['GC_to_AT_over_AT_to_GC']:>13.3f}"
              f"{r['GC3']:>7.3f}{r.get('CpG_frac_obs_q', np.nan):>10.2e}"
              f"{r.get('TiTv_q', np.nan):>10.2e}")

    print("\nGC3-matched comparison for the organs the tissue effect hinges on")
    print("   (BH across every organ x metric test below - without it, testing")
    print("    9 metrics on 5 organs yields ~2 hits at p<0.05 by chance alone)")
    focus_organs = ["Testis", "Brain", "Cerebellum", "Blood", "Spleen"]
    matched, flat = {}, []
    for organ in focus_organs:
        m = gc3_matched(d, rates, organ)
        if not m:
            continue
        matched[organ] = m
        for k, v in m.items():
            if k != "_n":
                flat.append((organ, k, v["p"]))
    if flat:
        qs = t5_confounders._bh(np.array([x[2] for x in flat]))
        qmap = {(o, k): q for (o, k, _), q in zip(flat, qs)}
        n_raw = sum(1 for _, _, pv in flat if pv < 0.05)
        n_bh = int((qs < 0.05).sum())
        print(f"   {len(flat)} tests: {n_raw} at raw p<0.05, "
              f"{n_bh} surviving BH q<0.05")
        for organ, m in matched.items():
            print(f"   --- {organ} (n={m['_n']} GC3-matched pairs)")
            for k, v in m.items():
                if k == "_n":
                    continue
                q = qmap[(organ, k)]
                flag = "*" if q < 0.05 else " "
                print(f"      {k:<24} organ {v['case']:.5f}  ctrl {v['ctrl']:.5f}"
                      f"  rb={v['rb']:+.3f}  p={v['p']:.2g}  q={q:.2f}{flag}")

    rates.to_csv(os.path.join(DERIVED, "t16_spectrum_rates.tsv.gz"), sep="\t",
                 compression={"method": "gzip", "mtime": 0})
    res.to_csv(os.path.join(DERIVED, "t16_organ_spectrum.tsv"), sep="\t")
    p = make_figure(res, rates, d)
    print(f"\nfigure written to: {p}")
    return counts, rates, res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--processes", type=int, default=4)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    report(processes=a.processes, force=a.force)
