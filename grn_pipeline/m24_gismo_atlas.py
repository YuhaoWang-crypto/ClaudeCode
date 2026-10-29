"""
Module 24 — GISMO annotation of every gene in the M14 pathway atlas.

For each of the 18 M14 pathways, every member gene gets deep-mammalian
constraint (GISMO gene loss, GISMO-mis fixed missense; Liao et al. Nat Genet
2026) next to human constraint (gnomAD v2.1.1 LOEUF), dosage sensitivity
(pHaplo), mouse knockout lethality and OMIM inheritance, plus the rule-based
tags from gismo.interpret.

Questions answered (all computed):
  1. Is each pathway more constrained than the genome? (one-sided
     Mann-Whitney vs all scored genes, Benjamini-Hochberg across pathways)
  2. Do switch-like ("high") pathways carry more constrained genes than
     feedback ("medium") or flux-state pathways?
  3. Fiber redundancy: are singleton nodes (layer with one gene - nothing can
     substitute) more constrained than members of paralogue fibers
     (layer with >= 2 genes)? This tests the M1/M11 premise that fiber
     members are interchangeable.

  4. Ascertainment control: are lethal pathway genes more constrained than
     other mouse-lethal genes? (M14 lists famous genes: 70% lethal vs 20%.)
  5. QC: genes whose GISMO "losses" are likely technical - GC-rich CDS
     (top 20%, see gismo.core.cds_gc) or assembly-quality-sensitive.

Rigor: ✅ scores, ranks and test statistics are computed; ⚠️ pathway gene
lists are the compact M14 cascades (4-15 genes each), not full pathway
databases, so per-pathway tests have low power; mouse lists are matched by
upper-cased symbol (approximate orthology).
"""
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from grn_pipeline.m14_atlas import PATHWAYS

NON_GENES = {"NICD"}            # M14 placeholders that are not gene symbols
OUT_TSV = "gismo/data/m14_pathway_gismo.tsv"


def _bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); prev = 1.0
    for rank, i in zip(range(n, 0, -1), o[::-1]):
        prev = min(prev, p[i] * n / rank); q[i] = prev
    return q


def _scores():
    """GISMO table joined with reference data when available."""
    try:
        from gismo.benchmark import _gene_list, merged
        d = merged()
        sym = d.symbol.astype(str)
        d["mouse_lethal"] = sym.isin(_gene_list("mouse_lethal.tsv", upper=True))
        d["omim_AR"] = sym.isin(_gene_list("all_ar.tsv"))
        d["omim_AD"] = sym.isin(_gene_list("all_ad.tsv"))
        full = True
    except Exception as e:                       # offline: GISMO columns only
        from gismo.core import load_scores
        print(f"  (reference data unavailable: {type(e).__name__}; GISMO-only annotation)")
        d, full = load_scores(), False
    return d.drop_duplicates("symbol"), full


def annotate():
    from gismo.interpret import classify
    d, full = _scores()
    idx = d.set_index("symbol")
    rows = []
    for pw, spec in PATHWAYS.items():
        L = len(spec["layers"])
        for li, (layer, members) in enumerate(spec["layers"]):
            genes = [g for g in members if g not in NON_GENES]
            for g in genes:
                r = {"pathway": pw, "switch": spec["switch"], "layer": layer,
                     "layer_index": li, "n_layers": L,
                     "fiber_size": len(genes),
                     "node_type": "singleton" if len(genes) == 1 else "paralogue",
                     "gene": g, "found": g in idx.index}
                if r["found"]:
                    s = idx.loc[g]
                    for c in ["gene_id", "chrom", "n_loss", "gismo", "gismo_decile",
                              "gismo_hq_decile", "loss_n50_shift", "gc_cds", "gc_cds_pct",
                              "gismo_mis", "gismo_mis_decile", "loeuf", "cds_length",
                              "pHaplo", "pTriplo", "mouse_lethal", "omim_AR", "omim_AD"]:
                        if c in s.index:
                            r[c] = s[c]
                    r["tags"] = "; ".join(classify(s.to_dict()))
                rows.append(r)
    ann = pd.DataFrame(rows)
    # within-fiber rank by GISMO (1 = most constrained member of its paralogue family)
    ann["fiber_rank_gismo"] = (ann.groupby(["pathway", "layer"]).gismo
                                  .rank(method="min"))
    return ann, d, full


def pathway_summary(ann, genome):
    bg = genome.gismo.dropna(); bgm = genome.gismo_mis.dropna()
    rows = []
    for pw, g in ann.groupby("pathway", sort=False):
        x, xm = g.gismo.dropna(), g.gismo_mis.dropna()
        r = {"pathway": pw, "switch": g.switch.iloc[0], "n_genes": len(g),
             "n_scored": len(x),
             "median_gismo_decile": float(g.gismo_decile.median()),
             "median_mis_decile": float(g.gismo_mis_decile.median()),
             "pct_top2_gismo": round(float((g.gismo_decile <= 2).mean() * 100), 1),
             "pct_top2_mis": round(float((g.gismo_mis_decile <= 2).mean() * 100), 1),
             "p_gismo": float(mannwhitneyu(x, bg, alternative="less").pvalue) if len(x) else np.nan,
             "p_mis": float(mannwhitneyu(xm, bgm, alternative="less").pvalue) if len(xm) else np.nan}
        if "mouse_lethal" in g:
            r["pct_mouse_lethal"] = round(float(g.mouse_lethal.mean() * 100), 1)
            r["median_loeuf"] = float(g.loeuf.median())
        gcr = g.tags.fillna("").str.contains("GC-rich CDS")
        r["median_gismo_decile_noGCrisk"] = float(g.gismo_decile[~gcr].median())
        r["n_gc_risk"] = int(gcr.sum())
        least = g.sort_values("gismo").dropna(subset=["gismo"])
        k = 2 if len(least) >= 4 else 1
        r["most_constrained"] = ",".join(least.gene.head(k))
        r["least_constrained"] = ",".join(least.gene.tail(k)[::-1])
        rows.append(r)
    s = pd.DataFrame(rows)
    s["q_gismo"] = _bh(s.p_gismo.fillna(1)); s["q_mis"] = _bh(s.p_mis.fillna(1))
    return s.sort_values("median_gismo_decile")


def contrasts(ann, genome):
    u = ann.drop_duplicates("gene")                  # each gene counted once
    out = {}
    bg = genome.gismo.dropna()
    out["all_pathway_genes"] = {
        "n": int(u.gismo.notna().sum()),
        "median_gismo_decile": float(u.gismo_decile.median()),
        "median_mis_decile": float(u.gismo_mis_decile.median()),
        "pct_top2_gismo": round(float((u.gismo_decile <= 2).mean() * 100), 1),
        "p_vs_genome_gismo": float(mannwhitneyu(u.gismo.dropna(), bg, alternative="less").pvalue),
        "p_vs_genome_mis": float(mannwhitneyu(u.gismo_mis.dropna(), genome.gismo_mis.dropna(),
                                              alternative="less").pvalue)}
    if "mouse_lethal" in u:
        out["all_pathway_genes"]["pct_mouse_lethal"] = round(float(u.mouse_lethal.mean() * 100), 1)
        out["all_pathway_genes"]["genome_pct_mouse_lethal"] = round(
            float(genome.mouse_lethal.mean() * 100), 1)

    def two(a, b, col):
        a, b = a[col].dropna(), b[col].dropna()
        return {"median_a": float(a.median()), "median_b": float(b.median()),
                "n_a": len(a), "n_b": len(b),
                "p_a_more_constrained": float(mannwhitneyu(a, b, alternative="less").pvalue)}

    # 1b. ascertainment control: M14 lists famous genes (70% mouse-lethal vs
    #     20% genome-wide). Compare lethal pathway genes with all other
    #     mouse-lethal genes instead of with the whole genome.
    if "mouse_lethal" in u and "mouse_lethal" in genome:
        pl = u[u.mouse_lethal == True]                          # noqa: E712
        ctrl = genome[genome.mouse_lethal & ~genome.symbol.isin(u.gene)]
        out["vs_other_mouse_lethal"] = {
            m: {"median_decile_pathway": float(pl[dc].median()),
                "median_decile_other_lethal": float(ctrl[dc].median()),
                "n_pathway": int(pl[m].notna().sum()), "n_other": int(ctrl[m].notna().sum()),
                "p_pathway_more_constrained": float(mannwhitneyu(
                    pl[m].dropna(), ctrl[m].dropna(), alternative="less").pvalue)}
            for m, dc in [("gismo", "gismo_decile"), ("gismo_mis", "gismo_mis_decile")]}

    # 2. switch class (gene-level, a gene in two classes counted in both)
    sw = {}
    for cls in ["high", "medium", "flux"]:
        g = ann[ann.switch == cls].drop_duplicates("gene")
        sw[cls] = {"n": int(g.gismo.notna().sum()),
                   "median_gismo_decile": float(g.gismo_decile.median()),
                   "median_mis_decile": float(g.gismo_mis_decile.median())}
        if "mouse_lethal" in g:
            sw[cls]["pct_mouse_lethal"] = round(float(g.mouse_lethal.mean() * 100), 1)
    hi = ann[ann.switch == "high"].drop_duplicates("gene")
    rest = ann[ann.switch != "high"].drop_duplicates("gene")
    rest = rest[~rest.gene.isin(hi.gene)]
    sw["high_vs_rest_gismo"] = two(hi, rest, "gismo")
    sw["high_vs_rest_mis"] = two(hi, rest, "gismo_mis")
    out["switch_class"] = sw

    # 3. singleton vs paralogue nodes (gene-level; a gene that is a singleton
    #    anywhere is called singleton)
    single = set(ann.loc[ann.node_type == "singleton", "gene"])
    s = u[u.gene.isin(single)]; p = u[~u.gene.isin(single)]
    out["singleton_vs_paralogue"] = {"gismo": two(s, p, "gismo"),
                                     "gismo_mis": two(s, p, "gismo_mis")}
    if "loeuf" in u:
        out["singleton_vs_paralogue"]["loeuf"] = two(s, p, "loeuf")
        out["singleton_vs_paralogue"]["pct_mouse_lethal"] = (
            round(float(s.mouse_lethal.mean() * 100), 1), round(float(p.mouse_lethal.mean() * 100), 1))
    out["singleton_vs_paralogue"]["genes_singleton"] = sorted(s.gene)

    # 4. evolutionary symmetry breaking inside paralogue fibers: decile spread
    fib = (ann[ann.node_type == "paralogue"].dropna(subset=["gismo_decile"])
           .groupby(["pathway", "layer"]))
    spread = fib.gismo_decile.agg(lambda x: x.max() - x.min())
    spread_m = fib.gismo_mis_decile.agg(lambda x: x.max() - x.min())
    wide = spread[spread >= 5].index
    out["fiber_spread"] = {
        "n_fibers": int(len(spread)),
        "median_gismo_range": float(spread.median()),
        "median_mis_range": float(spread_m.median()),
        "pct_gismo_range_ge5": round(float((spread >= 5).mean() * 100), 1),
        "wide_fibers": [f"{pw}:{ly}" for pw, ly in wide],
    }
    return out


def report(fig_path="figures/m24_gismo_atlas.png", out_tsv=OUT_TSV):
    print("=" * 68)
    print("MODULE 24 — GISMO annotation of the M14 pathway atlas")
    print("=" * 68)
    ann, genome, full = annotate()
    genome = genome.dropna(subset=["gismo"])
    missing = ann.loc[~ann.found, "gene"].tolist()
    print(f"{len(ann)} pathway-gene entries, {ann.gene.nunique()} unique genes, "
          f"{ann.found.sum()} scored; not in GISMO table: {missing or 'none'}")
    ann.to_csv(out_tsv, sep="\t", index=False, float_format="%.4g")
    print(f"per-gene annotation -> {out_tsv}")

    s = pathway_summary(ann, genome)
    cols = ["pathway", "switch", "n_scored", "median_gismo_decile",
            "median_gismo_decile_noGCrisk", "median_mis_decile",
            "pct_top2_gismo", "pct_mouse_lethal", "q_gismo", "q_mis",
            "most_constrained", "least_constrained"]
    cols = [c for c in cols if c in s]
    with pd.option_context("display.width", 220, "display.max_columns", 20):
        print(s[cols].round(3).to_string(index=False))

    c = contrasts(ann, genome)
    a = c["all_pathway_genes"]
    print(f"\nall {a['n']} pathway genes: median GISMO decile {a['median_gismo_decile']}, "
          f"GISMO-mis decile {a['median_mis_decile']} (genome = 5.5); "
          f"p(GISMO)={a['p_vs_genome_gismo']:.2g}, p(mis)={a['p_vs_genome_mis']:.2g}"
          + (f"; mouse-lethal {a['pct_mouse_lethal']}% vs genome {a['genome_pct_mouse_lethal']}%"
             if "pct_mouse_lethal" in a else ""))
    if "vs_other_mouse_lethal" in c:
        for m, v in c["vs_other_mouse_lethal"].items():
            print(f"  ascertainment control [{m:9s}]: mouse-lethal pathway genes d{v['median_decile_pathway']} "
                  f"(n={v['n_pathway']}) vs other mouse-lethal genes d{v['median_decile_other_lethal']} "
                  f"(n={v['n_other']}), p={v['p_pathway_more_constrained']:.2g}")
    for k in ["high", "medium", "flux"]:
        v = c["switch_class"][k]
        print(f"  switch={k:6s}: n={v['n']:3d}  GISMO d{v['median_gismo_decile']}  "
              f"mis d{v['median_mis_decile']}  lethal {v.get('pct_mouse_lethal', 'NA')}%")
    h = c["switch_class"]
    print(f"  high vs rest: p(GISMO)={h['high_vs_rest_gismo']['p_a_more_constrained']:.2g}, "
          f"p(mis)={h['high_vs_rest_mis']['p_a_more_constrained']:.2g}")
    sp = c["singleton_vs_paralogue"]
    for m in ["gismo", "gismo_mis", "loeuf"]:
        if m in sp:
            v = sp[m]
            print(f"  singleton vs paralogue [{m:9s}]: median {v['median_a']:.3g} (n={v['n_a']}) "
                  f"vs {v['median_b']:.3g} (n={v['n_b']}), p={v['p_a_more_constrained']:.2g}")
    if "pct_mouse_lethal" in sp:
        print(f"  mouse-lethal: singleton {sp['pct_mouse_lethal'][0]}% vs paralogue "
              f"{sp['pct_mouse_lethal'][1]}%")
    f = c["fiber_spread"]
    print(f"  paralogue fibers: {f['n_fibers']}; median within-fiber decile range "
          f"GISMO {f['median_gismo_range']}, mis {f['median_mis_range']}; "
          f"{f['pct_gismo_range_ge5']}% span >=5 GISMO deciles: {', '.join(f['wide_fibers'])}")
    sens = ann[ann.tags.fillna("").str.contains("assembly-sensitive")].drop_duplicates("gene")
    c["assembly_sensitive"] = sens[["gene", "gismo_decile", "gismo_hq_decile",
                                    "loss_n50_shift", "loeuf"]].to_dict("records") \
        if "loeuf" in sens else sens[["gene", "gismo_decile", "gismo_hq_decile"]].to_dict("records")
    print(f"  assembly-sensitive pathway genes (GISMO decile improves >=3 with high-quality "
          f"assemblies only): " + (", ".join(f"{r['gene']} d{int(r['gismo_decile'])}->d{int(r['gismo_hq_decile'])}"
                                            for r in c["assembly_sensitive"]) or "none"))
    gcf = ann[ann.tags.fillna("").str.contains("GC-rich CDS")].drop_duplicates("gene")
    c["gc_dropout_risk"] = gcf[["gene", "gismo_decile", "gc_cds"]].to_dict("records")
    print("  GC-dropout-risk pathway genes (GISMO d8-10 with top-20% CDS GC): "
          + (", ".join(f"{r['gene']} d{int(r['gismo_decile'])} GC {r['gc_cds']:.0%}"
                       for r in c["gc_dropout_risk"]) or "none"))
    _figure(ann, s, fig_path)
    print(f"figure written to: {fig_path}")
    return {"summary": s.to_dict("records"), "contrasts": c, "n_genes": int(ann.gene.nunique())}


def _figure(ann, s, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    order = s.pathway.tolist()[::-1]                 # most constrained on top
    cmap = plt.get_cmap("RdYlBu", 10)                # d1 red -> d10 blue
    fig = plt.figure(figsize=(17, 9))
    ax = fig.add_axes([0.13, 0.08, 0.55, 0.85])
    swc = {"high": "#c0392b", "medium": "#2980b9", "flux": "#27ae60"}
    for y, pw in enumerate(order):
        g = ann[ann.pathway == pw].reset_index(drop=True)
        for x, r in g.iterrows():
            dec = r.get("gismo_decile")
            col = "#dddddd" if pd.isna(dec) else cmap(int(dec) - 1)
            edge = "k" if r.get("mouse_lethal") is True else "none"
            marker = "s" if r.node_type == "singleton" else "o"
            ax.scatter(x, y, s=260, c=[col], marker=marker, edgecolors=edge, linewidths=1.4)
            gc_flag = "GC-rich CDS" in str(r.get("tags", ""))
            ax.text(x, y - 0.38, r.gene + ("†" if gc_flag else ""), ha="center", va="top",
                    fontsize=5.6, color="#7f4f00" if gc_flag else "k")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=8)
    for t, pw in zip(ax.get_yticklabels(), order):
        t.set_color(swc[ann.loc[ann.pathway == pw, "switch"].iloc[0]])
    ax.set_xlim(-0.7, ann.groupby("pathway").size().max() - 0.3)
    ax.set_ylim(-0.9, len(order) - 0.4)
    ax.set_xticks([]); ax.set_xlabel("genes in cascade order (upstream → downstream)")
    ax.set_title("GISMO decile per gene (red = most constrained, d1; blue = d10)\n"
                 "■ singleton node  ● paralogue-fiber member  black edge = mouse KO lethal; "
                 "label colour = switch class;  † GC-rich CDS: 'losses' may be assembly dropout",
                 fontsize=9)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0.5, 10.5))
    cb = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.01, ticks=range(1, 11))
    cb.set_label("GISMO decile")

    ax2 = fig.add_axes([0.76, 0.08, 0.21, 0.85])
    y = np.arange(len(order))
    ss = s.set_index("pathway").loc[order]
    ax2.barh(y + 0.18, ss.median_gismo_decile, 0.36, color="#b5838d", label="GISMO")
    ax2.barh(y - 0.18, ss.median_mis_decile, 0.36, color="#6a9fb5", label="GISMO-mis")
    ax2.axvline(5.5, color="k", ls="--", lw=0.8)
    ax2.text(5.6, len(order) - 0.6, "genome", fontsize=7)
    for yi, (qg, qm) in enumerate(zip(ss.q_gismo, ss.q_mis)):
        star = ("*" if min(qg, qm) < 0.05 else "")
        ax2.text(0.2, yi, star, va="center", fontsize=11)
    ax2.set_yticks(y); ax2.set_yticklabels([])
    ax2.set_xlim(0, 10); ax2.set_xlabel("median decile (lower = more constrained)")
    ax2.set_ylim(-0.9, len(order) - 0.4)
    ax2.legend(fontsize=7, loc="lower right"); ax2.set_title("* BH q<0.05 vs genome", fontsize=8)
    fig.savefig(path, dpi=125)
    plt.close(fig)


if __name__ == "__main__":
    report()
