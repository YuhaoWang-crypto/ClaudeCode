"""
Module 23 — Does a single-cell foundation model (Geneformer) see the M11
fibers?  In-silico deletion on real HEK293T cells vs simple baselines.

M1/M11 are SYMBOLIC claims: paralogues with isomorphic input trees (RAS
H/K/N, RAF A/B/1, MEK1/2, ERK1/2, RSK1-3, ...) collapse onto one fiber
representative. Geneformer-V2-104M (Theodoris lab, pretrained on ~104M
human single-cell transcriptomes) never saw our graph, so it is an
independent, data-driven test of that collapse.

GPU half: grn_pipeline/m23_geneformer_modal.py (Modal, H100, fp32) runs on
all 2,885 cells of the 10x Genomics HEK293T dataset and, per cell, deletes
every detected MAPK-graph gene singly and in pairs, reading the <cls> cell
embedding of the penultimate layer. Each real deletion is paired with a
deletion of a RANK-MATCHED null gene (random non-MAPK gene within +/-3 rank
positions in the same cell): Geneformer deletion effects scale with rank
position, so without this control any "effect" is mostly expression level.

This CPU half turns the cached cell-level records into statistics.

Pre-registered readouts (fixed before the full run):
  * D  single-gene deletion shift (1-cos) vs its rank-matched null.
  * A1 deletion-direction ALIGNMENT  cos(d_a, d_b) of the two displacement
       vectors, minus the same quantity for the rank-matched null pair.
       M11 prediction: fiber-mates are input-equivalent, so deleting either
       should move the cell the same way -> within-fiber pairs score higher
       than cross-fiber pairs.
  * A2 deletion EPISTASIS: projection of d_ab - d_a - d_b on the joint
       direction, minus null. Redundancy prediction: within-fiber > cross.
  * B  context-aware gene-embedding cosine of the two genes.
  Primary test: AUROC (within-fiber vs cross-fiber pairs) of A1, with a
  gene-label permutation null (fiber labels shuffled across genes, fiber
  sizes kept). Every Geneformer score is compared against two cheap
  baselines on the SAME pairs: co-expression (Spearman over all cells) and
  expression-level similarity (-|delta log mean expression|).
"""
import gzip
import json
import os

import numpy as np
from scipy.stats import rankdata, wilcoxon
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "data", "m23_geneformer_raw.json.gz")
FIG = os.path.join(HERE, "..", "figures", "m23_geneformer.png")
MIN_GENE_CELLS = 30          # gene testable if deleted in >= 30 cells
MIN_PAIR_CELLS = 20          # pair testable if co-detected in >= 20 cells
N_PERM = 20000
PARALOG_FIBERS = {"GTPase", "MAP3K", "MAP2K", "MAPK", "RSK", "GEF", "DUSP",
                  "ligand"}  # TF fiber = input-equivalent but NOT paralogues


def load(path=RAW):
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"{path} missing - run `modal run -m grn_pipeline.m23_geneformer_modal`")
    with gzip.open(path, "rt") as fh:
        return json.load(fh)


def auroc(score, label):
    """Mann-Whitney AUROC of score for label==1 vs label==0 (ties -> 0.5)."""
    score, label = np.asarray(score, float), np.asarray(label, bool)
    n1, n0 = label.sum(), (~label).sum()
    if n1 == 0 or n0 == 0:
        return np.nan
    r = rankdata(score)
    return (r[label].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = p[o] * len(p) / (np.arange(len(p)) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q); out[o] = np.minimum(q, 1)
    return out


def gene_table(r):
    rows = []
    for g in r["genes"]:
        real = np.array(r["single"][g]["real"])
        null = np.array([np.nan if x is None else x for x in r["single"][g]["null"]])
        ok = ~np.isnan(null)
        if ok.sum() < MIN_GENE_CELLS:
            continue
        real, null = real[ok], null[ok]
        rows.append({"gene": g, "fiber": r["fiber"][g], "n": int(ok.sum()),
                     "real": float(np.median(real)),
                     "null": float(np.median(null)),
                     "ratio": float(np.median(real) / np.median(null)),
                     "p": float(wilcoxon(real, null).pvalue)})
    for row, q in zip(rows, bh([x["p"] for x in rows])):
        row["q"] = float(q)
    return rows


def pair_table(r, genes):
    rows = []
    for key, v in r["pair"].items():
        a, b = key.split("|")
        if a not in genes or b not in genes or len(v["align_real"]) < MIN_PAIR_CELLS:
            continue
        ar, an = np.array(v["align_real"]), np.array(v["align_null"])
        er, en = np.array(v["epi_real"]), np.array(v["epi_null"])
        ma, mb = r["mean_expr"][a], r["mean_expr"][b]
        rows.append({"pair": key, "a": a, "b": b, "n": len(ar),
                     "align": float(np.mean(ar - an)),
                     "align_p": float(wilcoxon(ar, an).pvalue),
                     "epi": float(np.mean(er - en)),
                     "epi_p": float(wilcoxon(er, en).pvalue),
                     "gemb": float(np.mean(r["gene_emb_cos"][key])),
                     "coexpr": float(np.nan_to_num(r["coexpr"][key])),
                     "exprsim": -abs(np.log(ma) - np.log(mb))})
    return rows


SCORERS = [("align", "GF deletion alignment"),
           ("epi", "GF deletion epistasis"),
           ("gemb", "GF gene-embedding cos"),
           ("coexpr", "co-expression (baseline)"),
           ("exprsim", "expression level (baseline)")]


def fiber_auroc(pairs, fiber, genes, seed=0):
    """AUROC of each scorer for within- vs cross-fiber pairs, with a
    gene-label permutation p-value (one-sided, AUROC > observed)."""
    rng = np.random.default_rng(seed)
    a = np.array([genes.index(p["a"]) for p in pairs])
    b = np.array([genes.index(p["b"]) for p in pairs])
    lab = np.array([fiber[g] for g in genes])
    scores = {k: np.array([p[k] for p in pairs]) for k, _ in SCORERS}
    within = lab[a] == lab[b]
    obs = {k: auroc(s, within) for k, s in scores.items()}
    null = {k: [] for k in scores}
    for _ in range(N_PERM):
        pl = rng.permutation(lab)
        w = pl[a] == pl[b]
        if w.all() or not w.any():
            continue
        for k, s in scores.items():
            null[k].append(auroc(s, w))
    out = {}
    for k in scores:
        nk = np.array(null[k])
        out[k] = {"auroc": float(obs[k]),
                  "p": float((1 + np.sum(nk >= obs[k])) / (1 + len(nk))),
                  "null_hi": float(np.quantile(nk, 0.95))}
    return out, int(within.sum()), int(len(pairs))


def report(path=RAW, make_fig=True):
    print("=" * 68)
    print("MODULE 23 — Geneformer in-silico deletion vs the M11 fibers")
    print("=" * 68)
    r = load(path)
    m = r["meta"]
    print(f"model {m['model']} @ {m['revision'][:10]}  ({m['precision']}, "
          f"{m['layer']})")
    print(f"data  10x HEK293T: {m['n_cells']}/{m['n_cells_total']} cells, "
          f"median {m['median_genes_per_cell']:.0f} genes/cell; "
          f"{m['forward_passes']} forward passes")

    gt = gene_table(r)
    genes = [x["gene"] for x in gt]
    noise = np.array(r["noise"])
    sig = np.median([x["real"] for x in gt])
    print(f"\nnumerical noise floor (no-op re-run) 1-cos: median "
          f"{np.median(noise):.1e}, 99th pct {np.quantile(noise, .99):.1e}  "
          f"vs median single-gene shift {sig:.1e}")

    print(f"\n[D] single-gene deletion vs rank-matched null "
          f"({len(genes)} testable genes, >= {MIN_GENE_CELLS} cells)")
    print(f"  {'gene':8s} {'fiber':8s} {'n':>5s} {'real':>9s} {'null':>9s} "
          f"{'ratio':>6s} {'q':>8s}")
    for x in sorted(gt, key=lambda x: -x["ratio"]):
        print(f"  {x['gene']:8s} {x['fiber']:8s} {x['n']:5d} {x['real']:9.2e} "
              f"{x['null']:9.2e} {x['ratio']:6.2f} {x['q']:8.1e}")
    n_up = sum(x["q"] < 0.05 and x["ratio"] > 1 for x in gt)
    print(f"  -> {n_up}/{len(gt)} MAPK genes shift the cell MORE than a "
          f"rank-matched random gene (BH q<0.05)")

    pairs = pair_table(r, genes)
    res, n_within, n_pairs = fiber_auroc(pairs, r["fiber"], genes)
    par_genes = [g for g in genes if r["fiber"][g] in PARALOG_FIBERS]
    ppairs = [p for p in pairs if p["a"] in par_genes and p["b"] in par_genes]
    pres, pn_within, pn_pairs = fiber_auroc(ppairs, r["fiber"], par_genes)

    print(f"\n[A/B] within- vs cross-fiber pairs: {n_within} within / "
          f"{n_pairs} pairs (>= {MIN_PAIR_CELLS} co-detected cells)")
    print(f"  {'scorer':30s} {'AUROC':>6s} {'perm p':>7s}   "
          f"paralog fibers only ({pn_within}/{pn_pairs})")
    for k, name in SCORERS:
        print(f"  {name:30s} {res[k]['auroc']:6.3f} {res[k]['p']:7.3f}   "
              f"{pres[k]['auroc']:6.3f} {pres[k]['p']:7.3f}")
    print("\n  within-fiber pairs (excess over rank-matched null):")
    for p in sorted(pairs, key=lambda p: -p["align"]):
        if r["fiber"][p["a"]] == r["fiber"][p["b"]]:
            print(f"    {p['pair']:16s} [{r['fiber'][p['a']]:6s}] n={p['n']:4d} "
                  f"align {p['align']:+.3f} (p={p['align_p']:.1g})  "
                  f"epi {p['epi']:+.3f}  coexpr {p['coexpr']:+.3f}")

    gf, base = res["align"], max(res["coexpr"]["auroc"], res["exprsim"]["auroc"])
    supported = gf["auroc"] > 0.5 and gf["p"] < 0.05
    beats = gf["auroc"] > base
    verdict = (("SUPPORTED" if supported else "NOT SUPPORTED") +
               "; Geneformer " + ("beats" if beats else "does not beat") +
               " the best baseline")
    print(f"\nPRE-REGISTERED VERDICT (A1 alignment, all fibers): M11 fiber "
          f"collapse {verdict} (AUROC {gf['auroc']:.2f}, p={gf['p']:.3f}; "
          f"best baseline {base:.2f})")

    if make_fig:
        _figure(r, gt, pairs, res, noise)
    return {"genes": gt, "pairs": pairs, "auroc": res, "auroc_paralog": pres,
            "n_within": n_within, "n_pairs": n_pairs, "n_up": n_up,
            "noise_med": float(np.median(noise)), "signal_med": float(sig),
            "verdict": verdict, "supported": supported, "beats": beats,
            "meta": m}


def _figure(r, gt, pairs, res, noise):
    fig, ax = plt.subplots(2, 2, figsize=(13, 9.5))
    # (a) single-gene shift vs rank-matched null
    a = ax[0, 0]
    gt = sorted(gt, key=lambda x: -x["ratio"])
    x = np.arange(len(gt))
    a.bar(x - 0.2, [g["real"] for g in gt], 0.4, label="MAPK gene deleted",
          color="#2a6fbb")
    a.bar(x + 0.2, [g["null"] for g in gt], 0.4,
          label="rank-matched random gene", color="#bbbbbb")
    a.axhline(np.quantile(noise, .99), color="#c0392b", ls=":",
              label="numerical noise (99th pct)")
    a.set_yscale("log"); a.set_xticks(x)
    a.set_xticklabels([g["gene"] for g in gt], rotation=60, fontsize=8)
    a.set_ylabel("cell-embedding shift 1-cos (median)")
    a.set_title("[D] single-gene deletion vs rank-matched null")
    a.legend(fontsize=8)
    # (b, c) excess alignment / epistasis, within vs cross
    for axx, key, title in [(ax[0, 1], "align", "[A1] deletion-direction alignment"),
                            (ax[1, 0], "epi", "[A2] deletion epistasis")]:
        w = [p[key] for p in pairs if r["fiber"][p["a"]] == r["fiber"][p["b"]]]
        c = [p[key] for p in pairs if r["fiber"][p["a"]] != r["fiber"][p["b"]]]
        rng = np.random.default_rng(1)
        axx.scatter(rng.uniform(-.15, .15, len(c)), c, s=14, color="#999999",
                    label=f"cross-fiber ({len(c)})")
        axx.scatter(1 + rng.uniform(-.15, .15, len(w)), w, s=26,
                    color="#d35400", label=f"within-fiber ({len(w)})")
        for p in pairs:
            if r["fiber"][p["a"]] == r["fiber"][p["b"]] and r["fiber"][p["a"]] != "TF":
                axx.annotate(p["pair"], (1.18, p[key]), fontsize=6.5)
        axx.axhline(0, color="k", lw=0.6)
        axx.set_xticks([0, 1]); axx.set_xticklabels(["cross", "within"])
        axx.set_xlim(-0.5, 2.0)
        axx.set_ylabel(f"{key}: real - rank-matched null")
        axx.set_title(f"{title}  AUROC={res[key]['auroc']:.2f} "
                      f"(perm p={res[key]['p']:.2f})")
        axx.legend(fontsize=8, loc="lower left")
    # (d) AUROC of every scorer vs permutation null
    d = ax[1, 1]
    names = [n for _, n in SCORERS]
    vals = [res[k]["auroc"] for k, _ in SCORERS]
    hi = [res[k]["null_hi"] for k, _ in SCORERS]
    cols = ["#2a6fbb"] * 3 + ["#7f8c8d"] * 2
    y = np.arange(len(names))
    d.barh(y, vals, color=cols)
    d.scatter(hi, y, marker="|", s=300, color="#c0392b",
              label="95th pct of label-permutation null")
    d.axvline(0.5, color="k", lw=0.6)
    d.set_yticks(y); d.set_yticklabels(names, fontsize=9); d.invert_yaxis()
    d.set_xlim(0, 1); d.set_xlabel("AUROC: within- vs cross-fiber pairs")
    d.set_title("[A/B] who recovers the M11 fibers?")
    d.legend(fontsize=8, loc="upper right")
    m = r["meta"]
    fig.suptitle(f"M23 Geneformer-V2-104M in-silico deletion on {m['n_cells']} "
                 f"HEK293T cells — test of the M11 fibration", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIG, dpi=130)
    plt.close(fig)
    print(f"figure -> figures/m23_geneformer.png")


if __name__ == "__main__":
    report()
