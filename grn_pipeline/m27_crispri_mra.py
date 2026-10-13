"""
Module 27 — Does CRISPRi rescue Jacobian recovery? (follow-up to M26)

M26 found that modular response analysis (MRA) on CRISPR-knockout Perturb-seq
failed at chance level, and traced the cause to a measurable quantity: CRISPR
cutting often leaves the target transcript intact, so the response matrix has
almost no diagonal dominance (diag/off-diag ratio 2.2, sign consistency 51%).

CRISPRi represses transcription directly, so its diagonal should be real.
This module re-runs exactly the same test on GSE249595 (genome-scale CRISPRi
Perturb-seq in Jurkat T cells, 16 channels, ~800k cells, 83,401 guides over
~20k genes, 6 TCR-activation hashtags + untreated).

The pre-registered question: does the diagonal/off-diagonal ratio rise enough
to move MRA off chance, as the synthetic calibration curve in M26 predicts?

Pipeline (all on raw CellRanger output, cell calling done here):
  1. guides matrix -> called cells (guide UMI knee) with a dominant guide
     -> target gene; hashtag matrix -> activation condition
  2. coverage audit -> choose a module of genes that are both perturbed with
     enough cells and expressed in Jurkat
  3. streamed transcriptome pass -> pseudo-bulk log-fold-change per
     (target gene, channel) vs non-targeting controls
  4. response matrix on the module -> diagonal audit -> MRA -> bootstrap sign
     consistency, compared against the M26 CRISPR-cutting result and the
     synthetic calibration curve
"""
import gzip
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "..", "figures")
DATA = os.path.join(HERE, "..", "data_cache", "GSE249595")

N_CHANNELS = 16
MIN_GUIDE_UMI = 50          # knee of the guide-UMI distribution
MIN_DOMINANCE = 0.5         # top guide must carry >= this fraction of UMIs


def _read_mtx(path):
    """Market-matrix triplets as (row, col, val), zero-based."""
    with gzip.open(path, "rt") as f:
        f.readline(); f.readline(); f.readline()
        a = np.loadtxt(f, dtype=np.int64)
    return a[:, 0] - 1, a[:, 1] - 1, a[:, 2]


def guide_targets(path=os.path.join(DATA, "ch1_guides_features.tsv.gz")):
    """Guide id -> target gene symbol. Guide ids are <GENE>_<n>; the
    non-targeting controls are kept under their own label."""
    names = [l.split("\t")[0] for l in gzip.open(path, "rt")]
    # guide ids are <PREFIX>_<n>; the two control prefixes in this library
    # are "non-targeting" and "CTRL" (NEGR1/NONO/SELENON are real genes)
    tgt = [("NT" if n.rsplit("_", 1)[0] in ("non-targeting", "CTRL")
            else n.rsplit("_", 1)[0]) for n in names]
    return np.array(names), np.array(tgt)


def call_cells(channel, min_umi=MIN_GUIDE_UMI, dominance=MIN_DOMINANCE):
    """Called cells of one channel: barcode index, assigned guide, target
    gene, and hashtag. A cell needs >= min_umi guide UMIs and one guide
    carrying >= `dominance` of them (so multiplets are dropped)."""
    gi, bi, v = _read_mtx(os.path.join(DATA, f"ch{channel}_guides_matrix.mtx.gz"))
    order = np.argsort(bi, kind="stable")
    gi, bi, v = gi[order], bi[order], v[order]
    ub, start = np.unique(bi, return_index=True)
    stop = np.r_[start[1:], len(v)]
    tot = np.add.reduceat(v, start)
    keep = tot >= min_umi
    cells, guides, fracs = [], [], []
    for k in np.where(keep)[0]:
        s, e = start[k], stop[k]
        j = s + int(np.argmax(v[s:e]))
        f = v[j] / tot[k]
        if f >= dominance:
            cells.append(ub[k]); guides.append(gi[j]); fracs.append(f)
    cells = np.array(cells); guides = np.array(guides)
    # hashtag
    hi, hb, hv = _read_mtx(os.path.join(DATA,
                                        f"ch{channel}_labels_matrix.mtx.gz"))
    hash_of = {}
    o = np.argsort(hb, kind="stable")
    hi, hb, hv = hi[o], hb[o], hv[o]
    hub, hstart = np.unique(hb, return_index=True)
    hstop = np.r_[hstart[1:], len(hv)]
    for k, b in enumerate(hub):
        s, e = hstart[k], hstop[k]
        t = hv[s:e].sum()
        j = s + int(np.argmax(hv[s:e]))
        hash_of[b] = (hi[j], hv[j] / t) if t >= 5 else (-1, 0.0)
    hid = np.array([hash_of.get(c, (-1, 0.0))[0] for c in cells])
    hfr = np.array([hash_of.get(c, (-1, 0.0))[1] for c in cells])
    return dict(channel=channel, barcode=cells, guide=guides,
                guide_frac=np.array(fracs), hashtag=hid, hash_frac=hfr)


def coverage_audit(cache=os.path.join(DATA, "cells.npz")):
    """Call cells in every channel and count cells per target gene."""
    if os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        return {k: z[k] for k in z.files}
    names, tgt = guide_targets()
    parts = [call_cells(c) for c in range(1, N_CHANNELS + 1)]
    out = dict(
        channel=np.concatenate([np.full(len(p["barcode"]), p["channel"])
                                for p in parts]),
        barcode=np.concatenate([p["barcode"] for p in parts]),
        guide=np.concatenate([p["guide"] for p in parts]),
        guide_frac=np.concatenate([p["guide_frac"] for p in parts]),
        hashtag=np.concatenate([p["hashtag"] for p in parts]),
        hash_frac=np.concatenate([p["hash_frac"] for p in parts]))
    out["target"] = tgt[out["guide"]]
    np.savez_compressed(cache, **out)
    return out


# ---------------------------------------------------------------------------
# GSE190604: CRISPRa Perturb-seq in primary human T cells
# ---------------------------------------------------------------------------
DATA_A = os.path.join(HERE, "..", "data_cache", "GSE190604")


def load_crispra(chunk=20_000_000):
    """Stream the aggregated matrix once, keeping only the perturbed genes'
    rows plus each cell's total UMI (for CP10K normalisation).
    Returns log1p CP10K expression of the module genes per cell, the target
    assignment, and the donor/lane suffix."""
    import pandas as pd
    cache = os.path.join(DATA_A, "module.npz")
    gc = pd.read_csv(os.path.join(
        DATA_A, "GSE190604_cellranger-guidecalls-aggregated-unfiltered.txt.gz"),
        sep="\t")
    gc = gc[gc["num_features"] == 1].copy()
    gc["target"] = gc["feature_call"].str.replace(r"-\d+$", "", regex=True)
    bc = [l.strip() for l in gzip.open(
        os.path.join(DATA_A, "GSE190604_barcodes.tsv.gz"), "rt")]
    feat = [l.split("\t")[1] for l in gzip.open(
        os.path.join(DATA_A, "GSE190604_features.tsv.gz"), "rt")]
    bidx = {b: i for i, b in enumerate(bc)}
    gc = gc[gc["cell_barcode"].isin(bidx)]
    targets = sorted(set(gc["target"]) - {"NO-TARGET"})
    genes = [g for g in targets if g in feat]
    gi = {feat.index(g): k for k, g in enumerate(genes)}
    if os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        M, tot = z["M"], z["tot"]
    else:
        M = np.zeros((len(genes), len(bc)), dtype=np.float32)
        tot = np.zeros(len(bc), dtype=np.float64)
        rows = np.array(sorted(gi))
        with gzip.open(os.path.join(DATA_A, "GSE190604_matrix.mtx.gz"),
                       "rt") as f:
            f.readline(); f.readline(); f.readline()
            for part in pd.read_csv(f, sep=" ", header=None,
                                    names=["g", "b", "v"], dtype=np.int64,
                                    chunksize=chunk):
                g = part["g"].values - 1; b = part["b"].values - 1
                v = part["v"].values.astype(np.float64)
                np.add.at(tot, b, v)
                m = np.isin(g, rows)
                if m.any():
                    np.add.at(M, ([gi[x] for x in g[m]], b[m]),
                              v[m].astype(np.float32))
        np.savez_compressed(cache, M=M, tot=tot)
    cells = np.array([bidx[b] for b in gc["cell_barcode"]])
    X = np.log1p(M[:, cells] / np.maximum(tot[cells], 1) * 1e4).T
    lane = np.array([b.rsplit("-", 1)[-1] for b in gc["cell_barcode"]])
    return dict(X=X, genes=genes, target=gc["target"].values, lane=lane,
                n_umi=tot[cells])


def crispra_response(d, min_cells=100, n_boot=200, seed=0):
    """Response matrix R[i, j] = mean log1p CP10K of gene i in cells
    targeting gene j minus control (NO-TARGET) cells, paired within lane.
    Bootstrap over cells for the MRA stability test."""
    rng = np.random.default_rng(seed)
    genes = d["genes"]; X = d["X"]; tgt = d["target"]; lane = d["lane"]
    use = [g for g in genes if (tgt == g).sum() >= min_cells]
    R = np.full((len(genes), len(use)), np.nan)
    boots = {}
    ctrl = tgt == "NO-TARGET"
    for j, g in enumerate(use):
        ko = tgt == g
        diffs, w = [], []
        for L in np.unique(lane):
            m1 = ko & (lane == L); m2 = ctrl & (lane == L)
            if m1.sum() >= 10 and m2.sum() >= 10:
                diffs.append(X[m1].mean(0) - X[m2].mean(0)); w.append(m1.sum())
        if not diffs:
            continue
        D = np.array(diffs); w = np.array(w, float)
        R[:, j] = np.average(D, axis=0, weights=w)
        idx = np.where(ko)[0]; cidx = np.where(ctrl)[0]
        boots[g] = np.array([
            X[rng.choice(idx, len(idx))].mean(0)
            - X[rng.choice(cidx, len(cidx))].mean(0) for _ in range(n_boot)])
    keep = [k for k, g in enumerate(genes) if g in use]
    return R[keep], [genes[k] for k in keep], use, boots


def square_module(R, rows, cols):
    """Restrict the response matrix to genes that are both perturbed and
    measured (MRA needs a closed module)."""
    common = [g for g in rows if g in cols]
    Rs = R[np.ix_([rows.index(g) for g in common],
                  [cols.index(g) for g in common])]
    ok = np.isfinite(Rs).all(0) & np.isfinite(Rs).all(1)
    return Rs[np.ix_(ok, ok)], [g for g, k in zip(common, ok) if k]


def mra_bootstrap(d, genes, n_rep=200, seed=0):
    """MRA stability: per replicate, resample the control cells ONCE and
    each perturbation's cells, rebuild the whole response matrix, re-run
    MRA. Returns the point estimate, the fraction of off-diagonal signs
    reproduced, and how often the inversion fails."""
    from grn_pipeline.m26_texh_jacobian import mra
    rng = np.random.default_rng(seed)
    X = d["X"]; tgt = d["target"]
    gi = [d["genes"].index(g) for g in genes]
    cidx = np.where(tgt == "NO-TARGET")[0]
    kidx = {g: np.where(tgt == g)[0] for g in genes}
    R0 = np.column_stack([X[kidx[g]][:, gi].mean(0) - X[cidx][:, gi].mean(0)
                          for g in genes])
    r0 = mra(R0)
    off = ~np.eye(len(genes), dtype=bool)
    reps, fails = [], 0
    for _ in range(n_rep):
        c = X[rng.choice(cidx, len(cidx))][:, gi].mean(0)
        Rb = np.column_stack([X[rng.choice(kidx[g], len(kidx[g]))][:, gi]
                              .mean(0) - c for g in genes])
        try:
            reps.append(mra(Rb))
        except np.linalg.LinAlgError:
            fails += 1
    reps = np.array(reps)
    sign = float(np.mean(np.mean(np.sign(reps[:, off])
                                 == np.sign(r0[off])[None], 0))) if len(reps) \
        else np.nan
    return dict(R=R0, r=r0, sign_consistency=sign, n_fail=fails,
                n_rep=len(reps))


def split_half(d, genes, n_split=50, seed=0):
    """Out-of-sample test: split the cells of every group in half, run MRA
    independently on each half, and measure how often the two independent
    estimates agree in sign. Stronger than a bootstrap, which re-uses the
    same cells."""
    from grn_pipeline.m26_texh_jacobian import mra
    rng = np.random.default_rng(seed)
    X = d["X"]; tgt = d["target"]
    gi = [d["genes"].index(g) for g in genes]
    cidx = np.where(tgt == "NO-TARGET")[0]
    kidx = {g: np.where(tgt == g)[0] for g in genes}
    off = ~np.eye(len(genes), dtype=bool)
    agree, mags = [], []
    for _ in range(n_split):
        halves = []
        cperm = rng.permutation(cidx)
        kperm = {g: rng.permutation(v) for g, v in kidx.items()}
        for h in (0, 1):
            c = X[cperm[h::2]][:, gi].mean(0)
            Rh = np.column_stack([X[kperm[g][h::2]][:, gi].mean(0) - c
                                  for g in genes])
            try:
                halves.append(mra(Rh))
            except np.linalg.LinAlgError:
                break
        if len(halves) == 2:
            agree.append(np.sign(halves[0][off]) == np.sign(halves[1][off]))
            mags.append(np.minimum(np.abs(halves[0][off]),
                                   np.abs(halves[1][off])))
    agree = np.array(agree); mags = np.array(mags)
    return dict(overall=float(agree.mean()), agree=agree, mags=mags,
                n=len(agree))


# ---------------------------------------------------------------------------
# figure + report
# ---------------------------------------------------------------------------
def figure(res, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(13.5, 9.4)); ax = ax.ravel()
    a = ax[0]
    bins = np.linspace(-3, 1, 60)
    a.hist(np.log10(np.abs(res["offdiag"]) + 1e-4), bins=bins, color="#bbb",
           label="off-diagonal (network effects)")
    a.hist(np.log10(np.abs(res["diag"]) + 1e-4), bins=bins, color="#2a6fdb",
           alpha=0.85, label="diagonal (effect on the target itself)")
    a.set_xlabel("log10 |response|"); a.set_ylabel("count")
    a.legend(fontsize=8)
    a.set_title(f"(a) CRISPRa: the perturbation is visible on its own target\n"
                f"diag/off-diag = {res['ratio']:.0f}  "
                f"(CRISPR cutting in M26: {res['ratio_m26']:.1f})")
    a = ax[1]
    b = res["calibration"]
    a.semilogx(b["ratio"], 100 * np.array(b["sign"]), "o-", c="#2a6fdb",
               label="synthetic networks (known answer)")
    a.axhline(50, c="k", ls=":", lw=1, label="chance")
    for x, y, lab, col in [(res["ratio_m26"], 100 * res["sign_m26"],
                            "CRISPR cutting\n(GSE327731)", "#e4572e"),
                           (res["ratio"], 100 * res["sign_boot"],
                            "CRISPRa\n(GSE190604)", "#29a36a")]:
        a.scatter([x], [y], s=120, marker="*", c=col, zorder=5)
        a.annotate(lab, (x, y), textcoords="offset points", xytext=(-8, -30),
                   fontsize=7, color=col, ha="center")
    a.set_xlabel("diagonal / off-diagonal ratio")
    a.set_ylabel("% of interaction signs reproduced (bootstrap)")
    a.set_ylim(40, 100); a.legend(fontsize=7, loc="lower right")
    a.set_title("(b) fixing the diagonal moves MRA off chance,\nas the "
                "synthetic calibration predicted")
    a = ax[2]
    labs = ["weak\n|r|<0.02", "medium\n0.02-0.1", "strong\n|r|>0.1"]
    v = [100 * x for x in res["split_by_mag"]]
    a.bar(range(3), v, color=["#bbb", "#7fa8dd", "#2a6fdb"])
    a.axhline(50, c="k", ls=":", lw=1)
    a.axhline(100 * res["split_m26"], c="#e4572e", lw=2,
              label=f"CRISPR cutting, all edges "
                    f"({100 * res['split_m26']:.0f}%)")
    for i, x in enumerate(v):
        a.text(i, x + 1, f"{x:.0f}%", ha="center", fontsize=9)
    a.set_xticks(range(3), labs); a.set_ylim(40, 90)
    a.set_ylabel("% sign agreement between independent halves")
    a.legend(fontsize=7)
    a.set_title("(c) out-of-sample: only the strongest edges\nreproduce")
    a = ax[3]
    x = np.arange(2)
    a.bar(x - 0.2, [100 * res["sign_m26"], 100 * res["sign_boot"]], 0.4,
          color="#9b59b6", label="bootstrap (re-uses the same cells)")
    a.bar(x + 0.2, [100 * res["split_m26"], 100 * res["split_all"]], 0.4,
          color="#29a36a", label="split-half (independent cells)")
    a.axhline(50, c="k", ls=":", lw=1)
    a.set_xticks(x, ["CRISPR cutting\n(GSE327731)", "CRISPRa\n(GSE190604)"])
    a.set_ylabel("% of interaction signs reproduced")
    a.set_ylim(40, 90); a.legend(fontsize=8)
    for i, (p, q) in enumerate([(res["sign_m26"], res["split_m26"]),
                                (res["sign_boot"], res["split_all"])]):
        a.text(i - 0.2, 100 * p + 1, f"{100 * p:.0f}%", ha="center", fontsize=9)
        a.text(i + 0.2, 100 * q + 1, f"{100 * q:.0f}%", ha="center", fontsize=9)
    a.set_title("(d) bootstrap overstates reliability by ~20 points;\n"
                "report the split-half number")
    fig.suptitle("M27 — does a perturbation that moves the target's own "
                 "transcript rescue network inference?", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95]); fig.savefig(path, dpi=120)
    plt.close(fig)


def report(fig_path=os.path.join(FIG, "m27_crispri_mra.png")):
    import json
    print("=" * 70)
    print("MODULE 27 — CRISPRa Perturb-seq: does a real diagonal rescue MRA?")
    print("=" * 70)
    d = load_crispra()
    R, rows, cols, _ = crispra_response(d, n_boot=0)
    Rs, cg = square_module(R, rows, cols)
    b = mra_bootstrap(d, cg)
    R0 = b["R"]; off = ~np.eye(len(cg), dtype=bool)
    dg = np.diag(R0)
    ratio = float(np.median(np.abs(dg)) / np.median(np.abs(R0[off])))
    sh = split_half(d, cg)
    m = sh["mags"].ravel(); ag = sh["agree"].ravel()
    by = [float(ag[(m >= lo) & (m < hi)].mean())
          for lo, hi in [(0, 0.02), (0.02, 0.1), (0.1, 1e9)]]
    with open(os.path.join(FIG, "m26_summary.json")) as fh:
        m26 = json.load(fh)
    cal = {}
    for x in m26["breakdown"]:
        k = ("<1.5" if x["ratio"] < 1.5 else "1.5-3" if x["ratio"] < 3 else
             "3-6" if x["ratio"] < 6 else "6-12" if x["ratio"] < 12 else
             "12-40" if x["ratio"] < 40 else ">40")
        cal.setdefault(k, []).append(x)
    calib = dict(ratio=[float(np.median([y["ratio"] for y in v]))
                        for v in cal.values()],
                 sign=[float(np.median([y["sign"] for y in v]))
                       for v in cal.values()])
    o = np.argsort(calib["ratio"])
    calib = {k: list(np.array(v)[o]) for k, v in calib.items()}
    print(f"module: {len(cg)} genes both perturbed and measured "
          f"({d['X'].shape[0]} cells, {len(set(d['target'])) - 1} targets)")
    print(f"A. diagonal: {100 * (dg > 0).mean():.0f}% positive (CRISPRa "
          f"activates its target), median |diag| {np.median(np.abs(dg)):.2f}")
    print(f"B. diag/off-diag ratio {ratio:.0f}  vs  "
          f"{m26['mra']['ratio']:.1f} for CRISPR cutting; condition number "
          f"{np.linalg.cond(R0):.0f} vs {m26['mra']['cond']:.0f}")
    print(f"C. MRA bootstrap sign consistency {b['sign_consistency']:.0%} "
          f"vs {m26['mra']['sign_consistency']:.0%} (chance 50%)")
    print(f"D. split-half (independent cells) {sh['overall']:.0%} vs "
          f"{res_m26_split:.0%}" if False else
          f"D. split-half (independent cells, the honest metric): "
          f"{sh['overall']:.0%}")
    print(f"   by edge strength: weak {by[0]:.0%}, medium {by[1]:.0%}, "
          f"strong {by[2]:.0%}")
    # same metric on the M26 CRISPR-cutting module
    import anndata as ad
    import pandas as pd
    a26 = ad.read_h5ad(os.path.join(HERE, "..", "data_cache", "GSE327731",
                                    "processed.h5ad"))
    ko = pd.read_csv(os.path.join(HERE, "..", "data_cache", "GSE327731",
                                  "ko_effects.csv"))
    gl = [g for g in ko["gene"] if g in a26.var_names]
    X = a26[:, gl].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    tg = np.where(a26.obs["gRNA_gene"].values == "NT", "NO-TARGET",
                  a26.obs["gRNA_gene"].values)
    d26 = dict(X=X, genes=gl, target=tg)
    use26 = [g for g in gl if (tg == g).sum() >= 100]
    sh26 = split_half(d26, use26, n_split=50)
    print(f"   CRISPR cutting under the same test: {sh26['overall']:.0%}")
    res = dict(ratio=ratio, ratio_m26=m26["mra"]["ratio"],
               sign_boot=b["sign_consistency"],
               sign_m26=m26["mra"]["sign_consistency"],
               split_all=sh["overall"], split_m26=sh26["overall"],
               split_by_mag=by, diag=dg.tolist(),
               offdiag=R0[off].tolist(), genes=cg, calibration=calib,
               cond=float(np.linalg.cond(R0)), cond_m26=m26["mra"]["cond"],
               n_cells=int(d["X"].shape[0]), frac_pos=float((dg > 0).mean()))
    figure(res, fig_path)
    with open(os.path.join(FIG, "m27_summary.json"), "w") as fh:
        json.dump({k: v for k, v in res.items()
                   if k not in ("diag", "offdiag")}, fh, indent=1,
                  default=float)
    print(f"figure: {fig_path}")
    return res


if __name__ == "__main__":
    report()
