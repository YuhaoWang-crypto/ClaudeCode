"""Benchmark: gene-pair code panel vs NNLS scoring (numpy + matplotlib).

    python -m cell_state_library.benchmark

All targets are simulated from the same aggregate author references plus noise,
so this measures robustness of the scoring methods under stated assumptions.
It is NOT biological validation: no independent donors or studies are involved.

Task: given a target log2FC profile of a known cell type, call it
"none" (no IFN response), "IFNA" or "IFNG". Only stimuli that pass the
identifiability check count toward balanced accuracy; IFNG targets in
IFNa-only subsets are reported separately.

Methods (same decision protocol for all):
  nnls_full   v0.1/author NNLS on every reference gene of the cell type
  pair16      8-pair ternary code (16 genes)
  nnls16      NNLS restricted to the same 16 panel genes (gate bypassed)
The "signal vs none" threshold of every method is calibrated on pure-noise
targets to a 5% false-positive rate under the undistorted noise model; then
the class is the best-scoring identifiable stimulus.
"""
import json
import math
import sqlite3
from pathlib import Path

import numpy as np

from . import build, identifiability as ident, panel_select
from .pair_code import STRENGTH_GRID

ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures"
N = 300            # simulated targets per class per configuration
FPR = 0.05
SEED = 20261020


# ---------------------------------------------------------------- data
def load_db(db=build.DB):
    con = sqlite3.connect(db)
    out = {}
    for ct in build.ACTIVE:
        rows = con.execute(
            "SELECT v.gene, p.condition_id, v.log2fc FROM reference_values v JOIN reference_profiles p USING(profile_id) "
            "WHERE p.cell_type_id=?", (ct,)).fetchall()
        genes = sorted({r[0] for r in rows})
        gi = {g: i for i, g in enumerate(genes)}
        a, g = np.zeros(len(genes)), np.zeros(len(genes))
        for gene, cond, v in rows:
            (a if cond == "IFNA" else g)[gi[gene]] = v
        status = con.execute("SELECT status FROM identifiability WHERE cell_type_id=?", (ct,)).fetchone()[0]
        panel = [dict(zip(("gene_i", "gene_j", "d_ifna", "d_ifng"), r)) for r in con.execute(
            "SELECT gene_i, gene_j, d_ifna, d_ifng FROM pair_panels WHERE cell_type_id=? ORDER BY rank", (ct,))]
        out[ct] = {"genes": genes, "gi": gi, "a": a, "g": g, "status": status,
                   "classes": ident.identifiable_stimuli(status), "panel": panel}
    con.close()
    return out


# ---------------------------------------------------------------- methods
def nnls2(B, a, g):
    """Vectorised v0.1 two-reference NNLS on unit-normalised columns; NaN = missing gene.

    Returns (x_a, x_g, relative residual) per row of B.
    """
    M = ~np.isnan(B)
    Bz = np.where(M, B, 0.0)
    aa, gg = (M * a ** 2).sum(1), (M * g ** 2).sum(1)
    ag, ab, gb = (M * a * g).sum(1), (Bz * a).sum(1), (Bz * g).sum(1)
    bb = (Bz ** 2).sum(1)
    na, ng, nb = np.sqrt(aa), np.sqrt(gg), np.sqrt(np.maximum(bb, 1e-300))
    c, p, q = ag / (na * ng), ab / (na * nb), gb / (ng * nb)
    det = np.maximum(1 - c * c, 1e-12)
    xa_f, xg_f = (p - c * q) / det, (q - c * p) / det
    cand = [(np.zeros_like(p), np.zeros_like(p)), (np.maximum(p, 0), np.zeros_like(p)), (np.zeros_like(p), np.maximum(q, 0))]
    feas = (xa_f >= 0) & (xg_f >= 0)
    cand.append((np.where(feas, xa_f, 0), np.where(feas, xg_f, 0)))

    def loss(xa, xg):
        return 1 - 2 * (xa * p + xg * q) + xa ** 2 + xg ** 2 + 2 * c * xa * xg
    losses = np.stack([np.where((k == 3) & ~feas, np.inf, loss(*x)) for k, x in enumerate(cand)])
    best = losses.argmin(0)
    xa = np.choose(best, [x[0] for x in cand])
    xg = np.choose(best, [x[1] for x in cand])
    return xa, xg, np.sqrt(np.maximum(losses.min(0), 0))


def pair_scores(B, panel, gi, sigma, tau, classes):
    ll_none = panel_select.loglik_np(panel, B, gi, sigma, tau, "d_ifna", (0.0,))
    ll = {c: panel_select.loglik_np(panel, B, gi, sigma, tau, {"IFNA": "d_ifna", "IFNG": "d_ifng"}[c], STRENGTH_GRID)
          for c in classes}
    return ll_none, ll


def decide_nnls(B, a, g, classes, threshold):
    xa, xg, res = nnls2(B, a, g)
    signal = res < threshold
    if classes == ["IFNA", "IFNG"]:
        lab = np.where(xa >= xg, "IFNA", "IFNG")
    else:
        lab = np.full(len(res), classes[0] if classes else "none")
    return np.where(signal & bool(classes), lab, "none"), res


def decide_pair(B, panel, gi, sigma, tau, classes, threshold):
    ll_none, ll = pair_scores(B, panel, gi, sigma, tau, classes)
    if not ll:
        return np.full(B.shape[0], "none"), np.zeros(B.shape[0])
    stack = np.stack([ll[c] for c in classes])
    llr = stack.max(0) - ll_none
    lab = np.array(classes)[stack.argmax(0)]
    return np.where(llr > threshold, lab, "none"), llr


# ---------------------------------------------------------------- simulation
def simulate(rng, ref, n, sigma, strength, shift=0.0, saturate=False, dropout=0.0):
    signal = strength * ref
    if saturate:
        signal = 2.0 * np.tanh(signal / 2.0)   # monotone compression of large changes (platform ceiling)
    B = signal[None, :] + shift + rng.normal(0, sigma, (n, len(ref)))
    if dropout:
        B[rng.random(B.shape) < dropout] = np.nan
    return B


def run_config(d, rng, sigma, strength, panel=None, shift=0.0, saturate=False, dropout=0.0, methods=None):
    """Accuracy per true class for each method on one cell type."""
    a, g, gi, classes = d["a"], d["g"], d["gi"], d["classes"]
    panel = d["panel"] if panel is None else panel
    tau = sigma * math.sqrt(2)
    pidx = sorted({gi[p["gene_i"]] for p in panel} | {gi[p["gene_j"]] for p in panel})
    methods = methods or ["nnls_full", "pair16", "nnls16"]

    def cols(B, m):
        return B[:, pidx] if m == "nnls16" else B

    def ref_cols(v, m):
        return v[pidx] if m == "nnls16" else v

    # thresholds from the undistorted null (analyst does not know about the distortion)
    null = simulate(rng, np.zeros_like(a), 2000, sigma, 0.0, dropout=dropout)
    thr = {}
    for m in methods:
        if m == "pair16":
            _, s = decide_pair(null, panel, gi, sigma, tau, classes, np.inf)
            thr[m] = np.quantile(s, 1 - FPR)
        else:
            _, s = decide_nnls(cols(null, m), ref_cols(a, m), ref_cols(g, m), classes, -np.inf)
            thr[m] = np.quantile(s, FPR)
    out = {m: {} for m in methods}
    for truth, ref in (("none", np.zeros_like(a)), ("IFNA", a), ("IFNG", g)):
        B = simulate(rng, ref, N, sigma, 0.0 if truth == "none" else strength, shift, saturate, dropout)
        for m in methods:
            if m == "pair16":
                lab, _ = decide_pair(B, panel, gi, sigma, tau, classes, thr[m])
            else:
                lab, _ = decide_nnls(cols(B, m), ref_cols(a, m), ref_cols(g, m), classes, thr[m])
            out[m][truth] = {"acc": float(np.mean(lab == truth)),
                             "called": {k: float(np.mean(lab == k)) for k in ("none", "IFNA", "IFNG")}}
    for m in methods:
        scored = ["none"] + classes
        out[m]["balanced_acc"] = float(np.mean([out[m][c]["acc"] for c in scored]))
    return out


def summarize(per_ct):
    methods = next(iter(per_ct.values())).keys()
    return {m: float(np.mean([r[m]["balanced_acc"] for r in per_ct.values()])) for m in methods}


# ---------------------------------------------------------------- transfer (leave one subset out)
def transfer(data, rng, sigma=0.5, strength=0.5, k=6):
    """Universal panel on genes shared by all 13 subsets, learned without the held-out subset."""
    common = sorted(set.intersection(*[set(d["genes"]) for d in data.values()]))
    ci = {g: i for i, g in enumerate(common)}
    tau = sigma * math.sqrt(2)
    res = {}
    for ct, d in data.items():
        others = [o for o in data if o != ct]
        A = np.array([[data[o]["a"][data[o]["gi"][x]] for x in common] for o in others])
        G = np.array([[data[o]["g"][data[o]["gi"][x]] for x in common] for o in others])
        a_med, g_med = np.median(A, 0), np.median(G, 0)
        a_own = np.array([d["a"][d["gi"][x]] for x in common])
        g_own = np.array([d["g"][d["gi"][x]] for x in common])
        panel = panel_select.select_panel(common, a_med, g_med, ["IFNA", "IFNG"], k, sigma, tau)
        classes = d["classes"]
        null = simulate(rng, np.zeros(len(common)), 2000, sigma, 0.0)
        _, s = decide_pair(null, panel, ci, sigma, tau, classes, np.inf)
        thr_p = np.quantile(s, 1 - FPR)
        _, s = decide_nnls(null, a_med, g_med, classes, -np.inf)
        thr_t = np.quantile(s, FPR)
        _, s = decide_nnls(null, a_own, g_own, classes, -np.inf)
        thr_o = np.quantile(s, FPR)
        acc = {"pair_transfer": [], "nnls_transfer": [], "nnls_own_ref": []}
        for truth, ref in (("none", None), ("IFNA", a_own), ("IFNG", g_own)):
            if truth != "none" and truth not in classes:
                continue
            B = simulate(rng, np.zeros(len(common)) if ref is None else ref, N, sigma, 0.0 if ref is None else strength)
            acc["pair_transfer"].append(np.mean(decide_pair(B, panel, ci, sigma, tau, classes, thr_p)[0] == truth))
            acc["nnls_transfer"].append(np.mean(decide_nnls(B, a_med, g_med, classes, thr_t)[0] == truth))
            acc["nnls_own_ref"].append(np.mean(decide_nnls(B, a_own, g_own, classes, thr_o)[0] == truth))
        res[ct] = {m: float(np.mean(v)) for m, v in acc.items()}
        res[ct]["panel_genes"] = 2 * len(panel)
    return {"common_genes": len(common), "pairs": k, "per_cell_type": res,
            "mean": {m: float(np.mean([r[m] for r in res.values()])) for m in ("pair_transfer", "nnls_transfer", "nnls_own_ref")}}


# ---------------------------------------------------------------- main
def main():
    rng = np.random.default_rng(SEED)
    data = load_db()
    results = {"protocol": __doc__.strip(), "n_per_class": N, "fpr_calibration": FPR, "seed": SEED}

    base = {}
    for sigma in (0.25, 0.5, 1.0):
        for strength in (1.0, 0.5, 0.25):
            per_ct = {ct: run_config(d, rng, sigma, strength) for ct, d in data.items()}
            base[f"sigma={sigma},strength={strength}"] = {"mean_balanced_acc": summarize(per_ct), "per_cell_type": per_ct}
    results["noise_strength_grid"] = base

    robust = {}
    for name, kw in (("clean", {}), ("global_shift_+0.5", {"shift": 0.5}), ("saturation_2tanh(x/2)", {"saturate": True}),
                     ("dropout_30pct", {"dropout": 0.3})):
        per_ct = {ct: run_config(d, rng, 0.5, 0.5, **kw) for ct, d in data.items()}
        robust[name] = {"mean_balanced_acc": summarize(per_ct), "per_cell_type": per_ct}
    results["robustness_sigma0.5_strength0.5"] = robust

    size = {}
    for k in (1, 2, 4, 8, 16):
        per_ct = {}
        for ct, d in data.items():
            panel = panel_select.select_panel(d["genes"], d["a"], d["g"], d["classes"], k, 0.5, 0.5 * math.sqrt(2))
            r = run_config(d, rng, 0.5, 0.5, panel=panel, methods=["pair16", "nnls16"])
            per_ct[ct] = r
        size[2 * k] = summarize(per_ct)
    results["panel_size_sigma0.5_strength0.5"] = {str(k): {"pair": v["pair16"], "nnls_same_genes": v["nnls16"]} for k, v in size.items()}

    results["leave_one_subset_out_transfer"] = transfer(data, rng)

    # IFNG targets in IFNa-only subsets: what do the methods call them?
    weak = {}
    key = "sigma=0.5,strength=1.0"
    for ct, d in data.items():
        if "IFNG" not in d["classes"]:
            weak[ct] = {m: base[key]["per_cell_type"][ct][m]["IFNG"]["called"] for m in ("nnls_full", "pair16")}
    results["ifng_targets_in_ifna_only_subsets_sigma0.5_strength1"] = weak

    (ROOT / "reports" / "benchmark.json").write_text(json.dumps(results, indent=2))
    figures(data, results)
    print_summary(results)


def print_summary(r):
    print("mean balanced accuracy over 13 subsets (identifiable classes only)")
    for k, v in r["noise_strength_grid"].items():
        print(f"  {k:26s} " + "  ".join(f"{m}={x:.3f}" for m, x in v["mean_balanced_acc"].items()))
    print("robustness (sigma=0.5, strength=0.5)")
    for k, v in r["robustness_sigma0.5_strength0.5"].items():
        print(f"  {k:26s} " + "  ".join(f"{m}={x:.3f}" for m, x in v["mean_balanced_acc"].items()))
    print("panel size (genes): pair vs nnls on same genes")
    for k, v in r["panel_size_sigma0.5_strength0.5"].items():
        print(f"  {k:>3s} genes  pair={v['pair']:.3f}  nnls={v['nnls_same_genes']:.3f}")
    t = r["leave_one_subset_out_transfer"]
    print(f"LOSO transfer on {t['common_genes']} shared genes: " + "  ".join(f"{m}={x:.3f}" for m, x in t["mean"].items()))


def figures(data, r):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    FIG.mkdir(exist_ok=True)
    cts = list(data)
    ident_rows = json.loads((ROOT / "reports" / "identifiability.json").read_text())["cell_types"]

    fig, ax = plt.subplots(figsize=(10, 4))
    x = np.arange(len(cts))
    ax.bar(x - 0.2, [ident_rows[c]["z_ifna"] for c in cts], 0.4, label="IFN-α", color="#3b75af")
    ax.bar(x + 0.2, [ident_rows[c]["z_ifng"] for c in cts], 0.4, label="IFN-γ", color="#d9822b")
    ax.axhline(ident.Z_MIN, color="k", ls="--", lw=1, label=f"gate z={ident.Z_MIN:g}")
    for i, c in enumerate(cts):
        if ident_rows[c]["author_ifng_zero_rule"]:
            ax.text(i + 0.2, ident_rows[c]["z_ifng"] * 1.08, "A", ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax.set_yscale("log")
    ax.set_xticks(x, cts, rotation=45, ha="right")
    ax.set_ylabel("z of a full-strength response (σ=0.5)")
    ax.set_title("Identifiability per subset  (A = author code zeroes IFN-γ)")
    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1))
    fig.tight_layout()
    fig.savefig(FIG / "fig1_identifiability.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 4))
    per = r["noise_strength_grid"]["sigma=0.5,strength=0.5"]["per_cell_type"]
    for k, (m, col) in enumerate((("nnls_full", "#555555"), ("pair16", "#2a9d8f"), ("nnls16", "#e9c46a"))):
        ax.bar(x + (k - 1) * 0.27, [per[c][m]["balanced_acc"] for c in cts], 0.27, label=m, color=col)
    ax.set_xticks(x, cts, rotation=45, ha="right")
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("balanced accuracy")
    ax.set_title("σ=0.5, half-strength response: full-profile NNLS vs 16-gene pair code vs 16-gene NNLS")
    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1))
    fig.tight_layout()
    fig.savefig(FIG / "fig2_benchmark_per_subset.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 3.6))
    rob = r["robustness_sigma0.5_strength0.5"]
    names = list(rob)
    xx = np.arange(len(names))
    for k, (m, col) in enumerate((("nnls_full", "#555555"), ("pair16", "#2a9d8f"), ("nnls16", "#e9c46a"))):
        ax.bar(xx + (k - 1) * 0.27, [rob[n]["mean_balanced_acc"][m] for n in names], 0.27, label=m, color=col)
    ax.set_xticks(xx, names, fontsize=8)
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("mean balanced accuracy (13 subsets)")
    ax.set_title("Robustness to measurement distortions")
    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1))
    fig.tight_layout()
    fig.savefig(FIG / "fig3_robustness.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.5, 3.6))
    ps = r["panel_size_sigma0.5_strength0.5"]
    ks = [int(k) for k in ps]
    ax.plot(ks, [ps[str(k)]["pair"] for k in ks], "o-", color="#2a9d8f", label="pair code")
    ax.plot(ks, [ps[str(k)]["nnls_same_genes"] for k in ks], "s-", color="#e9c46a", label="NNLS, same genes")
    ref = r["noise_strength_grid"]["sigma=0.5,strength=0.5"]["mean_balanced_acc"]["nnls_full"]
    ax.axhline(ref, color="#555555", ls="--", lw=1, label="NNLS, all reference genes")
    ax.set_xscale("log", base=2)
    ax.set_xticks(ks, [str(k) for k in ks])
    ax.set_xlabel("genes measured")
    ax.set_ylabel("mean balanced accuracy")
    ax.set_ylim(0, 1.02)
    ax.set_title("How small can the panel be? (σ=0.5, half strength)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "fig4_panel_size.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
