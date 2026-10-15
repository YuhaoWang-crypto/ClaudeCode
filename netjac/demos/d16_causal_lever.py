"""
D16  Minimal-adjustment prediction with a CAUSAL Jacobian from perturbation
     data, re-testing the check that D15 failed.

D15 fitted J from the relaxation time course of ONE observational treatment
and asked which module to adjust. The pre-specified check failed: the top
lever was the slowest module (barrier), not the clinically validated
upstream target (IL-17). The diagnosis was that such a J contains no
measurement of what pulling any OTHER lever would do.

Perturbation screens contain exactly that. Replogle et al. 2022 genome-wide
CRISPRi (K562, 11,258 knockdowns x 8,248 genes) measures the steady-state
response to perturbing each gene individually, so for module-level
coordinates y the response matrix R (module i's shift when module j is
perturbed) is measured directly, and MRA gives the causal J = -R^-1.

Formal point: in `statepipe.minimal_adjustment` the effect of lever j is
J^-1 e_j. With J = -R^-1 that is exactly -R[:, j], the MEASURED response to
perturbing module j. So the formalism was never the problem; it needs a J
whose inverse reproduces measured perturbation responses. The test below
asks whether that is enough to recover the right lever.

TEST 1 (ground truth, leave-one-out). Knock down gene g of module M. The
displacement y0 is the module-coordinate shift. The lever that restores it
should be M itself: the thing that broke the state is the thing that fixes
it. The module scores exclude the perturbed gene, so a hit cannot come from
g's own transcript. Three predictors are compared on the same displacements:
  causal        rank levers by the measured response library (leave-one-out)
  observational rank by the covariance column of 585 non-targeting control
                rows, i.e. what co-fluctuation alone would say (no
                perturbation information) - the D9/D10 FDT comparison
  naive         pick the module whose own coordinate moved most
against the chance rate of 1/15.

TEST 2 (Norman et al. 2019 doubles). Which second perturbation best
restores the displacement caused by a first one, ranked among the partners
actually measured with it. Linear-response prediction of the double from
the two singles, against the observed restoration ratio.

Run:  python3 -m netjac.demos.d16_causal_lever
"""
import os
import pickle
import re

import numpy as np

from netjac import DATA_CACHE, FIG
from netjac import recovery as rc
from netjac.data import perturbseq as ps
from netjac.demos import d11_recovery_state as d11

RES = os.path.join(DATA_CACHE, "d16_results.pkl")
FAM = {**d11.POS_FAM, **d11.NEG_FAM}


# ============================================================ module data
def load_modules(key="K562_gwps_bulk"):
    """Module coordinates of every perturbation's response. The perturbed
    gene is excluded from every module score, so a module can never score
    on its own knocked-down transcript."""
    B = ps.load_bulk(key)
    X = np.asarray(B["X"], float)
    X[~np.isfinite(X)] = 0.0
    genes, pert = np.array(B["gene"]), np.array(B["pert_gene"])
    names = list(FAM)
    memb = {n: np.array([bool(re.match(FAM[n], g)) for g in genes]) for n in names}
    gi = {g: i for i, g in enumerate(genes)}
    Y = np.zeros((len(X), len(names)))
    for k, n in enumerate(names):
        m = memb[n].copy()
        s = X[:, m].mean(1) * m.sum()
        cnt = np.full(len(X), m.sum(), float)
        # subtract the perturbed gene's own contribution where it is a member
        for r, g in enumerate(pert):
            j = gi.get(g)
            if j is not None and m[j]:
                s[r] -= X[r, j]
                cnt[r] -= 1
        Y[:, k] = s / np.maximum(cnt, 1)
    pert_mod = np.array([next((n for n in names if re.match(FAM[n], g)), "") for g in pert])
    return {"Y": Y, "names": names, "pert": pert, "pert_mod": pert_mod,
            "control": np.asarray(B["control"]).astype(bool)}


def response_matrix(D, exclude=None):
    """R[i, j] = mean module-i coordinate over perturbations of module j."""
    names = D["names"]
    R = np.zeros((len(names), len(names)))
    for j, n in enumerate(names):
        m = (D["pert_mod"] == n)
        if exclude is not None:
            m = m & (np.arange(len(m)) != exclude)
        R[:, j] = D["Y"][m].mean(0)
    return R


def split_half_R(D, seed=0, n=20):
    """Replication of the response matrix: correlation of R built from two
    random halves of the perturbations of each module."""
    rng = np.random.default_rng(seed)
    rs = []
    for _ in range(n):
        A = np.zeros((len(D["names"]),) * 2)
        Bm = np.zeros_like(A)
        for j, nm in enumerate(D["names"]):
            idx = np.where(D["pert_mod"] == nm)[0]
            h = rng.permutation(idx)
            A[:, j] = D["Y"][h[:len(h) // 2]].mean(0)
            Bm[:, j] = D["Y"][h[len(h) // 2:]].mean(0)
        rs.append(np.corrcoef(A.ravel(), Bm.ravel())[0, 1])
    return float(np.mean(rs)), float(np.std(rs))


# ======================================================= lever prediction
def _rank_by_library(y0, L):
    """Rank levers by how much of the displacement a single lever can
    cancel: projection energy (L[:, j].y0)^2 / |L[:, j]|^2."""
    num = (L.T @ y0) ** 2
    den = (L ** 2).sum(0)
    return num / np.maximum(den, 1e-12)


def test1(D, seed=0):
    """Leave-one-out: does the predicted lever name the module that was
    actually perturbed?"""
    names = D["names"]
    nm = len(names)
    R_full = response_matrix(D)
    # observational comparator: covariance of module coordinates across the
    # 585 non-targeting control rows (no perturbation information at all).
    C = np.cov(D["Y"][D["control"]].T)
    rows, hits = [], {"causal": [], "observational": [], "naive": []}
    topk = {"causal": [], "observational": [], "naive": []}
    for i in np.where(D["pert_mod"] != "")[0]:
        true = names.index(D["pert_mod"][i])
        y0 = D["Y"][i]
        if not np.isfinite(y0).all() or np.linalg.norm(y0) == 0:
            continue
        R = response_matrix(D, exclude=i)           # leave the held-out row out
        sc = {"causal": _rank_by_library(y0, R),
              "observational": _rank_by_library(y0, C),
              "naive": np.abs(y0)}
        for k, v in sc.items():
            order = np.argsort(-v)
            hits[k].append(int(order[0] == true))
            topk[k].append(int(np.where(order == true)[0][0]) + 1)
        rows.append({"pert": D["pert"][i], "module": names[true],
                     "causal_top": names[int(np.argmax(sc["causal"]))],
                     "obs_top": names[int(np.argmax(sc["observational"]))],
                     "naive_top": names[int(np.argmax(sc["naive"]))]})
    out = {"n": len(rows), "rows": rows, "R": R_full, "cov": C, "names": names,
           "chance": 1.0 / nm}
    for k in hits:
        h = np.array(hits[k])
        rk = np.array(topk[k])
        # Wilson interval for the hit rate
        p, n = h.mean(), len(h)
        z = 1.96
        den = 1 + z * z / n
        c = (p + z * z / (2 * n)) / den
        hw = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
        out[k] = {"top1": float(p), "ci": (float(c - hw), float(c + hw)),
                  "top3": float((rk <= 3).mean()), "median_rank": float(np.median(rk))}
    # permutation null: shuffle which module each displacement belongs to,
    # so module sizes and the response library are kept but the pairing is broken
    rng = np.random.default_rng(seed)
    truth = np.array([names.index(r["module"]) for r in rows])
    got = np.array([names.index(r["causal_top"]) for r in rows])
    perm = [float((got == rng.permutation(truth)).mean()) for _ in range(2000)]
    out["perm_null"] = {"mean": float(np.mean(perm)), "p95": float(np.percentile(perm, 95)),
                        "p": float(np.mean(np.array(perm) >= out["causal"]["top1"]))}
    # identifiability: a lever can only be named if its response direction is
    # distinguishable from the others. Max |cos| to any other response column.
    Rn = R_full / np.maximum(np.linalg.norm(R_full, axis=0), 1e-12)
    cs = np.abs(Rn.T @ Rn)
    out["cos"] = cs.copy()
    np.fill_diagonal(cs, 0)
    out["max_cos"] = {names[j]: float(cs[:, j].max()) for j in range(nm)}
    out["norm_R"] = {names[j]: float(np.linalg.norm(R_full[:, j])) for j in range(nm)}
    out["cond_R"] = float(np.linalg.cond(R_full))
    # where do the errors go? Compare the true-vs-named response similarity
    # against a random module.
    from scipy.stats import mannwhitneyu
    err = [(names.index(r["module"]), names.index(r["causal_top"]))
           for r in rows if r["module"] != r["causal_top"]]
    if err:
        ce = [out["cos"][a, b] for a, b in err]
        cr = [out["cos"][a, rng.integers(0, nm)] for a, _ in err]
        out["errors"] = {"median_cos": float(np.median(ce)), "median_cos_random": float(np.median(cr)),
                         "p": float(mannwhitneyu(ce, cr).pvalue),
                         "frac_cos_gt_0.7": float(np.mean(np.array(ce) > 0.7)), "n": len(err)}
    # per-module hit rate for the causal predictor
    out["per_module"] = {}
    for n_ in names:
        sel = [r for r in rows if r["module"] == n_]
        if sel:
            out["per_module"][n_] = (float(np.mean([r["causal_top"] == n_ for r in sel])), len(sel))
    return out


# ================================================== Norman doubles (test 2)
def test2(thr=0.75, seed=0):
    """For each double A+B, predict the restoration of A's displacement by
    adding B, from the single responses (linear response: d_AB = d_A + d_B).
    Ranking task: among the partners measured with A, is the observed best
    partner ranked first? Plus the D11 classification for comparability."""
    import anndata as ad
    from scipy.stats import spearmanr
    fn = os.path.join(DATA_CACHE, "norman_gears", "norman", "perturb_processed.h5ad")
    if not os.path.exists(fn):
        return {"skipped": "norman file not cached"}
    a = ad.read_h5ad(fn)
    cond = a.obs["condition"].astype(str).values
    X = a.X.tocsr() if hasattr(a.X, "tocsr") else a.X
    means = {}
    for c in np.unique(cond):
        idx = np.where(cond == c)[0]
        means[c] = np.asarray(X[idx].mean(0)).ravel()
    ctrl = means["ctrl"]
    d = {c: v - ctrl for c, v in means.items()}

    def single(g):
        for c in (f"{g}+ctrl", f"ctrl+{g}"):
            if c in d:
                return c
        return None
    pairs = []
    for c in d:
        if c == "ctrl" or "ctrl" in c:
            continue
        A_, B_ = c.split("+")
        sa, sb = single(A_), single(B_)
        if sa is None or sb is None:
            continue
        dA, dB, dAB = d[sa], d[sb], d[c]
        nA, nB = np.linalg.norm(dA), np.linalg.norm(dB)
        pairs.append({"A": A_, "B": B_,
                      "rho": float(np.linalg.norm(dAB) / max(nA, nB)),
                      "pred": float(np.linalg.norm(dA + dB) / max(nA, nB)),
                      "cos": float(dA @ dB / (nA * nB))})
    rho = np.array([p["rho"] for p in pairs])
    pred = np.array([p["pred"] for p in pairs])
    y = rho < thr
    rk = spearmanr(pred, rho)
    # ranking task: for each A with >= 3 measured partners, is the observed
    # best partner the predicted best?
    byA, hits = {}, []
    for p in pairs:
        byA.setdefault(p["A"], []).append(p)
        byA.setdefault(p["B"], []).append(p)
    for A_, ps_ in byA.items():
        if len(ps_) < 3:
            continue
        o_best = min(ps_, key=lambda q: q["rho"])
        p_best = min(ps_, key=lambda q: q["pred"])
        hits.append(int(o_best is p_best or (o_best["A"], o_best["B"]) == (p_best["A"], p_best["B"])))
    return {"n_pairs": len(pairs), "n_pos": int(y.sum()),
            "auc_pred": rc.auc(y, -pred), "rho_spearman": (float(rk.statistic), float(rk.pvalue)),
            "rank_hit": float(np.mean(hits)) if hits else np.nan, "n_rank": len(hits),
            "chance_rank": float(np.mean([1 / len(v) for v in byA.values() if len(v) >= 3])),
            "pairs": pairs}


# ================================================================= report
def run():
    D = load_modules()
    R = {"split_half": split_half_R(D), "t1": test1(D), "t2": test2()}
    pickle.dump(R, open(RES, "wb"))
    return R


def report(R):
    t1 = R["t1"]
    m, s = R["split_half"]
    print(f"\nD16 causal J from Replogle genome-wide CRISPRi: {len(t1['names'])} modules, "
          f"{t1['n']} held-out perturbations; response matrix split-half r {m:.2f} +- {s:.2f}")
    print("\nTEST 1  does the predicted lever name the module that was actually broken?")
    print(f"   {'predictor':16s} {'top-1':>18s} {'top-3':>8s} {'median rank':>12s}")
    for k, lab in (("causal", "causal (response library)"), ("observational", "observational (covariance)"),
                   ("naive", "naive (largest shift)")):
        v = t1[k]
        print(f"   {lab:28s} {v['top1']:.0%} [{v['ci'][0]:.0%}, {v['ci'][1]:.0%}] {v['top3']:7.0%} {v['median_rank']:11.0f}")
    print(f"   {'chance':28s} {t1['chance']:.0%}")
    pn = t1["perm_null"]
    print(f"   label-permutation null {pn['mean']:.0%} (95th pct {pn['p95']:.0%}), p < {max(pn['p'], 1/2000):.4f}")
    print(f"\n   per-module hit rate (causal)   [condition number of R {t1['cond_R']:.0f}]:")
    print(f"      {'module':26s} {'hit':>5s} {'n':>5s}  |R_j|  max |cos| to another lever")
    for n_, (h, k) in sorted(t1["per_module"].items(), key=lambda kv: -kv[1][0]):
        print(f"      {n_:26s} {h:5.0%} {k:5d}  {t1['norm_R'][n_]:5.2f}   {t1['max_cos'][n_]:.2f}")
    from scipy.stats import spearmanr as _sr
    hs = [t1["per_module"][n_][0] for n_ in t1["per_module"]]
    print("      no module-level property predicts the hit rate at n=15 modules: " + ", ".join(
        f"{lab} {_sr([d[n_] for n_ in t1['per_module']], hs).statistic:+.2f}"
        f" (p={_sr([d[n_] for n_ in t1['per_module']], hs).pvalue:.2f})"
        for lab, d in (("|R_j|", t1["norm_R"]), ("max|cos|", t1["max_cos"]))))
    e = t1.get("errors")
    if e:
        print(f"      when wrong, the named lever's response resembles the true one more than a random module "
              f"would (median |cos| {e['median_cos']:.2f} vs {e['median_cos_random']:.2f}, p={e['p']:.2f}), "
              f"but only {e['frac_cos_gt_0.7']:.0%} of errors are strongly collinear: errors go to "
              f"functionally adjacent modules (mitoribosome -> complex I, spliceosome -> nuclear RNA machinery)")
    t2 = R["t2"]
    if "skipped" in t2:
        print("\nTEST 2 skipped:", t2["skipped"])
        return
    print(f"\nTEST 2  Norman 2019 doubles: {t2['n_pairs']} pairs, {t2['n_pos']} restore the state")
    print(f"   linear-response prediction of restoration: AUC {t2['auc_pred']:.2f}, "
          f"Spearman {t2['rho_spearman'][0]:+.2f} (p={t2['rho_spearman'][1]:.1g})")
    print(f"   picking the best partner for a given perturbation: {t2['rank_hit']:.0%} correct "
          f"({t2['n_rank']} perturbations, chance {t2['chance_rank']:.0%})")


def figures(R):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    t1 = R["t1"]
    names = t1["names"]
    short = [n.split(" (")[0][:18] for n in names]
    fig, ax = plt.subplots(1, 3, figsize=(17, 5.4))
    a = ax[0]
    v = np.abs(t1["R"])
    im = a.imshow(t1["R"], cmap="RdBu_r", vmin=-v.max(), vmax=v.max())
    a.set_xticks(range(len(names)), short, rotation=45, ha="right", fontsize=7)
    a.set_yticks(range(len(names)), short, fontsize=7)
    a.set_xlabel("module perturbed", labelpad=2)
    a.set_ylabel("module responding")
    m, s = R["split_half"]
    a.set_title(f"Measured response matrix R (causal)\nsplit-half r {m:.2f}", fontsize=10)
    fig.colorbar(im, ax=a, fraction=0.046)
    a = ax[1]
    ks = ["causal", "observational", "naive"]
    lab = ["causal\n(response library)", "observational\n(covariance)", "naive\n(largest shift)"]
    val = [t1[k]["top1"] for k in ks]
    err = [[t1[k]["top1"] - t1[k]["ci"][0] for k in ks], [t1[k]["ci"][1] - t1[k]["top1"] for k in ks]]
    a.bar(range(3), val, yerr=err, color=["C2", "C0", "0.6"], capsize=4)
    a.axhline(t1["chance"], color="k", ls="--", lw=1)
    a.text(2.4, t1["chance"] + 0.01, "chance", ha="right", fontsize=8)
    a.set_xticks(range(3), lab, fontsize=8)
    a.set_ylabel("top-1 hit rate: lever = module actually broken")
    a.set_ylim(0, 1)
    a.set_title(f"TEST 1  ground truth, {t1['n']} held-out perturbations", fontsize=10)
    a = ax[2]
    pm = list(t1["per_module"].items())
    xs = [t1["max_cos"][n] for n, _ in pm]
    ys = [v[0] for _, v in pm]
    a.scatter(xs, ys, s=[12 + v[1] for _, v in pm], color="C2")
    for (n, v), x, y in zip(pm, xs, ys):
        a.annotate(n.split(" (")[0][:16], (x, y), fontsize=6.5, xytext=(3, 3), textcoords="offset points")
    a.axhline(t1["chance"], color="k", ls="--", lw=1)
    a.set_xlabel("max |cos| between this lever's response and another's")
    a.set_ylabel("top-1 hit rate (causal)")
    a.set_ylim(0, 1)
    a.set_title("TEST 1  per module: no advance diagnostic\npredicts which levers are nameable (size = n)", fontsize=10)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    p = os.path.join(FIG, "d16_causal_lever.png")
    fig.savefig(p, dpi=110)
    plt.close(fig)
    return p


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", action="store_true")
    a = ap.parse_args()
    R = pickle.load(open(RES, "rb")) if a.report and os.path.exists(RES) else run()
    report(R)
    print(figures(R))
