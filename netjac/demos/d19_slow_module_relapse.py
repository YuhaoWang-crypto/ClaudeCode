"""
D19  Do unrestored SLOW modules mark relapse risk?

D14 found an ordering: under treatment the inflammatory and proliferative
modules restore first and the structural modules (barrier / lipid in skin,
epithelial metabolism in gut) restore last. That generated a falsifiable
prediction - the tissue that still carries an unrestored structural module
is the tissue that relapses - which nothing so far has tested.

PART A (main)  GSE30768, Suarez-Farinas / Krueger (PMID 22348003).
Four psoriasis patients on efalizumab (anti-CD11a), each sampled at
  LS       lesion before treatment
  Week 12  on treatment, clinically clear
  FLARE    relapse after the drug was stopped
plus non-lesional skin (NL) in two of them. All four relapsed, so the
original design of the test (relapser vs non-relapser) is NOT available here and
is not attempted. What this design does allow, within patient:
  Q1  at week 12, with the skin clinically clear, which modules have
      actually returned to the non-lesional state and which have not?
  Q2  does the flare return to the SAME point on the disease axis as the
      original lesion, i.e. is the lesion state an attractor?
  Q3  the sharp one: is what remains at week 12 already pointing along the
      direction the flare will take? cos(residual at wk12, rise to flare).
      Note the sign bias is conservative: the week-12 sample enters the
      residual with + and the rise with -, so its noise pushes this cosine
      DOWN, never up.
  Q4  n = 4, so any per-patient ordering is anecdote; the numbers are
      printed, no test is run.

PART B (secondary)  GSE229605: resolved vs non-lesional epidermis of
psoriasis patients, paired within patient, at baseline and after T-cell
stimulation (OKT3 vs IgG control). It asks whether skin that has been
clear for a long time still carries a displacement, and whether it
over-responds to a challenge. Caveat: this is a 116-gene targeted panel
that covers the inflammatory modules but only 2/12 barrier genes, so it
CANNOT test the barrier prediction - it tests the inflammatory side only.

Run:  python3 -m netjac.demos.d19_slow_module_relapse
"""
import os
import pickle

import numpy as np

from netjac import DATA_CACHE, FIG
from netjac.data import geo
from netjac.demos import d14_treatment_restoration as d14

RES = os.path.join(DATA_CACHE, "d19_results.pkl")
MOD = {k: v for k, v in d14.SKIN.items()}          # includes the housekeeping control
SLOW = "barrier / lipid (down in lesions)"


# ================================================================ part A
def load_a():
    X, M = geo.series_matrix("GSE30768")
    G = geo.collapse_to_genes(geo.quantile_normalise(X.astype(float)),
                              geo.gpl_symbols("GPL571"))
    pats = sorted(set(M.patient))
    ref = G[M.index[M.group == "NL"]].values.mean(1)       # only 2 patients have NL
    S = {}
    for p in pats:
        m = M[M.patient == p]
        S[p] = {g: G[m.index[m.group == g][0]].values for g in m.group}
    return {"genes": np.array(G.index), "ref": ref, "S": S, "pats": pats,
            "n_ref": int((M.group == "NL").sum())}


def _disease_genes(D, thr=1.0, min_consistent=3):
    """LS - reference, averaged over patients; keep genes above the fold
    threshold whose sign is consistent in at least `min_consistent` of 4."""
    d = np.array([D["S"][p]["LS"] - D["ref"] for p in D["pats"]])
    mean = d.mean(0)
    consistent = (np.sign(d) == np.sign(mean)).sum(0)
    keep = (np.abs(mean) > thr) & (consistent >= min_consistent)
    return np.where(keep)[0], mean


def _axis_loo(D, dis, p):
    """Cohort disease axis from the OTHER patients."""
    d = np.array([D["S"][q]["LS"] - D["ref"] for q in D["pats"] if q != p])
    return d.mean(0)


def _phi(x, ref, ls, u, idx):
    uu = u[idx] / np.linalg.norm(u[idx])
    return float((x - ref)[idx] @ uu) / float((ls - ref)[idx] @ uu)


def part_a(D, n_rand=2000, seed=0):
    dis, cohort = _disease_genes(D)
    gi = {g: i for i, g in enumerate(D["genes"])}
    midx = {}
    for name, genes in MOD.items():
        ii = [gi[g] for g in genes if g in gi and abs(cohort[gi[g]]) > 0.585]
        midx[name] = np.array(ii) if len(ii) >= 3 else None
    rng = np.random.default_rng(seed)
    rows = []
    for p in D["pats"]:
        s = D["S"][p]
        u = _axis_loo(D, dis, p)
        r = {"patient": p,
             "phi_wk12": _phi(s["Week 12"], D["ref"], s["LS"], u, dis),
             "phi_flare": _phi(s["FLARE"], D["ref"], s["LS"], u, dis)}
        for name, ii in midx.items():
            r[f"m:{name}"] = (_phi(s["Week 12"], D["ref"], s["LS"], u, ii)
                              if ii is not None else np.nan)
        # Q3: residual at week 12 vs the rise into the flare
        resid = (s["Week 12"] - D["ref"])[dis]
        rise = (s["FLARE"] - s["Week 12"])[dis]
        c = float(resid @ rise / (np.linalg.norm(resid) * np.linalg.norm(rise)))
        r["cos_resid_rise"] = c
        # null: random gene sets of the same size
        null = []
        for _ in range(n_rand):
            j = rng.integers(0, len(D["genes"]), len(dis))
            a = (s["Week 12"] - D["ref"])[j]
            b = (s["FLARE"] - s["Week 12"])[j]
            null.append(float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b))))
        r["cos_null_mean"] = float(np.mean(null))
        r["cos_null_p95"] = float(np.percentile(null, 95))
        r["p_vs_null"] = float(np.mean(np.array(null) >= c))
        # flare amplitude, for the anecdotal ordering
        r["flare_rise"] = float(np.linalg.norm(rise))
        rows.append(r)
    out = {"rows": rows, "n_dis": len(dis), "modules": [k for k, v in midx.items() if v is not None],
           "mod_n": {k: (0 if v is None else len(v)) for k, v in midx.items()}, "n_ref": D["n_ref"]}
    # module restoration summary across patients
    out["restoration"] = {}
    for name in MOD:
        v = np.array([r.get(f"m:{name}", np.nan) for r in rows], float)
        if np.isfinite(v).sum() >= 3:
            out["restoration"][name] = (float(1 - np.nanmedian(v)), float(np.nanstd(1 - v)))
    # Cleanest control for Q3. Both cos(resid_i, rise_i) and the random-gene
    # null carry a built-in negative bias, because the week-12 sample enters
    # the residual with + and the rise with -. Comparing a patient's residual
    # against their OWN rise vs OTHER patients' rises keeps that bias on the
    # own-pairing only, so own > others is conservative evidence that the
    # residual is specific to that patient's coming flare.
    resid = {p: (D["S"][p]["Week 12"] - D["ref"])[dis] for p in D["pats"]}
    rise = {p: (D["S"][p]["FLARE"] - D["S"][p]["Week 12"])[dis] for p in D["pats"]}

    def cosv(a, b):
        return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)))
    own, other = [], []
    for p in D["pats"]:
        own.append(cosv(resid[p], rise[p]))
        other += [cosv(resid[p], rise[q]) for q in D["pats"] if q != p]
    out["specificity"] = {"own": own, "own_mean": float(np.mean(own)),
                          "other_mean": float(np.mean(other)),
                          "delta": float(np.mean(own) - np.mean(other)),
                          "n_own_above_other_mean": int(sum(o > np.mean(other) for o in own))}
    return out


# ================================================================ part B
def part_b():
    X, M = geo.series_matrix("GSE229605")
    L = np.log2(X.astype(float).clip(lower=1) + 1)
    L = geo.quantile_normalise(L)
    pid = M["title"].str.extract(r"(2020-\d+-\d+)")[0]
    M = M.assign(pid=pid)
    genes = np.array(L.index.astype(str))
    gi = {g: i for i, g in enumerate(genes)}
    midx = {n: np.array([gi[g] for g in gs if g in gi]) for n, gs in MOD.items()}
    midx = {n: v for n, v in midx.items() if len(v) >= 3}
    out = {"modules": {n: len(v) for n, v in midx.items()}, "pairs": {}, "challenge": {}}
    from scipy.stats import wilcoxon
    for tp in ("Baseline", "IgG treated", "OKT3 treated"):
        sel = M[M["treatment/time point"] == tp]
        res, nl = [], []
        for p, m in sel.groupby("pid"):
            a = m.index[m.tissue == "Psoriasis resolved"]
            b = m.index[m.tissue == "Psoriasis non lesional"]
            if len(a) and len(b):
                res.append(L[a[0]].values)
                nl.append(L[b[0]].values)
        if len(res) < 5:
            continue
        R, N = np.array(res), np.array(nl)
        d = R - N
        rec = {}
        for n, ii in midx.items():
            v = d[:, ii].mean(1)
            rec[n] = (float(np.mean(v)), float(wilcoxon(v).pvalue) if len(v) > 5 else np.nan)
        out["pairs"][tp] = {"n": len(res), "modules": rec}
    # does resolved skin over-respond to T-cell stimulation?
    for tis in ("Psoriasis resolved", "Psoriasis non lesional"):
        okt, igg = [], []
        for p, m in M[M.tissue == tis].groupby("pid"):
            a = m.index[m["treatment/time point"] == "OKT3 treated"]
            b = m.index[m["treatment/time point"] == "IgG treated"]
            if len(a) and len(b):
                okt.append(L[a[0]].values)
                igg.append(L[b[0]].values)
        if len(okt) < 5:
            continue
        d = np.array(okt) - np.array(igg)
        out["challenge"][tis] = {"n": len(okt),
                                 "modules": {n: float(np.mean(d[:, ii].mean(1))) for n, ii in midx.items()}}
    return out


# ================================================================ report
def run():
    R = {"A": part_a(load_a())}
    try:
        R["B"] = part_b()
    except Exception as e:
        R["B"] = {"error": str(e)[:120]}
    pickle.dump(R, open(RES, "wb"))
    return R


def report(R):
    A = R["A"]
    print(f"\nD19 PART A  psoriasis efalizumab, 4 patients with lesion / week 12 / FLARE "
          f"(GSE30768); {A['n_dis']} disease genes; reference = {A['n_ref']} non-lesional samples")
    print("\n  Q1  module restoration at week 12 (1 = back to non-lesional, 0 = still the lesion):")
    for n, (med, sd) in sorted(A["restoration"].items(), key=lambda kv: -kv[1][0]):
        tag = "  <- the slow module" if n == SLOW else ("  [negative control]" if n.startswith("housekeeping") else "")
        print(f"      {n:34s} {med:+.2f} (sd {sd:.2f}, {A['mod_n'][n]} genes){tag}")
    print("\n  Q2  where the flare lands on the disease axis (1 = the original lesion):")
    for r in A["rows"]:
        print(f"      {r['patient']}  week 12 phi {r['phi_wk12']:+.2f}   FLARE phi {r['phi_flare']:+.2f}")
    ph = np.array([r["phi_flare"] for r in A["rows"]])
    print(f"      median FLARE phi {np.median(ph):+.2f}  (range {ph.min():+.2f} to {ph.max():+.2f})")
    print("\n  Q3  is the week-12 residual already pointing along the coming flare?")
    print("      NOTE this question is CONFOUNDED by design: with one post-treatment sample")
    print("      per patient, the residual and the rise share that sample (with opposite")
    print("      signs), so neither the raw cosine nor the own-vs-other contrast is clean.")
    for r in A["rows"]:
        print(f"      {r['patient']}  cos(residual, rise) {r['cos_resid_rise']:+.3f}   "
              f"random-gene null {r['cos_null_mean']:+.3f} (95th {r['cos_null_p95']:+.3f})  p={r['p_vs_null']:.3f}")
    cs = np.array([r["cos_resid_rise"] for r in A["rows"]])
    print(f"      all four positive: {bool((cs > 0).all())}; median {np.median(cs):+.3f}")
    sp = A["specificity"]
    print(f"      specificity control - own rise {sp['own_mean']:+.3f} vs OTHER patients' rises "
          f"{sp['other_mean']:+.3f}; difference {sp['delta']:+.3f}; "
          f"{sp['n_own_above_other_mean']}/4 patients above the cross-patient mean")
    print("      -> the own-pairing carries the shared-sample bias and the cross-pairing does")
    print("         not, and the gap is about the size of that bias, so this is NOT evidence")
    print("         either way. VERDICT: not answerable here; needs >= 2 post-clearance")
    print("         samples per patient, or a cohort where some patients do NOT relapse.")
    print("\n  Q4  per-patient ordering (n = 4, anecdote only, no test):")
    for r in sorted(A["rows"], key=lambda r: -r["phi_wk12"]):
        print(f"      {r['patient']}  residual at wk12 (phi) {r['phi_wk12']:+.2f}   "
              f"size of the flare rise {r['flare_rise']:6.1f}")
    B = R.get("B", {})
    if "error" in B:
        print("\nD19 PART B skipped:", B["error"])
        return
    print(f"\nD19 PART B  resolved vs non-lesional epidermis, paired within patient (GSE229605, "
          f"116-gene panel; barrier coverage too thin to test the slow module)")
    for tp, v in B["pairs"].items():
        print(f"   {tp} (n={v['n']} patients): resolved minus non-lesional, log2")
        for n, (m, p) in sorted(v["modules"].items(), key=lambda kv: -kv[1][0]):
            star = " *" if p == p and p < 0.05 else ""
            print(f"      {n:34s} {m:+.2f}  p={p:.3f}{star}")
    if B["challenge"]:
        print("   response to T-cell stimulation (OKT3 minus IgG, log2):")
        for tis, v in B["challenge"].items():
            print(f"      {tis} (n={v['n']}): " +
                  ", ".join(f"{n.split(' (')[0][:22]} {m:+.2f}" for n, m in
                            sorted(v["modules"].items(), key=lambda kv: -kv[1])[:4]))


def figures(R):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    A = R["A"]
    fig, ax = plt.subplots(1, 3, figsize=(16.5, 5.2))
    a = ax[0]
    items = sorted(A["restoration"].items(), key=lambda kv: kv[1][0])
    cols = ["C3" if n == SLOW else ("0.6" if n.startswith("housekeeping") else "C2") for n, _ in items]
    a.barh(range(len(items)), [v[0] for _, v in items], xerr=[v[1] for _, v in items],
           color=cols, capsize=3)
    a.axvline(1, color="k", ls=":", lw=1)
    a.axvline(0, color="k", lw=0.8)
    a.set_yticks(range(len(items)), [n.split(" (")[0][:26] for n, _ in items], fontsize=8)
    a.set_xlabel("restoration at week 12 (1 = non-lesional)")
    a.set_title("Q1  clinically clear at week 12 —\nbut which modules actually returned?", fontsize=10)
    a = ax[1]
    for r in A["rows"]:
        a.plot([0, 1, 2], [1, r["phi_wk12"], r["phi_flare"]], "o-", label=r["patient"])
    a.axhline(1, color="k", ls=":", lw=1)
    a.axhline(0, color="k", lw=0.8)
    a.set_xticks([0, 1, 2], ["lesion", "week 12\n(clear)", "FLARE"])
    a.set_ylabel("position on the disease axis (1 = own lesion)")
    a.set_title("Q2  the flare returns to the lesion state", fontsize=10)
    a.legend(fontsize=8)
    a = ax[2]
    xs = np.arange(len(A["rows"]))
    a.bar(xs - 0.2, [r["cos_resid_rise"] for r in A["rows"]], 0.4, color="C3",
          label="residual at wk 12 vs rise to flare")
    a.bar(xs + 0.2, [r["cos_null_p95"] for r in A["rows"]], 0.4, color="0.7",
          label="random genes, 95th pct")
    a.axhline(0, color="k", lw=0.8)
    a.set_xticks(xs, [r["patient"] for r in A["rows"]])
    a.set_ylabel("cosine")
    a.set_title("Q3  CONFOUNDED: residual and rise share the\nweek-12 sample — not answerable here", fontsize=10)
    a.legend(fontsize=7)
    fig.tight_layout()
    p = os.path.join(FIG, "d19_slow_module_relapse.png")
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
