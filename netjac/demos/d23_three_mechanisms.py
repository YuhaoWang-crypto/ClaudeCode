"""
D23  Three mechanisms, one disease: is anything about the molecular response
specific to the drug?

D20 found that the clinical response to secukinumab (anti-IL-17A) and to
etanercept (anti-TNF) are nearly the same vector (cos +0.67), which killed
the drug-level reading of a lever ranking. That was two drugs. With a third,
mechanistically different drug - tofacitinib (JAK inhibitor, with a placebo
arm) - the question can be asked properly:

  Q1  do the three cohorts even share a disease axis? (If not, nothing else
      is comparable.)
  Q2  does the restoration ordering - inflammation fast, barrier slow -
      repeat under all three mechanisms?
  Q3  is there ANY drug-specific component? Compare the cosine between two
      patients on the SAME drug with the cosine between two patients on
      DIFFERENT drugs. If between-drug is lower than within-drug, a
      drug-specific component exists, however small. This is the test D20
      could not run.
  Q4  does the module that each drug acts on most directly move most under
      that drug? JAK inhibition blocks interferon signalling directly, so
      the interferon module should fall furthest under tofacitinib
      relative to the other modules. A pre-specified, mechanism-driven
      prediction that can fail.
  Q5  placebo control (tofacitinib cohort): restoration must stay near zero.

Cohorts
  GSE11903   etanercept (anti-TNF), 15 patients, weeks 0/1/2/4/12, GPL571
  GSE102725  secukinumab (anti-IL-17A), 10 patients, weeks 0/16/32, GPL17586
  GSE69967   tofacitinib (JAK) vs placebo, per-patient day 0 to week 12, GPL570

Run:  python3 -m netjac.demos.d23_three_mechanisms
"""
import os
import pickle

import numpy as np

from netjac import DATA_CACHE, FIG
from netjac.data import geo
from netjac.demos import d14_treatment_restoration as d14
from netjac.demos import d20_secukinumab_closed_loop as d20

RES = os.path.join(DATA_CACHE, "d23_results.pkl")
MOD = dict(d14.SKIN)
SLOW = "barrier / lipid (down in lesions)"


# ------------------------------------------------------------------ loading
def load_tof():
    """GSE69967: tofacitinib vs placebo, lesional time course + own non-lesional."""
    X, M = geo.series_matrix("GSE69967")
    G = geo.collapse_to_genes(geo.quantile_normalise(X.astype(float)),
                              geo.gpl_symbols("GPL570"))
    S, arm, resp = {}, {}, {}
    for p in sorted(set(M.patient)):
        m = M[M.patient == p]
        nl = m.index[m.tissue == "non-lesional"]
        d0 = m.index[(m.tissue == "lesional") & (m.time == "Day 0")]
        w12 = m.index[(m.tissue == "lesional") & (m.time == "Week 12")]
        if not (len(nl) and len(d0) and len(w12)):
            continue
        S[p] = {"NL": G[nl[0]].values, 0: G[d0[0]].values, 12: G[w12[0]].values}
        arm[p] = m.treatment.iloc[0]
        resp[p] = m.pasi75_w12.iloc[0]
    return {"genes": np.array(G.index), "S": S,
            "pats": [p for p in S if arm[p] == "Tofacitinib"],
            "placebo": [p for p in S if arm[p] == "Placebo"],
            "resp": resp, "arm": arm}


def load_etn():
    D = d20.load_etn()
    D["S"] = {p: {"NL": v["NL"], 0: v[0], 12: v[12]} for p, v in D["S"].items()}
    return D


def load_sec():
    D = d20.load_sec()
    # week 16 is secukinumab's clinical-response visit (everyone at PASI 98)
    D["S"] = {p: {"NL": v["NL"], 0: v[0], 12: v[16]} for p, v in D["S"].items()}
    D["pats"] = D["pats"]
    return D


def _disease(D, pats):
    return np.array([D["S"][p][0] - D["S"][p]["NL"] for p in pats]).mean(0)


def _effect(D, p):
    return D["S"][p][12] - D["S"][p][0]


def _cos(u, v):
    u = u - u.mean()
    v = v - v.mean()
    return float(u @ v / max(np.linalg.norm(u) * np.linalg.norm(v), 1e-12))


# ------------------------------------------------------------------ analysis
def run(seed=0):
    from scipy import stats
    C = {"etanercept (anti-TNF)": (load_etn(), "resp"),
         "secukinumab (anti-IL-17A)": (load_sec(), "pats"),
         "tofacitinib (JAK)": (load_tof(), "pats")}
    out = {"n": {k: len(v[0][v[1]]) for k, v in C.items()}}

    # common gene space
    shared = None
    for D, _ in C.values():
        shared = set(D["genes"]) if shared is None else shared & set(D["genes"])
    shared = sorted(shared)
    out["n_genes"] = len(shared)
    idx = {k: np.array([{g: i for i, g in enumerate(D["genes"])}[g] for g in shared])
           for k, (D, _) in C.items()}

    # Q1 one disease axis?
    dis = {k: _disease(D, D[key])[idx[k]] for k, (D, key) in C.items()}
    names = list(C)
    out["q1"] = {f"{a} vs {b}": _cos(dis[a], dis[b])
                 for i, a in enumerate(names) for b in names[i + 1:]}
    # per-patient spread within each cohort, as the yardstick
    within = {}
    for k, (D, key) in C.items():
        v = [(D["S"][p][0] - D["S"][p]["NL"])[idx[k]] for p in D[key]]
        within[k] = float(np.median([_cos(v[i], v[j]) for i in range(len(v))
                                     for j in range(i + 1, len(v))]))
    out["q1_within_cohort"] = within

    # disease genes shared by all three cohorts
    keep = np.where(np.all([np.abs(dis[k]) > 1.0 for k in names], 0))[0]
    if len(keep) < 50:
        keep = np.where(np.all([np.abs(dis[k]) > 0.585 for k in names], 0))[0]
    out["n_disease_genes"] = len(keep)

    # Q2 restoration ordering under each mechanism
    rest = {}
    for k, (D, key) in C.items():
        gi = {g: i for i, g in enumerate(D["genes"])}
        d = _disease(D, D[key])
        rows = {}
        for mn, genes in MOD.items():
            mi = np.array([gi[g] for g in genes if g in gi and abs(d[gi[g]]) > 0.585],
                          dtype=int)
            if len(mi) < 3:
                rows[mn] = {"defined": False, "n": len(mi)}
                continue
            u = d[mi] / np.linalg.norm(d[mi])
            base = [float((D["S"][p][0] - D["S"][p]["NL"])[mi] @ u) for p in D[key]]
            now = [float((D["S"][p][12] - D["S"][p]["NL"])[mi] @ u) for p in D[key]]
            b = float(np.median(base))
            vals = [1 - n / b for n in now]
            rows[mn] = {"defined": True, "n": len(mi), "med": float(np.median(vals)),
                        "sd": float(np.std(vals))}
        rest[k] = rows
    out["q2"] = rest

    # Q3 is any of the response drug-specific?
    eff = {k: {p: _effect(D, p)[idx[k]][keep] for p in D[key]}
           for k, (D, key) in C.items()}
    win, bet = [], []
    for a in names:
        va = list(eff[a].values())
        win += [_cos(va[i], va[j]) for i in range(len(va)) for j in range(i + 1, len(va))]
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            bet += [_cos(x, y) for x in eff[a].values() for y in eff[b].values()]
    u_p = stats.mannwhitneyu(win, bet, alternative="greater").pvalue
    # DRUG AND STUDY ARE PERFECTLY CONFOUNDED: each drug comes from exactly one
    # cohort on one platform, so "between drugs" is also "between studies and
    # between platforms". The one pair that shares a platform family
    # (etanercept GPL571 and tofacitinib GPL570 are both Affymetrix 3' arrays,
    # both from the same group) is the only semi-comparable contrast.
    PLATFORM = {"etanercept (anti-TNF)": "affy3", "tofacitinib (JAK)": "affy3",
                "secukinumab (anti-IL-17A)": "hta"}
    same_pf = [_cos(x, y) for x in eff["etanercept (anti-TNF)"].values()
               for y in eff["tofacitinib (JAK)"].values()]
    diff_pf = [_cos(x, y) for a in names for b in names
               if a < b and PLATFORM[a] != PLATFORM[b]
               for x in eff[a].values() for y in eff[b].values()]
    win_pf = []
    for a in ("etanercept (anti-TNF)", "tofacitinib (JAK)"):
        va = list(eff[a].values())
        win_pf += [_cos(va[i], va[j]) for i in range(len(va))
                   for j in range(i + 1, len(va))]
    out["q3"] = {"within_median": float(np.median(win)), "n_within": len(win),
                 "between_median": float(np.median(bet)), "n_between": len(bet),
                 "p": float(u_p),
                 "pairs": {f"{a} vs {b}":
                           float(np.median([_cos(x, y) for x in eff[a].values()
                                            for y in eff[b].values()]))
                           for i, a in enumerate(names) for b in names[i + 1:]},
                 "same_platform_between": float(np.median(same_pf)),
                 "diff_platform_between": float(np.median(diff_pf)),
                 "within_same_platform_drugs": float(np.median(win_pf)),
                 "p_same_platform": float(stats.mannwhitneyu(
                     win_pf, same_pf, alternative="greater").pvalue),
                 "within_by_drug": {a: float(np.median(
                     [_cos(list(eff[a].values())[i], list(eff[a].values())[j])
                      for i in range(len(eff[a])) for j in range(i + 1, len(eff[a]))]))
                     for a in names}}

    # Q4 mechanism-specific module: interferon should fall furthest under JAK
    rel = {}
    for k in names:
        r = {mn: v["med"] for mn, v in rest[k].items() if v["defined"]}
        if "interferon" in r and len(r) > 1:
            others = [v for mn, v in r.items()
                      if mn != "interferon" and not mn.startswith("housekeeping")]
            rel[k] = {"interferon": r["interferon"],
                      "other_median": float(np.median(others)),
                      "excess": float(r["interferon"] - np.median(others))}
    out["q4"] = rel

    # Q5 placebo control
    D = C["tofacitinib (JAK)"][0]
    gi = {g: i for i, g in enumerate(D["genes"])}
    d = _disease(D, D["pats"])
    pl = {}
    for mn, genes in MOD.items():
        mi = np.array([gi[g] for g in genes if g in gi and abs(d[gi[g]]) > 0.585],
                      dtype=int)
        if len(mi) < 3:
            continue
        u = d[mi] / np.linalg.norm(d[mi])
        base = [float((D["S"][p][0] - D["S"][p]["NL"])[mi] @ u) for p in D["placebo"]]
        now = [float((D["S"][p][12] - D["S"][p]["NL"])[mi] @ u) for p in D["placebo"]]
        b = float(np.median(base))
        pl[mn] = float(np.median([1 - n / b for n in now]))
    out["q5_placebo"] = {"n": len(D["placebo"]), "restoration": pl}
    out["_names"] = names
    return out


# ------------------------------------------------------------------ report
def report(R):
    names = R["_names"]
    print("=" * 78)
    print("D23  three mechanisms, one disease")
    print("=" * 78)
    print(f"  cohorts: " + ", ".join(f"{k} n={R['n'][k]}" for k in names))
    print(f"  {R['n_genes']} genes shared by all three platforms, "
          f"{R['n_disease_genes']} of them disease genes in all three cohorts")

    print("\n  Q1  is it the same disease axis?")
    for k, v in R["q1"].items():
        print(f"    {k:52s} cos {v:+.3f}")
    print("    for scale, median cosine between two PATIENTS inside one cohort:")
    for k, v in R["q1_within_cohort"].items():
        print(f"    {k:52s} cos {v:+.3f}")

    print("\n  Q2  restoration at the clinical-response visit (1 = back to non-lesional)")
    mods = [m for m in MOD if any(R["q2"][k][m]["defined"] for k in names)]
    print(f"    {'module':34s}" + "".join(f"{k.split(' (')[0][:14]:>16s}" for k in names))
    for m in mods:
        row = f"    {m.split(' (')[0][:34]:34s}"
        for k in names:
            v = R["q2"][k][m]
            cell = ("%+.2f (n%d)" % (v["med"], v["n"])) if v["defined"] else "undef"
            row += "%16s" % cell
        print(row)

    q = R["q3"]
    print("\n  Q3  is any of the response drug-specific?")
    print(f"    two patients on the SAME drug:      median cos {q['within_median']:+.3f} "
          f"({q['n_within']} pairs)")
    print(f"    two patients on DIFFERENT drugs:    median cos {q['between_median']:+.3f} "
          f"({q['n_between']} pairs)")
    print(f"    within > between, Mann-Whitney p = {q['p']:.2e}")
    for k, v in q["within_by_drug"].items():
        print(f"      within {k:30s} {v:+.3f}")
    for k, v in q["pairs"].items():
        print(f"      between {k:48s} {v:+.3f}")
    print("    CONFOUND: each drug comes from one cohort on one platform, so "
          "'different drug' also means\n      'different study and platform'. "
          "Etanercept and tofacitinib share a platform family and a group:")
    print(f"      within those two drugs            {q['within_same_platform_drugs']:+.3f}")
    print(f"      between them (same platform)      {q['same_platform_between']:+.3f}"
          f"   Mann-Whitney p = {q['p_same_platform']:.3f}")
    print(f"      between drugs across platforms    {q['diff_platform_between']:+.3f}")

    print("\n  Q4  does the interferon module fall furthest under the JAK inhibitor?")
    for k, v in R["q4"].items():
        print(f"    {k:28s} interferon {v['interferon']:+.2f} vs other modules "
              f"{v['other_median']:+.2f}  (excess {v['excess']:+.2f})")

    p = R["q5_placebo"]
    print(f"\n  Q5  placebo arm of the tofacitinib trial (n = {p['n']}), restoration at week 12")
    for mn, v in sorted(p["restoration"].items(), key=lambda kv: -kv[1]):
        print(f"    {mn:42s} {v:+.2f}")


def figures(R):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    names = R["_names"]
    fig, ax = plt.subplots(1, 3, figsize=(17, 5.4))

    a = ax[0]
    mods = [m for m in MOD if all(R["q2"][k][m]["defined"] for k in names)]
    xs = np.arange(len(mods))
    w = 0.26
    for i, k in enumerate(names):
        a.bar(xs + (i - 1) * w, [R["q2"][k][m]["med"] for m in mods], w,
              label=k.split(" (")[0])
    a.axhline(0, color="k", lw=0.8)
    a.axhline(1, color="k", ls=":", lw=1)
    a.set_xticks(xs, [m.split(" (")[0][:20] for m in mods], rotation=30, ha="right",
                 fontsize=8)
    a.set_ylabel("restoration (1 = non-lesional)")
    a.set_title("Q2  same ordering under three mechanisms", fontsize=10)
    a.legend(fontsize=8)

    a = ax[1]
    q = R["q3"]
    vals = [q["within_median"], q["same_platform_between"], q["diff_platform_between"]]
    a.bar(range(3), vals, 0.55, color=["C2", "C0", "C3"])
    a.set_xticks(range(3), ["same drug\n(same study)",
                            "different drug\nsame platform",
                            "different drug\nother platform"], fontsize=8.5)
    a.set_ylabel("median cos between two patients' drug effects")
    a.set_title("Q3  drug and study are confounded:\nmost of the 'specificity' is "
                "the platform", fontsize=10)
    for i, v in enumerate(vals):
        a.text(i, v + 0.01, f"{v:+.2f}", ha="center", fontsize=9)

    a = ax[2]
    ks = list(R["q4"])
    a.bar(np.arange(len(ks)) - 0.2, [R["q4"][k]["interferon"] for k in ks], 0.38,
          label="interferon module", color="C0")
    a.bar(np.arange(len(ks)) + 0.2, [R["q4"][k]["other_median"] for k in ks], 0.38,
          label="median of the other modules", color="0.7")
    a.set_xticks(range(len(ks)), [k.split(" (")[0][:14] for k in ks], fontsize=9)
    a.axhline(0, color="k", lw=0.8)
    a.set_ylabel("restoration")
    a.set_title("Q4  mechanism-specific prediction:\nJAK should clear interferon first",
                fontsize=10)
    a.legend(fontsize=8)

    fig.tight_layout()
    fn = os.path.join(FIG, "d23_three_mechanisms.png")
    fig.savefig(fn, dpi=140)
    print(f"\nfigure: {fn}")


def main(force=False):
    if os.path.exists(RES) and not force:
        R = pickle.load(open(RES, "rb"))
    else:
        R = run()
        pickle.dump(R, open(RES, "wb"))
    report(R)
    figures(R)
    return R


if __name__ == "__main__":
    main(force=bool(os.environ.get("D23_FORCE")))
