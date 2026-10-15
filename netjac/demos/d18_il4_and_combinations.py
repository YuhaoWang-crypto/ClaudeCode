"""
D18  Closing two gaps left by D17: IL-4, and measured COMBINATION responses.

D17 built a keratinocyte lever library from GSE305982 and got the psoriasis
answer right (IL-17A) but left two holes: the panel had no IL-4 (dupilumab
blocks IL-4Ralpha, which carries IL-4 and IL-13), and it had no combinations,
so nothing tested what D16 proved a linear J cannot do - predict the effect
of pulling two levers at once.

PART A  IL-4 and TNF in the library (GSE36287, Chiricozzi / Krueger).
Primary human keratinocytes from 3 donors, untreated or treated with
IFN-alpha, IFN-gamma, IL-4, IL-13, IL-17A or TNF. Now both indications have
their full validated lever set, so the ground truth is a SET, not a single
label:
  psoriasis         IL-17A (secukinumab, ixekizumab) and TNF (etanercept,
                    adalimumab) are validated; IL-4 / IL-13 are not.
  atopic dermatitis IL-4 and IL-13 (dupilumab via IL-4Ralpha, tralokinumab,
                    lebrikizumab) are validated; IL-17A is not.
GSE36287 is on GPL570, the same platform as the atopic dermatitis cohort.

PART B  Measured combinations (GSE24767, same group). Primary keratinocytes
from the same donors treated with IL-17, TNF (1 and 10 ng/ml), and
IL-17 + TNF. With singles AND the combination measured in the same donors,
three things can be separated:
  additive prediction   d_IL17 + d_TNF        (all a linear J can produce)
  measured combination  d_combo
  synergy               s = d_combo - d_IL17 - d_TNF
and the question D16 could not ask with ground truth becomes direct: is the
synergy - the part no linear model can generate - itself aligned with the
psoriasis lesion?

Run:  python3 -m netjac.demos.d18_il4_and_combinations
"""
import os
import pickle

import numpy as np

from netjac import DATA_CACHE, FIG
from netjac.data import geo
from netjac.demos import d14_treatment_restoration as d14
from netjac.demos import d17_keratinocyte_lever as d17

RES = os.path.join(DATA_CACHE, "d18_results.pkl")

# clinically validated lever sets, independent of every dataset used here
VALID = {"psoriasis": {"IL-17A", "TNF"},
         "atopic dermatitis": {"IL-4", "IL-13"}}
LEVERS_A = {"IFN-alpha": "IFNa-treated", "IFN-gamma": "IFNg-treated", "IL-4": "IL4-treated",
            "IL-13": "IL13-treated", "IL-17A": "IL17A-treated", "TNF": "TNF-treated"}


def _genes(gse, gpl, log=False):
    X, M = geo.series_matrix(gse)
    X = X.astype(float)
    if log:
        X = np.log2(X.clip(lower=16))
    G = geo.collapse_to_genes(geo.quantile_normalise(X), geo.gpl_symbols(gpl))
    return G, M


# ============================================================ part A
def library_il4():
    """Lever library with IL-4 and TNF: mean treated - mean untreated."""
    G, M = _genes("GSE36287", "GPL570")
    ctrl = G[M.index[M.treatment == "untreated control"]].values.mean(1)
    cols, names = [], []
    for name, tag in LEVERS_A.items():
        idx = M.index[M.treatment == tag]
        if len(idx):
            cols.append(G[idx].values.mean(1) - ctrl)
            names.append(name)
    return {"genes": np.array(G.index), "R": np.column_stack(cols), "names": names,
            "n_ctrl": int((M.treatment == "untreated control").sum())}


def part_a(DIS, n_boot=2000, seed=0):
    lib = library_il4()
    out = {"names": lib["names"], "res": {}, "crossed": {}}
    for ind, dv in DIS.items():
        out["res"][ind] = d17.rank(dv["d"], dv["genes"], lib["R"], lib["genes"], lib["names"])
    # crossed contrast per lever: aligned with its own indication more than the other?
    gi = {g: i for i, g in enumerate(lib["genes"])}
    shared = sorted(set(gi) & set(DIS["psoriasis"]["genes"]) & set(DIS["atopic dermatitis"]["genes"]))
    A = lib["R"][[gi[g] for g in shared]]
    A = A - A.mean(0)
    Y = {}
    for ind, dv in DIS.items():
        di = {g: i for i, g in enumerate(dv["genes"])}
        y = np.asarray(dv["d"], float)[[di[g] for g in shared]]
        Y[ind] = y - y.mean()
    rng = np.random.default_rng(seed)

    def cosv(a, y):
        return float(a @ y / max(np.linalg.norm(a) * np.linalg.norm(y), 1e-12))
    for lever in lib["names"]:
        j = lib["names"].index(lever)
        pro = "psoriasis" if lever in VALID["psoriasis"] else (
            "atopic dermatitis" if lever in VALID["atopic dermatitis"] else None)
        if pro is None:
            continue
        other = "atopic dermatitis" if pro == "psoriasis" else "psoriasis"
        obs = cosv(A[:, j], Y[pro]) - cosv(A[:, j], Y[other])
        bs = [cosv(A[b, j], Y[pro][b]) - cosv(A[b, j], Y[other][b])
              for b in (rng.integers(0, len(shared), len(shared)) for _ in range(n_boot))]
        out["crossed"][lever] = {"for": pro, "delta": obs,
                                 "ci": tuple(np.percentile(bs, [2.5, 97.5])),
                                 "p": float(np.mean(np.array(bs) <= 0))}
    out["n_genes"] = len(shared)
    return out


# ============================================================ part B
def part_b(DIS, dose="10 ng/ml", n_boot=2000, seed=0):
    """Paired within-donor single and combination responses."""
    G, M = _genes("GSE24767", "GPL6947")
    tag = {"ctrl": "none", "IL17": "IL17", "TNF": f"{dose} TNFα",
           "combo": f"IL17 + {dose} TNFα"}
    # the treatment strings use 'TNFα'; match loosely
    def pick(did, key):
        want = {"ctrl": lambda s: s == "none",
                "IL17": lambda s: s == "IL17",
                "TNF": lambda s: s.startswith(dose.split()[0] + " ng/ml"),
                "combo": lambda s: s.startswith("IL17 +") and s.split("+")[1].strip().startswith(dose.split()[0])}[key]
        idx = [i for i in M.index if M.loc[i, "ID"] == did and want(M.loc[i, "treatment"])]
        return G[idx].values.mean(1) if idx else None
    donors = [d for d in sorted(set(M["ID"]))
              if all(pick(d, k) is not None for k in ("ctrl", "IL17", "TNF", "combo"))]
    d17v, dtnf, dcom = [], [], []
    for d in donors:
        c = pick(d, "ctrl")
        d17v.append(pick(d, "IL17") - c)
        dtnf.append(pick(d, "TNF") - c)
        dcom.append(pick(d, "combo") - c)
    # is the synergy reproducible across donors, or is 'non-additivity' noise?
    syn_d = [c - a - b for a, b, c in zip(d17v, dtnf, dcom)]
    com_d = list(dcom)
    pair_syn = [float(np.corrcoef(syn_d[i], syn_d[j])[0, 1])
                for i in range(len(syn_d)) for j in range(i + 1, len(syn_d))]
    pair_com = [float(np.corrcoef(com_d[i], com_d[j])[0, 1])
                for i in range(len(com_d)) for j in range(i + 1, len(com_d))]
    d17v, dtnf, dcom = map(lambda a: np.mean(a, 0), (d17v, dtnf, dcom))
    add = d17v + dtnf
    syn = dcom - add
    genes = np.array(G.index)
    out = {"donors": donors, "dose": dose,
           "synergy_donor_r": (float(np.mean(pair_syn)), [round(x, 2) for x in pair_syn]),
           "combo_donor_r": (float(np.mean(pair_com)), [round(x, 2) for x in pair_com]),
           "nonadditivity": float(np.linalg.norm(syn) / np.linalg.norm(dcom)),
           "cos_add_combo": float(add @ dcom / (np.linalg.norm(add) * np.linalg.norm(dcom))),
           "gain": float(np.linalg.norm(dcom) / np.linalg.norm(add))}
    # alignment with the psoriasis lesion
    dv = DIS["psoriasis"]
    di = {g: i for i, g in enumerate(dv["genes"])}
    gi = {g: i for i, g in enumerate(genes)}
    shared = sorted(set(gi) & set(di))
    y = np.asarray(dv["d"], float)[[di[g] for g in shared]]
    y = y - y.mean()
    vecs = {"IL-17 alone": d17v, "TNF alone": dtnf, "additive IL-17 + TNF": add,
            "measured IL-17 + TNF": dcom, "synergy (combo - additive)": syn}
    rng = np.random.default_rng(seed)

    def cosv(a, b):
        return float(a @ b / max(np.linalg.norm(a) * np.linalg.norm(b), 1e-12))
    out["cos"], out["ci"] = {}, {}
    V = {k: (v[[gi[g] for g in shared]] - v[[gi[g] for g in shared]].mean()) for k, v in vecs.items()}
    for k, v in V.items():
        out["cos"][k] = cosv(v, y)
        bs = [cosv(v[b], y[b]) for b in (rng.integers(0, len(shared), len(shared)) for _ in range(n_boot))]
        out["ci"][k] = tuple(np.percentile(bs, [2.5, 97.5]))
    # does the measured combination beat the additive prediction? paired bootstrap
    bs = [cosv(V["measured IL-17 + TNF"][b], y[b]) - cosv(V["additive IL-17 + TNF"][b], y[b])
          for b in (rng.integers(0, len(shared), len(shared)) for _ in range(n_boot))]
    out["delta_measured_minus_additive"] = (float(np.mean(bs)), *np.percentile(bs, [2.5, 97.5]),
                                            float(np.mean(np.array(bs) <= 0)))
    out["n_genes"] = len(shared)
    # Is the synergy real at SPECIFIC genes even though the global vector is
    # donor-noisy? Leave-one-donor-out: pick the top synergy genes in the other
    # donors, then measure the held-out donor's synergy on them vs random genes.
    loo = []
    for h in range(len(syn_d)):
        others = np.mean([syn_d[i] for i in range(len(syn_d)) if i != h], 0)
        top = np.argsort(-others)[:20]
        rngl = np.random.default_rng(h)
        rand = [float(np.mean(syn_d[h][rngl.integers(0, len(syn_d[h]), 20)])) for _ in range(500)]
        loo.append({"held_out": donors[h], "mean_syn_top": float(np.mean(syn_d[h][top])),
                    "mean_syn_random": float(np.mean(rand)),
                    "p": float(np.mean(np.array(rand) >= np.mean(syn_d[h][top])))})
    out["loo_synergy"] = loo
    # top synergy genes
    o = np.argsort(-syn)
    out["top_synergy"] = [(genes[i], float(syn[i]), float(dcom[i]), float(add[i])) for i in o[:15]]
    return out


# ================================================================= report
def run():
    DIS = d17.disease_vectors()
    R = {"A": part_a(DIS), "B": {d: part_b(DIS, dose=d) for d in ("10 ng/ml", "1 ng/ml")},
         "valid": VALID}
    pickle.dump(R, open(RES, "wb"))
    return R


def report(R):
    A = R["A"]
    print(f"\nD18 part A: keratinocyte lever library WITH IL-4 and TNF (GSE36287, 3 donors, "
          f"{A['n_genes']} shared genes)")
    print(f"   validated lever sets: " + "; ".join(f"{k} -> {sorted(v)}" for k, v in R["valid"].items()))
    for ind, v in A["res"].items():
        top = v["top"]
        ok = "HIT" if top in R["valid"][ind] else "miss"
        print(f"\n   {ind}: top lever {top} [{ok}]")
        for j in np.argsort(-v["energy"]):
            n = v["names"][j]
            mark = " <- validated" if n in R["valid"][ind] else ""
            print(f"      {n:10s} cos {v['cos'][j]:+.3f}  coef {v['coef'][j]:+.3f}  p={v['p'][j]:.4f}{mark}")
    print("\n   crossed specificity (lever vs its own indication minus the other):")
    for lever, v in A["crossed"].items():
        ok = "OK" if v["ci"][0] > 0 else "n.s."
        print(f"      {lever:10s} for {v['for']:18s} {v['delta']:+.3f} "
              f"[{v['ci'][0]:+.3f}, {v['ci'][1]:+.3f}] p={v['p']:.4f} [{ok}]")
    for dose, B in R["B"].items():
        print(f"\nD18 part B: measured IL-17 + TNF ({dose}) in the same donors "
              f"(GSE24767, donors {B['donors']}, {B['n_genes']} genes)")
        print(f"   non-additivity |combo - (IL17+TNF)| / |combo| = {B['nonadditivity']:.2f}; "
              f"|combo| / |additive| = {B['gain']:.2f}; cos(additive, combo) = {B['cos_add_combo']:.2f}")
        print(f"   reproducible across donors? synergy vector between-donor r {B['synergy_donor_r'][0]:+.2f} "
              f"{B['synergy_donor_r'][1]}; combination response r {B['combo_donor_r'][0]:+.2f} "
              f"{B['combo_donor_r'][1]} (synergy is noise if the first is ~0)")
        print("   alignment with the psoriasis lesion:")
        for k in ("IL-17 alone", "TNF alone", "additive IL-17 + TNF", "measured IL-17 + TNF",
                  "synergy (combo - additive)"):
            c, ci = B["cos"][k], B["ci"][k]
            print(f"      {k:28s} cos {c:+.3f} [{ci[0]:+.3f}, {ci[1]:+.3f}]")
        d = B["delta_measured_minus_additive"]
        print(f"   measured - additive: {d[0]:+.3f} [{d[1]:+.3f}, {d[2]:+.3f}] p={d[3]:.4f}")
        print("   leave-one-donor-out: synergy of the held-out donor at the top-20 synergy genes "
              "found in the others, vs random genes:")
        for l in B["loo_synergy"]:
            print(f"      donor {l['held_out']}: {l['mean_syn_top']:+.2f} log2 vs random "
                  f"{l['mean_syn_random']:+.2f}, p={l['p']:.3f}")
        if dose == "10 ng/ml":
            print("   strongest synergy genes (combo vs additive, log2): " +
                  ", ".join(f"{g} {s:+.1f}" for g, s, _, _ in B["top_synergy"][:8]))


def figures(R):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    fig, ax = plt.subplots(1, 4, figsize=(20.5, 5.2))
    A = R["A"]
    a = ax[0]
    names = A["names"]
    w = 0.38
    for k, ind in enumerate(("psoriasis", "atopic dermatitis")):
        v = A["res"][ind]
        col = ["C3" if n in R["valid"][ind] else "0.7" for n in v["names"]]
        a.bar(np.arange(len(names)) + (k - 0.5) * w, v["cos"], w, color=col,
              edgecolor="k" if k else None, linewidth=0.6, hatch="" if k == 0 else "//",
              label=f"{ind} (validated: {', '.join(sorted(R['valid'][ind]))})")
    a.axhline(0, color="k", lw=0.8)
    a.set_xticks(range(len(names)), names, rotation=30, ha="right", fontsize=8)
    a.set_ylabel("cosine with lesional displacement")
    a.set_title("A  library with IL-4 and TNF (GSE36287)\nred = clinically validated for that indication", fontsize=9.5)
    a.legend(fontsize=7)
    a = ax[1]
    cr = A["crossed"]
    ks = list(cr)
    vals = [cr[k]["delta"] for k in ks]
    err = [[cr[k]["delta"] - cr[k]["ci"][0] for k in ks], [cr[k]["ci"][1] - cr[k]["delta"] for k in ks]]
    a.bar(range(len(ks)), vals, yerr=err, capsize=4,
          color=["C2" if cr[k]["ci"][0] > 0 else "C1" for k in ks])
    a.axhline(0, color="k", lw=0.8)
    a.set_xticks(range(len(ks)), [f"{k}\nfor {cr[k]['for'][:12]}" for k in ks], fontsize=7)
    a.set_ylabel("cos(own indication) - cos(other)")
    a.set_title("A  crossed specificity\n(green = CI excludes 0)", fontsize=9.5)
    a = ax[2]
    B = R["B"]["10 ng/ml"]
    ks = ["IL-17 alone", "TNF alone", "additive IL-17 + TNF", "measured IL-17 + TNF",
          "synergy (combo - additive)"]
    vals = [B["cos"][k] for k in ks]
    err = [[B["cos"][k] - B["ci"][k][0] for k in ks], [B["ci"][k][1] - B["cos"][k] for k in ks]]
    cols = ["0.7", "0.7", "C0", "C2", "C3"]
    a.bar(range(len(ks)), vals, yerr=err, color=cols, capsize=4)
    a.axhline(0, color="k", lw=0.8)
    a.set_xticks(range(len(ks)), [k.replace(" (combo - additive)", "\n(combo - additive)").replace(" IL-17 + TNF", "\nIL-17 + TNF")
                                  for k in ks], fontsize=7)
    a.set_ylabel("cosine with the psoriasis lesion")
    a.set_title(f"B  measured combination vs linear superposition\n"
                f"measured - additive {B['delta_measured_minus_additive'][0]:+.3f} "
                f"[{B['delta_measured_minus_additive'][1]:+.3f}, {B['delta_measured_minus_additive'][2]:+.3f}]",
                fontsize=9.5)
    a = ax[3]
    loo = B["loo_synergy"]
    x = np.arange(len(loo))
    a.bar(x - 0.2, [l["mean_syn_top"] for l in loo], 0.4, color="C3", label="top-20 synergy genes\n(chosen in the OTHER donors)")
    a.bar(x + 0.2, [l["mean_syn_random"] for l in loo], 0.4, color="0.7", label="random genes")
    a.axhline(0, color="k", lw=0.8)
    a.set_xticks(x, [f"held-out\ndonor {l['held_out']}" for l in loo], fontsize=8)
    a.set_ylabel("synergy in the held-out donor (log2)")
    a.set_title(f"B  synergy is gene-specific, not global\n"
                f"whole-vector donor r {B['synergy_donor_r'][0]:+.2f} vs combination {B['combo_donor_r'][0]:+.2f}",
                fontsize=9.5)
    a.legend(fontsize=7)
    fig.tight_layout()
    p = os.path.join(FIG, "d18_il4_and_combinations.png")
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
