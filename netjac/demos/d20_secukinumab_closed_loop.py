"""
D20  Closing the loop on IL-17A: does the drug that blocks the predicted
lever move the tissue the way the prediction says?

D17/D18 ranked nine measured keratinocyte cytokine responses against the
psoriasis displacement and put IL-17A first. That is a prediction about a
drug that did not enter the ranking: blocking IL-17A should move lesional
skin along MINUS the measured IL-17A response direction. Secukinumab is
exactly that antibody, and GSE102725 is exactly that experiment.

PART A (main)  GSE102725 - PSORITUS (NCT02362789, PMID 39565228).
Ten patients, lesional skin at baseline / week 16 / week 32 plus their own
non-lesional skin at baseline. Everyone receives open-label secukinumab
300 mg to week 16; only patients reaching PASI 98 are then randomised to
continue (5) or to withdraw to placebo (5) until week 32.
  Q1  does the observed drug effect (week 16 - baseline) align with MINUS
      the measured IL-17A keratinocyte response, and does IL-17A win the
      ranking over the other eight levers?
  Q2  the boring explanation: treatment reverses the disease, and every
      psoriasis-related cytokine looks like the disease axis. So repeat Q1
      after projecting the disease axis out of both sides. What survives is
      lever information, not "this looks like psoriasis".
  Q3  magnitude: least squares within one source (GSE36287, six levers,
      magnitude-comparable) - is the IL-17A coefficient the one that carries
      the drug effect?
  Q4  CROSSED control with GSE11903 (etanercept, anti-TNF, 15 patients,
      weeks 0/1/2/4/12). Two drugs x the same library. If the library carries
      drug-specific information, IL-17A should beat TNF more for secukinumab
      than for etanercept. This is the test that can fail, and it is the one
      that says whether lever ranking is about the DRUG or only the DISEASE.

PART B  what full IL-17A blockade does not restore, and the withdrawal arm.
  Q5  module restoration at week 16 in PASI-98 responders: a fifth drug, the
      strictest clinical endpoint so far. Does the barrier module lag again?
  Q6  week 32, randomised: placebo-withdrawal vs continued secukinumab on the
      disease axis. This is the control arm D19 did not have.
  Q7  D19's open question, now with that control arm: inside the placebo arm,
      does the week-16 residual point along the relapse direction? n = 5, and
      the secukinumab arm gives the same statistic where no relapse happens,
      which is the bias estimate D19 was missing.

Run:  python3 -m netjac.demos.d20_secukinumab_closed_loop
"""
import os
import pickle
import re

import numpy as np

from netjac import DATA_CACHE, FIG
from netjac.data import geo
from netjac.data import keratinocyte as kc
from netjac.demos import d14_treatment_restoration as d14

RES = os.path.join(DATA_CACHE, "d20_results.pkl")
MOD = dict(d14.SKIN)
SLOW = "barrier / lipid (down in lesions)"
# Visit -> week. The paper takes lesional biopsies at baseline, week 16 and
# week 32; the three lesional visits are V2 < V9 < V14 in that order.
# Checked against the data in report(): V9 must be clear and the placebo arm
# must rebound by V14, which is what the trial reports.
WEEK = {"V2": 0, "V9": 16, "V14": 32}


# ================================================================== loading
def load_sec():
    """GSE102725, secukinumab withdrawal trial."""
    X, M = geo.series_matrix("GSE102725")
    G = geo.collapse_to_genes(geo.quantile_normalise(X.astype(float)),
                              geo.gpl_symbols("GPL17586"))
    rows = {}
    for gsm, t in zip(M.index, M.title):
        pat = re.search(r"Skin_(\S+)\s", t).group(1)
        visit = re.search(r"\s(V\d+)\b", t).group(1)
        arm = re.search(r"\((V|P)\)", t).group(1)
        rows[gsm] = {"patient": pat, "week": WEEK[visit],
                     "lesional": "non-les" not in t, "arm": arm}
    pats = sorted({r["patient"] for r in rows.values()})
    S, arm = {}, {}
    for p in pats:
        g = {k: v for k, v in rows.items() if v["patient"] == p}
        arm[p] = next(iter(g.values()))["arm"]
        S[p] = {"NL": G[[k for k, v in g.items() if not v["lesional"]][0]].values}
        for w in (0, 16, 32):
            k = [k for k, v in g.items() if v["lesional"] and v["week"] == w]
            if k:
                S[p][w] = G[k[0]].values
    return {"genes": np.array(G.index), "S": S, "pats": pats, "arm": arm,
            "sec": [p for p in pats if arm[p] == "V"],
            "pbo": [p for p in pats if arm[p] == "P"]}


def load_etn():
    """GSE11903, etanercept time course (responders used for the week-12 effect)."""
    X, M = geo.series_matrix("GSE11903")
    G = geo.collapse_to_genes(geo.quantile_normalise(X.astype(float)),
                              geo.gpl_symbols("GPL571"))
    S, pats = {}, []
    for p in sorted(set(M.Patient)):
        m = M[M.Patient == p]
        nl = m.index[(m.Condition == "non-lesional")]
        ls0 = m.index[(m.Condition == "lesional") & (m.Time == "0")]
        ls12 = m.index[(m.Condition == "lesional") & (m.Time == "12 weeks")]
        if len(nl) and len(ls0) and len(ls12):
            S[p] = {"NL": G[nl[0]].values, 0: G[ls0[0]].values, 12: G[ls12[0]].values,
                    "group": m.Group.iloc[0]}
            pats.append(p)
    return {"genes": np.array(G.index), "S": S, "pats": pats,
            "resp": [p for p in pats if S[p]["group"] == "responder"]}


# ================================================================== helpers
def _paired(D, pats, a, b):
    """Mean within-patient difference of states a and b (keys of D['S'][p])."""
    return np.array([D["S"][p][a] - D["S"][p][b] for p in pats]).mean(0)


def _disease_genes(D, pats, thr=1.0, frac=0.7):
    d = np.array([D["S"][p][0] - D["S"][p]["NL"] for p in pats])
    mean = d.mean(0)
    consistent = (np.sign(d) == np.sign(mean)).mean(0)
    return np.where((np.abs(mean) > thr) & (consistent >= frac))[0], mean


def _cos(u, v):
    u = u - u.mean()
    v = v - v.mean()
    return float(u @ v / max(np.linalg.norm(u) * np.linalg.norm(v), 1e-12))


def _boot_cos(u, v, n=2000, seed=0):
    rng = np.random.default_rng(seed)
    b = [_cos(u[i], v[i]) for i in (rng.integers(0, len(u), len(u)) for _ in range(n))]
    return tuple(np.percentile(b, [2.5, 97.5]))


def _resid(x, d):
    """Remove the component along the disease axis d."""
    dh = (d - d.mean())
    dh = dh / max(np.linalg.norm(dh), 1e-12)
    x = x - x.mean()
    return x - (x @ dh) * dh


def _align(lib, genes, idx):
    """Shared genes between a cohort's selected genes and the lever library."""
    gi = {g: i for i, g in enumerate(lib["genes"])}
    keep = [(j, gi[genes[j]]) for j in idx if genes[j] in gi]
    return np.array([k[0] for k in keep]), np.array([k[1] for k in keep])


# ================================================================== part A
def part_a(SEC, ETN, lib, n_rand=2000, seed=0):
    out = {}
    idx, dis = _disease_genes(SEC, SEC["pats"])
    out["n_disease_genes"] = len(idx)
    ci, li = _align(lib, SEC["genes"], idx)
    out["n_shared"] = len(ci)

    drug = _paired(SEC, SEC["pats"], 16, 0)[ci]      # observed secukinumab effect
    d = dis[ci]                                      # disease displacement
    R = lib["R"][li]                                 # measured lever responses
    out["cos_drug_vs_minus_disease"] = _cos(drug, -d)   # the boring baseline

    rows = []
    for j, name in enumerate(lib["names"]):
        r = -R[:, j]                                 # blocking lever j
        c = _cos(drug, r)
        lo, hi = _boot_cos(drug, r, seed=seed)
        cr = _cos(_resid(drug, d), _resid(r, d))
        lo2, hi2 = _boot_cos(_resid(drug, d), _resid(r, d), seed=seed)
        rows.append({"lever": name, "source": lib["source"][j], "cos": c, "ci": (lo, hi),
                     "cos_vs_disease": _cos(-r, d), "cos_resid": cr, "ci_resid": (lo2, hi2)})
    rows.sort(key=lambda r: -r["cos"])
    out["levers"] = rows
    out["top"] = rows[0]["lever"]
    out["top_resid"] = max(rows, key=lambda r: r["cos_resid"])["lever"]

    # null: shuffle the gene labels of the IL-17A response
    k = lib["names"].index("IL-17A")
    rng = np.random.default_rng(seed)
    r17 = -R[:, k]
    null = [_cos(drug, r17[rng.permutation(len(r17))]) for _ in range(n_rand)]
    out["il17a_vs_shuffle"] = {"cos": _cos(drug, r17), "null_p95": float(np.percentile(null, 95)),
                              "p": float(np.mean(np.array(null) >= _cos(drug, r17)))}

    # Q3 magnitude: least squares within GSE36287 only (comparable magnitudes)
    same = [j for j, s in enumerate(lib["source"]) if s == "GSE36287"]
    A = np.column_stack([-R[:, j] - (-R[:, j]).mean() for j in same])
    y = drug - drug.mean()
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    out["lstsq"] = {"levers": [lib["names"][j] for j in same], "coef": coef.tolist(),
                    "r2": float(1 - np.sum((y - A @ coef) ** 2) / np.sum(y ** 2))}
    k1 = lib["names"].index("IL-17A")
    a1 = (-R[:, k1] - (-R[:, k1]).mean())[:, None]
    c1, *_ = np.linalg.lstsq(a1, y, rcond=None)
    out["lstsq_il17a_only"] = {"coef": float(c1[0]),
                               "r2": float(1 - np.sum((y - a1 @ c1) ** 2) / np.sum(y ** 2))}

    # Q4 crossed: two drugs x the same library, on genes that are disease
    # genes in BOTH cohorts, so the contrast cannot come from gene selection.
    eidx, edis = _disease_genes(ETN, ETN["resp"])
    sg = set(SEC["genes"][idx]) & set(ETN["genes"][eidx]) & set(lib["genes"])
    gi = {g: i for i, g in enumerate(lib["genes"])}
    si = {g: i for i, g in enumerate(SEC["genes"])}
    ei = {g: i for i, g in enumerate(ETN["genes"])}
    sh = sorted(sg)
    out["n_crossed_genes"] = len(sh)
    dr = {"secukinumab (anti-IL-17A)": _paired(SEC, SEC["pats"], 16, 0)[[si[g] for g in sh]],
          "etanercept (anti-TNF)": _paired(ETN, ETN["resp"], 12, 0)[[ei[g] for g in sh]]}
    Rs = lib["R"][[gi[g] for g in sh]]
    cross = {}
    for dn, dv in dr.items():
        cross[dn] = {n: _cos(dv, -Rs[:, j]) for j, n in enumerate(lib["names"])}
    out["crossed"] = cross
    # difference-of-differences with a gene bootstrap
    j17, jtnf = lib["names"].index("IL-17A"), lib["names"].index("TNF")
    rng = np.random.default_rng(seed)

    def dod(sel):
        a = _cos(dr["secukinumab (anti-IL-17A)"][sel], -Rs[sel, j17]) \
            - _cos(dr["secukinumab (anti-IL-17A)"][sel], -Rs[sel, jtnf])
        b = _cos(dr["etanercept (anti-TNF)"][sel], -Rs[sel, j17]) \
            - _cos(dr["etanercept (anti-TNF)"][sel], -Rs[sel, jtnf])
        return a - b
    allg = np.arange(len(sh))
    obs = dod(allg)
    bs = [dod(rng.integers(0, len(sh), len(sh))) for _ in range(n_rand)]
    out["crossed_dod"] = {"obs": obs, "ci": tuple(np.percentile(bs, [2.5, 97.5])),
                          "p": float(np.mean(np.array(bs) <= 0))}
    # How much of the ranking is just "this lever looks like psoriasis"?
    from scipy.stats import spearmanr
    rho, pv = spearmanr([r["cos"] for r in rows], [r["cos_vs_disease"] for r in rows])
    out["rank_vs_disease_rho"] = {"rho": float(rho), "p": float(pv)}
    a = [cross["secukinumab (anti-IL-17A)"][n] for n in lib["names"]]
    b = [cross["etanercept (anti-TNF)"][n] for n in lib["names"]]
    rho2, pv2 = spearmanr(a, b)
    out["drug_ranking_rho"] = {"rho": float(rho2), "p": float(pv2)}
    out["drug_vs_drug_cos"] = _cos(dr["secukinumab (anti-IL-17A)"], dr["etanercept (anti-TNF)"])
    return out


# ================================================================== part B
def part_b(SEC):
    out = {}
    idx, dis = _disease_genes(SEC, SEC["pats"])
    gi = {g: i for i, g in enumerate(SEC["genes"])}
    gl = SEC["genes"]

    def module_idx(names):
        return np.array([gi[g] for g in names if g in gi
                         and abs(dis[gi[g]]) > 0.585], dtype=int)

    # Q5 module restoration at week 16 (all 10 patients, all on secukinumab).
    # Normalised by the COHORT-MEDIAN baseline, not each patient's own: a
    # patient whose own module baseline is near zero otherwise divides by
    # noise (that gave a restoration sd of 46 in the first run).
    rest = {}
    for mn, genes in MOD.items():
        present = [g for g in genes if g in gi]
        mi = module_idx(genes)
        if len(mi) < 3:
            rest[mn] = {"med": np.nan, "sd": np.nan, "n": len(mi),
                        "n_present": len(present), "n_total": len(genes)}
            continue
        u = dis[mi] / np.linalg.norm(dis[mi])
        base = [float((SEC["S"][p][0] - SEC["S"][p]["NL"])[mi] @ u) for p in SEC["pats"]]
        now = [float((SEC["S"][p][16] - SEC["S"][p]["NL"])[mi] @ u) for p in SEC["pats"]]
        b = float(np.median(base))
        vals = [1 - n / b for n in now]
        rest[mn] = {"med": float(np.median(vals)), "sd": float(np.std(vals)), "n": len(mi),
                    "n_present": len(present), "n_total": len(genes)}
    out["restoration_wk16"] = rest

    # Q6 week 32 by arm, on the cohort disease axis (leave-one-out per patient)
    u_all = {}
    for p in SEC["pats"]:
        others = [q for q in SEC["pats"] if q != p]
        u_all[p] = np.array([SEC["S"][q][0] - SEC["S"][q]["NL"] for q in others]).mean(0)
    phi = {}
    for p in SEC["pats"]:
        u = u_all[p][idx] / np.linalg.norm(u_all[p][idx])
        den = float((SEC["S"][p][0] - SEC["S"][p]["NL"])[idx] @ u)
        phi[p] = {w: float((SEC["S"][p][w] - SEC["S"][p]["NL"])[idx] @ u) / den
                  for w in (0, 16, 32) if w in SEC["S"][p]}
    out["phi"] = phi
    out["phi_by_arm"] = {a: {w: float(np.median([phi[p][w] for p in SEC[k] if w in phi[p]]))
                             for w in (0, 16, 32)}
                         for a, k in (("secukinumab", "sec"), ("placebo", "pbo"))}

    # Q7 residual at week 16 vs the rise from week 16 to week 32, by arm.
    # The same statistic in the secukinumab arm is the bias / noise estimate:
    # there the week-32 change is not a relapse.
    cosr = {}
    for p in SEC["pats"]:
        if 32 not in SEC["S"][p]:
            continue
        resid = (SEC["S"][p][16] - SEC["S"][p]["NL"])[idx]
        rise = (SEC["S"][p][32] - SEC["S"][p][16])[idx]
        cosr[p] = {"arm": SEC["arm"][p], "cos": _cos(resid, rise),
                   "resid_norm": float(np.linalg.norm(resid)),
                   "rise_norm": float(np.linalg.norm(rise)),
                   "phi16": phi[p][16], "phi32": phi[p][32]}
    out["resid_vs_rise"] = cosr
    for a, k in (("secukinumab", "sec"), ("placebo", "pbo")):
        v = [cosr[p]["cos"] for p in SEC[k] if p in cosr]
        out[f"cos_{k}"] = (float(np.median(v)), len(v))

    # The arm difference is NOT evidence on its own: the week-16 sample enters
    # the residual with + and the rise with -, so its noise biases the cosine
    # down by an amount that shrinks as |rise| grows. The placebo arm has the
    # larger |rise| by construction, so it gets the smaller bias for free.
    out["rise_norm_by_arm"] = {
        a: float(np.median([cosr[p]["rise_norm"] for p in SEC[k] if p in cosr]))
        for a, k in (("secukinumab", "sec"), ("placebo", "pbo"))}

    # Bias-free version: pair each patient's week-16 residual with OTHER
    # patients' relapse directions. Only the own-pairing carries the shared
    # sample, so own > others would be conservative evidence.
    pbo = [p for p in SEC["pbo"] if p in cosr]
    own, other = [], []
    for p in pbo:
        resid = (SEC["S"][p][16] - SEC["S"][p]["NL"])[idx]
        own.append(_cos(resid, (SEC["S"][p][32] - SEC["S"][p][16])[idx]))
        other += [_cos(resid, (SEC["S"][q][32] - SEC["S"][q][16])[idx])
                  for q in pbo if q != p]
    out["specificity_pbo"] = {"own": float(np.mean(own)), "other": float(np.mean(other)),
                              "delta": float(np.mean(own) - np.mean(other)),
                              "n_own": len(own), "n_other": len(other)}

    # Does a larger week-16 residual go with a larger relapse? n = 5, anecdote.
    from scipy.stats import spearmanr
    r = [cosr[p]["resid_norm"] for p in pbo]
    f = [cosr[p]["phi32"] for p in pbo]
    rho, pv = spearmanr(r, f)
    out["residual_vs_relapse"] = {"rho": float(rho), "p": float(pv), "n": len(pbo),
                                  "resid_norm": r, "phi32": f,
                                  "patients": pbo}
    return out


# ================================================================== report
def report(R):
    A, B = R["A"], R["B"]
    print("=" * 78)
    print("D20  closing the loop on IL-17A (secukinumab, GSE102725)")
    print("=" * 78)
    print(f"\nPART A  {A['n_disease_genes']} disease genes, {A['n_shared']} of them in the "
          f"keratinocyte lever library")
    print(f"  boring baseline: cos(drug effect, -disease axis) = {A['cos_drug_vs_minus_disease']:+.3f}"
          "   (treatment mostly just reverses the disease)")
    print("\n  Q1/Q2  which lever does the observed drug effect look like?")
    print(f"    {'lever':10s} {'source':10s} {'cos':>8s} {'95% CI':>18s} "
          f"{'cos|disease':>12s} {'resid cos':>10s}")
    for r in A["levers"]:
        mark = "  <- the drug's target" if r["lever"] == "IL-17A" else ""
        print(f"    {r['lever']:10s} {r['source']:10s} {r['cos']:+8.3f} "
              f"[{r['ci'][0]:+.3f}, {r['ci'][1]:+.3f}] {r['cos_vs_disease']:+12.3f} "
              f"{r['cos_resid']:+10.3f}{mark}")
    print(f"    top by raw cosine: {A['top']}   |   top after removing the disease axis: "
          f"{A['top_resid']}")
    s = A["il17a_vs_shuffle"]
    print(f"    IL-17A vs gene-shuffled IL-17A: {s['cos']:+.3f} vs 95th pct {s['null_p95']:+.3f}"
          f"  p={s['p']:.4f}")
    print("\n  Q3  least squares within GSE36287 (magnitude-comparable levers)")
    for n, c in zip(A["lstsq"]["levers"], A["lstsq"]["coef"]):
        print(f"    {n:10s} coef {c:+.3f}")
    print(f"    R2 = {A['lstsq']['r2']:.3f}   |   IL-17A alone: coef "
          f"{A['lstsq_il17a_only']['coef']:+.3f}, R2 = {A['lstsq_il17a_only']['r2']:.3f}")
    print(f"\n  Q4  crossed control, {A['n_crossed_genes']} genes that are disease genes in both cohorts")
    for dn, cs in A["crossed"].items():
        order = sorted(cs, key=lambda n: -cs[n])
        print(f"    {dn:28s} ranking: " + ", ".join(f"{n} {cs[n]:+.2f}" for n in order[:4]))
        print(f"    {'':28s} own target rank: IL-17A #{order.index('IL-17A') + 1}, "
              f"TNF #{order.index('TNF') + 1}")
    dd = A["crossed_dod"]
    print(f"    difference-of-differences (IL-17A minus TNF, secukinumab minus etanercept): "
          f"{dd['obs']:+.3f} [{dd['ci'][0]:+.3f}, {dd['ci'][1]:+.3f}] p={dd['p']:.4f}")
    print(f"    cos(secukinumab effect, etanercept effect) = {A['drug_vs_drug_cos']:+.3f}")
    print(f"    lever ranking vs 'looks like the disease': Spearman "
          f"{A['rank_vs_disease_rho']['rho']:+.2f} (p={A['rank_vs_disease_rho']['p']:.4f})")
    print(f"    the two drugs rank the nine levers the same way: Spearman "
          f"{A['drug_ranking_rho']['rho']:+.2f} (p={A['drug_ranking_rho']['p']:.4f})")

    print("\nPART B")
    print("  Q5  module restoration at week 16 (PASI 98 responders, 1 = back to non-lesional)")
    for mn, v in sorted(B["restoration_wk16"].items(),
                        key=lambda kv: -(kv[1]["med"] if kv[1]["med"] == kv[1]["med"] else -9)):
        tag = f"{v['n']} disrupted of {v['n_present']}/{v['n_total']} on the array"
        if v["med"] != v["med"]:
            print(f"    {mn:40s} undefined ({tag})")
        else:
            print(f"    {mn:40s} {v['med']:+.2f}  (sd {v['sd']:.2f}, {tag})")
    print("  Q6  position on the disease axis by arm (median)")
    for a, v in B["phi_by_arm"].items():
        print(f"    {a:12s} week 0 {v[0]:+.2f}   week 16 {v[16]:+.2f}   week 32 {v[32]:+.2f}")
    print("  Q7  cos(residual at week 16, change from week 16 to week 32)")
    for p, v in sorted(B["resid_vs_rise"].items(), key=lambda kv: kv[1]["arm"]):
        print(f"    {p:10s} arm {v['arm']}  cos {v['cos']:+.3f}  phi16 {v['phi16']:+.2f} "
              f"-> phi32 {v['phi32']:+.2f}")
    print(f"    median placebo {B['cos_pbo'][0]:+.3f} (n={B['cos_pbo'][1]})   "
          f"median secukinumab {B['cos_sec'][0]:+.3f} (n={B['cos_sec'][1]})")
    print(f"    but |rise| differs by arm (placebo {B['rise_norm_by_arm']['placebo']:.1f} vs "
          f"secukinumab {B['rise_norm_by_arm']['secukinumab']:.1f}), and the shared-sample bias "
          "shrinks as |rise| grows, so the arm gap is not evidence")
    sp = B["specificity_pbo"]
    print(f"    bias-free pairing inside the placebo arm: own {sp['own']:+.3f} vs "
          f"other patients' relapses {sp['other']:+.3f}  (delta {sp['delta']:+.3f})")
    rr = B["residual_vs_relapse"]
    print(f"    residual size vs relapse size (n={rr['n']}, anecdote): "
          f"Spearman rho {rr['rho']:+.2f}, p={rr['p']:.2f}")


def figures(R):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    A, B = R["A"], R["B"]
    fig, ax = plt.subplots(1, 3, figsize=(16.5, 5.2))

    a = ax[0]
    rows = A["levers"]
    y = np.arange(len(rows))[::-1]
    cols = ["C3" if r["lever"] == "IL-17A" else "C0" for r in rows]
    a.barh(y, [r["cos"] for r in rows], 0.38, color=cols)
    a.barh(y - 0.4, [r["cos_resid"] for r in rows], 0.38,
           color=["C1" if r["lever"] == "IL-17A" else "0.7" for r in rows])
    a.set_yticks(y - 0.2, [r["lever"] for r in rows], fontsize=9)
    a.axvline(0, color="k", lw=0.8)
    a.set_xlabel("cos(observed secukinumab effect, blocking this lever)")
    a.set_title("Q1/Q2  the drug effect vs each measured lever\n"
                "(the drug's own target is highlighted)", fontsize=10)
    import matplotlib.patches as mpatches
    a.legend(handles=[mpatches.Patch(color="C0", label="raw cosine"),
                      mpatches.Patch(color="0.7", label="after removing the disease axis"),
                      mpatches.Patch(color="C3", label="IL-17A, raw"),
                      mpatches.Patch(color="C1", label="IL-17A, disease axis removed")],
             fontsize=7.5, loc="lower right")

    a = ax[1]
    names = [r["lever"] for r in rows]
    w = 0.38
    xs = np.arange(len(names))
    for k, (dn, cs) in enumerate(A["crossed"].items()):
        a.bar(xs + (k - 0.5) * w, [cs[n] for n in names], w, label=dn)
    a.set_xticks(xs, names, rotation=45, ha="right", fontsize=8)
    a.axhline(0, color="k", lw=0.8)
    a.set_ylabel("cos(drug effect, blocking lever)")
    a.set_title("Q4  two drugs, same library:\ndoes the ranking know which drug it was?",
                fontsize=10)
    a.legend(fontsize=8)

    a = ax[2]
    for p, v in B["phi"].items():
        ws = sorted(v)
        a.plot(ws, [v[w] for w in ws], "o-",
               color="C0" if R["SEC_arm"][p] == "V" else "C3", alpha=0.8,
               label=None)
    a.plot([], [], "o-", color="C0", label="secukinumab continued")
    a.plot([], [], "o-", color="C3", label="placebo withdrawal")
    a.axhline(1, color="k", ls=":", lw=1)
    a.axhline(0, color="k", lw=0.8)
    a.set_xticks([0, 16, 32], ["week 0\n(lesion)", "week 16\n(PASI 98)", "week 32"])
    a.set_ylabel("position on the disease axis (1 = own lesion)")
    a.set_title("Q6  randomised withdrawal", fontsize=10)
    a.legend(fontsize=8)

    fig.tight_layout()
    fn = os.path.join(FIG, "d20_secukinumab_closed_loop.png")
    fig.savefig(fn, dpi=140)
    print(f"\nfigure: {fn}")


def main(force=False):
    if os.path.exists(RES) and not force:
        R = pickle.load(open(RES, "rb"))
    else:
        SEC, ETN = load_sec(), load_etn()
        lib = kc.lever_library()
        R = {"A": part_a(SEC, ETN, lib), "B": part_b(SEC), "SEC_arm": SEC["arm"],
             "lib": {"names": lib["names"], "source": lib["source"]}}
        pickle.dump(R, open(RES, "wb"))
    report(R)
    figures(R)
    return R


if __name__ == "__main__":
    main(force=bool(os.environ.get("D20_FORCE")))
