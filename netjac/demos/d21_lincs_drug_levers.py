"""
D21  From cytokine levers to real drugs: does the LINCS L1000 library rank
psoriasis therapies above everything else?

Steps 1-4 of the pipeline describe the state; step 5 asks what moves it back.
D17/D18/D20 used nine measured cytokine responses in primary keratinocytes.
This demo swaps that library for ~2,100 compounds (GSE70138, LINCS Phase 2)
and scores each one the standard Connectivity-Map way:

    connectivity(signature) = -cos(signature, disease displacement)

positive = the compound moves cells AGAINST the lesional displacement.

This is not a new score - it is the CMap logic. What is new here is the
control structure the earlier demos forced on us:

  Q1  POSITIVE CONTROL. 23 compounds in this library are psoriasis therapies
      or share their mechanism (topical steroids, vitamin D, calcineurin, JAK,
      PDE4, retinoids, fumarate), pre-specified from pharmacology, never from
      the ranking. Do they rank above a pre-specified set of 28 drugs with no
      psoriasis indication?
  Q2  MEASURED vs INFERRED. L1000 measures 978 genes and infers 11,350. The
      modules this project works with (barrier, IL-17 response, T cell) have
      ZERO landmark genes. So the whole analysis is run in both spaces and the
      answers are compared; if they disagree, the all-gene answer is a model
      output, not a measurement.
  Q3  CELL TYPE. D17: a lever library from the wrong cell type scored 0/4.
      LINCS has no keratinocyte - it has MCF7, A375, PC3, HT29, HA1E, YAPC,
      HELA. Does the ranking survive from one line to the next?
  Q4  MAGNITUDE CONFOUND. Compounds with violent signatures (HDAC, HSP90,
      proteasome) float to the top of any connectivity ranking. Correlate the
      score with signature norm, and re-run after removing the effect.
  Q5  SPECIFICITY. Rank against atopic dermatitis with the same library. A
      compound that reverses every disease displacement equally is reversing
      "being a sick cell", not this disease.
  Q6  THE SLOW MODULE. Nothing in D14-D20 restored the barrier module. Does
      any compound raise it? Honest answer first: how many barrier genes does
      L1000 even measure?

Run:  python3 -m netjac.demos.d21_lincs_drug_levers
"""
import os
import pickle

import numpy as np

from netjac import DATA_CACHE, FIG
from netjac.data import lincs
from netjac.demos import d14_treatment_restoration as d14
from netjac.demos import d17_keratinocyte_lever as d17

RES = os.path.join(DATA_CACHE, "d21_results.pkl")
MOD = dict(d14.SKIN)
SLOW = "barrier / lipid (down in lesions)"

# ---------------------------------------------------------------- drug sets
# Pre-specified from pharmacology (approved psoriasis therapies and their
# mechanisms), NOT from any ranking in this demo.
POSITIVE = {
    "glucocorticoid": ["clobetasol", "dexamethasone", "dexamethasone-acetate",
                       "triamcinolone", "triamcinolone-acetonide",
                       "fluocinolone-acetonide", "hydrocortisone",
                       "methylprednisolone", "prednisone"],
    "vitamin D receptor": ["calcitriol", "ercalcitriol"],
    "calcineurin": ["cyclosporin-a", "tacrolimus"],
    "JAK": ["tofacitinib", "ruxolitinib", "S-ruxolitinib", "baricitinib"],
    "PDE4": ["apremilast", "roflumilast", "rolipram"],
    "retinoid": ["tazarotene", "tretinoin"],
    "fumarate": ["dimethyl-fumarate"],
}
# Common drugs with no psoriasis indication.
NEGATIVE = ["metformin", "atorvastatin", "simvastatin", "amlodipine", "losartan",
            "lisinopril", "omeprazole", "ciprofloxacin", "sildenafil", "loratadine",
            "warfarin", "furosemide", "levothyroxine", "haloperidol", "fluoxetine",
            "diazepam", "caffeine", "glipizide", "metoprolol", "ranitidine",
            "cimetidine", "digoxin", "verapamil", "pioglitazone", "rosiglitazone",
            "sitagliptin", "allopurinol", "aspirin"]
# Drugs reported to trigger or worsen psoriasis: should not be top reversers.
AGGRAVATOR = ["chloroquine", "propranolol"]
# Immunosuppressants without a first-line psoriasis role: in neither set.
NEITHER = ["azathioprine", "mycophenolic-acid", "sirolimus", "everolimus",
           "thalidomide", "colchicine", "dapsone", "ibuprofen"]
POS_ALL = [d for v in POSITIVE.values() for d in v]
# Machinery check (Q0): drugs that must lower the proliferation module, which
# is the ONLY module here with measured (landmark) genes. If these are not
# recovered, the scoring code or the data is broken and nothing else counts.
ANTIPROLIFERATIVE = ["alisertib", "alvocidib", "barasertib", "camptothecin",
                     "cytarabine", "dinaciclib", "docetaxel", "doxorubicin",
                     "etoposide", "fluorouracil", "gemcitabine", "hydroxyurea",
                     "irinotecan", "paclitaxel", "palbociclib", "roscovitine"]


# ---------------------------------------------------------------- scoring
def _centre(x):
    return x - x.mean()


def connectivity(Z, genes, disease, dgenes):
    """-cos(signature, disease displacement) for every signature."""
    di = {g: i for i, g in enumerate(dgenes)}
    keep = [(j, di[g]) for j, g in enumerate(genes) if g in di]
    zj = np.array([k[0] for k in keep])
    d = _centre(np.asarray(disease, float)[[k[1] for k in keep]])
    S = Z[:, zj].astype(np.float64)
    S = S - S.mean(1, keepdims=True)
    num = S @ d
    den = np.linalg.norm(S, axis=1) * np.linalg.norm(d)
    return -num / np.maximum(den, 1e-12), len(keep)


def updown_score(Z, genes, disease, dgenes, k=100):
    """Rank-based alternative to the cosine, closer to the classic CMap score:
    for each signature, mean rank of the disease-DOWN genes minus mean rank of
    the disease-UP genes, scaled to [-1, 1]. Positive = reverses the disease.
    Run so that the conclusion does not rest on one choice of metric."""
    di = {g: i for i, g in enumerate(dgenes)}
    keep = [(j, di[g]) for j, g in enumerate(genes) if g in di]
    zj = np.array([x[0] for x in keep])
    d = np.asarray(disease, float)[[x[1] for x in keep]]
    k = min(k, len(d) // 3)
    up = zj[np.argsort(-d)[:k]]
    dn = zj[np.argsort(d)[:k]]
    R = np.argsort(np.argsort(Z, axis=1), axis=1).astype(np.float32)
    n = Z.shape[1]
    return (R[:, dn].mean(1) - R[:, up].mean(1)) / (n / 2.0), k


def per_compound(score, sig, by_cell=False):
    """Median score per compound (median within cell line first)."""
    import pandas as pd
    T = pd.DataFrame({"pert": sig.pert_iname.values, "cell": sig.cell_id.values,
                      "s": score})
    cell = T.groupby(["pert", "cell"]).s.median()
    if by_cell:
        return cell
    return cell.groupby("pert").median()


# ---------------------------------------------------------------- statistics
def auc(pos, neg):
    """P(a random positive scores above a random negative)."""
    pos, neg = np.asarray(pos), np.asarray(neg)
    gt = (pos[:, None] > neg[None, :]).sum()
    eq = (pos[:, None] == neg[None, :]).sum()
    return float((gt + 0.5 * eq) / (len(pos) * len(neg)))


def auc_perm_p(values, is_pos, n=20000, seed=0):
    """Label-permutation p for the AUC of the labelled subset."""
    rng = np.random.default_rng(seed)
    v, m = np.asarray(values), np.asarray(is_pos, bool)
    obs = auc(v[m], v[~m])
    null = []
    for _ in range(n):
        p = rng.permutation(len(v))
        null.append(auc(v[p][m], v[p][~m]))
    return obs, float(np.mean(np.array(null) >= obs))


# ---------------------------------------------------------------- analysis
def _disease_gene_symbols(dv, thr=1.0):
    """Genes the diseases actually move, plus every module gene."""
    out = set()
    for v in dv.values():
        g = np.asarray(v["genes"])
        out |= set(g[np.abs(np.asarray(v["d"], float)) > thr])
    for gl in MOD.values():
        out |= set(gl)
    return out


def run(space="landmark", seed=0):
    """space "landmark" = the 978 MEASURED genes, weighted by the displacement.
    space "disease"  = the genes psoriasis / AD actually move (|log2FC| > 1)
    plus the module genes, almost all of which L1000 INFERS rather than
    measures. These are the two choices a practitioner would make; running
    both shows whether the conclusion depends on the choice."""
    sig = lincs.sig_info()
    dv = d17.disease_vectors()
    if space == "landmark":
        M = lincs.signature_matrix(space="landmark")
    else:
        M = lincs.signature_matrix(space="all", symbols=_disease_gene_symbols(dv))
    order = {s: i for i, s in enumerate(M["sig_id"])}
    sig = sig[sig.sig_id.isin(order) & (sig.pert_type == "trt_cp")].copy()
    rows = np.array([order[s] for s in sig.sig_id])
    Z = M["Z"][rows]
    out = {"space": space, "n_sig": len(sig), "n_genes": len(M["genes"]),
           "cells": sorted(set(sig.cell_id))}

    scores, cov = {}, {}
    for ind, v in dv.items():
        scores[ind], cov[ind] = connectivity(Z, M["genes"], v["d"], v["genes"])
    out["n_scored_genes"] = cov

    comp = {ind: per_compound(s, sig) for ind, s in scores.items()}
    psor = comp["psoriasis"]
    out["n_compounds"] = int(len(psor))

    pos = [d for d in POS_ALL if d in psor.index]
    neg = [d for d in NEGATIVE if d in psor.index]
    agg = [d for d in AGGRAVATOR if d in psor.index]
    out["sets"] = {"positive": pos, "negative": neg, "aggravator": agg}

    # Q1 positive control
    lab = psor.loc[pos + neg]
    a, p = auc_perm_p(lab.values, [True] * len(pos) + [False] * len(neg), seed=seed)
    out["q1"] = {"auc_pos_vs_neg": a, "p": p,
                 "auc_pos_vs_all": auc(psor.loc[pos].values,
                                       psor.drop(index=pos).values),
                 "median_pos": float(psor.loc[pos].median()),
                 "median_neg": float(psor.loc[neg].median()),
                 "median_all": float(psor.median()),
                 "median_aggravator": float(psor.loc[agg].median()) if agg else np.nan}
    ud, k = updown_score(Z, M["genes"], dv["psoriasis"]["d"], dv["psoriasis"]["genes"])
    udc = per_compound(ud, sig)
    a2, p2 = auc_perm_p(udc.loc[pos + neg].values,
                        [True] * len(pos) + [False] * len(neg), seed=seed)
    from scipy.stats import spearmanr as _sp
    out["q1_updown"] = {"k": k, "auc": a2, "p": p2,
                        "rho_with_cosine": float(_sp(psor.loc[udc.index].values,
                                                     udc.values)[0]),
                        "top10": [str(i) for i in
                                  udc.sort_values(ascending=False).head(10).index]}
    out["q1_by_class"] = {k: float(psor.loc[[d for d in v if d in psor.index]].median())
                          for k, v in POSITIVE.items()}
    out["top20"] = [(str(i), float(x)) for i, x in psor.sort_values(ascending=False)
                    .head(20).items()]
    out["positive_ranks"] = {d: int((psor > psor[d]).sum()) + 1 for d in pos}

    # Q0 machinery check: does the proliferation module behave?
    gi = {g: i for i, g in enumerate(M["genes"])}
    prol = [g for g in MOD["proliferation"] if g in gi]
    out["q0"] = {"genes": prol}
    if len(prol) >= 3:
        pm = per_compound(Z[:, [gi[g] for g in prol]].mean(1), sig)
        cyt = [d for d in ANTIPROLIFERATIVE if d in pm.index]
        out["q0"].update({
            "n_cytotoxic": len(cyt),
            "auc_cytotoxic_lowers_proliferation":
                auc(-pm.loc[cyt].values, -pm.loc[neg].values),
            "median_cytotoxic": float(pm.loc[cyt].median()),
            "median_negative": float(pm.loc[neg].median()),
            "median_all": float(pm.median())})

    # Q3 cell lines
    cellwise = per_compound(scores["psoriasis"], sig, by_cell=True).unstack()
    big = [c for c in cellwise.columns if cellwise[c].notna().sum() >= 500]
    out["q3"] = {"cells": big,
                 "auc": {c: auc(cellwise[c].loc[pos].dropna().values,
                                cellwise[c].loc[neg].dropna().values) for c in big},
                 "between_cell_spearman": float(
                     cellwise[big].corr(method="spearman").values[
                         np.triu_indices(len(big), 1)].mean())}

    # Q4 magnitude confound
    import pandas as pd
    norm = pd.Series(np.linalg.norm(Z.astype(np.float64), axis=1),
                     index=sig.pert_iname.values).groupby(level=0).median()
    from scipy.stats import spearmanr
    j = psor.index.intersection(norm.index)
    rho, pv = spearmanr(psor.loc[j].values, norm.loc[j].values)
    out["q4"] = {"rho_score_vs_norm": float(rho), "p": float(pv),
                 "median_norm_pos": float(norm.loc[pos].median()),
                 "median_norm_all": float(norm.median()),
                 # does the drug do anything at all in an unstimulated line?
                 "norm_by_class": {k: float(norm.loc[[d for d in v if d in norm.index]]
                                            .median()) for k, v in POSITIVE.items()},
                 "norm_cytotoxic": float(norm.loc[[d for d in ANTIPROLIFERATIVE
                                                   if d in norm.index]].median())}

    # Q5 specificity: psoriasis vs atopic dermatitis
    ad = comp["atopic dermatitis"]
    j = psor.index.intersection(ad.index)
    rho2, pv2 = spearmanr(psor.loc[j].values, ad.loc[j].values)
    out["q5"] = {"rho_psoriasis_vs_AD": float(rho2), "p": float(pv2),
                 "pos_delta": float((psor.loc[pos] - ad.loc[pos]).median()),
                 "all_delta": float((psor.loc[j] - ad.loc[j]).median()),
                 "top20_overlap": len(set(x for x, _ in out["top20"]) &
                                      set(map(str, ad.sort_values(ascending=False)
                                              .head(20).index)))}

    # Q6 the slow module
    G = lincs.gene_info()
    lm = set(G[G.pr_is_lm == "1"].pr_gene_symbol)
    allg = set(G.pr_gene_symbol)
    out["q6_coverage"] = {mn: {"measured": sum(g in lm for g in gl),
                               "measured_or_inferred": sum(g in allg for g in gl),
                               "total": len(gl)} for mn, gl in MOD.items()}
    present = [g for g in MOD[SLOW] if g in set(M["genes"])]
    if present:
        gi = {g: i for i, g in enumerate(M["genes"])}
        bar = Z[:, [gi[g] for g in present]].mean(1)
        bs = per_compound(bar, sig)
        out["q6"] = {"genes_used": present,
                     "top": [(str(i), float(x)) for i, x in
                             bs.sort_values(ascending=False).head(10).items()],
                     "positive_median": float(bs.loc[pos].median()),
                     "all_median": float(bs.median())}
    else:
        out["q6"] = {"genes_used": [], "note": "no barrier gene in this space"}
    out["_psor"] = {str(k): float(v) for k, v in psor.items()}
    out["_ad"] = {str(k): float(v) for k, v in ad.items()}
    return out


# ---------------------------------------------------------------- report
def report(R):
    for space in ("landmark", "disease"):
        A = R[space]
        print("=" * 78)
        print(f"D21  LINCS drug ranking for psoriasis - gene space: {space} "
              f"({A['n_genes']} genes, {A['n_scored_genes']['psoriasis']} shared "
              f"with the disease vector)")
        print("=" * 78)
        print(f"  {A['n_sig']} compound signatures, {A['n_compounds']} compounds, "
              f"{len(A['cells'])} cell lines")
        z = A.get("q0", {})
        if "auc_cytotoxic_lowers_proliferation" in z:
            print(f"\n  Q0  machinery check on the proliferation module "
                  f"({len(z['genes'])} genes: {', '.join(z['genes'])})")
            print(f"      {z['n_cytotoxic']} antiproliferatives vs the unrelated-drug set: "
                  f"AUC {z['auc_cytotoxic_lowers_proliferation']:.3f} "
                  f"(median {z['median_cytotoxic']:+.3f} vs {z['median_negative']:+.3f}, "
                  f"all {z['median_all']:+.3f})")
        q = A["q1"]
        print(f"\n  Q1  positive control ({len(A['sets']['positive'])} psoriasis drugs vs "
              f"{len(A['sets']['negative'])} unrelated drugs)")
        print(f"      AUC {q['auc_pos_vs_neg']:.3f}  (permutation p = {q['p']:.4f})"
              f"   AUC vs all compounds {q['auc_pos_vs_all']:.3f}")
        print(f"      median connectivity: positives {q['median_pos']:+.3f}, "
              f"negatives {q['median_neg']:+.3f}, all {q['median_all']:+.3f}, "
              f"aggravators {q['median_aggravator']:+.3f}")
        print("      by mechanism class: " +
              ", ".join(f"{k} {v:+.3f}" for k, v in A["q1_by_class"].items()))
        print("      ranks of the psoriasis drugs (of "
              f"{A['n_compounds']}): " +
              ", ".join(f"{d} #{r}" for d, r in
                        sorted(A["positive_ranks"].items(), key=lambda kv: kv[1])[:8]))
        u = A["q1_updown"]
        print(f"      rank-based (CMap-style, top {u['k']} up / {u['k']} down) instead of "
              f"cosine: AUC {u['auc']:.3f} (p = {u['p']:.4f}); the two metrics agree "
              f"at Spearman {u['rho_with_cosine']:+.2f}")
        print("      top 20 by connectivity: " +
              ", ".join(n for n, _ in A["top20"]))
        print(f"\n  Q3  cell lines: AUC " +
              ", ".join(f"{c} {v:.2f}" for c, v in A["q3"]["auc"].items()))
        print(f"      mean between-cell-line Spearman of the whole ranking: "
              f"{A['q3']['between_cell_spearman']:+.3f}")
        print("      median |z| by mechanism class: " +
              ", ".join(f"{k} {v:.0f}" for k, v in A["q4"]["norm_by_class"].items()) +
              f", antiproliferatives {A['q4']['norm_cytotoxic']:.0f}, "
              f"library {A['q4']['median_norm_all']:.0f}")
        print(f"\n  Q4  score vs signature magnitude: Spearman "
              f"{A['q4']['rho_score_vs_norm']:+.3f} (p={A['q4']['p']:.2g}); "
              f"median |z| positives {A['q4']['median_norm_pos']:.1f} vs all "
              f"{A['q4']['median_norm_all']:.1f}")
        print(f"\n  Q5  psoriasis vs atopic dermatitis ranking: Spearman "
              f"{A['q5']['rho_psoriasis_vs_AD']:+.3f}; top-20 overlap "
              f"{A['q5']['top20_overlap']}/20; positives' psoriasis-minus-AD "
              f"{A['q5']['pos_delta']:+.3f} vs all compounds {A['q5']['all_delta']:+.3f}")
        print()
    A = R["landmark"]
    print("  Q6  module coverage in L1000 (measured / measured+inferred / total)")
    for mn, c in A["q6_coverage"].items():
        print(f"      {mn:42s} {c['measured']:2d} / {c['measured_or_inferred']:2d} "
              f"/ {c['total']}")
    B = R["disease"]["q6"]
    if B.get("genes_used"):
        print(f"      barrier genes available (inferred): {', '.join(B['genes_used'])}")
        print(f"      compounds raising the barrier module most: " +
              ", ".join(n for n, _ in B["top"][:8]))
        print(f"      psoriasis drugs median {B['positive_median']:+.3f} vs all "
              f"{B['all_median']:+.3f}")


def figures(R):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    A = R["landmark"]
    fig, ax = plt.subplots(1, 4, figsize=(21, 5.2))

    a = ax[3]
    z = A["q0"]
    a.bar(range(3), [z["median_cytotoxic"], z["median_negative"], z["median_all"]],
          color=["C3", "C0", "0.7"])
    a.set_xticks(range(3), ["16 anti-\nproliferatives", "28 unrelated\ndrugs",
                            "all 1768"], fontsize=9)
    a.axhline(0, color="k", lw=0.8)
    a.set_ylabel("proliferation module (mean z)")
    a.set_title(f"Q0  machinery check, measured genes\nAUC "
                f"{z['auc_cytotoxic_lowers_proliferation']:.2f} — the code works",
                fontsize=10)

    a = ax[0]
    psor = A["_psor"]
    allv = np.array(list(psor.values()))
    pos = np.array([psor[d] for d in A["sets"]["positive"]])
    neg = np.array([psor[d] for d in A["sets"]["negative"]])
    a.hist(allv, bins=60, color="0.8", label=f"all {len(allv)} compounds")
    for v, c, lab in ((pos, "C2", "psoriasis drugs"), (neg, "C3", "unrelated drugs")):
        a.plot(v, np.full(len(v), a.get_ylim()[1] * 0.06), "|", color=c, ms=14,
               mew=2, label=lab)
    a.axvline(0, color="k", lw=0.8)
    a.set_xlabel("connectivity  (−cos with the psoriasis displacement)")
    a.set_ylabel("compounds")
    a.set_title(f"Q1  positive control, measured genes only\n"
                f"AUC {A['q1']['auc_pos_vs_neg']:.2f} (p={A['q1']['p']:.3f}); "
                f"rank-based metric {A['q1_updown']['auc']:.2f}", fontsize=10)
    a.legend(fontsize=8)

    a = ax[1]
    cells = A["q3"]["cells"]
    vals = [A["q3"]["auc"][c] for c in cells]
    a.bar(range(len(cells)), vals, color="C0")
    a.axhline(0.5, color="k", ls=":", lw=1)
    a.set_xticks(range(len(cells)), cells, rotation=45, ha="right", fontsize=9)
    a.set_ylabel("AUC of psoriasis drugs vs unrelated drugs")
    a.set_title("Q3  the answer depends on the cell line\n(none is a keratinocyte)",
                fontsize=10)

    a = ax[2]
    ad = R["landmark"]["_ad"]
    shared = sorted(set(psor) & set(ad))
    a.plot([psor[k] for k in shared], [ad[k] for k in shared], ".", ms=3,
           color="0.7", alpha=0.6)
    px = [psor[d] for d in A["sets"]["positive"]]
    py = [ad[d] for d in A["sets"]["positive"]]
    a.plot(px, py, "o", color="C2", ms=5, label="psoriasis drugs")
    xs = [psor[k] for k in shared] + [ad[k] for k in shared]
    lim = [float(np.percentile(xs, 0.2)), float(np.percentile(xs, 99.8))]
    a.plot(lim, lim, "k:", lw=1)
    a.set_xlim(lim)
    a.set_ylim(lim)
    a.set_xlabel("connectivity vs psoriasis")
    a.set_ylabel("connectivity vs atopic dermatitis")
    a.set_title(f"Q5  specificity: Spearman "
                f"{A['q5']['rho_psoriasis_vs_AD']:+.2f}", fontsize=10)
    a.legend(fontsize=8)

    fig.tight_layout()
    fn = os.path.join(FIG, "d21_lincs_drug_levers.png")
    fig.savefig(fn, dpi=140)
    print(f"\nfigure: {fn}")


def main(force=False):
    if os.path.exists(RES) and not force:
        R = pickle.load(open(RES, "rb"))
    else:
        R = {s: run(space=s) for s in ("landmark", "disease")}
        pickle.dump(R, open(RES, "wb"))
    report(R)
    figures(R)
    return R


if __name__ == "__main__":
    main(force=bool(os.environ.get("D21_FORCE")))
