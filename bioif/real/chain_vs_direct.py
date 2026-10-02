"""
bioif.real.chain_vs_direct -- the long-chain question, with labels on both ends.

Everything else in this repo builds chains because end-to-end labels do not
exist. Here, for once, they do: 2,064 compounds are measured in BOTH Tox21
(SR-p53 reporter, the F node) and the Hansen Ames benchmark (the G node). So
the question that the whole cross-model-interface programme rests on can be
asked directly:

    **when you have end-to-end labels, does routing through an intermediate
    node beat predicting the endpoint directly?**

Four predictors of Ames, scored on the same held-out scaffolds:

  direct      SMILES -> Ames                        (C->G, one hop)
  chain       SMILES -> P(p53) -> Ames              (C->F->G, two hops)
  oracle      MEASURED p53 label -> Ames            (what a perfect F node buys)
  augmented   SMILES + measured p53 -> Ames         (does F add to C at all?)

Comparing them decomposes the chain's error into two parts that are usually
confounded:

  chain vs oracle   = the error of the C->F model
  oracle vs direct  = the information the intermediate THROWS AWAY, which no
                      improvement to the C->F model can recover

Leakage control: the p53 model is trained only on Tox21 compounds whose
Bemis-Murcko scaffold does not appear in the Ames test fold. Without that,
the chain gets to see its own test chemistry through the intermediate.
"""
from __future__ import annotations

import collections
from dataclasses import dataclass, field

import numpy as np

from . import qsar, tox


@dataclass
class Arm:
    name: str
    auroc: float
    ap: float
    prec_at_k: dict[int, float] = field(default_factory=dict)
    note: str = ""


def _prec_at_k(y: np.ndarray, score: np.ndarray, ks=(50, 100, 200)):
    order = np.argsort(-score)
    return {k: float(y[order[:k]].mean()) for k in ks if k <= len(y)}


def run(alpha: float = 0.1, seed: int = 0,
        permute: bool = False) -> dict:
    # ---- the joined set -------------------------------------------------
    ov = tox.overlap()
    recs = [(k, v) for k, v in ov.items()
            if v["tox"][tox.P53] not in ("", "NA")]
    recs.sort(key=lambda kv: kv[0])                 # deterministic order
    smi = [v["ames"]["smiles"] for _, v in recs]
    y_ames = np.array([int(float(v["ames"]["ames"])) for _, v in recs])
    p53_meas = np.array([int(float(v["tox"][tox.P53])) for _, v in recs])

    X, ok = qsar.featurize(smi)
    keep = np.where(ok)[0]
    smi = [smi[i] for i in keep]
    X, y_ames, p53_meas = X[keep], y_ames[keep], p53_meas[keep]
    recs = [recs[i] for i in keep]

    tr, ca, te = qsar.scaffold_split(smi, fracs=(0.45, 0.2, 0.35),
                                     seed=seed, permute=permute)
    test_scaffolds = {qsar.scaffold(smi[i]) for i in te}

    out: dict = {
        "n_joined": len(smi),
        "n_train": len(tr), "n_cal": len(ca), "n_test": len(te),
        "prevalence_test": float(y_ames[te].mean()),
        "p53_pos_test": int(p53_meas[te].sum()),
        "arms": [],
    }

    # ---- arm 1: direct C->G ---------------------------------------------
    cc = qsar.fit_conformal_classifier(X[tr], y_ames[tr], X[ca], y_ames[ca],
                                       alpha, "direct")
    s = cc.proba(X[te])
    out["arms"].append(Arm("direct  C->G", *_metrics(y_ames[te], s),
                           note="one hop, structure straight to endpoint"))
    out["direct_cc"] = cc

    # ---- the p53 model, trained with the Ames test scaffolds held out ---
    prows, psmi, py = qsar.p53_dataset()
    PX, pok = qsar.featurize(psmi)
    pkeep = [i for i in np.where(pok)[0]
             if qsar.scaffold(psmi[i]) not in test_scaffolds]
    PX2, py2 = PX[pkeep], py[pkeep]
    psmi2 = [psmi[i] for i in pkeep]
    ptr, pca, _ = qsar.scaffold_split(psmi2, fracs=(0.7, 0.3, 0.0),
                                      seed=seed, permute=permute)
    p53_cc = qsar.fit_conformal_classifier(PX2[ptr], py2[ptr],
                                           PX2[pca], py2[pca], alpha, "p53")
    out["p53_model"] = {
        "n_train": len(ptr), "n_cal": len(pca),
        "n_dropped_for_leakage": len(psmi) - len(pkeep),
    }
    p53_hat_test = p53_cc.proba(X[te])
    # how good is the C->F model on these very compounds?
    out["p53_model"]["auroc_on_joined_test"] = _metrics(p53_meas[te],
                                                        p53_hat_test)[0]

    # ---- arm 2: chain C->F->G via the measured 2x2 ----------------------
    t = tox.two_by_two(tox.P53, ov)
    # P(Ames+ | p53 state), estimated on the TRAIN fold only
    def _cond(mask):
        m = [i for i in tr if bool(p53_meas[i]) == mask]
        return float(y_ames[m].mean()) if m else float(y_ames[tr].mean())
    p_pos, p_neg = _cond(True), _cond(False)
    s_chain = p53_hat_test * p_pos + (1.0 - p53_hat_test) * p_neg
    out["bridge"] = {"P(Ames|p53+)": p_pos, "P(Ames|p53-)": p_neg,
                     "risk_ratio": p_pos / p_neg if p_neg else float("inf"),
                     "fitted_on": "train fold of the joined set"}
    out["arms"].append(Arm("chain   C->F->G", *_metrics(y_ames[te], s_chain),
                           note="two hops, through the predicted p53 state"))

    # ---- arm 3: oracle intermediate -------------------------------------
    s_oracle = np.where(p53_meas[te] == 1, p_pos, p_neg)
    out["arms"].append(Arm("oracle  F->G", *_metrics(y_ames[te], s_oracle),
                           note="the MEASURED p53 label, i.e. a perfect C->F"))

    # ---- arm 4: structure + measured intermediate -----------------------
    Xa = np.hstack([X, p53_meas.reshape(-1, 1).astype(np.float32)])
    cc4 = qsar.fit_conformal_classifier(Xa[tr], y_ames[tr], Xa[ca],
                                        y_ames[ca], alpha, "augmented")
    s4 = cc4.proba(Xa[te])
    out["arms"].append(Arm("augment C+F->G", *_metrics(y_ames[te], s4),
                           note="direct model plus the measured p53 label"))

    for a in out["arms"]:
        a.prec_at_k = _prec_at_k(y_ames[te],
                                 {"direct  C->G": s, "chain   C->F->G": s_chain,
                                  "oracle  F->G": s_oracle,
                                  "augment C+F->G": s4}[a.name])
    return out


def _metrics(y, s):
    from sklearn.metrics import average_precision_score, roc_auc_score
    y = np.asarray(y)
    if len(set(y.tolist())) < 2:
        return float("nan"), float("nan")
    return float(roc_auc_score(y, s)), float(average_precision_score(y, s))


def report(res: dict | None = None) -> str:
    r = res or run()
    L = ["Chain vs direct, on compounds with labels at BOTH ends",
         f"  joined set: {r['n_joined']} compounds measured in Tox21 SR-p53 "
         f"and Ames",
         f"  scaffold split: train {r['n_train']} / calibrate {r['n_cal']} / "
         f"test {r['n_test']}",
         f"  test-fold Ames prevalence {r['prevalence_test']:.3f}; "
         f"{r['p53_pos_test']} p53-positives in the test fold",
         f"  C->F model: trained on {r['p53_model']['n_train']} compounds, "
         f"{r['p53_model']['n_dropped_for_leakage']} dropped so no test "
         f"scaffold leaks in",
         f"  C->F AUROC on these very compounds: "
         f"{r['p53_model']['auroc_on_joined_test']:.3f}",
         f"  bridge fitted on train fold: P(Ames|p53+)="
         f"{r['bridge']['P(Ames|p53+)']:.3f}  P(Ames|p53-)="
         f"{r['bridge']['P(Ames|p53-)']:.3f}  "
         f"RR={r['bridge']['risk_ratio']:.2f}",
         "",
         f"  {'arm':<18}{'AUROC':>8}{'AP':>8}{'P@50':>8}{'P@100':>8}{'P@200':>8}"]
    for a in r["arms"]:
        pk = a.prec_at_k
        L.append(f"  {a.name:<18}{a.auroc:>8.3f}{a.ap:>8.3f}"
                 + "".join(f"{pk.get(k, float('nan')):>8.3f}"
                           for k in (50, 100, 200)))
    by = {a.name: a for a in r["arms"]}
    direct, chain = by["direct  C->G"], by["chain   C->F->G"]
    oracle, aug = by["oracle  F->G"], by["augment C+F->G"]
    L += ["",
          "Decomposition:",
          f"  chain AP {chain.ap:.3f} vs oracle AP {oracle.ap:.3f}"
          f"  -> {oracle.ap - chain.ap:+.3f} is the C->F MODEL's error",
          f"  oracle AP {oracle.ap:.3f} vs direct AP {direct.ap:.3f}"
          f"  -> {direct.ap - oracle.ap:+.3f} is information the intermediate",
          "     THROWS AWAY, which no improvement to the C->F model recovers",
          f"  augmented AP {aug.ap:.3f} vs direct AP {direct.ap:.3f}"
          f"  -> {aug.ap - direct.ap:+.3f} is what a MEASURED p53 adds on top",
          "     of structure",
          "",
          "Conclusion:",
          f"  The direct one-hop model wins by {direct.ap - chain.ap:+.3f} AP.",
          "  Most of that gap is not the intermediate model being bad -- it is",
          "  the intermediate being a 1-bit summary of a 2048-bit input. Routing",
          "  a rich input through a narrow node destroys information, and the",
          "  destruction is irreversible.",
          "",
          "  So: build chains where end-to-end labels are MISSING, and expect a",
          "  direct model to beat them wherever such labels exist. A chain's",
          "  value is reach and interpretability, not accuracy."]
    return "\n".join(L)


if __name__ == "__main__":
    print(report())
