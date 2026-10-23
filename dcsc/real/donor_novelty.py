"""
Real-data test of the "no similar diseased donor in the library" call
(dcsc/donor_match.py) on lab-centred donor profiles from the 1M Census library
plus 3 external labs.

Leave-one-collection-out: for every Census collection that has query donors,
the matcher is fitted and calibrated on donors of all OTHER collections only,
then each donor of the held-out collection is queried. Ground truth per query:

    healthy      disease == "normal"
    matchable    its disease occurs in the fitted library (another collection)
    novel        its disease occurs nowhere else -> the right answer is
                 "no similar diseased donor"

Compared with the previous behaviour (always return the majority disease of
the nearest diseased donors), which by construction mislabels every novel
query it does not call healthy.

    python3 -m dcsc.real.donor_novelty
"""
import json

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from ..donor_match import DonorMatcher
from .donor_search import add_collection, build_profiles, center, features

ALPHAS_MATCH = [0.01, 0.02, 0.05, 0.10, 0.20, 0.30]


def run(profile="program", alpha_abnormal=0.10):
    P = add_collection(build_profiles())
    X, keep = center(features(P, profile), P, "lab-centred")
    dis, coll, ext = P["disease"].astype(str), P["collection"], P["external"]
    recs = []
    for C in np.unique(coll[keep]):
        q = np.where(keep & (coll == C))[0]
        libm = keep & ~ext & (coll != C)
        m = DonorMatcher(X[libm], dis[libm], coll[libm], alpha_abnormal=alpha_abnormal, alpha_match=0.10)
        lib_diseases = set(dis[libm][dis[libm] != "normal"])
        for i in q:
            r = m.query(X[i])
            truth = "healthy" if dis[i] == "normal" else ("matchable" if dis[i] in lib_diseases else "novel")
            recs.append(dict(donor=P["donor"][i], collection=C, external=bool(ext[i]), disease=dis[i], truth=truth,
                             p_abnormal=r["p_abnormal"], p_match=r["p_match"], nearest_disease=r["nearest_disease"],
                             best_similarity=r["best_similarity"]))
    return pd.DataFrame(recs)


def outcomes(df, alpha_abnormal, alpha_match, always_match=False):
    abn = df.p_abnormal < alpha_abnormal
    if always_match:
        call = np.where(abn, "match", "healthy")
    else:
        call = np.where(~abn, "healthy", np.where(df.p_match < alpha_match, "no_similar_disease", "match"))
    correct_label = df.nearest_disease == df.disease
    out = {}
    for t in ["healthy", "matchable", "novel"]:
        s = df.truth == t
        c = call[s]
        out[t] = dict(n=int(s.sum()),
                      called_healthy=float(np.mean(c == "healthy")),
                      match_correct=float(np.mean((c == "match") & correct_label[s])),
                      match_wrong=float(np.mean((c == "match") & ~correct_label[s])),
                      no_similar=float(np.mean(c == "no_similar_disease")))
    return out


def operating_points(df):
    """Outcome shares at several (alpha_abnormal, alpha_match), split into
    'flagged abnormal' and what the novelty layer does among those."""
    corr = df.nearest_disease == df.disease
    rows = []
    for aa in [0.05, 0.10, 0.20]:
        abn = df.p_abnormal < aa
        for am in [None, 0.10, 0.20, 0.30]:
            nov = (df.p_match < am) if am else np.zeros(len(df), bool)
            r = dict(alpha_abnormal=aa, alpha_match=am if am else "off (previous behaviour)",
                     healthy_false_alarm=float(abn[df.truth == "healthy"].mean()))
            for t in ["matchable", "novel"]:
                s = (df.truth == t) & abn
                r[f"{t}_flagged"] = float(abn[df.truth == t].mean())
                r[f"{t}|flagged_match_correct"] = float(((~nov) & corr)[s].mean())
                r[f"{t}|flagged_match_wrong"] = float(((~nov) & ~corr)[s].mean())
                r[f"{t}|flagged_no_similar"] = float(nov[s].mean())
            rows.append(r)
    return rows


def main(out="dcsc/real/results_donor_novelty.json"):
    res = {}
    for profile in ["program", "comp+prog"]:
        df = run(profile)
        r = dict(n_queries=len(df), truth_counts=df.truth.value_counts().to_dict(),
                 novel_diseases=df[df.truth == "novel"].disease.value_counts().to_dict(),
                 matchable_diseases=df[df.truth == "matchable"].disease.value_counts().to_dict())
        # how well does p_match separate matchable from novel among abnormal sick queries?
        sick = df[df.truth != "healthy"]
        r["auroc_p_match_matchable_vs_novel"] = float(roc_auc_score(sick.truth == "matchable", sick.p_match))
        r["auroc_p_abnormal_sick_vs_healthy"] = float(roc_auc_score(df.truth != "healthy", -df.p_abnormal))
        r["operating_points"] = operating_points(df)
        r["previous_behaviour"] = outcomes(df, 0.10, None, always_match=True)
        r["sweep"] = {f"alpha_abnormal=0.05, alpha_match={a}": outcomes(df, 0.05, a) for a in ALPHAS_MATCH}
        r["sweep_abnormal_0.10"] = {f"alpha_abnormal=0.10, alpha_match={a}": outcomes(df, 0.10, a)
                                    for a in [0.05, 0.10, 0.20]}
        r["by_disease_default"] = (df.assign(call=np.where(df.p_abnormal >= 0.10, "healthy",
                                                           np.where(df.p_match < 0.20, "no_similar_disease", "match")))
                                   .groupby(["truth", "disease"]).call.value_counts(normalize=True)
                                   .unstack(fill_value=0).round(3).reset_index().to_dict("records"))
        res[profile] = r
        df.to_csv(f"data_cache/c1m/donor_novelty_{profile.replace('+', '_')}.csv", index=False)
        print(profile, json.dumps({k: r[k] for k in ["truth_counts", "auroc_p_match_matchable_vs_novel",
                                                      "auroc_p_abnormal_sick_vs_healthy"]}))
        print(" previous:", r["previous_behaviour"])
        print(" default (0.10 / 0.20):", outcomes(df, 0.10, 0.20))
    with open(out, "w") as f:
        json.dump(res, f, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))


if __name__ == "__main__":
    main()
