"""
Real-data test of two-stage retrieval (dcsc/twostage.py).

A. Norman 2019  — "cocktail" pursuit: name both genes of 131 double
                  perturbations from a 105-entry library of singles.
B. sci-Plex3    — single drug retrieval, same cell line and across lines.
C. Census       — cell-level k-NN over ~80k real cells, querying with cells
                  from a held-out dataset (other lab / chemistry). This is the
                  setting codes exist for: a big library of cells.

Every stage-1 variant uses codes only; stage 2 rescores the M shortlisted
items on their full vectors. Ground truth: A/B the known perturbation; C the
exact top-10 neighbours in the full 2000-gene space, plus author cell-type and
disease labels.

    python3 -m dcsc.real.twostage_eval            # all parts (~15 min)
    python3 -m dcsc.real.twostage_eval norman     # one part
"""
import json
import sys
import time

import numpy as np
import pandas as pd

from .. import delta as DL
from ..twostage import TwoStageIndex, exact_knn

OUT = "dcsc/real/results_twostage.json"


def _match(got, truth, rel):
    return len(got) == len(truth) and all(any(rel[g, t] for g in got) for t in truth)


# ---------------------------------------------------------------------------
def part_norman():
    from .norman import compute_deltas
    d = compute_deltas()
    singles, doubles = list(d["singles"]), list(d["doubles"])
    L, Q = d["Ls"].astype(np.float32), d["Dd"].astype(np.float32)
    Ln = L / np.linalg.norm(L, axis=1, keepdims=True)
    dup = ((Ln @ Ln.T) > 0.9) | np.eye(len(L), dtype=bool)
    truth = [{singles.index(g) for g in p.split("_")} for p in doubles]

    codes = {"ternary K32 (8 B)": dict(code="ternary", K=32),
             "int8 K32 (32 B)": dict(code="int8", K=32),
             "int8 K64 (64 B)": dict(code="int8", K=64)}
    rows = []

    def score(name, gots, touched):
        rows.append(dict(method=name,
                         pair=np.mean([g == t for g, t in zip(gots, truth)]),
                         pair_dup=np.mean([_match(g, t, dup) for g, t in zip(gots, truth)]),
                         at_least_one=np.mean([len(g & t) > 0 for g, t in zip(gots, truth)]),
                         full_vectors_touched=float(np.mean(touched)) / len(L)))

    score("exact: pursuit over all full vectors",
          [set(DL.residual_pursuit(q, L, None)["stimuli"]) for q in Q], [len(L) * 2] * len(Q))
    for cname, kw in codes.items():
        ix = TwoStageIndex(L, **kw)
        # code only: pursuit on decoded codes in program space, nothing full touched
        if kw["code"] == "int8":
            Zdec = ix.Zhat
        else:
            mag = np.array([np.abs((L @ ix.W)[:, j][ix.T[:, j] != 0]).mean() if (ix.T[:, j] != 0).any() else 0
                            for j in range(ix.W.shape[1])])
            Zdec = ix.T * mag
        score(f"code only, {cname}",
              [set(DL.residual_pursuit(q @ ix.W, Zdec, None)["stimuli"]) for q in Q], [0] * len(Q))
        for M in (4, 8, 16):
            for refresh in (False, True):
                for sym in ([True, False] if kw["code"] == "ternary" else [False]):
                    gots, touched = [], []
                    for q in Q:
                        cnt = [0]
                        orig = ix.shortlist

                        def counted(*a, **k):
                            r = orig(*a, **k)
                            cnt[0] += len(r)
                            return r
                        ix.shortlist = counted
                        gots.append(set(ix.pursuit(q, M=M, refresh=refresh, symmetric=sym)["stimuli"]))
                        ix.shortlist = orig
                        touched.append(cnt[0])
                    score(f"two-stage {cname}, M={M}, {'refresh' if refresh else 'static'}"
                          f"{', symmetric' if sym else ''}", gots, touched)
        print(rows[-1])
    return rows


# ---------------------------------------------------------------------------
def part_sciplex():
    from .sciplex import load_deltas
    dd = load_deltas()
    D1, D2, hv, active, lines, target = dd["D1"], dd["D2"], dd["hv"], dd["active"], dd["lines"], dd["target"]
    rows = []
    for la in lines:
        lib = sorted(active[la])
        L = np.array([D1[la][x][hv] for x in lib], np.float32)
        Ln = L / np.linalg.norm(L, axis=1, keepdims=True)
        idx = {"ternary K32 (8 B)": TwoStageIndex(L, K=32, code="ternary"),
               "int8 K32 (32 B)": TwoStageIndex(L, K=32, code="int8")}
        for lb in lines:
            qs = [x for x in lib if x in active[lb] and x in D2[lb]]
            for q_name in qs:
                q = D2[lb][q_name][hv].astype(np.float32)
                t = lib.index(q_name)
                full_rank = np.argsort(-(Ln @ q))
                rec = dict(cross=la != lb, exact_top1=full_rank[0] == t,
                           exact_target=target[lib[full_rank[0]]] == target[q_name])
                for cname, ix in idx.items():
                    for sym in ([True, False] if "ternary" in cname else [False]):
                        tag = f"{cname}{' sym' if sym else ''}"
                        c1 = ix.code_only(q, 1, symmetric=sym)[0]
                        rec[f"{tag} | code only"] = c1 == t
                        rec[f"{tag} | code only target"] = target[lib[c1]] == target[q_name]
                        for M in (3, 5, 10):
                            sl = ix.shortlist(q, M, symmetric=sym)
                            top, _ = ix.search(q, 1, M=M, symmetric=sym)
                            rec[f"{tag} | M={M} recall"] = t in sl
                            rec[f"{tag} | M={M}"] = top[0] == t
                            rec[f"{tag} | M={M} target"] = target[lib[top[0]]] == target[q_name]
                rows.append(rec)
    df = pd.DataFrame(rows)
    out = df.groupby("cross").mean().T
    out.columns = ["same line", "other line"]
    out["n_queries"] = None
    print(out.round(3).to_string())
    return dict(table=out.drop(columns="n_queries").round(4).to_dict(),
                n_queries=df.groupby("cross").size().to_dict())


# ---------------------------------------------------------------------------
def part_census(held_out=("COMBAT (10x 5' v1)", "SLE_1.2M (10x 3' v2)", "Wilk (Seq-Well)"),
                n_q=300, n_hvg=2000):
    from .census import load
    X, obs, genes, _ = load()
    rng = np.random.default_rng(0)
    results, bio_rows, timing = [], [], []
    for test in held_out:
        lib = np.where((obs.dataset != test).values)[0]
        qry = rng.choice(np.where((obs.dataset == test).values)[0], n_q, replace=False)
        Xl = X[lib]
        mu = np.asarray(Xl.mean(0)).ravel()
        var = np.asarray(Xl.multiply(Xl).mean(0)).ravel() - mu ** 2
        hv = np.argsort(-var)[:n_hvg]
        sd = np.sqrt(var[hv]) + 1e-3
        F = np.clip((X[:, hv].toarray() - mu[hv]) / sd, -10, 10).astype(np.float32)
        Fl, Fq = F[lib], F[qry]
        fit = rng.choice(len(lib), 20000, replace=False)
        truth = [exact_knn(Fl, q, 10) for q in Fq]
        leaf, sick = obs.leaf.values[lib], (obs.disease.values[lib] != "normal")
        q_leaf, q_sick = obs.leaf.values[qry], obs.disease.values[qry] != "normal"
        t_exact = _time(lambda: exact_knn(Fl, Fq[0], 10))
        timing.append(dict(held_out=test, method="exact full scan", ms=t_exact * 1e3, n_lib=len(lib)))
        configs = {"ternary K32 sym (8 B)": (dict(code="ternary", K=32), True),
                   "ternary K32 asym (8 B)": (dict(code="ternary", K=32), False),
                   "ternary K128 asym (32 B)": (dict(code="ternary", K=128), False),
                   "int8 K32 (32 B)": (dict(code="int8", K=32), False),
                   "int8 K64 (64 B)": (dict(code="int8", K=64), False),
                   "faiss PQ m=8 (8 B)": (dict(code="pq", pq_m=8), False),
                   "faiss PQ m=32 (32 B)": (dict(code="pq", pq_m=32), False),
                   "faiss PQ m=64 (64 B)": (dict(code="pq", pq_m=64), False)}
        for name, (kw, sym) in configs.items():
            t0 = time.time()
            ix = TwoStageIndex(Fl, fit_rows=fit, **kw)
            build = time.time() - t0
            rec = {"code only": []}
            for M in (100, 1000, 5000):
                rec[f"M={M}"] = []
            nb = {"code only": [], "M=1000": []}
            for q, tr in zip(Fq, truth):
                tr = set(tr)
                c = ix.code_only(q, 10, symmetric=sym)
                rec["code only"].append(len(tr & set(c)) / 10)
                nb["code only"].append(c)
                for M in (100, 1000, 5000):
                    got, _ = ix.search(q, 10, M=M, symmetric=sym)
                    rec[f"M={M}"].append(len(tr & set(got)) / 10)
                    if M == 1000:
                        nb["M=1000"].append(got)
            results.append(dict(held_out=test, code=name, bytes=ix.bytes, build_s=round(build, 1),
                                **{f"recall@10 {k}": float(np.mean(v)) for k, v in rec.items()}))
            for k, lst in nb.items():
                bio_rows.append(dict(held_out=test, code=name, method=k,
                                     same_leaf=float(np.mean([np.mean(leaf[n] == ql) for n, ql in zip(lst, q_leaf)])),
                                     same_status=float(np.mean([np.mean(sick[n] == qs) for n, qs in zip(lst, q_sick)]))))
            if kw["code"] != "pq":
                t1 = _time(lambda: ix.coarse_scores(Fq[0], symmetric=sym))
            else:
                import faiss
                faiss.omp_set_num_threads(1)
                qn = np.concatenate([Fq[0] / np.linalg.norm(Fq[0]), np.zeros(ix.pad)])[None].astype(np.float32)
                t1 = _time(lambda: ix.pq.search(qn, 1000))
            t2 = _time(lambda: ix.search(Fq[0], 10, M=1000, symmetric=sym)) - (
                _time(lambda: ix.coarse_scores(Fq[0], symmetric=sym)))
            timing.append(dict(held_out=test, method=f"{name}: stage 1", ms=t1 * 1e3, n_lib=len(lib)))
            timing.append(dict(held_out=test, method=f"{name}: stage 2 (M=1000)", ms=max(t2, 0) * 1e3,
                               n_lib=len(lib)))
            print(results[-1])
            del ix
        for k, lst in {"exact": truth}.items():
            bio_rows.append(dict(held_out=test, code="-", method="exact full scan",
                                 same_leaf=float(np.mean([np.mean(leaf[n] == ql) for n, ql in zip(lst, q_leaf)])),
                                 same_status=float(np.mean([np.mean(sick[n] == qs) for n, qs in zip(lst, q_sick)]))))
        bio_rows.append(dict(held_out=test, code="-", method="random library cell",
                             same_leaf=float(np.mean([np.mean(leaf == ql) for ql in q_leaf[:50]])),
                             same_status=float(np.mean([np.mean(sick == qs) for qs in q_sick[:50]]))))
        del F, Fl
    return dict(knn=results, biology=bio_rows, timing=timing)


def part_census_space(held_out=("COMBAT (10x 5' v1)", "SLE_1.2M (10x 3' v2)", "Wilk (Seq-Well)"),
                      n_q=300, n_hvg=2000):
    """Which space should stage 2 rescore in? Exact k-NN in raw z-scored genes
    vs in the top-K principal components (denoised), scored by biology."""
    from .census import load
    X, obs, genes, _ = load()
    rng = np.random.default_rng(0)
    rows = []
    for test in held_out:
        lib = np.where((obs.dataset != test).values)[0]
        qry = rng.choice(np.where((obs.dataset == test).values)[0], n_q, replace=False)
        Xl = X[lib]
        mu = np.asarray(Xl.mean(0)).ravel()
        var = np.asarray(Xl.multiply(Xl).mean(0)).ravel() - mu ** 2
        hv = np.argsort(-var)[:n_hvg]
        sd = np.sqrt(var[hv]) + 1e-3
        F = np.clip((X[:, hv].toarray() - mu[hv]) / sd, -10, 10).astype(np.float32)
        Fl, Fq = F[lib], F[qry]
        fit = rng.choice(len(lib), 20000, replace=False)
        _, _, vt = np.linalg.svd(Fl[fit], full_matrices=False)
        leaf, sick = obs.leaf.values[lib], (obs.disease.values[lib] != "normal")
        q_leaf, q_sick = obs.leaf.values[qry], obs.disease.values[qry] != "normal"
        spaces = {"raw 2000 genes": (Fl, Fq)}
        for K in (16, 32, 64, 200):
            W = vt[:K].T
            spaces[f"top {K} PCs (float)"] = (Fl @ W, Fq @ W)
        for name, (A, B) in spaces.items():
            nb = [exact_knn(A, q, 10) for q in B]
            rows.append(dict(held_out=test, space=name,
                             same_leaf=float(np.mean([np.mean(leaf[n] == ql) for n, ql in zip(nb, q_leaf)])),
                             same_status=float(np.mean([np.mean(sick[n] == qs) for n, qs in zip(nb, q_sick)]))))
            print(rows[-1])
        del F, Fl
    return rows


def _time(fn, reps=10):
    fn()
    t = time.perf_counter()
    for _ in range(reps):
        fn()
    return (time.perf_counter() - t) / reps


def main():
    parts = sys.argv[1:] or ["norman", "sciplex", "census"]
    try:
        res = json.load(open(OUT))
    except FileNotFoundError:
        res = {}
    for p in parts:
        res[p] = {"norman": part_norman, "sciplex": part_sciplex, "census": part_census,
                  "census_space": part_census_space}[p]()
        with open(OUT, "w") as f:
            json.dump(res, f, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))


if __name__ == "__main__":
    main()
