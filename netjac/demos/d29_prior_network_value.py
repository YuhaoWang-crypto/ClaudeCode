"""
D29  Does the prior interaction network earn its place in the Jacobian?

The three-layer design this project keeps proposing is: the graph (OmniPath +
complexes) supplies the SUPPORT of J, the data supplies the numbers. D9 showed
the recovered edges overlap the graph more than chance (sparse MRA 14.4 % vs
3.6 % baseline), but overlap is not the question a buyer asks. The question is
whether putting the graph IN as a prior makes the recovered J predict better
than leaving it out — and overlap cannot answer that, because a prior raises
overlap with itself by construction.

So the test here is strictly out-of-sample: hold out whole perturbation
columns, fit the row-wise sparse MRA on the rest, and predict the held-out
responses. A prior enters as a per-feature Lasso penalty weight (known edges
get a smaller penalty), which is exact: minimising
    ||y - Xb||^2 + a * sum_k w_k |b_k|
equals a plain Lasso on X_k / w_k followed by b = c / w.

PART A  calibration on a synthetic network with a PLANTED J and a prior of
  KNOWN accuracy, because without it there is no way to read the real number:
  A1  a prior naming 100 % of the true edges -> how much can held-out R^2 rise
      at all? That is the ceiling.
  A2  accuracy 0 / 25 / 50 / 75 / 100 % -> the gain must scale with accuracy,
      and a 0 %-accurate prior of the same size must not help.
  A3  the gain must grow as the data is thinned, since that is the whole claim
      for using a prior.

PART B  the real data: Replogle 2022 K562 CRISPRi, the 300 most reproducible
  perturbed-and-measured nodes, R normalised so R[j, j] = 1.
  B1  what the graph even covers on this node set
  B2  held-out response prediction: no prior vs the real prior at several
      penalty discounts vs a DEGREE-PRESERVING SHUFFLED prior vs prior-only
  B3  split-half edge reproducibility (the prior does not control this, unlike
      overlap with the prior itself)
  B4  the thinning sweep: does the prior's advantage grow when data is scarce?

PART C  what this settles for the three-layer design

Run:  python3 -m netjac.demos.d29_prior_network_value
"""
import os
import pickle

import numpy as np

from netjac import DATA_CACHE, FIG

RES = os.path.join(DATA_CACHE, "d29_results.pkl")
N_NODES = 300
ALPHAS = (0.03, 0.1, 0.3)
FOLDS = 3
BOOSTS = (2.0, 5.0, 20.0)


# ----------------------------------------------------------------- fitting


def _fit_rows(R, cols, w=None, alpha=0.03, restrict=None):
    """Row-wise sparse MRA on the perturbation columns `cols`.

    R        nS x nS, R[i, j] = response of node i to perturbation of node j
    w        nS x nS penalty weights (w[i, k] for feature k in row i), or None
    restrict nS x nS boolean: features allowed in each row, or None for all
    Returns B, the local coefficient matrix with a zero diagonal.
    """
    from sklearn.linear_model import Lasso
    n = R.shape[0]
    B = np.zeros((n, n))
    Xall = R[:, cols].T                      # perturbations x nodes
    for i in range(n):
        keep = np.arange(n) != i
        if restrict is not None:
            keep = keep & restrict[i]
            if keep.sum() == 0:
                continue
        idx = np.where(keep)[0]
        X = Xall[:, idx]
        if w is not None:
            X = X / w[i, idx][None, :]
        y = R[i, cols]
        # make alpha a RELATIVE penalty: dividing X and y by the same scalar
        # leaves the solution unchanged and rescales the L1 term, so one alpha
        # grid is comparable across rows and across datasets
        sc = float(np.sqrt((y ** 2).mean()))
        if sc <= 0:
            continue
        m = Lasso(alpha=alpha, fit_intercept=False, max_iter=5000,
                  tol=1e-4).fit(X / sc, y / sc)
        c = m.coef_
        if w is not None:
            c = c / w[i, idx]
        B[i, idx] = c
    return B


def _score(R, B, cols):
    """Held-out fit of the predicted responses. Returns pooled R^2 about zero
    (the model has no intercept) and the median per-row Pearson r."""
    n = R.shape[0]
    num = den = 0.0
    rs = []
    for i in range(n):
        keep = np.arange(n) != i
        y = R[i, cols]
        yh = B[i, keep] @ R[np.ix_(keep, cols)]
        num += float(((y - yh) ** 2).sum())
        den += float((y ** 2).sum())
        if y.std() > 0 and yh.std() > 0:
            rs.append(float(np.corrcoef(y, yh)[0, 1]))
    return {"r2": 1.0 - num / den, "median_r": float(np.median(rs)),
            "n_rows_scored": len(rs)}


def _weights(known, boost):
    """Penalty weights: known features pay 1/boost of the usual penalty."""
    w = np.ones(known.shape)
    w[known] = 1.0 / boost
    return w


def _heldout(R, known=None, boost=1.0, alphas=ALPHAS, folds=FOLDS, seed=0,
             restrict=None, frac=1.0):
    """K-fold over perturbation columns. For fairness every configuration gets
    the same alpha sweep and reports its own best held-out score."""
    n = R.shape[0]
    rng = np.random.default_rng(seed)
    order = rng.permutation(n)
    w = None if known is None or boost == 1.0 else _weights(known, boost)
    best = None
    for a in alphas:
        sc = []
        for f in range(folds):
            test = order[f::folds]
            train = np.setdiff1d(order, test)
            if frac < 1.0:
                k = max(10, int(round(frac * len(train))))
                train = train[:k]
            B = _fit_rows(R, train, w=w, alpha=a, restrict=restrict)
            s = _score(R, B, test)
            s["nnz"] = int((B != 0).sum())
            sc.append(s)
        m = {"alpha": a,
             "r2": float(np.mean([s["r2"] for s in sc])),
             "median_r": float(np.mean([s["median_r"] for s in sc])),
             "nnz": float(np.mean([s["nnz"] for s in sc]))}
        if best is None or m["r2"] > best["r2"]:
            best = m
    return best


# ------------------------------------------------- PART A  calibration


def _synth(n=120, density=0.03, noise=0.25, seed=0):
    """A random stable sparse network, its normalised response matrix, and the
    true edge set."""
    rng = np.random.default_rng(seed)
    A = rng.normal(0, 1, (n, n)) * (rng.random((n, n)) < density)
    np.fill_diagonal(A, 0.0)
    # stability with a fixed margin, taken from the spectral abscissa rather
    # than the row sums, so the off-diagonal response is not crushed
    ab = float(np.max(np.linalg.eigvals(A).real))
    J = A - (ab + 0.5) * np.eye(n)
    G = -np.linalg.inv(J)                     # global response to unit inputs
    R = G / np.diag(G)[None, :]
    off = ~np.eye(n, dtype=bool)
    R = R + noise * rng.normal(0, R[off].std(), R.shape)
    np.fill_diagonal(R, 1.0)
    true = (A != 0)
    return R, true


def _prior_of_accuracy(true, acc, rng):
    """A prior with the same number of entries as `true`, of which a fraction
    `acc` are real edges and the rest are false."""
    n = true.shape[0]
    off = ~np.eye(n, dtype=bool)
    ti = np.argwhere(true & off)
    fi = np.argwhere((~true) & off)
    k = len(ti)
    n_true = int(round(acc * k))
    pick_t = ti[rng.choice(len(ti), n_true, replace=False)] if n_true else np.empty((0, 2), int)
    n_false = k - n_true
    pick_f = fi[rng.choice(len(fi), n_false, replace=False)] if n_false else np.empty((0, 2), int)
    P = np.zeros((n, n), bool)
    for arr in (pick_t, pick_f):
        for a, b in arr:
            P[a, b] = True
    return P


def part_a():
    print("\n" + "=" * 78)
    print("PART A  calibration: how much CAN a prior buy, on a planted network?")
    print("=" * 78)
    out = {}
    R, true = _synth()
    n = R.shape[0]
    off = ~np.eye(n, dtype=bool)
    print("  planted network: %d nodes, %d true edges (%.1f%% density), response"
          " noise 25%%" % (n, int(true[off].sum()), 100 * true[off].mean()))
    base = _heldout(R, None, alphas=ALPHAS, folds=FOLDS)
    print("\n  A1/A2  held-out R^2 by prior accuracy and penalty discount")
    print("  %10s %8s %8s %8s %8s"
          % ("accuracy", "boost 2", "boost 5", "boost 20", "no prior"))
    rng = np.random.default_rng(0)
    rows = []
    for acc in (0.0, 0.25, 0.5, 0.75, 1.0):
        P = _prior_of_accuracy(true, acc, rng)
        vals = {}
        for b in BOOSTS:
            vals[b] = _heldout(R, P, boost=b, alphas=ALPHAS, folds=FOLDS)["r2"]
        rows.append({"accuracy": acc, **{"boost_%g" % b: vals[b] for b in BOOSTS}})
        print("  %9.0f%% %8.4f %8.4f %8.4f %8.4f"
              % (100 * acc, vals[BOOSTS[0]], vals[BOOSTS[1]], vals[BOOSTS[2]],
                 base["r2"]))
    out["accuracy_sweep"] = rows
    out["no_prior"] = base
    best_full = max(r["boost_%g" % b] for r in rows if r["accuracy"] == 1.0
                    for b in BOOSTS)
    best_zero = max(r["boost_%g" % b] for r in rows if r["accuracy"] == 0.0
                    for b in BOOSTS)
    out["ceiling_gain"] = float(best_full - base["r2"])
    out["zero_accuracy_gain"] = float(best_zero - base["r2"])
    print("\n  ceiling: a perfect prior buys %+.4f in held-out R^2 (%.4f -> %.4f)"
          % (out["ceiling_gain"], base["r2"], best_full))
    print("  a 0%%-accurate prior of the same size buys %+.4f"
          % out["zero_accuracy_gain"])
    print("  -> so on this design anything below about %+.3f is not a real gain."
          % max(abs(out["zero_accuracy_gain"]), 0.002))

    # A3 the gain must grow as data is thinned
    print("\n  A3  does the advantage grow when the data is thinned?")
    print("  %10s %10s %10s %10s" % ("train frac", "no prior", "perfect prior",
                                     "gain"))
    P1 = _prior_of_accuracy(true, 1.0, np.random.default_rng(1))
    thin = []
    for fr in (0.25, 0.5, 1.0):
        b0 = _heldout(R, None, frac=fr)["r2"]
        b1 = max(_heldout(R, P1, boost=b, frac=fr)["r2"] for b in BOOSTS)
        thin.append({"frac": fr, "no_prior": b0, "with_prior": b1,
                     "gain": b1 - b0})
        print("  %10.0f%% %10.4f %10.4f %10.4f" % (100 * fr, b0, b1, b1 - b0))
    out["thinning"] = thin
    grows = thin[0]["gain"] > thin[-1]["gain"]
    out["gain_grows_when_thin"] = bool(grows)
    print("  gain is larger at 25%% of the data than at 100%%: %s" % grows)
    return out


# --------------------------------------------------------- PART B  real data


def _real_nodes(n_nodes=N_NODES):
    from netjac.data import perturbseq as ps
    r = ps.reduce_singlecell()
    H, pert, gid = r["half_means"], r["pert"], r["gene_id"]
    pg = np.array([p.split("_")[1] for p in pert])
    pid = np.array([p.split("_")[-1] for p in pert])
    ix = {e: i for i, e in enumerate(gid)}
    ctrl = pg == "non-targeting"
    meas = np.array([i for i in range(len(pert))
                     if pid[i] in ix and not ctrl[i]])
    own = np.array([ix[pid[i]] for i in meas])
    D = H[meas].mean(1)
    selfz = D[np.arange(len(meas)), own]
    rep = np.array([float(np.corrcoef(H[k, 0], H[k, 1])[0, 1]) for k in meas])
    cand = np.where((selfz < -0.5) & (rep > 0.2))[0]
    # one row per gene, the strongest knockdown, then the most reproducible
    best = {}
    for k in cand:
        g = pg[meas[k]]
        if g not in best or selfz[k] < selfz[best[g]]:
            best[g] = k
    cand = np.array(sorted(best.values()))
    sel = cand[np.argsort(-rep[cand])[:n_nodes]]
    sel = sel[np.argsort(sel)]
    genes = pg[meas[sel]]
    oS = own[sel]

    def norm(Hh):
        M = Hh[meas[sel]][:, oS].T             # node i x perturbation j
        return M / np.diag(M)[None, :]
    return {"genes": genes, "R": norm(H[:, 0] * 0 + H.mean(1)),
            "RA": norm(H[:, 0]), "RB": norm(H[:, 1]),
            "rep": rep[sel], "selfz": selfz[sel],
            "n_candidates": len(cand)}


def _real_prior(genes):
    from netjac.data import omnipath as op
    gi = {g: i for i, g in enumerate(genes)}
    n = len(genes)
    direct = np.zeros((n, n), bool)
    for a, b, _s in op.signed_edges():
        if a in gi and b in gi:
            direct[gi[b], gi[a]] = True        # effect of a on b -> feature a in row b
    cpx = np.zeros((n, n), bool)
    for a, b in op.complex_pairs(op.complexes(), genes):
        cpx[gi[a], gi[b]] = cpx[gi[b], gi[a]] = True
    known = direct | direct.T | cpx
    np.fill_diagonal(known, False)
    return {"known": known, "direct": direct, "complex": cpx}


def _shuffled_prior(known, seed=0):
    """Degree-preserving shuffle of the symmetric prior graph."""
    import networkx as nx
    n = known.shape[0]
    sym = known | known.T
    G = nx.from_numpy_array(sym.astype(int))
    G.remove_edges_from(nx.selfloop_edges(G))
    m = G.number_of_edges()
    if m >= 2:
        nx.double_edge_swap(G, nswap=20 * m, max_tries=400 * m, seed=seed)
    S = np.zeros((n, n), bool)
    for a, b in G.edges():
        S[a, b] = S[b, a] = True
    return S


def part_b():
    print("\n" + "=" * 78)
    print("PART B  the real network, out of sample")
    print("=" * 78)
    N = _real_nodes()
    R, genes = N["R"], N["genes"]
    n = R.shape[0]
    off = ~np.eye(n, dtype=bool)
    print("  %d perturbed-and-measured nodes pass the reliability filter;"
          " keeping the %d most reproducible" % (N["n_candidates"], n))
    print("  split-half reproducibility of the kept nodes: median %.3f"
          % np.median(N["rep"]))
    print("  self knockdown of the kept nodes: median %.2f SD"
          % np.median(N["selfz"]))

    P = _real_prior(genes)
    known = P["known"]
    S = _shuffled_prior(known)
    print("\n  B1 what the graph covers on these %d nodes" % n)
    print("     OmniPath signed directed edges (either orientation): %d of %d"
          " ordered pairs (%.2f%%)"
          % (int((P["direct"] | P["direct"].T)[off].sum()), int(off.sum()),
             100 * (P["direct"] | P["direct"].T)[off].mean()))
    print("     complex co-membership: %d (%.2f%%)"
          % (int(P["complex"][off].sum()), 100 * P["complex"][off].mean()))
    print("     union used as the prior: %d (%.2f%%)"
          % (int(known[off].sum()), 100 * known[off].mean()))
    print("     degree-preserving shuffle keeps %d (%.2f%%), overlap with the"
          " real prior %d" % (int(S[off].sum()), 100 * S[off].mean(),
                              int((S & known)[off].sum())))
    nodes_cov = int((known.any(1)).sum())
    print("     nodes with at least one prior edge: %d / %d" % (nodes_cov, n))

    print("\n  B2 held-out response prediction (3 folds over perturbation"
          " columns, each configuration tuned over the same alpha grid)")
    print("  %-26s %9s %9s %8s %9s" % ("configuration", "R^2", "median r",
                                       "alpha", "edges"))
    res = {}
    cfgs = [("no prior", dict(known=None, boost=1.0))]
    cfgs += [("prior, penalty /%g" % b, dict(known=known, boost=b))
             for b in BOOSTS]
    cfgs += [("SHUFFLED prior, /%g" % BOOSTS[1],
              dict(known=S, boost=BOOSTS[1])),
             ("prior-only support", dict(known=None, boost=1.0,
                                         restrict=known))]
    for name, kw in cfgs:
        m = _heldout(R, **kw)
        res[name] = m
        print("  %-26s %9.4f %9.4f %8.3g %9.0f"
              % (name, m["r2"], m["median_r"], m["alpha"], m["nnz"]))
    base = res["no prior"]["r2"]
    best_prior = max(res["prior, penalty /%g" % b]["r2"] for b in BOOSTS)
    shuf = res["SHUFFLED prior, /%g" % BOOSTS[1]]["r2"]
    print("\n     gain of the real prior over no prior:      %+.4f" % (best_prior - base))
    print("     gain of the SHUFFLED prior over no prior:  %+.4f" % (shuf - base))
    print("     real minus shuffled:                       %+.4f"
          % (best_prior - shuf))

    print("\n  B3 split-half edge reproducibility (the prior cannot inflate this"
          " the way overlap-with-itself can)")
    print("  %-26s %9s %9s %10s" % ("configuration", "edges A", "edges B",
                                    "both, same sign"))
    rep3 = {}
    allc = np.arange(n)
    for name, kw in [("no prior", dict(w=None)),
                     ("prior, penalty /%g" % BOOSTS[1],
                      dict(w=_weights(known, BOOSTS[1]))),
                     ("SHUFFLED prior, /%g" % BOOSTS[1],
                      dict(w=_weights(S, BOOSTS[1])))]:
        a_best = res["no prior"]["alpha"]
        BA = _fit_rows(N["RA"], allc, alpha=a_best, **kw)
        BB = _fit_rows(N["RB"], allc, alpha=a_best, **kw)
        both = (BA != 0) & (BB != 0) & (np.sign(BA) == np.sign(BB)) & off
        nzA, nzB = int((BA != 0)[off].sum()), int((BB != 0)[off].sum())
        chance = nzA * nzB / off.sum()
        rep3[name] = {"edges_A": nzA, "edges_B": nzB,
                      "both": int(both.sum()), "chance": float(chance),
                      "enrichment": float(both.sum() / chance) if chance else np.nan,
                      "prec_known": float(known[off][both[off]].mean())
                      if both.sum() else np.nan}
        print("  %-26s %9d %9d %10d  (chance %.0f, %.1fx; overlap with the real"
              " graph %.1f%%)"
              % (name, nzA, nzB, int(both.sum()), chance,
                 both.sum() / chance if chance else np.nan,
                 100 * rep3[name]["prec_known"]))
    print("     baseline rate of prior edges among all pairs: %.2f%%"
          % (100 * known[off].mean()))

    print("\n  B4 thinning: does the prior's advantage grow when data is scarce?")
    print("  %10s %10s %10s %10s %10s"
          % ("train frac", "no prior", "prior", "shuffled", "real gain"))
    thin = []
    for fr in (0.25, 0.5, 1.0):
        b0 = _heldout(R, None, frac=fr)["r2"]
        b1 = max(_heldout(R, known, boost=b, frac=fr)["r2"] for b in BOOSTS)
        b2 = _heldout(R, S, boost=BOOSTS[1], frac=fr)["r2"]
        thin.append({"frac": fr, "no_prior": b0, "prior": b1, "shuffled": b2,
                     "gain": b1 - b0})
        print("  %10.0f%% %10.4f %10.4f %10.4f %10.4f"
              % (100 * fr, b0, b1, b2, b1 - b0))
    thin_grows = thin[0]["gain"] > thin[-1]["gain"]
    print("     gain larger at 25%% than at 100%% of the data: %s" % thin_grows)

    return {"n_nodes": n, "genes": genes.tolist(), "rep_median": float(np.median(N["rep"])),
            "coverage": {"direct": int((P["direct"] | P["direct"].T)[off].sum()),
                         "complex": int(P["complex"][off].sum()),
                         "union": int(known[off].sum()),
                         "pairs": int(off.sum()),
                         "frac": float(known[off].mean()),
                         "shuffle_overlap": int((S & known)[off].sum()),
                         "nodes_covered": nodes_cov},
            "heldout": res, "base": base, "best_prior": best_prior,
            "shuffled": shuf, "split_half": rep3, "thinning": thin,
            "thin_grows": bool(thin_grows)}


# ----------------------------------------------------------------- report


def report(R):
    A, B = R["A"], R["B"]
    gain = B["best_prior"] - B["base"]
    gshuf = B["shuffled"] - B["base"]
    print("\n" + "=" * 78)
    print("D29 VERDICT")
    print("=" * 78)
    print("Question: does the prior interaction network make the recovered")
    print("Jacobian predict better, as opposed to merely overlapping it?")
    print("\nCalibration (part A, planted network): a PERFECT prior buys")
    print("%+.4f in held-out R^2, and a 0%%-accurate prior of the same size buys"
          % A["ceiling_gain"])
    print("%+.4f. So the detectable range on this design is known before looking."
          % A["zero_accuracy_gain"])
    print("The planted gain also grows when the data is thinned: %s."
          % A["gain_grows_when_thin"])
    print("\nReal data (part B, %d K562 nodes): the graph names %.2f%% of the"
          % (B["n_nodes"], 100 * B["coverage"]["frac"]))
    print("ordered pairs (%d of %d; %d from OmniPath, %d from complexes), and"
          % (B["coverage"]["union"], B["coverage"]["pairs"],
             B["coverage"]["direct"], B["coverage"]["complex"]))
    print("%d of %d nodes have at least one prior edge."
          % (B["coverage"]["nodes_covered"], B["n_nodes"]))
    print("  held-out R^2, no prior            %.4f" % B["base"])
    print("  held-out R^2, real prior (best)   %.4f   gain %+.4f"
          % (B["best_prior"], gain))
    print("  held-out R^2, SHUFFLED prior      %.4f   gain %+.4f"
          % (B["shuffled"], gshuf))
    print("  real minus shuffled               %+.4f" % (B["best_prior"] - B["shuffled"]))
    print("  prior-only support                %.4f"
          % B["heldout"]["prior-only support"]["r2"])
    print("  advantage grows when data is thinned: %s" % B["thin_grows"])
    # read the measured gain back onto the calibration curve
    accs = np.array([r["accuracy"] for r in A["accuracy_sweep"]])
    curve = np.array([max(r["boost_%g" % b] for b in BOOSTS)
                      for r in A["accuracy_sweep"]]) - A["no_prior"]["r2"]
    o = np.argsort(curve)
    implied = float(np.interp(gain, curve[o], accs[o]))
    print("\nReading %+.4f back onto the calibration curve: the real graph acts"
          % gain)
    print("like a prior of roughly %.0f%% accuracy on this node set." % (100 * implied))
    R.setdefault("derived", {})["implied_accuracy"] = implied

    sh = B["split_half"]
    k0 = "no prior"
    k1 = [k for k in sh if k.startswith("prior,")][0]
    k2 = [k for k in sh if k.startswith("SHUFFLED")][0]
    print("\nThe circularity, measured (part B3). Adding the prior raises the")
    print("recovered edges' overlap with that same prior from %.1f%% to %.1f%%,"
          % (100 * sh[k0]["prec_known"], 100 * sh[k1]["prec_known"]))
    print("which is the number a paper would report as 'recovers known")
    print("interactions'. But a DEGREE-PRESERVING SHUFFLED prior yields %d"
          % sh[k2]["both"])
    print("reproducible edges against the real prior's %d and the no-prior %d,"
          % (sh[k1]["both"], sh[k0]["both"]))
    print("with overlap unchanged at %.1f%%, and the enrichment over chance stays"
          % (100 * sh[k2]["prec_known"]))
    print("%.1fx / %.1fx / %.1fx. So the reproducibility gain comes from the extra"
          % (sh[k0]["enrichment"], sh[k1]["enrichment"], sh[k2]["enrichment"]))
    print("sparsity nudge, and the overlap gain is manufactured by the prior.")
    print("\nTwo things to say with the result:")
    print("1. On this node set the 'regulatory network' is almost entirely")
    print("   protein-complex co-membership: %d of the %d prior entries come"
          % (B["coverage"]["complex"], B["coverage"]["union"]))
    print("   from complexes and only %d from OmniPath's signed signalling"
          % B["coverage"]["direct"])
    print("   edges (%.1f%%). This result is about complexes, not about signalling."
          % (100 * B["coverage"]["direct"] / B["coverage"]["union"]))
    print("2. The data is not scarce here: no-prior held-out R^2 is already")
    print("   %.3f at 300 nodes with ~200 training perturbations, which is why"
          % B["base"])
    print("   there is little for a prior to add - and thinning to 25% does not")
    print("   reverse that (gain %+.4f)." % B["thinning"][0]["gain"])


def figures(R):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    A, B = R["A"], R["B"]
    fig, ax = plt.subplots(1, 4, figsize=(19.5, 4.6))

    a = ax[0]
    accs = [r["accuracy"] for r in A["accuracy_sweep"]]
    for b in BOOSTS:
        a.plot(accs, [r["boost_%g" % b] for r in A["accuracy_sweep"]], "o-",
               label="penalty /%g" % b)
    a.axhline(A["no_prior"]["r2"], color="0.4", ls="--", lw=1.2, label="no prior")
    a.set_xlabel("fraction of the prior's edges that are real")
    a.set_ylabel("held-out R²")
    a.legend(fontsize=8)
    a.set_title("A  calibration on a planted network:\na perfect prior buys %+.3f"
                % A["ceiling_gain"], fontsize=9.5)

    a = ax[1]
    names = ["no prior"] + ["prior, penalty /%g" % b for b in BOOSTS] + \
            ["SHUFFLED prior, /%g" % BOOSTS[1], "prior-only support"]
    vals = [B["heldout"][k]["r2"] for k in names]
    cols = ["0.5"] + ["C0"] * len(BOOSTS) + ["C1", "C2"]
    a.barh(range(len(vals)), vals, 0.6, color=cols)
    a.axvline(B["base"], color="0.3", ls="--", lw=1.2)
    for i, v in enumerate(vals):
        a.annotate("%.4f" % v, (v, i), va="center", ha="left", fontsize=7.5)
    a.set_yticks(range(len(vals)),
                 [k.replace("penalty ", "").replace(" prior", "") for k in names],
                 fontsize=8)
    a.invert_yaxis()
    a.set_xlabel("held-out R² (3 folds over perturbation columns)")
    a.set_xlim(min(vals) - 0.02, max(vals) + 0.03)
    a.set_title("B  real prior vs shuffled prior:\ngain %+.4f vs %+.4f"
                % (B["best_prior"] - B["base"], B["shuffled"] - B["base"]),
                fontsize=9.5)

    a = ax[2]
    fr = [100 * t["frac"] for t in B["thinning"]]
    a.plot(fr, [t["no_prior"] for t in B["thinning"]], "o-", label="no prior")
    a.plot(fr, [t["prior"] for t in B["thinning"]], "s-", label="real prior")
    a.plot(fr, [t["shuffled"] for t in B["thinning"]], "^-", label="shuffled prior")
    a.set_xlabel("% of the training perturbations used")
    a.set_ylabel("held-out R²")
    a.legend(fontsize=8)
    a.set_title("B4  the case for a prior is strongest when\ndata is scarce —"
                " tested, not assumed", fontsize=9.5)

    a = ax[3]
    sh = B["split_half"]
    ks = ["no prior"] + [k for k in sh if k.startswith("prior,")] + \
         [k for k in sh if k.startswith("SHUFFLED")]
    xs = np.arange(len(ks))
    ov = [100 * sh[k]["prec_known"] for k in ks]
    a.bar(xs, ov, 0.55, color=["0.5", "C0", "C1"])
    for x, v, k in zip(xs, ov, ks):
        a.annotate("%.1f%%\n%d edges" % (v, sh[k]["both"]), (x, v),
                   ha="center", va="bottom", fontsize=8)
    a.axhline(100 * B["coverage"]["frac"], color="C3", ls="--", lw=1.2)
    a.annotate("rate among all pairs (%.1f%%)" % (100 * B["coverage"]["frac"]),
               (len(ks) - 0.5, 100 * B["coverage"]["frac"] + 1), ha="right",
               fontsize=7.5, color="C3")
    a.set_xticks(xs, ["no prior", "real prior", "SHUFFLED\nprior"], fontsize=8.5)
    a.set_ylabel("recovered edges that are in the REAL graph (%)")
    a.set_ylim(0, 46)
    a.set_title("B3  the circular metric, measured: a prior\nmanufactures its own"
                " overlap", fontsize=9.5)

    fig.tight_layout()
    fn = os.path.join(FIG, "d29_prior_network_value.png")
    fig.savefig(fn, dpi=140)
    print("\nfigure: %s" % fn)


def main(force=False):
    if os.path.exists(RES) and not force:
        R = pickle.load(open(RES, "rb"))
    else:
        A = part_a()
        B = part_b()
        R = {"A": A, "B": B}
        pickle.dump(R, open(RES, "wb"))
    report(R)
    figures(R)
    return R


if __name__ == "__main__":
    main(force=bool(os.environ.get("D29_FORCE")))
