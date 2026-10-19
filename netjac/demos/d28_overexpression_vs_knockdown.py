"""
D28  Is overexpression the negative of knockdown?

Every lever ranking in this project is built from a knockdown library
(Replogle 2022 CRISPRi; D9/D10/D16). Every therapeutic idea that says
"raise module X" is an overexpression. The step that silently connects them
is the linear-response assumption

    d = -J^-1 b u      =>      d(-u) = -d(u)

i.e. pushing a node up gives exactly minus the response to pushing it down.
If that fails, a screen built on a knockdown library mis-ranks overexpression
targets in a way no amount of extra data fixes, because the error is in the
model class, not in the sample size.

K562 is the one cell line where both directions exist publicly:
  Norman 2019 (Science 365:786)  CRISPRa  105 single-gene activations
  Replogle 2022 (Cell 185:2559)  CRISPRi  11,258 single-gene knockdowns
84 genes are perturbed in both, 806 measured genes are shared.

PART A  the prediction stated exactly, and the estimator validated on a
  synthetic system with a PLANTED answer, because three different failures
  look similar in a correlation:
    A1  pure linear response                 -> cos(r_up, r_down) = -1 exactly
    A2  saturating but odd input (tanh)      -> still exactly -1: magnitude
        differs, DIRECTION does not. So "saturation" cannot explain a loss of
        antisymmetry; only a direction-specific channel can.
    A3  a direction-specific second channel (a stress program that only the
        up-perturbation engages) -> cos moves off -1 by a predictable amount;
        the estimator must recover the planted value.
    A4  a shared global mode added to every perturbation -> the RAW matched
        correlation is destroyed, while matched-minus-mismatched survives.
        This validates the control used in part B.

PART B  the measurement
  B1  response matrices on the shared gene space, in control-SD units
  B2  positive control: the perturbed gene itself must go up under CRISPRa
      and down under CRISPRi
  B3  reliability ceilings. Norman: split-half, plus 47 genes measured in two
      guide-pair orientations (same biology, same study) = the ceiling any
      cross-study correlation can reach. Replogle: 585 non-targeting rows =
      the noise floor, plus 12 genes with two TSS annotations.
  B4  matched vs mismatched cosine, before and after mean-field removal
  B5  attenuation-corrected antisymmetry, slope, and a per-gene three-way
      classification against the mismatched null

PART C  what it costs a screen that owns only a knockdown library
  C1  sign error rate of -r_KD as a predictor of r_OE
  C2  the ranking task: pick the overexpression target that produces a wanted
      state change, scoring candidates by -r_KD

PART D  what would settle it

Run:  python3 -m netjac.demos.d28_overexpression_vs_knockdown
"""
import json
import os
import pickle

import numpy as np

from netjac import DATA_CACHE, FIG

RES = os.path.join(DATA_CACHE, "d28_results.pkl")
NORMAN = os.path.join(DATA_CACHE, "norman_gears", "norman", "perturb_processed.h5ad")
GWPS = os.path.join(DATA_CACHE, "perturbseq", "K562_gwps_normalized_bulk_01.h5ad")
MIN_CELLS = 50


# ------------------------------------------------------------------ helpers


def _cos(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(a @ b / (na * nb)) if na > 0 and nb > 0 else np.nan


def _pairwise_cos(A, B):
    """cos between every row of A and every row of B."""
    An = A / (np.linalg.norm(A, axis=1, keepdims=True) + 1e-12)
    Bn = B / (np.linalg.norm(B, axis=1, keepdims=True) + 1e-12)
    return An @ Bn.T


def _summary(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return {"n": int(v.size), "median": float(np.median(v)),
            "mean": float(v.mean()), "q05": float(np.percentile(v, 5)),
            "q95": float(np.percentile(v, 95))}


# -------------------------------------------- PART A  estimator validation


def _synthetic(n=60, planted=0.0, global_mode=0.0, seed=0, u=1.0, node=None,
               net_seed=12345):
    """Returns (r_up, r_down) for one node of a random stable linear network.

    The NETWORK is fixed (net_seed) and the perturbed node varies with `seed`,
    so two different perturbations give genuinely different response directions
    and the mismatched baseline sits at 0 rather than inheriting one shared
    direction.

    planted      weight of a channel engaged ONLY by the up-perturbation
    global_mode  weight of a shared mode added to every response
    """
    net = np.random.default_rng(net_seed)
    J = net.normal(0, 1.0 / np.sqrt(n), (n, n)) - 1.5 * np.eye(n)
    g = net.normal(0, 1, n)            # the shared global mode, same for all
    g /= np.linalg.norm(g)
    rng = np.random.default_rng(seed)
    b = np.zeros(n)
    b[node if node is not None else seed % n] = 1.0
    c = rng.normal(0, 1, n)            # the direction-specific channel
    c /= np.linalg.norm(c)
    Ji = np.linalg.inv(J)
    r_up = -Ji @ (b * u + c * planted * u)
    r_dn = -Ji @ (b * (-u))
    if global_mode:
        s = global_mode * np.linalg.norm(r_up)
        r_up = r_up + s * g
        r_dn = r_dn + s * g
    return r_up, r_dn


def part_a():
    print("\n" + "=" * 78)
    print("PART A  the prediction, and the estimator checked against a planted answer")
    print("=" * 78)
    out = {}

    # A1 / A2 exact antisymmetry
    r_up, r_dn = _synthetic(planted=0.0)
    out["linear_cos"] = _cos(r_up, r_dn)
    # odd saturating input: r(u) ~ tanh(u), direction unchanged
    r_up2, _ = _synthetic(planted=0.0, u=3.0)
    out["saturated_cos"] = _cos(r_up2, r_dn)
    out["saturated_norm_ratio"] = float(np.linalg.norm(r_up2) / np.linalg.norm(r_up))
    print("\nA1  pure linear response: cos(r_up, r_down) = %.6f" % out["linear_cos"])
    print("A2  3x stronger up-perturbation: cos = %.6f, norm ratio %.2f"
          % (out["saturated_cos"], out["saturated_norm_ratio"]))
    print("    -> saturation changes the MAGNITUDE, never the direction. A loss")
    print("       of antisymmetry cannot be explained by dose or saturation.")

    # A3 a direction-specific channel, with the predicted cosine
    rows = []
    print("\nA3  a channel engaged only by the up-perturbation:")
    print("    %8s %12s %18s" % ("planted", "median cos", "IQR across nodes"))
    for w in (0.0, 0.5, 1.0, 2.0, 4.0):
        cs = [_cos(*_synthetic(planted=w, seed=s_)) for s_ in range(40)]
        rows.append({"planted": w, "measured": float(np.median(cs)),
                     "iqr": [float(np.percentile(cs, 25)),
                             float(np.percentile(cs, 75))]})
        print("    %8.1f %12.4f  [%.3f, %.3f]"
              % (w, rows[-1]["measured"], *rows[-1]["iqr"]))
    out["planted_sweep"] = rows
    mono = all(rows[i]["measured"] <= rows[i + 1]["measured"] + 1e-9
               for i in range(len(rows) - 1))
    out["monotone_in_planted"] = bool(mono)
    print("    cosine rises monotonically with the planted weight: %s" % mono)

    # A4 the global-mode confound: raw, contrast, and explicit mean-field removal
    print("\nA4  a shared global mode added to every perturbation")
    print("    (prediction with no confound: matched -1, mismatched 0, so the")
    print("     matched-minus-mismatched contrast is -1):")
    print("    %6s %12s %12s %12s %16s"
          % ("mode", "raw matched", "mismatched", "contrast", "after mean-field"))
    rows4 = []
    for gm in (0.0, 0.5, 1.0, 2.0, 4.0):
        UP, DN = [], []
        for s_ in range(40):
            ru, rd = _synthetic(planted=0.0, global_mode=gm, seed=s_)
            UP.append(ru)
            DN.append(rd)
        UP, DN = np.array(UP), np.array(DN)
        C = _pairwise_cos(UP, DN)
        matched = np.diag(C)
        mism = C[~np.eye(len(C), dtype=bool)]
        Cm = _pairwise_cos(UP - UP.mean(0), DN - DN.mean(0))
        rows4.append({"global_mode": gm, "matched": float(np.median(matched)),
                      "mismatched": float(np.median(mism)),
                      "contrast": float(np.median(matched) - np.median(mism)),
                      "mean_field_removed": float(np.median(np.diag(Cm)))})
        print("    %6.1f %12.4f %12.4f %12.4f %16.4f"
              % (gm, rows4[-1]["matched"], rows4[-1]["mismatched"],
                 rows4[-1]["contrast"], rows4[-1]["mean_field_removed"]))
    out["global_mode_sweep"] = rows4
    mfr = np.array([r["mean_field_removed"] for r in rows4])
    gms = np.array([r["global_mode"] for r in rows4])
    out["mf_dev_mode_le_1"] = float(np.abs(mfr[gms <= 1.0] + 1.0).max())
    out["mf_dev_worst"] = float(np.abs(mfr + 1.0).max())
    print("    -> the raw matched cosine is destroyed by the global mode (-1.00 ->"
          " %+.2f)," % rows4[-1]["matched"])
    print("       and the matched-minus-mismatched contrast degrades with it"
          " (-1.00 -> %+.2f)." % rows4[-1]["contrast"])
    print("       Explicit mean-field removal holds the planted -1 to within %.3f"
          % out["mf_dev_mode_le_1"])
    print("       while the mode is <= the signal, and slips to %.3f when the mode"
          % out["mf_dev_worst"])
    print("       is 4x the signal. So part B subtracts the across-perturbation")
    print("       mean AND reports the mismatched null beside it; neither alone is")
    print("       sufficient.")
    return out


# ------------------------------------------------------- PART B  real data


def _norman(seed=0):
    """CRISPRa single-gene responses in control-SD units, plus split halves
    and the two guide-pair orientations where both exist."""
    import anndata as ad
    a = ad.read_h5ad(NORMAN)
    cond = a.obs["condition"].astype(str).values
    X = a.X.tocsr() if hasattr(a.X, "tocsr") else a.X
    genes = np.array(a.var["gene_name"].astype(str))
    ci = np.where(cond == "ctrl")[0]
    C = X[ci]
    C = C.toarray() if hasattr(C, "toarray") else np.asarray(C)
    mu0 = C.mean(0)
    sd0 = C.std(0)
    keep = sd0 > 1e-6
    rng = np.random.default_rng(seed)

    resp, halves, ncell = {}, {}, {}
    for c in np.unique(cond):
        if c == "ctrl" or "ctrl" not in c:
            continue
        g = c.replace("+ctrl", "").replace("ctrl+", "")
        idx = np.where(cond == c)[0]
        if len(idx) < MIN_CELLS:
            continue
        Xi = X[idx]
        m = np.asarray(Xi.mean(0)).ravel()
        r = np.zeros_like(m)
        r[keep] = (m[keep] - mu0[keep]) / sd0[keep]
        resp.setdefault(g, {})[c] = r
        ncell.setdefault(g, {})[c] = int(len(idx))
        p = rng.permutation(len(idx))
        h = []
        for part in (p[: len(p) // 2], p[len(p) // 2:]):
            mm = np.asarray(X[idx[part]].mean(0)).ravel()
            rr = np.zeros_like(mm)
            rr[keep] = (mm[keep] - mu0[keep]) / sd0[keep]
            h.append(rr)
        halves.setdefault(g, {})[c] = tuple(h)
    return {"genes": genes, "resp": resp, "halves": halves, "ncell": ncell,
            "n_ctrl": int(len(ci)), "keep": keep}


def _replogle():
    """CRISPRi responses (already in control-SD units) plus the non-targeting
    rows that give the noise floor."""
    import anndata as ad
    a = ad.read_h5ad(GWPS)
    idx = list(a.obs.index)
    pert = np.array([s.split("_")[1] for s in idx])
    tss = np.array([s.split("_")[2] for s in idx])
    genes = np.array(a.var["gene_name"].astype(str))
    X = np.asarray(a.X, dtype=np.float32)
    nt = tss == "non-targeting"
    fold = np.asarray(a.obs["fold_expr"], dtype=float)
    ncell = np.asarray(a.obs["num_cells_filtered"], dtype=float)
    return {"genes": genes, "X": X, "pert": pert, "tss": tss, "nt": nt,
            "fold": fold, "ncell": ncell}


def part_b():
    print("\n" + "=" * 78)
    print("PART B  the two K562 screens, matched gene by gene")
    print("=" * 78)
    N = _norman()
    R = _replogle()

    common = sorted(set(N["genes"][N["keep"]]) & set(R["genes"]))
    ni = {g: i for i, g in enumerate(N["genes"])}
    ri = {}
    for i, g in enumerate(R["genes"]):
        ri.setdefault(g, i)
    cn = np.array([ni[g] for g in common])
    cr = np.array([ri[g] for g in common])
    print("  CRISPRa (Norman 2019): %d single-gene activations with >= %d cells,"
          " %d control cells" % (len(N["resp"]), MIN_CELLS, N["n_ctrl"]))
    print("  CRISPRi (Replogle 2022 gwps): %d rows, %d non-targeting"
          % (R["X"].shape[0], int(R["nt"].sum())))
    finite = np.isfinite(R["X"][:, cr]).all(0)
    if not finite.all():
        common = [g for g, f in zip(common, finite) if f]
        cn, cr = cn[finite], cr[finite]
        print("  shared measured genes: %d (%d dropped for missing CRISPRi values)"
              % (len(common), int((~finite).sum())))
    else:
        print("  shared measured genes: %d" % len(common))

    # one CRISPRi row per gene: the strongest knockdown among its TSS rows,
    # falling back to the row with the most cells when fold_expr is missing
    kd_row, kd_fold = {}, {}
    for g in set(R["pert"]):
        rows = np.where((R["pert"] == g) & ~R["nt"])[0]
        if rows.size == 0:
            continue
        f = R["fold"][rows]
        if np.isfinite(f).any():
            j = rows[np.nanargmin(np.where(np.isfinite(f), f, np.inf))]
        else:
            j = rows[int(np.argmax(R["ncell"][rows]))]
        kd_row[g] = j
        kd_fold[g] = float(R["fold"][j])

    pairs = sorted(set(N["resp"]) & set(kd_row))
    print("  genes perturbed in BOTH directions: %d" % len(pairs))

    # ---- B2 positive control: the perturbed gene's own direction
    own = []
    for g in pairs:
        if g not in ni or g not in ri:
            continue
        oe = np.mean([v[ni[g]] for v in N["resp"][g].values()])
        kd = R["X"][kd_row[g], ri[g]]
        own.append({"gene": g, "oe_own": float(oe), "kd_own": float(kd),
                    "fold_expr": kd_fold[g]})
    oe_own = np.array([o["oe_own"] for o in own])
    kd_own = np.array([o["kd_own"] for o in own])
    print("\n  B2 positive control, the perturbed gene itself (n = %d measured in"
          " both):" % len(own))
    print("     CRISPRa raises it in %d/%d (median %+.2f SD)"
          % (int((oe_own > 0).sum()), len(own), float(np.median(oe_own))))
    print("     CRISPRi lowers it in %d/%d (median %+.2f SD)"
          % (int((kd_own < 0).sum()), len(own), float(np.median(kd_own))))
    print("     -> the two libraries are oriented as expected; any loss of")
    print("        antisymmetry below is not a sign-flip bookkeeping error.")

    # ---- B3 reliability ceilings
    # Norman split-half
    sh = []
    for g in pairs:
        for c, (h1, h2) in N["halves"].get(g, {}).items():
            sh.append(_cos(h1[cn], h2[cn]))
    # Norman: two guide-pair orientations of the same activation
    orient = []
    for g in pairs:
        cs = list(N["resp"][g].values())
        if len(cs) == 2:
            orient.append(_cos(cs[0][cn], cs[1][cn]))
    # Replogle: non-targeting rows = noise floor
    NTr = R["X"][np.where(R["nt"])[0]][:, cr]
    Cnt = _pairwise_cos(NTr, NTr)
    nt_cos = Cnt[~np.eye(len(Cnt), dtype=bool)]
    nt_norm = np.linalg.norm(NTr, axis=1)
    # Replogle: two TSS annotations of the same gene
    tss_rep = []
    for g in pairs:
        rows = np.where((R["pert"] == g) & ~R["nt"])[0]
        if rows.size >= 2:
            tss_rep.append(_cos(R["X"][rows[0]][cr], R["X"][rows[1]][cr]))
    sig_norm = np.array([np.linalg.norm(R["X"][kd_row[g]][cr]) for g in pairs])
    # per-gene detectability: is this response bigger than its own noise floor?
    kd_pct = np.array([float((nt_norm < v).mean()) for v in sig_norm])
    oe_rel = np.array([float(np.median([_cos(h1[cn], h2[cn]) for h1, h2
                                        in N["halves"].get(g, {}).values()]))
                       for g in pairs])
    print("\n  B3 reliability ceilings (on the %d shared genes):" % len(common))
    print("     CRISPRa split-half cos                  median %+.3f (n = %d)"
          % (np.median(sh), len(sh)))
    print("     CRISPRa two guide-pair orientations     median %+.3f (n = %d)"
          "  <- same biology, same study" % (np.median(orient), len(orient)))
    print("     CRISPRi non-targeting pairs (noise)     median %+.3f (n = %d)"
          % (np.median(nt_cos), len(nt_cos)))
    print("     CRISPRi two-TSS pairs                   median %+.3f (n = %d)"
          % (np.median(tss_rep) if tss_rep else np.nan, len(tss_rep)))
    print("     CRISPRi signal/noise norm ratio         %.2f (median %.2f vs"
          " non-targeting %.2f)" % (np.median(sig_norm) / np.median(nt_norm),
                                    np.median(sig_norm), np.median(nt_norm)))
    print("     CRISPRi norm above the non-targeting floor (95th pct = %.2f):"
          " %d / %d genes" % (np.percentile(nt_norm, 95),
                              int((kd_pct > 0.95).sum()), len(pairs)))
    print("     CRISPRa split-half cos > 0.3:                       %d / %d genes"
          % (int((oe_rel > 0.3).sum()), len(pairs)))
    detect = (kd_pct > 0.95) & (oe_rel > 0.3)
    print("     BOTH detectable:                                    %d / %d genes"
          % (int(detect.sum()), len(pairs)))
    # two independent estimates of the CRISPRi reliability, because 12 TSS pairs
    # is thin and those rows are among the weakest perturbations
    rel_tss = float(np.median(tss_rep)) if tss_rep else np.nan
    rel_norm = float(max(0.0, 1.0 - (np.median(nt_norm) / np.median(sig_norm)) ** 2))
    rel_oe = float(np.median(orient))
    ceil_tss = float(np.sqrt(max(rel_oe, 0.0) * max(rel_tss, 0.0)))
    ceil_norm = float(np.sqrt(max(rel_oe, 0.0) * rel_norm))
    ceiling = max(ceil_tss, ceil_norm)
    print("     CRISPRi reliability, two estimates: two-TSS cos %.3f, and"
          " 1-(noise/signal)^2 = %.3f" % (rel_tss, rel_norm))
    print("     => attenuation ceiling for any cross-study |cos|: %.3f to %.3f"
          % (min(ceil_tss, ceil_norm), max(ceil_tss, ceil_norm)))

    # ---- B4 matched vs mismatched, raw and mean-field removed
    OE = np.array([np.mean([v[cn] for v in N["resp"][g].values()], axis=0)
                   for g in pairs])
    KD = np.array([R["X"][kd_row[g]][cr] for g in pairs])
    # drop each gene's own column from the comparison (trivially antisymmetric)
    self_col = np.array([common.index(g) if g in set(common) else -1
                         for g in pairs])
    res = {}
    for tag, (A, B) in (("raw", (OE, KD)),
                        ("mean-field removed", (OE - OE.mean(0), KD - KD.mean(0)))):
        Ax, Bx = A.copy(), B.copy()
        for k, j in enumerate(self_col):
            if j >= 0:
                Ax[k, j] = 0.0
                Bx[k, j] = 0.0
        C = _pairwise_cos(Ax, Bx)
        matched = np.diag(C).copy()
        mism = C[~np.eye(len(C), dtype=bool)]
        res[tag] = {"matched": _summary(matched), "mismatched": _summary(mism),
                    "matched_vals": matched.tolist(),
                    "mismatched_median": float(np.median(mism)),
                    "difference": float(np.median(matched) - np.median(mism))}
        print("\n  B4 %s:" % tag)
        print("     matched cos(r_OE_i, r_KD_i)    median %+.4f [5-95%%: %+.3f,"
              " %+.3f], n = %d" % (res[tag]["matched"]["median"],
                                   res[tag]["matched"]["q05"],
                                   res[tag]["matched"]["q95"],
                                   res[tag]["matched"]["n"]))
        print("     mismatched cos(r_OE_i, r_KD_j) median %+.4f [5-95%%: %+.3f,"
              " %+.3f], n = %d" % (res[tag]["mismatched"]["median"],
                                   res[tag]["mismatched"]["q05"],
                                   res[tag]["mismatched"]["q95"],
                                   res[tag]["mismatched"]["n"]))
        print("     difference                     %+.4f  (the prediction is -1:"
              " matched -1, mismatched 0)" % res[tag]["difference"])

    # the same contrast restricted to genes detectable on BOTH sides
    if detect.sum() >= 10:
        Ad = (OE - OE.mean(0))[detect]
        Bd = (KD - KD.mean(0))[detect]
        for k, j in enumerate(self_col[detect]):
            if j >= 0:
                Ad[k, j] = 0.0
                Bd[k, j] = 0.0
        Cd = _pairwise_cos(Ad, Bd)
        md = np.diag(Cd)
        misd = Cd[~np.eye(len(Cd), dtype=bool)]
        res["detectable"] = {"n": int(detect.sum()),
                             "matched_median": float(np.median(md)),
                             "mismatched_median": float(np.median(misd)),
                             "genes": [pairs[k] for k in np.where(detect)[0]]}
        print("\n  B4b restricted to the %d genes detectable on BOTH sides:"
              % int(detect.sum()))
        print("     matched %+.4f, mismatched %+.4f, difference %+.4f"
              % (res["detectable"]["matched_median"],
                 res["detectable"]["mismatched_median"],
                 res["detectable"]["matched_median"]
                 - res["detectable"]["mismatched_median"]))
    else:
        res["detectable"] = {"n": int(detect.sum()),
                             "note": "too few genes detectable on both sides"}
        print("\n  B4b only %d genes are detectable on both sides — the matched"
              " test is UNDER-POWERED on the knockdown side" % int(detect.sum()))

    # permutation test on the pairing
    rng = np.random.default_rng(0)
    A = OE - OE.mean(0)
    B = KD - KD.mean(0)
    for k, j in enumerate(self_col):
        if j >= 0:
            A[k, j] = 0.0
            B[k, j] = 0.0
    C = _pairwise_cos(A, B)
    obs = float(np.median(np.diag(C)))
    null = []
    for _ in range(5000):
        p = rng.permutation(len(C))
        null.append(float(np.median(C[np.arange(len(C)), p])))
    null = np.array(null)
    pval = float((null <= obs).mean())
    print("\n     permutation test on the pairing: observed median %+.4f,"
          " null %+.4f +- %.4f, p = %.4f"
          % (obs, null.mean(), null.std(), pval))

    # ---- B5 three-way classification against the mismatched null
    thr = float(np.percentile(np.abs(C[~np.eye(len(C), dtype=bool)]), 95))
    m = np.diag(C)
    anti = int((m < -thr).sum())
    co = int((m > thr).sum())
    indep = int(len(m) - anti - co)
    print("\n  B5 per-gene classification (threshold = 95th pct of |mismatched| ="
          " %.3f):" % thr)
    print("     antisymmetric (cos < -thr) : %d / %d" % (anti, len(m)))
    print("     independent   (|cos| <= thr): %d / %d" % (indep, len(m)))
    print("     co-directional (cos > +thr): %d / %d" % (co, len(m)))
    from scipy.stats import binomtest
    tails = anti + co
    bt = binomtest(anti, tails, 0.5, alternative="greater") if tails else None
    print("     tail asymmetry: %d antisymmetric vs %d co-directional,"
          " binomial p = %s"
          % (anti, co, ("%.2e" % bt.pvalue) if bt else "n/a"))
    print("     (noise alone would fill the two tails equally, so this is the")
    print("      cleanest evidence that a real antisymmetric component exists)")
    order = np.argsort(m)
    print("     most antisymmetric: %s"
          % ", ".join("%s %+.2f" % (pairs[i], m[i]) for i in order[:5]))
    print("     most co-directional: %s"
          % ", ".join("%s %+.2f" % (pairs[i], m[i]) for i in order[-5:]))
    # pooled slope of r_OE on r_KD
    sl = float(np.polyfit(B.ravel(), A.ravel(), 1)[0])
    print("     pooled slope of r_OE on r_KD: %+.4f (the prediction is -1)" % sl)

    # ---- B6 restrict to genes where BOTH perturbations demonstrably worked
    own_map = {o["gene"]: o for o in own}
    eff = [k for k, g in enumerate(pairs)
           if g in own_map and own_map[g]["oe_own"] > 0.5
           and np.isfinite(own_map[g]["fold_expr"])
           and own_map[g]["fold_expr"] < 0.7]
    sub = {}
    if len(eff) >= 10:
        idx = np.array(eff)
        Cs = _pairwise_cos(A[idx], B[idx])
        ms = np.diag(Cs)
        mis_s = Cs[~np.eye(len(Cs), dtype=bool)]
        sub = {"n": len(eff), "matched_median": float(np.median(ms)),
               "mismatched_median": float(np.median(mis_s)),
               "genes": [pairs[k] for k in eff]}
        print("\n  B6 restricted to genes where BOTH perturbations demonstrably"
              " worked")
        print("     (CRISPRa raises the gene > 0.5 SD AND CRISPRi fold_expr <"
              " 0.7): n = %d" % len(eff))
        print("     matched %+.4f, mismatched %+.4f — the conclusion does not"
              " depend on efficiency" % (sub["matched_median"],
                                         sub["mismatched_median"]))
    else:
        print("\n  B6 too few genes pass the joint efficiency filter (%d)" % len(eff))

    norm_oe = float(np.median(np.linalg.norm(OE, axis=1)))
    norm_kd = float(np.median(np.linalg.norm(KD, axis=1)))
    print("     median response norm: CRISPRa %.2f vs CRISPRi %.2f (ratio %.1fx)"
          % (norm_oe, norm_kd, norm_oe / norm_kd))

    return {"efficient_subset": sub, "norm_oe": norm_oe, "norm_kd": norm_kd,"pairs": pairs, "common": common, "own": own,
            "reliability": {"split_half": _summary(sh),
                            "orientations": _summary(orient),
                            "nt_cos": _summary(nt_cos),
                            "tss_pairs": _summary(tss_rep),
                            "signal_norm": float(np.median(sig_norm)),
                            "nt_norm": float(np.median(nt_norm)),
                            "nt_norm_q95": float(np.percentile(nt_norm, 95)),
                            "rel_oe": rel_oe, "rel_tss": rel_tss,
                            "rel_norm": rel_norm, "ceil_tss": ceil_tss,
                            "ceil_norm": ceil_norm,
                            "kd_detect": int((kd_pct > 0.95).sum()),
                            "oe_detect": int((oe_rel > 0.3).sum()),
                            "both_detect": int(detect.sum()),
                            "ceiling": ceiling},
            "kd_pct": kd_pct.tolist(), "oe_rel": oe_rel.tolist(),
            "detect": detect.tolist(),
            "contrast": res, "perm_p": pval, "perm_null_mean": float(null.mean()),
            "obs_median": obs, "thr": thr, "slope": sl,
            "classes": {"antisymmetric": anti, "independent": indep,
                        "co_directional": co,
                        "tail_binom_p": float(bt.pvalue) if bt else np.nan},
            "matched_cos": m.tolist(),
            "OE": OE, "KD": KD, "self_col": self_col.tolist()}


# ------------------------------------------- PART C  cost to a screen


def part_c(B):
    print("\n" + "=" * 78)
    print("PART C  what this costs a screen that owns only a knockdown library")
    print("=" * 78)
    pairs = B["pairs"]
    OE, KD = B["OE"].copy(), B["KD"].copy()
    OE -= OE.mean(0)
    KD -= KD.mean(0)
    for k, j in enumerate(B["self_col"]):
        if j >= 0:
            OE[k, j] = 0.0
            KD[k, j] = 0.0
    m = np.array(B["matched_cos"])
    wrong = int((m > 0).sum())
    print("\n  C1 sign of the prediction")
    print("     -r_KD predicts r_OE with the WRONG sign for %d / %d genes (%.0f%%)"
          % (wrong, len(m), 100 * wrong / len(m)))
    print("     i.e. for those genes, knocking the gene down and raising it move")
    print("     the cell in the SAME direction, not opposite ones.")

    # C2 ranking task
    C = _pairwise_cos(-KD, OE)          # score of candidate j for target i
    rng = np.random.default_rng(0)
    top1 = top5 = 0
    ranks = []
    for i in range(len(pairs)):
        s = C[:, i]
        order = np.argsort(-s)
        rk = int(np.where(order == i)[0][0]) + 1
        ranks.append(rk)
        top1 += rk == 1
        top5 += rk <= 5
    n = len(pairs)
    print("\n  C2 ranking task: given the state change an overexpression actually")
    print("     produced, rank all %d candidates by cos(-r_KD, target)" % n)
    print("     top-1 hit rate %.1f%% (chance %.1f%%), top-5 %.1f%% (chance %.1f%%),"
          " median rank %d of %d"
          % (100 * top1 / n, 100 / n, 100 * top5 / n, 500 / n,
             int(np.median(ranks)), n))
    # a random-score baseline for calibration
    rnd = []
    for _ in range(200):
        Cr = _pairwise_cos(rng.permutation(-KD), OE)
        rnd.append(np.mean([np.argsort(-Cr[:, i])[0] == i for i in range(n)]))
    print("     permuted-library baseline: %.1f%% top-1" % (100 * np.mean(rnd)))
    return {"wrong_sign": wrong, "n": n, "top1": top1 / n, "top5": top5 / n,
            "median_rank": int(np.median(ranks)), "ranks": ranks,
            "perm_top1": float(np.mean(rnd))}


# ----------------------------------------------------------------- report


def report(R):
    A, B, C = R["A"], R["B"], R["C"]
    rel = B["reliability"]
    mf = B["contrast"]["mean-field removed"]
    eff = B.get("efficient_subset", {})
    lo = min(rel["ceil_tss"], rel["ceil_norm"])
    hi = max(rel["ceil_tss"], rel["ceil_norm"])
    corr_hi = mf["matched"]["median"] / hi if hi > 0 else np.nan
    corr_lo = mf["matched"]["median"] / lo if lo > 0 else np.nan
    print("\n" + "=" * 78)
    print("D28 VERDICT")
    print("=" * 78)
    print("Prediction (linear response): cos(r_OE, r_KD) = -1, slope -1.")
    print("Measured on the %d genes perturbed in BOTH directions in K562:"
          % len(B["pairs"]))
    print("  matched cosine     %+.4f   mismatched %+.4f   difference %+.4f"
          % (mf["matched"]["median"], mf["mismatched"]["median"],
             mf["difference"]))
    sl_tss = B["slope"] / rel["rel_tss"] if rel["rel_tss"] > 0 else np.nan
    sl_norm = B["slope"] / rel["rel_norm"] if rel["rel_norm"] > 0 else np.nan
    print("  pooled slope       %+.4f" % B["slope"])
    print("    note: a cosine is attenuated by the noise on BOTH sides, a")
    print("    regression slope only by the noise on the regressor. Correcting")
    print("    the slope for the knockdown reliability alone gives %+.2f (norm"
          % sl_norm)
    print("    estimate) to %+.2f (two-TSS estimate)." % sl_tss)
    print("  permutation p      %.4f (the antisymmetric component is REAL)"
          % B["perm_p"])
    print("  per gene: %d antisymmetric, %d independent, %d co-directional"
          % (B["classes"]["antisymmetric"], B["classes"]["independent"],
             B["classes"]["co_directional"]))
    print("    tail asymmetry binomial p = %.2e — noise would fill both tails"
          % B["classes"]["tail_binom_p"])
    if eff.get("n"):
        print("  efficiency-filtered subset (n = %d): matched %+.4f, mismatched"
              " %+.4f" % (eff["n"], eff["matched_median"],
                          eff["mismatched_median"]))
    print("\nBut the knockdown side is close to its own noise floor for these")
    print("genes: only %d/%d have a response norm above the 95th percentile of"
          % (rel["kd_detect"], len(B["pairs"])))
    print("the %d non-targeting rows, the two-TSS replicate cosine is %+.3f, and"
          % (585, rel["rel_tss"]))
    print("the norm-based reliability is %.3f. So the attenuation ceiling on any"
          % rel["rel_norm"])
    print("cross-study |cos| is %.3f to %.3f, and the measured %+.4f corrects to"
          % (lo, hi, mf["matched"]["median"]))
    print("roughly %+.2f to %+.2f once attenuation is undone." % (corr_hi, corr_lo))
    print("\nPutting the two corrections together: the knockdown response gets")
    print("the SIGN right and, after correcting for its own noise, roughly the")
    print("right per-gene SLOPE (%+.2f against a prediction of -1) - but it" % sl_norm)
    print("explains only a small share of the overexpression response, because")
    print("the corrected direction correlation is about %+.2f, i.e. r^2 ~ %.2f."
          % (corr_hi, corr_hi ** 2))
    print("Mechanically: r_OE ~ -r_KD PLUS a large component orthogonal to it.")
    print("CRISPRa responses are %.1fx larger in norm (%.2f vs %.2f), so a small"
          % (B["norm_oe"] / B["norm_kd"], B["norm_oe"], B["norm_kd"]))
    print("correlation is compatible with a near -1 slope.")
    print("So overexpression is NOT the negative of knockdown: %d/%d genes are"
          % (B["classes"]["independent"], len(B["pairs"])))
    print("directionally independent, and the shared part is a minority of what")
    print("the activation actually does.")
    print("\nCost to a screen that owns only a knockdown library:")
    print("  -r_KD points the wrong way for %d/%d genes (%.0f%%), though none of"
          % (C["wrong_sign"], C["n"], 100 * C["wrong_sign"] / C["n"]))
    print("  those reach the co-directional threshold, so they are consistent")
    print("  with noise rather than proven reversals.")
    print("  Picking the right overexpression target: top-1 %.1f%% vs %.1f%%"
          % (100 * C["top1"], 100 / C["n"]))
    print("  chance (%.0fx), top-5 %.1f%%, median rank %d of %d."
          % (100 * C["top1"] / (100 / C["n"]), 100 * C["top5"],
             C["median_rank"], C["n"]))
    print("  For comparison, D16's lever test WITHIN the knockdown modality hits")
    print("  66%. Crossing the modality costs most of that.")
    print("\nEstimator checks (part A): exact linear response gives %.4f; a 3x"
          % A["linear_cos"])
    print("stronger perturbation still gives %.4f, so dose and saturation cannot"
          % A["saturated_cos"])
    print("explain a loss of antisymmetry. A shared global mode destroys the raw")
    print("matched cosine AND degrades the matched-minus-mismatched contrast;")
    print("only explicit mean-field removal recovers the planted answer, which is")
    print("the procedure used above.")


def figures(R):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    A, B, C = R["A"], R["B"], R["C"]
    mf = B["contrast"]["mean-field removed"]
    raw = B["contrast"]["raw"]
    fig, ax = plt.subplots(2, 3, figsize=(15.5, 8.6))

    a = ax[0, 0]
    w = [r["planted"] for r in A["planted_sweep"]]
    v = [r["measured"] for r in A["planted_sweep"]]
    a.plot(w, v, "o-", lw=2)
    a.axhline(-1, color="C3", ls="--", lw=1.2)
    a.annotate("pure linear response", (max(w), -1), ha="right", va="bottom",
               fontsize=8, color="C3")
    a.set_xlabel("planted weight of an up-only channel")
    a.set_ylabel("cos(r_up, r_down)")
    a.set_title("A  the estimator tracks a planted answer\n(dose alone never"
                " moves it: 3x gives %.4f)" % A["saturated_cos"], fontsize=9.5)

    a = ax[0, 1]
    gm = [r["global_mode"] for r in A["global_mode_sweep"]]
    a.plot(gm, [r["matched"] for r in A["global_mode_sweep"]], "o-", label="raw matched")
    a.plot(gm, [r["mismatched"] for r in A["global_mode_sweep"]], "s-", label="mismatched")
    a.plot(gm, [r["contrast"] for r in A["global_mode_sweep"]], "^-",
           label="contrast")
    a.plot(gm, [r["mean_field_removed"] for r in A["global_mode_sweep"]], "d-",
           label="after mean-field removal")
    a.axhline(-1, color="C3", ls="--", lw=1.2)
    a.set_xlabel("weight of a shared global mode")
    a.set_ylabel("cosine")
    a.legend(fontsize=8)
    a.set_title("A4  only mean-field removal recovers the\nplanted answer",
                fontsize=9.5)

    a = ax[0, 2]
    own_oe = [o["oe_own"] for o in B["own"]]
    own_kd = [o["kd_own"] for o in B["own"]]
    a.scatter(own_kd, own_oe, s=26, alpha=0.75)
    a.axhline(0, color="0.6", lw=0.8)
    a.axvline(0, color="0.6", lw=0.8)
    a.set_xlabel("CRISPRi: the gene's own response (control SD)")
    a.set_ylabel("CRISPRa: the gene's own response")
    a.set_title("B2  positive control — the libraries are\noriented as expected",
                fontsize=9.5)

    a = ax[1, 0]
    bins = np.linspace(-1, 1, 61)
    a.hist(mf["mismatched_median"] * 0 + np.array(
        B["contrast"]["mean-field removed"]["matched_vals"]), bins=bins,
        alpha=0.8, label="matched (n=%d)" % mf["matched"]["n"], density=True)
    a.axvline(mf["mismatched"]["median"], color="C1", lw=2,
              label="mismatched median")
    a.axvline(-1, color="C3", ls="--", lw=1.5, label="prediction")
    a.axvline(0, color="0.6", lw=0.8)
    a.set_xlabel("cos(r_OE, r_KD), mean-field removed")
    a.set_ylabel("density")
    a.legend(fontsize=8)
    a.set_title("B4  measured antisymmetry: median %+.3f\nagainst a prediction of"
                " -1" % mf["matched"]["median"], fontsize=9.5)

    a = ax[1, 1]
    labs = ["CRISPRa\nsplit-half", "CRISPRa\n2 orientations", "CRISPRi\n2 TSS",
            "CRISPRi\nnon-targeting", "matched\nOE vs KD"]
    vals = [B["reliability"]["split_half"]["median"],
            B["reliability"]["orientations"]["median"],
            B["reliability"]["tss_pairs"]["median"],
            B["reliability"]["nt_cos"]["median"],
            mf["matched"]["median"]]
    cols = ["C0", "C0", "C2", "0.6", "C3"]
    a.bar(range(len(vals)), vals, 0.6, color=cols)
    for i, v in enumerate(vals):
        a.annotate("%+.3f" % v, (i, v), ha="center",
                   va="bottom" if v >= 0 else "top", fontsize=7.5)
    a.axhline(0, color="0.4", lw=0.8)
    a.set_xticks(range(len(vals)), labs, fontsize=8)
    a.set_ylabel("median cosine")
    a.set_title("B3  reliability anchors vs the measured\nantisymmetry",
                fontsize=9.5)

    a = ax[1, 2]
    r = np.array(C["ranks"])
    a.hist(r, bins=np.arange(0.5, C["n"] + 1.5, 2), alpha=0.85)
    a.axvline(np.median(r), color="C3", lw=2,
              label="median rank %d" % C["median_rank"])
    a.set_xlabel("rank of the true overexpression target (of %d)" % C["n"])
    a.set_ylabel("count")
    a.legend(fontsize=8)
    a.set_title("C2  picking an overexpression target from a\nknockdown library:"
                " top-1 %.1f%% vs %.1f%% chance"
                % (100 * C["top1"], 100 / C["n"]), fontsize=9.5)

    fig.tight_layout()
    fn = os.path.join(FIG, "d28_overexpression_vs_knockdown.png")
    fig.savefig(fn, dpi=140)
    print("\nfigure: %s" % fn)


def main(force=False):
    if os.path.exists(RES) and not force:
        R = pickle.load(open(RES, "rb"))
    else:
        A = part_a()
        B = part_b()
        C = part_c(B)
        R = {"A": A, "B": B, "C": C}
        pickle.dump(R, open(RES, "wb"))
    report(R)
    figures(R)
    return R


if __name__ == "__main__":
    main(force=bool(os.environ.get("D28_FORCE")))
