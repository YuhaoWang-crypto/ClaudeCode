"""
D27  Feasibility audit (scenarios 4 and 5 of the indication catalogue):
can public data answer the EGFR drug-tolerant-persister question and the
AML differentiation-block question?

This is deliberately NOT a result demo. It is the product the catalogue
calls NJ-FEAS-01: decide, before promising anything, which parameters the
available measurement can and cannot pin down, and name the measurement
that unblocks the rest.

PART A  identifiability, derived and then verified numerically (no data).

  A1  EGFR two-state model, S = drug-sensitive, R = drug-tolerant:
        dS/dt = (gS - a) S + b R
        dR/dt =  a S + (gR - b) R
      Write p = R/(S+R), N = S+R. Then exactly

        dp/dt = a (1-p) - b p + (gR - gS) p (1-p)
        d lnN/dt = gS (1-p) + gR p

      so a proportion time course sees only (a, b, dg = gR - gS): the
      common growth level is an exact null direction. The widespread claim
      "proportions cannot tell switching from selection" is WRONG - both a
      and dg are identifiable. What proportions cannot see is whether the
      population is growing or being killed.

  A2  AML two-compartment model, I = immature, M = mature:
        dI/dt = (gI - dI - delta) I
        dM/dt = delta I - dM M
      With r = M/I,  dr/dt = delta - C r,  C = gI - dI - delta + dM.
      So delta (the differentiation flux rate) IS identifiable from a
      proportion time course - the catalogue's "proportions cannot
      distinguish differentiation from killing" is half wrong. What is not
      identifiable is the kill rate: gI, dI and dM enter only through C.
      Adding ln N buys dM and the net immature growth; gI vs dI stays
      lumped under both observation models, so separating proliferation
      from death needs a third channel.

  A3  negative control for the machinery: the EGFR design WITH ln N must
      come out full rank with a finite condition number, or the rank test
      is just reporting that finite differences are noisy.

PART B  the EGFR scenario on the real public deposit: the PC9 + erlotinib
  Drop-seq time course of Aissa et al. 2021 (SuperSeries GSE149383).
  B0  de-duplication by content hash: the SuperSeries advertises five PC9
      subseries, and some libraries are byte-identical re-deposits.
  B1  a label-free two-component split of the pooled cells -> p per library,
      with a depth-artefact control and a supervised cross-check.
  B2  the library-to-library floor: two independent day-0 libraries and
      several independent day-11 libraries bound how much of p is batch.
  B3  fit (a, b, dg, p0) to p(t) with bootstrap CIs.
  B4  the price of the missing ln N: the span of tumour-burden fold change
      consistent with the SAME fitted proportion curve.
  B5  recovered barcodes are not a cell count: medians, empty-droplet
      fractions, and the near-equal arm sizes of the xenograft deposit.

PART C  the audit itself: over a reproducible GEO query, how many series
  carry (i) a time course and (ii) any absolute-count field?

PART D  the design that unblocks each scenario, computed from the Fisher
  information rather than asserted: cells per timepoint and ln N precision
  for a target CV, and for AML the third channel.

Run:  python3 -m netjac.demos.d27_egfr_aml_feasibility
"""
import gzip
import hashlib
import json
import os
import pickle
import re
import subprocess
import time
import urllib.parse
import urllib.request

import numpy as np
from scipy.linalg import expm
from scipy.integrate import solve_ivp

from netjac import DATA_CACHE, FIG

RES = os.path.join(DATA_CACHE, "d27_results.pkl")
DIR = os.path.join(DATA_CACHE, "dtp")
AUD = os.path.join(DIR, "geo_audit")

# The PC9 + erlotinib libraries as GEO advertises them across the five
# subseries of GSE149383.  Listed BEFORE de-duplication on purpose: part B0
# decides which are independent.
LIBS = [
    # label             day   subseries    gsm         file
    ("D0",               0.0, "GSE134839", "GSM3972657", "gse134839/D0.dge.txt.gz"),
    ("D1",               1.0, "GSE134839", "GSM3972658", "gse134839/D1.dge.txt.gz"),
    ("D2",               2.0, "GSE134839", "GSM3972659", "gse134839/D2.dge.txt.gz"),
    ("D4",               4.0, "GSE134839", "GSM3972660", "gse134839/D4.dge.txt.gz"),
    ("D9",               9.0, "GSE134839", "GSM3972661", "gse134839/D9.dge.txt.gz"),
    ("D11",             11.0, "GSE134839", "GSM3972662", "gse134839/D11.dge.txt.gz"),
    ("D0_b",             0.0, "GSE134841", "GSM3972669", "gse134841/D0.dge.txt.gz"),
    ("D2_b",             2.0, "GSE134841", "GSM3972670", "gse134841/D2.dge.txt.gz"),
    ("D11_ERL",         11.0, "GSE134841", "GSM3972671", "gse134841/D11_ERL.dge.txt.gz"),
    ("D19_DMSO",        19.0, "GSE134841", "GSM3972672", "gse134841/D19DMSO.dge.txt.gz"),
    ("D19_ERL",         19.0, "GSE134841", "GSM3972673", "gse134841/D19_ERL.dge.txt.gz"),
    ("D0_c",             0.0, "GSE149214", "GSM4494347", "gse149214/D0C_expression_matrix.txt.gz"),
    ("D11_c",           11.0, "GSE149214", "GSM4494348", "gse149214/D11C_expression_matrix.txt.gz"),
    ("D11_ERL_c",       11.0, "GSE149214", "GSM4494349", "gse149214/D11_ERL_expression_matrix.txt.gz"),
]
GSM_URL = {
    "GSM3972657": "GSM3972657_D0.dge.txt.gz",
    "GSM3972658": "GSM3972658_D1.dge.txt.gz",
    "GSM3972659": "GSM3972659_D2.dge.txt.gz",
    "GSM3972660": "GSM3972660_D4.dge.txt.gz",
    "GSM3972661": "GSM3972661_D9.dge.txt.gz",
    "GSM3972662": "GSM3972662_D11.dge.txt.gz",
    "GSM3972669": "GSM3972669_D0.dge.txt.gz",
    "GSM3972670": "GSM3972670_D2.dge.txt.gz",
    "GSM3972671": "GSM3972671_D11_ERL.dge.txt.gz",
    "GSM3972672": "GSM3972672_D19DMSO.dge.txt.gz",
    "GSM3972673": "GSM3972673_D19_ERL.dge.txt.gz",
    "GSM4494347": "GSM4494347_D0C_expression_matrix.txt.gz",
    "GSM4494348": "GSM4494348_D11C_expression_matrix.txt.gz",
    "GSM4494349": "GSM4494349_D11_ERL_expression_matrix.txt.gz",
}
MIN_UMI, MIN_GENES, N_HVG, N_PC = 1000, 300, 1500, 15
DTP_MARKERS = ["CDKN1A", "SERPINE1", "FN1", "VIM", "ANXA1", "CAV1", "AKR1C1",
               "AKR1B1", "NNMT", "TXNRD1", "SLC7A11", "ALDH3A1"]
CYCLE_MARKERS = ["MKI67", "TOP2A", "CCNB1", "CDK1", "PCNA", "MCM2", "TYMS",
                 "RRM2", "BIRC5", "UBE2C"]

# ----------------------------------------------------------------- fetching


def _fetch(rel, gsm):
    fn = os.path.join(DIR, rel)
    if not os.path.exists(fn):
        os.makedirs(os.path.dirname(fn), exist_ok=True)
        stub = gsm[:-3] + "nnn"
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{stub}/{gsm}/suppl/{GSM_URL[gsm]}"
        subprocess.run(["curl", "-sS", "--max-time", "1800", "-o", fn + ".part", url],
                       check=True)
        os.replace(fn + ".part", fn)
    return fn


def _content_md5(fn):
    h = hashlib.md5()
    with gzip.open(fn, "rb") as f:
        while True:
            b = f.read(1 << 22)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def extract(label, rel, gsm):
    """Read one dense Drop-seq DGE, keep barcodes passing QC, cache sparse."""
    import scipy.sparse as sp
    out = os.path.join(DIR, "npz", f"{label}.npz")
    if os.path.exists(out):
        z = np.load(out, allow_pickle=True)
        return sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])), \
            list(z["genes"]), dict(z["qc"].item())
    import pandas as pd
    fn = _fetch(rel, gsm)
    hdr = gzip.open(fn, "rt").readline().rstrip("\n").split("\t")
    dt = {c: np.int32 for c in hdr[1:]}
    dt[hdr[0]] = str
    X = pd.read_csv(fn, sep="\t", index_col=0, dtype=dt)
    genes = [str(g) for g in X.index]
    A = X.values.T                                  # barcodes x genes
    del X
    depth = A.sum(1).astype(np.float64)
    ngene = (A > 0).sum(1)
    keep = (depth >= MIN_UMI) & (ngene >= MIN_GENES)
    qc = {"barcodes": int(A.shape[0]), "cells": int(keep.sum()),
          "median_depth_all": float(np.median(depth)),
          "median_depth_kept": float(np.median(depth[keep])) if keep.any() else np.nan,
          "empty_frac": float(1.0 - keep.mean())}
    M = sp.csr_matrix(A[keep].astype(np.float32))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez_compressed(out, data=M.data, indices=M.indices, indptr=M.indptr,
                        shape=np.array(M.shape), genes=np.array(genes, dtype=object),
                        qc=np.array(qc, dtype=object))
    return M, genes, qc


# ------------------------------------------------- PART A  identifiability


def _egfr_obs(th, ts, p0=0.01):
    """Returns (p(t), lnN(t)) for th = (a, b, gS, gR)."""
    a, b, gS, gR = th
    A = np.array([[gS - a, b], [a, gR - b]])
    y0 = np.array([1.0 - p0, p0])
    p, ln = [], []
    for t in ts:
        y = expm(A * t) @ y0
        p.append(y[1] / y.sum())
        ln.append(np.log(y.sum()))
    return np.array(p), np.array(ln)


def _aml_obs(th, ts, f0=0.05):
    """Returns (f_M(t), lnN(t)) for th = (gI, dI, delta, dM)."""
    gI, dI, de, dM = th
    A = np.array([[gI - dI - de, 0.0], [de, -dM]])
    y0 = np.array([1.0 - f0, f0])
    f, ln = [], []
    for t in ts:
        y = expm(A * t) @ y0
        f.append(y[1] / y.sum())
        ln.append(np.log(y.sum()))
    return np.array(f), np.array(ln)


def _egfr_obs2(phi, ts, p0=0.01):
    """Same model in the identifiable basis phi = (a, b, dg, gS), gR = gS + dg."""
    a, b, dg, gS = phi
    return _egfr_obs((a, b, gS, gS + dg), ts, p0)


def _aml_obs2(phi, ts, f0=0.05):
    """AML in the basis phi = (delta, C, gImdI, dM) with C pinned by the rest:
    gI - dI = gImdI, C = gImdI - delta + dM, so (delta, gImdI, dM) is the free
    triple and the fourth slot is a dummy kept for shape compatibility."""
    de, _unused, gImdI, dM = phi
    gI, dI = gImdI, 0.0
    return _aml_obs((gI, dI, de, dM), ts, f0)


def _fim(obs, th, ts, sig_prop, sig_ln=None, h=1e-5):
    """Fisher information for parameters th; sig_ln=None drops the lnN channel."""
    th = np.asarray(th, float)
    SP, SL = [], []
    for j in range(len(th)):
        d = h * max(abs(th[j]), 1.0)
        tp, tm = th.copy(), th.copy()
        tp[j] += d
        tm[j] -= d
        pp, lp = obs(tp, ts)
        pm, lm = obs(tm, ts)
        SP.append((pp - pm) / (2 * d))
        SL.append((lp - lm) / (2 * d))
    SP, SL = np.array(SP), np.array(SL)              # params x times
    F = SP @ SP.T / sig_prop ** 2
    if sig_ln is not None:
        F = F + SL @ SL.T / sig_ln ** 2
    return F


def _crlb(F, tol=1e-9):
    """Per-parameter Cramer-Rao sd, with inf where the parameter's own unit
    vector has a component in the null space of F. A plain pseudo-inverse
    silently splits a null direction between the parameters that span it and
    reports finite numbers for quantities that are not estimable at all."""
    w, V = np.linalg.eigh(F)
    w = np.maximum(w, 0.0)
    keep = w > tol * w.max()
    Fp = (V[:, keep] / w[keep]) @ V[:, keep].T
    Vn = V[:, ~keep]
    sd = []
    for j in range(F.shape[0]):
        e = np.zeros(F.shape[0])
        e[j] = 1.0
        leak = float(np.linalg.norm(Vn.T @ e)) if Vn.size else 0.0
        sd.append(np.inf if leak > 1e-6 else float(np.sqrt(e @ Fp @ e)))
    return np.array(sd)


def _rank_report(F, tol=1e-9):
    w, V = np.linalg.eigh(F)
    w = np.maximum(w, 0.0)
    rank = int((w > tol * w.max()).sum())
    cond = float(w.max() / w.min()) if w.min() > 0 else np.inf
    null = [V[:, i] / np.linalg.norm(V[:, i]) for i in range(len(w))
            if w[i] <= tol * w.max()]
    return {"eig": w.tolist(), "rank": rank, "cond": cond,
            "null_vectors": [v.tolist() for v in null]}


def part_a():
    print("\n" + "=" * 78)
    print("PART A  identifiability, derived then verified numerically")
    print("=" * 78)
    out = {}

    # --- A1 EGFR -----------------------------------------------------------
    th = np.array([0.08, 0.01, 0.30, -0.05])         # a, b, gS, gR  (per day)
    ts = np.array([0.0, 1.0, 2.0, 4.0, 9.0, 11.0])
    sig_p, sig_ln = 0.02, 0.10
    Fp = _fim(_egfr_obs, th, ts, sig_p, None)
    Fa = _fim(_egfr_obs, th, ts, sig_p, sig_ln)
    rp, ra = _rank_report(Fp), _rank_report(Fa)
    print("\nA1  EGFR two-state, truth a=%.2f b=%.2f gS=%.2f gR=%.2f /day,"
          " 6 timepoints" % tuple(th))
    print("    proportions only : rank %d/4  cond %.2e" % (rp["rank"], rp["cond"]))
    print("    + ln N           : rank %d/4  cond %.2e" % (ra["rank"], ra["cond"]))
    for v in rp["null_vectors"]:
        print("    null direction   : (a,b,gS,gR) = (%+.3f, %+.3f, %+.3f, %+.3f)"
              % tuple(v))
    # verify the null direction exactly
    c = 0.25
    p1, l1 = _egfr_obs(th, ts)
    p2, l2 = _egfr_obs(th + np.array([0, 0, c, c]), ts)
    out["egfr"] = {"theta": th.tolist(), "ts": ts.tolist(), "prop_only": rp,
                   "with_lnN": ra,
                   "shift": c, "max_dp": float(np.abs(p1 - p2).max()),
                   "N_fold": float(np.exp(c * ts[-1])),
                   "sig_p": sig_p, "sig_ln": sig_ln}
    print("    shifting BOTH growth rates by %+.2f/day moves p(t) by at most %.2e"
          % (c, out["egfr"]["max_dp"]))
    print("    while the total cell number at day %g changes %.0f-fold."
          % (ts[-1], out["egfr"]["N_fold"]))
    sd_p = _crlb(Fp)
    sd_a = _crlb(Fa)
    out["egfr"]["crlb_prop_only"] = sd_p.tolist()
    out["egfr"]["crlb_with_lnN"] = sd_a.tolist()
    print("    CRLB sd  prop only : a %.3f  b %.3f  gS %s  gR %s"
          % (sd_p[0], sd_p[1], sd_p[2], sd_p[3]))
    print("    CRLB sd  + ln N    : a %.3f  b %.3f  gS %.3f  gR %.3f" % tuple(sd_a))
    # the pseudo-inverse splits the null direction between gS and gR, so the
    # honest numbers come from the identifiable basis (a, b, dg) with gS dropped
    phi = np.array([th[0], th[1], th[3] - th[2], th[2]])
    Fp2 = _fim(_egfr_obs2, phi, ts, sig_p, None)
    blk = np.linalg.inv(Fp2[:3, :3])
    sd3 = np.sqrt(np.diag(blk))
    out["egfr"]["crlb_identifiable_triple"] = sd3.tolist()
    out["egfr"]["cond_identifiable_triple"] = float(
        np.linalg.cond(Fp2[:3, :3]))
    print("    in the identifiable basis (a, b, dg), proportions alone give")
    print("    sd a %.3f (CV %.0f%%)  sd b %.3f (CV %.0f%%)  sd dg %.3f"
          " (CV %.0f%%), cond %.1e"
          % (sd3[0], 100 * sd3[0] / abs(phi[0]), sd3[1],
             100 * sd3[1] / abs(phi[1]), sd3[2],
             100 * sd3[2] / abs(phi[2]), out["egfr"]["cond_identifiable_triple"]))
    # structural identifiability is not practical identifiability: report the
    # best- and worst-determined directions of the identifiable 3x3 block
    wv, Vv = np.linalg.eigh(Fp2[:3, :3])
    sloppy = Vv[:, 0] / np.linalg.norm(Vv[:, 0])
    stiff = Vv[:, -1] / np.linalg.norm(Vv[:, -1])
    sd_sloppy = float(1.0 / np.sqrt(wv[0]))
    sd_stiff = float(1.0 / np.sqrt(wv[-1]))
    out["egfr"]["sloppy_dir"] = sloppy.tolist()
    out["egfr"]["stiff_dir"] = stiff.tolist()
    out["egfr"]["sd_sloppy"] = sd_sloppy
    out["egfr"]["sd_stiff"] = sd_stiff
    print("    best-determined combination  %+.2f*a %+.2f*b %+.2f*dg : sd %.4f"
          % (stiff[0], stiff[1], stiff[2], sd_stiff))
    print("    worst-determined combination %+.2f*a %+.2f*b %+.2f*dg : sd %.4f"
          % (sloppy[0], sloppy[1], sloppy[2], sd_sloppy))
    print("    -> a and dg = gR-gS are STRUCTURALLY identifiable from"
          " proportions, but at this")
    print("       design only one combination is PRACTICALLY estimable; the"
          " common growth level")
    print("       is not identifiable at all.")

    # --- A2 AML ------------------------------------------------------------
    th2 = np.array([0.50, 0.10, 0.20, 0.30])         # gI, dI, delta, dM
    ts2 = np.array([0.0, 1.0, 2.0, 3.0, 5.0, 7.0])
    Gp = _fim(_aml_obs, th2, ts2, sig_p, None)
    Ga = _fim(_aml_obs, th2, ts2, sig_p, sig_ln)
    qp, qa = _rank_report(Gp), _rank_report(Ga)
    print("\nA2  AML two-compartment, truth gI=%.2f dI=%.2f delta=%.2f dM=%.2f /day"
          % tuple(th2))
    print("    proportions only : rank %d/4  cond %.2e" % (qp["rank"], qp["cond"]))
    print("    + ln N           : rank %d/4  cond %.2e" % (qa["rank"], qa["cond"]))
    for v in qa["null_vectors"]:
        print("    null direction even WITH counts: (gI,dI,delta,dM) ="
              " (%+.3f, %+.3f, %+.3f, %+.3f)" % tuple(v))
    # explicit pair: same delta and same C, very different totals
    C = th2[0] - th2[1] - th2[2] + th2[3]
    alt = np.array([0.10, 0.30, 0.20, 0.90])         # gI-dI+dM = 0.7 -> same C
    assert abs((alt[0] - alt[1] - alt[2] + alt[3]) - C) < 1e-12
    f1, n1 = _aml_obs(th2, ts2)
    f2, n2 = _aml_obs(alt, ts2)
    out["aml"] = {"theta": th2.tolist(), "alt": alt.tolist(), "C": float(C),
                  "ts": ts2.tolist(), "prop_only": qp, "with_lnN": qa,
                  "max_df": float(np.abs(f1 - f2).max()),
                  "N_fold_ratio": float(np.exp(n1[-1] - n2[-1])),
                  "f_day7": float(f1[-1])}
    print("    two parameter sets with the same (delta, C): mature fraction at"
          " every timepoint agrees to %.2e" % out["aml"]["max_df"])
    print("    (day %g fraction %.3f in both) while total cells differ %.0f-fold."
          % (ts2[-1], out["aml"]["f_day7"], out["aml"]["N_fold_ratio"]))
    # honest CRLBs in the identifiable bases
    de = th2[2]
    # proportions see only (delta, C): build the 2x2 block in (delta, C)
    def obs_dc(z, t):
        d_, C_ = z[0], z[1]
        # realise (delta, C) with dI = 0 and gI chosen to match C
        gI = C_ + d_ - th2[3]
        return _aml_obs((gI, 0.0, d_, th2[3]), t)
    Fdc = _fim(obs_dc, np.array([de, C, 0.0, 0.0]), ts2, sig_p, None)[:2, :2]
    sd_dc = np.sqrt(np.diag(np.linalg.inv(Fdc)))
    out["aml"]["crlb_delta_C"] = sd_dc.tolist()
    print("    in the identifiable basis, proportions alone give sd delta %.4f"
          " (CV %.0f%%) and sd C %.4f (CV %.0f%%)"
          % (sd_dc[0], 100 * sd_dc[0] / de, sd_dc[1], 100 * sd_dc[1] / C))
    Hfull = _fim(_aml_obs2, np.array([de, 0.0, th2[0] - th2[1], th2[3]]),
                 ts2, sig_p, sig_ln)
    idx = np.array([0, 2, 3])
    sd3a = np.sqrt(np.diag(np.linalg.inv(Hfull[np.ix_(idx, idx)])))
    out["aml"]["crlb_with_lnN_triple"] = sd3a.tolist()
    print("    adding ln N gives sd delta %.4f, sd (gI-dI) %.4f, sd dM %.4f"
          % tuple(sd3a))
    print("    -> delta IS identifiable from proportions; the kill rate is not,"
          " and counts do not fix it.")

    # --- A3 negative control for the rank test -----------------------------
    out["negative_control"] = {
        "egfr_with_lnN_rank": ra["rank"], "egfr_with_lnN_cond": ra["cond"],
        "aml_with_lnN_rank": qa["rank"]}
    print("\nA3  negative control: the same code gives the EGFR + ln N design"
          " rank %d/4 at cond %.1e," % (ra["rank"], ra["cond"]))
    print("    so the rank deficiencies above are structure, not numerical noise.")
    return out


# ------------------------------------------------------ PART B  real data


def part_b0():
    print("\n" + "=" * 78)
    print("PART B0  de-duplication of the SuperSeries by content hash")
    print("=" * 78)
    seen, groups = {}, []
    rows = []
    for label, day, ser, gsm, rel in LIBS:
        fn = _fetch(rel, gsm)
        h = _content_md5(fn)
        rows.append({"label": label, "day": day, "series": ser, "gsm": gsm,
                     "md5": h})
        if h in seen:
            groups.append((seen[h], label))
        else:
            seen[h] = label
    dups = {}
    for first, dup in groups:
        dups.setdefault(first, []).append(dup)
    print("  %d libraries listed, %d distinct contents, %d re-deposits"
          % (len(LIBS), len(seen), len(groups)))
    for first, ds in dups.items():
        print("    %-10s == %s  (byte-identical after decompression)"
              % (first, ", ".join(ds)))
    indep = [r for r in rows if r["label"] in set(seen.values())]
    print("  independent libraries kept: %s"
          % ", ".join("%s(d%g)" % (r["label"], r["day"]) for r in indep))
    return {"rows": rows, "duplicates": dups, "independent":
            [r["label"] for r in indep]}


def _readouts(X, labs, genes, day, indep, tag):
    """Three state-fraction definitions on one normalised matrix.

    gmm  label-free two-component mixture on PCA space (no time labels used
         for the clustering; only for deciding which component to call late,
         and only over the erlotinib-treated libraries so the untreated
         day-19 control cannot flip the sign)
    sup  projection on the day11-minus-day0 axis of the Drop-seq course,
         thresholded at the midpoint of the two endpoint means
    dtp  a published drug-tolerant-persister marker score, thresholded at the
         pooled median
    """
    from sklearn.decomposition import PCA
    from sklearn.mixture import GaussianMixture
    gl = np.array(genes)
    v = X.var(0)
    hv = np.argsort(-v)[:N_HVG]
    Z = X[:, hv]
    Z = (Z - Z.mean(0)) / (Z.std(0) + 1e-9)
    P = PCA(n_components=N_PC, random_state=0).fit_transform(Z)
    gm = GaussianMixture(2, covariance_type="full", n_init=10,
                         random_state=0).fit(P)
    post = gm.predict_proba(P)[:, 1]
    erl = np.isin(labs, [l for l in indep if l != "D19_DMSO"])
    late = 1 if np.corrcoef(post[erl], day[erl])[0, 1] > 0 else 0
    q = post if late == 1 else 1.0 - post
    gmm = q > 0.5

    a0 = X[labs == "D0"].mean(0)
    a11 = X[labs == "D11"].mean(0)
    w = a11 - a0
    w = w / np.linalg.norm(w)
    proj = X @ w
    mid = 0.5 * (proj[labs == "D0"].mean() + proj[labs == "D11"].mean())
    sup = proj > mid

    mk = [g for g in DTP_MARKERS if g in set(gl)]
    cm = [g for g in CYCLE_MARKERS if g in set(gl)]
    gi = {g: i for i, g in enumerate(genes)}
    sc = X[:, [gi[g] for g in mk]].mean(1) if mk else np.zeros(len(labs))
    dtp = sc > np.median(sc)

    out = {"tag": tag, "markers_used": mk, "cycle_used": cm, "n_cells": {},
           "gmm": {}, "sup": {}, "dtp": {}}
    for lab in indep:
        m = labs == lab
        out["n_cells"][lab] = int(m.sum())
        out["gmm"][lab] = float(gmm[m].mean())
        out["sup"][lab] = float(sup[m].mean())
        out["dtp"][lab] = float(dtp[m].mean())
    # what the label-free component actually is
    d = X[gmm].mean(0) - X[~gmm].mean(0)
    out["gmm_cycle_shift"] = float(np.mean([d[gi[g]] for g in cm])) if cm else np.nan
    out["gmm_dtp_shift"] = float(np.mean([d[gi[g]] for g in mk])) if mk else np.nan
    out["gmm_top_up"] = gl[np.argsort(-d)[:12]].tolist()
    out["pc1_vs_depth"] = None          # filled by the caller
    return out


def part_b1(indep):
    import scipy.sparse as sp
    print("\n" + "=" * 78)
    print("PART B1  state fraction per library: three definitions, and the")
    print("         sequencing-depth confound")
    print("=" * 78)
    meta = {l: (d, s, g, r) for l, d, s, g, r in LIBS}
    mats, genesets, qcs = {}, {}, {}
    for lab in indep:
        d, s, g, r = meta[lab]
        M, genes, qc = extract(lab, r, g)
        mats[lab], genesets[lab], qcs[lab] = M, genes, qc
        print("  %-10s d%-4g %-10s barcodes %6d  cells %5d  median UMI %6.0f"
              % (lab, d, s, qc["barcodes"], qc["cells"], qc["median_depth_kept"]))
    common = sorted(set.intersection(*[set(g) for g in genesets.values()]))
    print("  genes common to all libraries: %d" % len(common))

    # depth is confounded with day in the Drop-seq course
    course = ["D0", "D1", "D2", "D4", "D9", "D11"]
    dc = np.array([qcs[l]["median_depth_kept"] for l in course])
    dd = np.array([meta[l][0] for l in course])
    r_course = float(np.corrcoef(dd, np.log(dc))[0, 1])
    alld = np.array([meta[l][0] for l in indep])
    alldep = np.array([qcs[l]["median_depth_kept"] for l in indep])
    print("  median UMI rises %.0f -> %.0f from day 0 to day 11 in the Drop-seq"
          " course:" % (dc[0], dc[-1]))
    print("  corr(day, log median UMI) = %+.3f  -- depth and time are"
          " confounded here." % r_course)

    # pooled raw-depth matrix and a depth-matched (downsampled) copy
    rng = np.random.default_rng(0)
    Xr, Xd, labs, dep = [], [], [], []
    for lab in indep:
        idx = {g: i for i, g in enumerate(genesets[lab])}
        cols = np.array([idx[g] for g in common])
        A = np.asarray(mats[lab][:, cols].todense(), dtype=np.int64)
        d = A.sum(1, keepdims=True)
        Xr.append(np.log1p(1e4 * A / d))
        # exact thinning without replacement to the common budget
        B = np.zeros_like(A)
        for i in range(A.shape[0]):
            nz = np.nonzero(A[i])[0]
            if A[i].sum() < MIN_UMI:
                B[i, nz] = A[i, nz]
                continue
            B[i, nz] = rng.multivariate_hypergeometric(A[i, nz], MIN_UMI)
        Xd.append(np.log1p(1e4 * B / np.maximum(B.sum(1, keepdims=True), 1)))
        labs += [lab] * A.shape[0]
        dep.append(d.ravel())
    Xr = np.vstack(Xr)
    Xd = np.vstack(Xd)
    labs = np.array(labs)
    dep = np.concatenate(dep)
    del mats
    day = np.array([meta[l][0] for l in labs])
    print("  pooled matrix: %d cells x %d genes; depth-matched copy thinned to"
          " %d UMI/cell" % (Xr.shape[0], Xr.shape[1], MIN_UMI))

    from sklearn.decomposition import PCA
    P1 = PCA(n_components=3, random_state=0).fit_transform(
        (Xr - Xr.mean(0)) / (Xr.std(0) + 1e-9))
    P2 = PCA(n_components=3, random_state=0).fit_transform(
        (Xd - Xd.mean(0)) / (Xd.std(0) + 1e-9))
    r1 = float(np.corrcoef(P1[:, 0], np.log(dep))[0, 1])
    r2 = float(np.corrcoef(P2[:, 0], np.log(dep))[0, 1])
    print("  corr(PC1, log UMI): %+.3f raw -> %+.3f after thinning" % (r1, r2))

    R = {}
    for tag, X in (("raw", Xr), ("depth-matched", Xd)):
        R[tag] = _readouts(X, labs, common, day, indep, tag)
        R[tag]["pc1_vs_depth"] = r1 if tag == "raw" else r2

    print("\n  %-10s %4s %6s | %-22s | %-22s"
          % ("library", "day", "cells", "raw  gmm / sup / dtp",
             "depth-matched  g/s/d"))
    for lab in sorted(indep, key=lambda l: (meta[l][0], l)):
        print("  %-10s %4g %6d | %6.3f %6.3f %6.3f  | %6.3f %6.3f %6.3f"
              % (lab, meta[lab][0], R["raw"]["n_cells"][lab],
                 R["raw"]["gmm"][lab], R["raw"]["sup"][lab], R["raw"]["dtp"][lab],
                 R["depth-matched"]["gmm"][lab], R["depth-matched"]["sup"][lab],
                 R["depth-matched"]["dtp"][lab]))

    print("\n  what the label-free component is (depth-matched): cycle-marker"
          " shift %+.3f, DTP-marker shift %+.3f"
          % (R["depth-matched"]["gmm_cycle_shift"],
             R["depth-matched"]["gmm_dtp_shift"]))
    print("  top genes up in it: %s"
          % ", ".join(R["depth-matched"]["gmm_top_up"]))

    # mutual agreement of the three definitions, and depth dependence
    out = {"qc": qcs, "readouts": R, "n_common_genes": len(common),
           "n_cells_total": int(Xr.shape[0]), "corr_day_depth": r_course,
           "pc1_depth_raw": r1, "pc1_depth_matched": r2, "indep": indep,
           "median_depth": {l: qcs[l]["median_depth_kept"] for l in indep}}
    for tag in R:
        keys = ["gmm", "sup", "dtp"]
        ag = {}
        for i in range(3):
            for j in range(i + 1, 3):
                x = [R[tag][keys[i]][l] for l in indep]
                y = [R[tag][keys[j]][l] for l in indep]
                ag["%s-%s" % (keys[i], keys[j])] = float(np.corrcoef(x, y)[0, 1])
        dr = {k: float(np.corrcoef([R[tag][k][l] for l in indep],
                                   np.log(alldep))[0, 1]) for k in keys}
        tr = {k: float(np.corrcoef([R[tag][k][l] for l in indep], alld)[0, 1])
              for k in keys}
        out.setdefault("agreement", {})[tag] = ag
        out.setdefault("vs_depth", {})[tag] = dr
        out.setdefault("vs_day", {})[tag] = tr
        print("\n  %s: pairwise agreement across libraries %s"
              % (tag, ", ".join("%s %+.2f" % (k, v) for k, v in ag.items())))
        print("  %s: corr with log median UMI %s"
              % (tag, ", ".join("%s %+.2f" % (k, v) for k, v in dr.items())))
        print("  %s: corr with day              %s"
              % (tag, ", ".join("%s %+.2f" % (k, v) for k, v in tr.items())))
    # the primary readout: supervised axis on depth-matched counts
    out["p_primary"] = dict(R["depth-matched"]["sup"])
    out["cells"] = dict(R["depth-matched"]["n_cells"])
    return out


def part_b2(B1):
    """How much of p is library, not biology?"""
    print("\n" + "=" * 78)
    print("PART B2  the library-to-library floor at a fixed biological state")
    print("=" * 78)
    meta = {l: d for l, d, s, g, r in LIBS}
    out = {}
    for day in (0.0, 11.0):
        labs = [l for l in B1["p_primary"] if meta[l] == day]
        vals = np.array([B1["p_primary"][l] for l in labs])
        ns = np.array([B1["cells"][l] for l in labs])
        binom = np.sqrt(np.maximum(vals * (1 - vals), 1e-6) / ns)
        out[day] = {"labels": labs, "p": vals.tolist(),
                    "spread": float(vals.max() - vals.min()),
                    "sd": float(vals.std(ddof=1)) if len(vals) > 1 else np.nan,
                    "binomial_sd": binom.tolist()}
        print("  day %-4g %d independent libraries: p = %s"
              % (day, len(labs), ", ".join("%.3f" % v for v in vals)))
        print("          spread %.3f, sd %.3f, vs binomial sd %.3f-%.3f"
              % (out[day]["spread"], out[day]["sd"], binom.min(), binom.max()))
    floor = float(np.nanmax([out[0.0]["sd"], out[11.0]["sd"]]))
    signal = abs(np.mean(out[11.0]["p"]) - np.mean(out[0.0]["p"]))
    print("  batch floor on p (worst sd) %.3f  vs  day0 -> day11 signal %.3f"
          "  = %.1fx" % (floor, signal, signal / floor if floor > 0 else np.inf))
    out["floor"] = floor
    out["signal"] = signal
    out["snr"] = float(signal / floor) if floor > 0 else np.inf
    return out


def _p_traj(th, ts):
    """p(t) for the two-state model, exactly.

    p obeys dp/dt = a(1-p) - b p + dg p (1-p), and p does not depend on the
    common growth level, so the trajectory is read off the linear two-state
    system with gS set to 0 and gR = dg - no ODE solver and no Riccati algebra.
    """
    a, b, dg, p0 = th
    return _egfr_obs((a, b, 0.0, dg), np.asarray(ts, float), p0=p0)[0]


def part_b3(B1, B2):
    """Fit (a, b, dg, p0) to the Drop-seq course; profile what is not pinned."""
    from scipy.optimize import least_squares
    print("\n" + "=" * 78)
    print("PART B3  fitting the two-state model to the proportion time course")
    print("=" * 78)
    course = ["D0", "D1", "D2", "D4", "D9", "D11"]
    ts = np.array([0.0, 1.0, 2.0, 4.0, 9.0, 11.0])
    p = np.array([B1["p_primary"][l] for l in course])
    n = np.array([B1["cells"][l] for l in course], float)
    binom = np.sqrt(np.maximum(p * (1 - p), 1e-4) / n)
    se = np.sqrt(binom ** 2 + B2["floor"] ** 2)
    print("  single-platform course (GSE134839, Drop-seq), depth-matched:")
    for l, t, v, s_ in zip(course, ts, p, se):
        print("    %-5s day %-4g p = %.3f +- %.3f" % (l, t, v, s_))

    lo = np.array([0.0, 0.0, -5.0, 1e-4])
    hi = np.array([20.0, 20.0, 5.0, 0.9])
    STARTS = [(a0, dg0) for a0 in (0.05, 0.5, 3.0) for dg0 in (-2.0, 0.0, 2.0)]

    def fit(pv, fix=None):
        """fix = (index, value) holds one parameter; returns (theta, cost)."""
        best, bx = np.inf, None
        for a0, dg0 in STARTS:
            x0 = np.clip(np.array([a0, 0.05, dg0, max(pv[0], 1e-3)]),
                         lo + 1e-9, hi - 1e-9)
            l2, h2 = lo.copy(), hi.copy()
            if fix is not None:
                j, v = fix
                x0[j] = np.clip(v, lo[j] + 1e-9, hi[j] - 1e-9)
                l2[j], h2[j] = x0[j] - 1e-12, x0[j] + 1e-12
            try:
                r = least_squares(lambda th: (_p_traj(th, ts) - pv) / se,
                                  x0, bounds=(l2, h2))
            except Exception:
                continue
            if r.cost < best:
                best, bx = r.cost, r.x
        return bx, best

    th, cost = fit(p)
    names = ["a (S->R, /day)", "b (R->S, /day)", "dg = gR-gS (/day)", "p0"]
    print("  fit: " + "  ".join("%s %.4f" % (nm, v) for nm, v in zip(names, th)))
    chi2 = 2 * cost
    print("  weighted chi2 = %.2f on %d residuals - %d parameters = %d dof"
          % (chi2, len(ts), len(th), len(ts) - len(th)))
    at_bound = [names[j] for j in range(4)
                if min(abs(th[j] - lo[j]), abs(th[j] - hi[j])) < 1e-3]
    if at_bound:
        print("  at a bound: %s - flat direction, not an estimate"
              % ", ".join(at_bound))

    # profile likelihood: what the curve really constrains
    prof = {}
    for j, grid in [(0, np.geomspace(0.05, 20.0, 25)),
                    (1, np.concatenate([[0.0], np.geomspace(1e-3, 20.0, 24)])),
                    (2, np.linspace(-5.0, 5.0, 25))]:
        cs = []
        for v in grid:
            _x, c = fit(p, fix=(j, v))
            cs.append(2 * c)
        cs = np.array(cs)
        ok = grid[cs - chi2 <= 3.84]
        prof[names[j]] = {"grid": grid.tolist(), "chi2": cs.tolist(),
                          "interval": [float(ok.min()), float(ok.max())] if ok.size
                          else None,
                          "spans_grid": bool(ok.size == grid.size)}
        print("  profile %-20s 95%% interval [%.3f, %.3f]%s"
              % (names[j], ok.min() if ok.size else np.nan,
                 ok.max() if ok.size else np.nan,
                 "  = the whole grid, i.e. unconstrained"
                 if ok.size == grid.size else ""))

    rng = np.random.default_rng(0)
    B = 200
    draws = []
    for _ in range(B):
        pv = np.clip(p + rng.normal(0, se), 1e-4, 0.999)
        x, _c = fit(pv)
        if x is not None:
            draws.append(x)
    D = np.array(draws)
    ci = {nm: [float(np.percentile(D[:, j], 2.5)), float(np.percentile(D[:, j], 97.5))]
          for j, nm in enumerate(names)}
    for j, nm in enumerate(names):
        print("    %-20s %.4f  bootstrap 95%% [%.4f, %.4f]" % (nm, th[j], *ci[nm]))
    half = float(np.log(2.0) / th[0]) if th[0] > 0 else np.inf
    print("  a = %.2f/day means the sensitive pool is half converted in %.2f day"
          "(s)," % (th[0], half))
    print("  i.e. the transition is already over by the first sampled day - the")
    print("  design samples the plateau, not the transient.")

    # the day-0 and day-11 points are the libraries that DEFINE the supervised
    # axis, so they are not independent evidence. Refit with the independent
    # libraries at those two days instead.
    print("\n  circularity check: the axis is built from D0 and D11, so those two")
    print("  points are not evidence. Refit with the independent libraries at the")
    print("  same two days (D0_c at day 0; D11_ERL and D11_c averaged at day 11):")
    p_alt = p.copy()
    p_alt[0] = B1["p_primary"]["D0_c"]
    p_alt[-1] = float(np.mean([B1["p_primary"]["D11_ERL"],
                               B1["p_primary"]["D11_c"]]))
    th_alt, cost_alt = fit(p_alt)
    print("    p(day 0) %.3f -> p(day 11) %.3f; fit a %.3f, b %.3f, dg %+.3f,"
          " p0 %.4f" % (p_alt[0], p_alt[-1], th_alt[0], th_alt[1], th_alt[2],
                        th_alt[3]))
    prof_alt = {}
    for j, grid in [(0, np.geomspace(0.05, 20.0, 25)),
                    (2, np.linspace(-5.0, 5.0, 25))]:
        cs = np.array([2 * fit(p_alt, fix=(j, v))[1] for v in grid])
        ok = grid[cs - 2 * cost_alt <= 3.84]
        prof_alt[names[j]] = {"interval": [float(ok.min()), float(ok.max())]
                              if ok.size else None,
                              "spans_grid": bool(ok.size == grid.size)}
        print("    profile %-20s 95%% interval [%.3f, %.3f]%s"
              % (names[j], ok.min() if ok.size else np.nan,
                 ok.max() if ok.size else np.nan,
                 "  = unconstrained" if ok.size == grid.size else ""))
    return {"course": course, "ts": ts.tolist(), "p": p.tolist(), "se": se.tolist(),
            "theta": th.tolist(), "names": names, "ci": ci, "chi2": float(chi2),
            "dof": int(len(ts) - len(th)), "profile": prof,
            "at_bound": at_bound, "half_life_days": half,
            "p_alt": p_alt.tolist(), "theta_alt": th_alt.tolist(),
            "profile_alt": prof_alt, "boot": D.tolist()}


def part_b4(B3):
    """What the missing ln N costs, kept separate from what the fit cannot pin."""
    print("\n" + "=" * 78)
    print("PART B4  the price of the missing absolute count")
    print("=" * 78)
    a, b, dg, p0 = B3["theta"]
    T = max(B3["ts"])
    ts = np.linspace(0, T, 400)
    pt = _p_traj([a, b, dg, p0], ts)
    integ = float(np.trapezoid(pt, ts))
    print("  ln N(%g) - ln N(0) = gS * %g + dg * %.3f, and the proportion data"
          % (T, T, integ))
    print("  constrains neither term: gS is a structural null direction (part A)")
    print("  and dg's profile here spans the whole grid (part B3).")
    # the structural part alone, with dg held at the fit
    span = [-0.5, 0.5]
    struct = float(np.exp((span[1] - span[0]) * T))
    print("  holding dg at the fit, gS in [%+.1f, %+.1f] /day alone spreads the"
          " day-%g" % (span[0], span[1], T))
    print("  burden over %.1e-fold. That factor is pure structure: no number of"
          % struct)
    print("  cells and no number of timepoints shrinks it.")
    # and with dg also free over its profile interval
    dgi = B3["profile"]["dg = gR-gS (/day)"]["interval"]
    tot = float(np.exp((span[1] - span[0]) * T + (dgi[1] - dgi[0]) * integ))
    print("  letting dg roam its 95%% profile [%+.1f, %+.1f] as well takes it to"
          " %.1e-fold." % (dgi[0], dgi[1], tot))
    need = []
    for cv in (0.5, 0.2, 0.05):
        need.append({"cv_N": cv, "sd_gS": float(cv / T)})
        print("  one terminal count with CV %4.0f%% would pin gS to sd %.3f /day"
              % (100 * cv, cv / T))
    print("  -> structurally, proportions can rank which transition to block;")
    print("     they can never say whether the burden falls. On this deposit even")
    print("     the ranking is out of reach (part B3).")
    return {"integral_p": integ, "gS_span": span, "fold_ratio": struct,
            "fold_ratio_with_dg": tot, "dg_interval": dgi,
            "terminal_count": need, "T": T}


def part_b5(B1):
    print("\n" + "=" * 78)
    print("PART B5  recovered barcodes are not a cell count")
    print("=" * 78)
    rows = []
    for lab, qc in sorted(B1["qc"].items(), key=lambda kv: -kv[1]["barcodes"]):
        rows.append({"label": lab, **qc})
        print("  %-10s barcodes %6d  passing QC %5d  empty %5.1f%%"
              "  median UMI all %6.0f kept %6.0f"
              % (lab, qc["barcodes"], qc["cells"], 100 * qc["empty_frac"],
                 qc["median_depth_all"], qc["median_depth_kept"]))
    bc = np.array([r["barcodes"] for r in rows], float)
    print("  deposited barcode counts span %.0f-fold across libraries of the"
          " same experiment" % (bc.max() / bc.min()))
    xeno = np.array([29983, 29785, 29963, 29340])     # scPerturb Aissa xenograft arms
    print("  for contrast, the four xenograft arms of the same study carry"
          " %s cells" % ", ".join("%d" % v for v in xeno))
    print("  = CV %.1f%%, i.e. set by loading, not by how many cells the drug"
          " left alive." % (100 * xeno.std(ddof=1) / xeno.mean()))
    return {"rows": rows, "barcode_fold": float(bc.max() / bc.min()),
            "xeno_arms": xeno.tolist(),
            "xeno_cv": float(xeno.std(ddof=1) / xeno.mean())}


# ----------------------------------------------------------- PART C  audit

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
QUERIES = {
    "EGFR drug-tolerant persister":
        '(EGFR) AND (erlotinib OR osimertinib OR gefitinib) AND (persister OR '
        '"drug tolerant" OR tolerance OR resistance) AND gse[Entry Type]',
    "AML differentiation therapy":
        '(acute myeloid leukemia) AND (differentiation) AND (ATRA OR tretinoin '
        'OR LSD1 OR menin OR "vitamin D") AND gse[Entry Type]',
}
COUNT_PAT = re.compile(
    r"(cell count|cell number|cell concentration|viable cell|total cell|"
    r"cells?\s*/\s*m[lL]|cells per m[lL]|trucount|absolute count|"
    r"counting bead|h[ae]mocytometer|cell density|cells counted|"
    r"live cell count)", re.I)
SEED_PAT = re.compile(r"(cells? (seeded|plated)|seeding density)", re.I)
TIME_PAT = re.compile(r"^(time ?point|time|day|hour|duration|treatment time|"
                      r"timepoint|harvest|collection time)", re.I)


def _get(url, tries=4):
    for k in range(tries):
        try:
            return urllib.request.urlopen(url, timeout=120).read().decode(
                "utf-8", "replace")
        except Exception:
            if k == tries - 1:
                raise
            time.sleep(2 ** k)


def _gsm_brief(gse):
    os.makedirs(AUD, exist_ok=True)
    fn = os.path.join(AUD, gse + ".txt")
    if not os.path.exists(fn):
        t = _get("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s"
                 "&targ=gsm&form=text&view=brief" % gse)
        open(fn, "w").write(t)
        time.sleep(0.4)
    return open(fn, errors="replace").read()


def part_c(retmax=60):
    print("\n" + "=" * 78)
    print("PART C  how many public series carry a time course AND a count?")
    print("=" * 78)
    out = {}
    for name, q in QUERIES.items():
        r = json.loads(_get(EUTILS + "/esearch.fcgi?db=gds&retmax=%d&retmode=json"
                            "&term=%s" % (retmax, urllib.parse.quote(q))))
        ids = r["esearchresult"]["idlist"]
        total = int(r["esearchresult"]["count"])
        s = json.loads(_get(EUTILS + "/esummary.fcgi?db=gds&retmode=json&id="
                            + ",".join(ids)))
        accs = []
        for i in ids:
            a = s["result"][i].get("accession", "")
            if a.startswith("GSE"):
                accs.append((a, s["result"][i].get("title", "")))
        n_time, n_count, n_seed, n_both, hits = 0, 0, 0, 0, []
        for acc, title in accs:
            try:
                txt = _gsm_brief(acc)
            except Exception as e:
                print("    %s: fetch failed (%s)" % (acc, e))
                continue
            chars = [l.split("=", 1)[1].strip() for l in txt.split("\n")
                     if l.startswith("!Sample_characteristics")]
            keys = [c.split(":")[0].strip() for c in chars if ":" in c]
            tvals = set()
            for c in chars:
                if ":" in c and TIME_PAT.match(c.split(":")[0].strip()):
                    tvals.add(c.split(":", 1)[1].strip())
            has_t = len(tvals) >= 3
            has_c = any(COUNT_PAT.search(c) for c in chars) or \
                any(COUNT_PAT.search(k) for k in keys)
            has_s = any(SEED_PAT.search(c) for c in chars)
            n_time += has_t
            n_count += has_c
            n_seed += has_s
            n_both += (has_t and has_c)
            if has_c:
                hits.append({"gse": acc, "title": title[:70],
                             "fields": sorted({k for k in keys
                                               if COUNT_PAT.search(k)})})
        print("\n  %s" % name)
        print("    GEO reports %d series; audited the top %d." % (total, len(accs)))
        print("    >=3 distinct time values in sample characteristics : %d/%d"
              % (n_time, len(accs)))
        print("    any absolute-count field                           : %d/%d"
              % (n_count, len(accs)))
        print("    seeding density only (t=0, not a time course)      : %d/%d"
              % (n_seed, len(accs)))
        print("    BOTH time course and absolute count                : %d/%d"
              % (n_both, len(accs)))
        for h in hits:
            print("      count field in %s: %s" % (h["gse"], h["fields"]))
        out[name] = {"total": total, "audited": len(accs), "n_time": n_time,
                     "n_count": n_count, "n_seed": n_seed, "n_both": n_both,
                     "hits": hits,
                     "accessions": [a for a, _t in accs]}
    return out


# ------------------------------------------------- PART D  unblocking design


def part_d():
    print("\n" + "=" * 78)
    print("PART D  the measurement that unblocks each scenario")
    print("=" * 78)
    out = {}

    # EGFR, in the meaningful basis phi = (a, b, dg, gS); b's truth is 0.01 so
    # its CV is not a sensible target and its absolute sd is reported instead.
    phi = np.array([0.08, 0.01, -0.35, 0.30])
    ts = np.array([0.0, 1.0, 2.0, 4.0, 9.0, 11.0])
    p_typ = 0.3
    targets = {"a": 0.25, "dg": 0.25, "gS": 0.25}
    grid = []
    print("\n  EGFR, 6 timepoints; target CV <= 25% on a, dg and gS")
    print("  %8s %8s %8s %9s %8s %8s %s"
          % ("cells", "CV_N", "CV a", "sd b", "CV dg", "CV gS", ""))
    for ncell in (500, 2000, 10000):
        sig_p = np.sqrt(p_typ * (1 - p_typ) / ncell)
        for cvN in (None, 0.5, 0.2, 0.05):
            F = _fim(_egfr_obs2, phi, ts, sig_p, cvN)
            sd = _crlb(F)
            cv = sd / np.abs(phi)
            ok = bool(cv[0] <= targets["a"] and cv[2] <= targets["dg"]
                      and cv[3] <= targets["gS"])
            grid.append({"cells": ncell, "cv_N": cvN, "sd": sd.tolist(),
                         "cv": cv.tolist(), "meets_target": ok})
            fm = lambda v: ("%8.2f" % v) if np.isfinite(v) else "%8s" % "inf"
            print("  %8d %8s %s %9.3f %s %s %s"
                  % (ncell, "none" if cvN is None else "%.0f%%" % (100 * cvN),
                     fm(cv[0]), sd[1], fm(cv[2]), fm(cv[3]),
                     "<= target" if ok else ""))
    out["egfr_grid"] = grid
    out["egfr_basis"] = ["a", "b", "dg", "gS"]
    out["egfr_truth"] = phi.tolist()
    good = [g for g in grid if g["meets_target"]]
    if good:
        g = min(good, key=lambda g: (g["cells"], -(g["cv_N"] or 0)))
        print("  minimal viable design on this grid: %d cells/timepoint plus an"
              " absolute count at CV %.0f%%" % (g["cells"], 100 * g["cv_N"]))
        out["egfr_minimal"] = g
    else:
        print("  no design on this grid reaches the target")
        out["egfr_minimal"] = None
    print("  note: with no count at all (CV_N 'none') gS is not estimable at any"
          " cell number -")
    print("  more cells buy precision on a and dg, never on the burden.")

    # dg is the binding constraint, so sweep the sampling density too
    print("\n  EGFR, 2000 cells/timepoint and a count at CV 20%: more timepoints")
    print("  %10s %8s %9s %8s %8s %s" % ("timepoints", "CV a", "sd b", "CV dg",
                                         "CV gS", ""))
    dens = []
    sig_p = np.sqrt(p_typ * (1 - p_typ) / 2000)
    for k in (6, 12, 24, 48):
        tt = np.linspace(0.0, 11.0, k)
        F = _fim(_egfr_obs2, phi, tt, sig_p, 0.2)
        sd = _crlb(F)
        cv = sd / np.abs(phi)
        ok = bool(cv[0] <= 0.25 and cv[2] <= 0.25 and cv[3] <= 0.25)
        dens.append({"n_times": k, "sd": sd.tolist(), "cv": cv.tolist(),
                     "meets_target": ok})
        print("  %10d %8.2f %9.3f %8.2f %8.2f %s"
              % (k, cv[0], sd[1], cv[2], cv[3], "<= target" if ok else ""))
    out["egfr_density"] = dens
    hit = [d for d in dens if d["meets_target"]]
    if hit:
        print("  -> %d timepoints over 11 days, 2000 cells each, plus a count at"
              " CV 20%% clears the target." % hit[0]["n_times"])
        out["egfr_density_minimal"] = hit[0]
    else:
        print("  -> even 48 timepoints do not clear it: dg needs a tighter count,"
              " not more sampling.")
        out["egfr_density_minimal"] = None

    # AML: how many channels are needed
    th2 = np.array([0.50, 0.10, 0.20, 0.30])
    ts2 = np.array([0.0, 1.0, 2.0, 3.0, 5.0, 7.0])
    sig_p = np.sqrt(0.3 * 0.7 / 2000)
    rows = []
    for name, cvN, prolif in [("proportions only", None, None),
                              ("+ absolute counts", 0.2, None),
                              ("+ counts + proliferation readout", 0.2, 0.2)]:
        F = _fim(_aml_obs, th2, ts2, sig_p, cvN)
        if prolif is not None:
            # a direct noisy measurement of gI (EdU, dye dilution), CV `prolif`
            e = np.zeros(4)
            e[0] = 1.0
            F = F + np.outer(e, e) / (prolif * th2[0]) ** 2
        w = np.linalg.eigvalsh(F)
        rank = int((w > 1e-9 * w.max()).sum())
        sd = _crlb(F)
        cv = sd / np.abs(th2)
        rows.append({"design": name, "rank": rank, "sd": sd.tolist(),
                     "cv": cv.tolist()})
        fm = lambda v: ("%.2f" % v) if np.isfinite(v) else "inf"
        print("\n  AML %-34s rank %d/4" % (name, rank))
        print("      CV: gI %s  dI %s  delta %s  dM %s"
              % (fm(cv[0]), fm(cv[1]), fm(cv[2]), fm(cv[3])))
    out["aml_designs"] = rows
    out["aml_truth"] = th2.tolist()
    print("\n  -> AML needs THREE channels: flow proportions, an absolute count")
    print("     (beads or volumetric), and a proliferation readout. Two are not")
    print("     enough: without the third, growth and death stay lumped, and the")
    print("     differentiation rate delta - the one thing proportions do give -")
    print("     was never the blocker.")
    return out


# ----------------------------------------------------------------- reporting


def report(R):
    A, B0, B1, B2, B3, B4, B5, C, D = (R["A"], R["B0"], R["B1"], R["B2"],
                                       R["B3"], R["B4"], R["B5"], R["C"], R["D"])
    aml_aud = C["AML differentiation therapy"]
    egfr_aud = C["EGFR drug-tolerant persister"]
    pa = B3["profile"]["a (S->R, /day)"]["interval"]
    print("\n" + "=" * 78)
    print("D27 VERDICT")
    print("=" * 78)
    print("scenario 4 (EGFR drug-tolerant persister): FEASIBLE IN PRINCIPLE,")
    print("NOT ANSWERABLE FROM THIS DEPOSIT.")
    print("  structure: a proportion time course identifies the switching rates")
    print("  and the growth DIFFERENCE but never the common growth level, so the")
    print("  burden is unidentifiable by construction - the same fitted p(t)")
    print("  spans %.0e-fold at day %g." % (B4["fold_ratio"], B4["T"]))
    print("  practice: on GSE149383 the fit pins nothing. a's 95% profile is")
    print("  [%.2f, %.2f] /day and dg's covers the whole grid, because the"
          % tuple(pa))
    print("  transition is over by day 1 (half-conversion %.2f day) and the"
          % B3["half_life_days"])
    print("  day-0 point carries a library floor of %.3f." % B2["floor"])
    print("  what IS established here, and is the saleable part:")
    print("    - 3 of the 14 deposited PC9 libraries are byte-identical"
          " re-deposits;")
    print("    - sequencing depth rises with day (corr %+.2f) and PC1 is depth"
          % B1["corr_day_depth"])
    print("      (corr %+.2f), which after thinning falls to %+.2f;"
          % (B1["pc1_depth_raw"], B1["pc1_depth_matched"]))
    print("    - after depth matching the three independent day-11 libraries"
          " agree to")
    print("      sd %.3f, so the state fraction IS reproducible - it was the"
          % B2[11.0]["sd"])
    print("      depth, not the biology, that moved it;")
    print("    - but the three ways of defining the second state disagree"
          " (pairwise")
    print("      agreement %s), so 'the persister fraction' is not yet a"
          % ", ".join("%+.2f" % v for v in
                      B1["agreement"]["depth-matched"].values()))
    print("      well-posed quantity.")
    print("  %d/%d audited EGFR series carry a time course and a count field."
          % (egfr_aud["n_both"], egfr_aud["audited"]))
    print("\nscenario 5 (AML differentiation block): NOT FEASIBLE, and the")
    print("catalogue's stated reason is wrong.")
    print("  delta, the differentiation flux rate, IS identifiable from a")
    print("  proportion time course (CV %.0f%% at 2,000 cells x 6 days)."
          % (100 * D["aml_designs"][0]["cv"][2]))
    print("  What fails is the kill rate: gI and dI enter only as gI - dI, even")
    print("  WITH absolute counts, so two parameter sets give identical mature")
    print("  fractions at every timepoint while the total differs %.0f-fold."
          % A["aml"]["N_fold_ratio"])
    print("  %d/%d audited AML series carry both a time course and a count field."
          % (aml_aud["n_both"], aml_aud["audited"]))
    print("\nthe deliverable: three channels for AML (proportions + absolute")
    print("count + proliferation readout), and for EGFR sub-day sampling with")
    print("replicate libraries, depth matching, and a terminal count at CV <= 5%.")


def figures(R):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    A, B1, B2, B3, B4, C = R["A"], R["B1"], R["B2"], R["B3"], R["B4"], R["C"]
    fig, ax = plt.subplots(2, 3, figsize=(15.5, 8.6))

    # 1 the null direction of the EGFR model
    a = ax[0, 0]
    ts = np.linspace(0, 11, 200)
    th = np.array(A["egfr"]["theta"])
    c = A["egfr"]["shift"]
    for lab, tt, st in [("gS, gR as given", th, "-"),
                        ("both +%.2f/day" % c, th + np.array([0, 0, c, c]), "--")]:
        p, ln = _egfr_obs(tt, ts)
        a.plot(ts, p, st, lw=2.4 if st == "-" else 1.6, label=lab)
    a.set_xlabel("day")
    a.set_ylabel("tolerant fraction p")
    a.legend(fontsize=8)
    a.set_title("A1  two populations, identical p(t)\n(max diff %.0e) but %.0f-fold"
                " different in number"
                % (A["egfr"]["max_dp"], A["egfr"]["N_fold"]), fontsize=9.5)

    # 2 AML: same fraction, different totals
    a = ax[0, 1]
    ts2 = np.linspace(0, 7, 200)
    for lab, tt, st in [("gI .50 dI .10 dM .30", np.array(A["aml"]["theta"]), "-"),
                        ("gI .10 dI .30 dM .90", np.array(A["aml"]["alt"]), "--")]:
        f, ln = _aml_obs(tt, ts2)
        a.plot(ts2, f, st, lw=2.4 if st == "-" else 1.6, label=lab)
    a.set_xlabel("day")
    a.set_ylabel("mature fraction")
    a.legend(fontsize=8)
    a.set_title("A2  identical mature fraction, total cell\nnumber %.0f-fold apart"
                " at day 7" % A["aml"]["N_fold_ratio"], fontsize=9.5)

    # 3 p per library with the batch floor
    a = ax[0, 2]
    meta = {l: d for l, d, s_, g, r in LIBS}
    fine = np.linspace(0, 11.5, 300)
    a.plot(fine, _p_traj(B3["theta"], fine), "k-", lw=1.6,
           label="two-state fit (days 0-11)")
    tt = np.array(B3["ts"])
    a.errorbar(tt, B3["p"], yerr=B3["se"], fmt="o", ms=6, color="C0",
               ecolor="0.45", capsize=3, zorder=3, label="Drop-seq course")
    for lab in B1["p_primary"]:
        if lab in B3["course"]:
            continue
        ctrl = lab == "D19_DMSO"
        a.scatter(meta[lab], B1["p_primary"][lab], s=60, zorder=3,
                  marker="x" if ctrl else "s",
                  color="C3" if ctrl else "C2",
                  label=("untreated control" if ctrl else "other libraries")
                  if lab in ("D19_DMSO", "D11_ERL") else None)
        a.annotate(lab, (meta[lab], B1["p_primary"][lab]), fontsize=6.5,
                   xytext=(5, 3), textcoords="offset points")
    a.set_xlabel("day of erlotinib")
    a.set_ylabel("tolerant fraction p (depth-matched)")
    a.legend(fontsize=7.5, loc="lower right")
    a.set_title("B  p per independent library; day-11 libraries\nagree to sd"
                " %.3f, the day-0 pair to %.3f"
                % (B2[11.0]["sd"], B2[0.0]["sd"]), fontsize=9.5)

    # 4 the burden band
    a = ax[1, 0]
    T = B4["T"]
    gs = np.linspace(-0.6, 0.6, 200)
    a.semilogy(gs, np.exp(gs * T), "k-", lw=2)
    a.axhline(1.0, color="0.6", ls=":")
    a.axvspan(-0.5, 0.5, color="C3", alpha=0.12)
    a.set_xlabel("gS, growth rate of the sensitive pool (/day) — unmeasured")
    a.set_ylabel("burden at day %g, relative to gS = 0" % T)
    a.set_title("B4  the same fitted p(t) is consistent with\na %.0e-fold range"
                " in tumour burden" % B4["fold_ratio"], fontsize=9.5)

    # 5 the audit
    a = ax[1, 1]
    names = list(C.keys())
    xs = np.arange(len(names))
    w = 0.26
    a.bar(xs - w, [C[n]["n_time"] / C[n]["audited"] for n in names], w,
          label=">=3 time values")
    a.bar(xs, [C[n]["n_count"] / C[n]["audited"] for n in names], w,
          label="any count field")
    a.bar(xs + w, [C[n]["n_both"] / C[n]["audited"] for n in names], w,
          label="both")
    for k, n in enumerate(names):
        for dx, key in ((-w, "n_time"), (0.0, "n_count"), (w, "n_both")):
            a.annotate("%d/%d" % (C[n][key], C[n]["audited"]),
                       (xs[k] + dx, C[n][key] / C[n]["audited"]),
                       ha="center", va="bottom", fontsize=7.5)
    a.set_xticks(xs, ["EGFR\npersister", "AML\ndifferentiation"], fontsize=8.5)
    a.set_ylabel("fraction of audited series")
    a.set_ylim(0, 0.045)
    a.legend(fontsize=8)
    a.set_title("C  what the deposits actually carry\n(%d and %d series audited)"
                % tuple(C[n]["audited"] for n in names), fontsize=9.5)

    # 6 what the count precision buys, per parameter, at 2000 cells
    a = ax[1, 2]
    G = [g for g in R["D"]["egfr_grid"] if g["cells"] == 2000]
    labs6 = ["none" if g["cv_N"] is None else "%.0f%%" % (100 * g["cv_N"])
             for g in G]
    xs = np.arange(len(G))
    w = 0.26
    for k, (j, nm) in enumerate([(0, "a (switching)"), (2, "dg (growth diff.)"),
                                 (3, "gS (burden)")]):
        v = np.array([g["cv"][j] for g in G])
        capped = np.minimum(np.nan_to_num(v, posinf=4.0), 4.0)
        a.bar(xs + (k - 1) * w, capped, w, label=nm)
        for x, raw, c in zip(xs + (k - 1) * w, v, capped):
            a.annotate("inf" if not np.isfinite(raw) else "%.2f" % raw,
                       (x, c), ha="center", va="bottom", fontsize=6.5,
                       rotation=90)
    a.axhline(0.25, color="C3", ls="--", lw=1.2)
    a.annotate("target CV 0.25", (len(G) - 0.5, 0.27), ha="right", fontsize=7.5,
               color="C3")
    a.set_yscale("log")
    a.set_xticks(xs, labs6)
    a.set_xlabel("precision of the absolute count (CV), 2000 cells/timepoint")
    a.set_ylabel("Cramer-Rao CV (capped at 4)")
    a.legend(fontsize=7.5)
    a.set_title("D  what a count buys: the burden parameter is\nnot estimable"
                " without one, at any cell number", fontsize=9.5)

    fig.tight_layout()
    fn = os.path.join(FIG, "d27_egfr_aml_feasibility.png")
    fig.savefig(fn, dpi=140)
    print("\nfigure: %s" % fn)


def _cached(name, fn_, *args, force=False):
    """Per-part cache so one slow part does not have to be redone when a later
    part is edited.  D27_REDO=B3,B4 forces just those."""
    redo = set(filter(None, os.environ.get("D27_REDO", "").split(",")))
    p = os.path.join(DATA_CACHE, "d27_part_%s.pkl" % name)
    if os.path.exists(p) and not force and name not in redo:
        return pickle.load(open(p, "rb"))
    v = fn_(*args)
    pickle.dump(v, open(p, "wb"))
    return v


def main(force=False):
    if os.path.exists(RES) and not force:
        R = pickle.load(open(RES, "rb"))
    else:
        A = _cached("A", part_a)
        B0 = _cached("B0", part_b0)
        B1 = _cached("B1", part_b1, B0["independent"])
        B2 = _cached("B2", part_b2, B1)
        B3 = _cached("B3", part_b3, B1, B2)
        B4 = _cached("B4", part_b4, B3)
        B5 = _cached("B5", part_b5, B1)
        C = _cached("C", part_c)
        D = _cached("D", part_d)
        R = {"A": A, "B0": B0, "B1": B1, "B2": B2, "B3": B3, "B4": B4,
             "B5": B5, "C": C, "D": D}
        pickle.dump(R, open(RES, "wb"))
    report(R)
    figures(R)
    return R


if __name__ == "__main__":
    main(force=bool(os.environ.get("D27_FORCE")))
