"""M24 - measure the RNA -> secreted-cytokine link instead of assuming it.

Any "peptide -> pathway -> RNA -> cytokine" model contains one step that is
almost always assumed rather than measured: that a cell's cytokine mRNA tells
you how much cytokine that cell secretes. TRAPS-seq (Wu et al., Nat Methods
2023; GEO GSE200690) is the dataset that can test it, because secreted IFN-g,
IL-2 and TNF are captured on the surface of the SAME cell whose transcriptome
is sequenced.

Two samples are used:
  EndTimePoint  4,894 cells, one capture window, 4 hashtags
  Tracing       3,223 cells, an unstimulated Ctrl arm and a 10 h stimulated arm,
                and TWO successive capture windows per cytokine (_1st, _2nd)

The Tracing design is what makes this more than a correlation: the Ctrl arm
gives the assay's own background, and the two windows let us ask whether
secretion is driven by current transcription or is a persistent cell property.

HONEST FINDING (see REPORT/VIRTUAL_TCELL_REPORT.md for the full argument):
  * mRNA and secretion are coupled, but weakly and unequally by cytokine.
    In EndTimePoint, own-mRNA explains R^2 = 0.25 (IFN-g), 0.12 (TNF) and
    0.02 (IL-2) of the cell-to-cell spread in secretion.
  * Secretion is far more predictable from the SAME cell's earlier secretion
    (rho 0.55-0.80) than from its mRNA, and adding mRNA on top of the earlier
    window buys +0.002 to +0.033 R^2.
  * The trend is nonetheless monotone: top-quartile mRNA cells secrete ~3-4x
    the median of mRNA-negative cells.
So the RNA -> cytokine arrow supports a RANKING, not a per-cell quantity.

Needs network on first run (~86 MB from NCBI FTP), cached under figures/_traps_cache.
Run:  python3 -m grn_pipeline.m24_rna_secretion_coupling
"""

from __future__ import annotations

import gzip
import os
import subprocess
import tarfile

import numpy as np
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "figures")
CACHE = os.path.join(FIGDIR, "_traps_cache")
URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE200nnn/GSE200690/suppl/"
       "GSE200690_RAW.tar")

END = {"rna": "GSM6042322_EndTimePoint_ProcessedData_EndTimePoint",
       "hto": "GSM6042323_EndTimePoint_ProcessedData_EndTimePoint",
       "srf": "GSM6042324_EndTimePoint_ProcessedData_EndTimePoint",
       "cyt": "GSM6042325_EndTimePoint_ProcessedData_EndTimePoint"}
TRACE = {"rna": "GSM6042326_Tracing_ProcessedData_Tracing",
         "hto": "GSM6042327_Tracing_ProcessedData_Tracing",
         "srf": "GSM6042328_Tracing_ProcessedData_Tracing",
         "cyt": "GSM6042329_Tracing_ProcessedData_Tracing"}
PAIRS_END = [("IL2", "IL-2"), ("IFNG", "IFN-g"), ("TNF", "TNF-a")]
PAIRS_TR = [("IL2", "IL-2"), ("IFNG", "IFN-y"), ("TNF", "TNF-a")]   # spelling differs


# ------------------------------------------------------------------ data access
def ensure_data() -> str:
    os.makedirs(CACHE, exist_ok=True)
    probe = os.path.join(CACHE, END["cyt"] + "_Cytokine_features.tsv.gz")
    if os.path.exists(probe):
        return CACHE
    tar = os.path.join(CACHE, "GSE200690_RAW.tar")
    if not os.path.exists(tar):
        subprocess.run(["curl", "-sL", "--max-time", "900", URL, "-o", tar], check=True)
    with tarfile.open(tar) as t:
        t.extractall(CACHE)
    os.remove(tar)
    return CACHE


def _lines(p):
    with gzip.open(p, "rt") as fh:
        return [l.rstrip("\n") for l in fh]


def _barcodes(p):
    # the RNA matrix carries a '-1' lane suffix that the antibody matrices lack
    return [b.split("-")[0] for b in _lines(p)]


def _dense(prefix, tag):
    f = _lines(f"{CACHE}/{prefix}_{tag}_features.tsv.gz")
    b = _barcodes(f"{CACHE}/{prefix}_{tag}_barcodes.tsv.gz")
    M = np.zeros((len(f), len(b)))
    with gzip.open(f"{CACHE}/{prefix}_{tag}_matrix.mtx.gz", "rt") as fh:
        for line in fh:
            if not line.startswith("%"):
                break
        for line in fh:
            r, c, v = line.split()
            M[int(r) - 1, int(c) - 1] = float(v)
    return f, b, M


def _rna(prefix, want):
    feats = _lines(f"{CACHE}/{prefix}_RNA_features.tsv.gz")
    b = _barcodes(f"{CACHE}/{prefix}_RNA_barcodes.tsv.gz")
    row2sym = {}
    for i, l in enumerate(feats):
        p = l.split("\t")
        s = p[1] if len(p) > 1 else p[0]
        if s in want:
            row2sym[i + 1] = s
    tot = np.zeros(len(b))
    out = {s: np.zeros(len(b)) for s in want}
    with gzip.open(f"{CACHE}/{prefix}_RNA_matrix.mtx.gz", "rt") as fh:
        for line in fh:
            if not line.startswith("%"):
                break
        for line in fh:
            r, c, v = line.split()
            r = int(r); c = int(c) - 1; v = float(v)
            tot[c] += v
            s = row2sym.get(r)
            if s is not None:
                out[s][c] += v
    return b, tot, out


def _align(P, pairs):
    """Returns per-cell arrays on the barcodes shared by all four matrices.

    Secretion is normalised by the cell's TOTAL SURFACE-ANTIBODY counts, not by
    its total cytokine counts: with only three cytokine features the latter makes
    them compositional and manufactures negative correlations between them.
    """
    want = {s for s, _ in pairs}
    rbc, rtot, rmat = _rna(P["rna"], want)
    cf, cbc, C = _dense(P["cyt"], "Cytokine")
    sf, sbc, S = _dense(P["srf"], "SurfaceMarker")
    hf, hbc, H = _dense(P["hto"], "HTO")
    shared = set(cbc) & set(sbc) & set(hbc)
    common = [b for b in rbc if b in shared]
    pos = {b: i for i, b in enumerate(rbc)}
    ri = np.array([pos[b] for b in common])
    ci = np.array([{b: i for i, b in enumerate(cbc)}[b] for b in common])
    si = np.array([{b: i for i, b in enumerate(sbc)}[b] for b in common])
    hi = np.array([{b: i for i, b in enumerate(hbc)}[b] for b in common])
    srf = S[:-1, si].sum(axis=0)
    rtt = rtot[ri]
    ok = (srf > 0) & (rtt > 0)
    return dict(ri=ri[ok], ci=ci[ok], si=si[ok], hi=hi[ok], srf=srf[ok], rtt=rtt[ok],
                C=C, S=S, H=H, cf=cf, sf=sf, hf=hf, rmat=rmat, n=int(ok.sum()))


# --------------------------------------------------------------------- analyses
def end_timepoint(A):
    cidx = {f.split("-")[0] + "-" + f.split("-")[1]: i for i, f in enumerate(A["cf"][:-1])}
    out = {}
    print("\n[EndTimePoint] own-mRNA vs own secreted protein, same cell")
    for sym, ab in PAIRS_END:
        rna = np.log1p(A["rmat"][sym][A["ri"]] / A["rtt"] * 1e4)
        prot = np.log1p(A["C"][cidx[ab], A["ci"]] / A["srf"] * 1e3)
        rho, p = stats.spearmanr(rna, prot)
        r, _ = stats.pearsonr(rna, prot)
        out[sym] = dict(rho=float(rho), r2=float(r * r), p=float(p))
        print(f"   {sym:5s} -> {ab:6s}  rho {rho:+.3f}  R2 {r*r:.3f}  (p={p:.1g})"
              f"   [{100*(1-r*r):.0f}% of the spread unexplained]")

    print("\n[EndTimePoint] specificity (rho, mRNA row -> secreted column)")
    print("            " + "".join(f"{ab:>9s}" for _, ab in PAIRS_END))
    for sx, _ in PAIRS_END:
        rna = np.log1p(A["rmat"][sx][A["ri"]] / A["rtt"] * 1e4)
        cells = []
        for _, ab in PAIRS_END:
            prot = np.log1p(A["C"][cidx[ab], A["ci"]] / A["srf"] * 1e3)
            cells.append(f"{stats.spearmanr(rna, prot)[0]:+9.3f}")
        print(f"   {sx:8s} " + "".join(cells))

    k69 = [i for i, f in enumerate(A["sf"]) if f.startswith("CD69")][0]
    cd69 = np.log1p(A["S"][k69, A["si"]] / A["srf"] * 1e3)
    Z = np.column_stack([cd69, np.log(A["rtt"]), np.log(A["srf"]), np.ones(A["n"])])
    print("\n[EndTimePoint] partial rho, controlling CD69 protein and both depths")
    for sym, ab in PAIRS_END:
        x = np.log1p(A["rmat"][sym][A["ri"]] / A["rtt"] * 1e4)
        y = np.log1p(A["C"][cidx[ab], A["ci"]] / A["srf"] * 1e3)
        rx = x - Z @ np.linalg.lstsq(Z, x, rcond=None)[0]
        ry = y - Z @ np.linalg.lstsq(Z, y, rcond=None)[0]
        rho, p = stats.spearmanr(rx, ry)
        out[sym]["partial_rho"] = float(rho)
        print(f"   {sym:5s} -> {ab:6s}  partial rho {rho:+.3f} (raw {out[sym]['rho']:+.3f})")
    return out


def tracing(A):
    cidx = {f.rsplit("-", 1)[0]: i for i, f in enumerate(A["cf"][:-1])}
    hs = A["H"][:-1, A["hi"]]
    lab = np.array([f.split("-")[0] for f in A["hf"][:-1]])[hs.argmax(axis=0)]
    pur = hs.max(axis=0) / np.maximum(hs.sum(axis=0), 1)
    pure = pur > 0.6
    ctrl = pure & np.char.endswith(lab.astype(str), "Ctrl")
    stim = pure & np.char.endswith(lab.astype(str), "10hr")
    print(f"\n[Tracing] hashtag arms: Ctrl {ctrl.sum()} cells, 10hr {stim.sum()} cells")

    def P(name, m):
        return A["C"][cidx[name], A["ci"]][m] / A["srf"][m] * 1e3

    def R(sym, m):
        return A["rmat"][sym][A["ri"]][m] / A["rtt"][m] * 1e4

    print("[Tracing] stimulation effect (median, counts per 1,000 antibody counts)")
    out = {"stim_effect": {}}
    for sym, ab in PAIRS_TR:
        for w in ("1st", "2nd"):
            k = f"{ab}_{w}"
            m0, m1 = np.median(P(k, ctrl)), np.median(P(k, stim))
            pv = stats.mannwhitneyu(P(k, stim), P(k, ctrl), alternative="greater").pvalue
            out["stim_effect"][k] = (float(m0), float(m1))
            print(f"   {k:10s} Ctrl {m0:7.2f}  10hr {m1:7.2f}  ratio {m1/max(m0,1e-9):5.2f}"
                  f"  p={pv:.1g}")
    for sym, _ in PAIRS_TR:
        print(f"   {sym+' mRNA':10s} detected in {100*(R(sym,ctrl)>0).mean():5.1f}% of Ctrl"
              f" vs {100*(R(sym,stim)>0).mean():5.1f}% of 10hr cells")

    print("\n[Tracing] mRNA at harvest vs secretion in each window (10hr cells)")
    out["windows"] = {}
    for sym, ab in PAIRS_TR:
        x = np.log1p(R(sym, stim))
        vals = []
        for w in ("1st", "2nd"):
            y = np.log1p(P(f"{ab}_{w}", stim))
            rho, _ = stats.spearmanr(x, y)
            vals.append(float(rho))
        out["windows"][sym] = vals
        print(f"   {sym:5s}  earlier window rho {vals[0]:+.3f}   later window rho {vals[1]:+.3f}")

    print("\n[Tracing] secretion as a persistent cell property, and what mRNA adds")
    out["persistence"] = {}
    for sym, ab in PAIRS_TR:
        a = np.log1p(P(f"{ab}_1st", stim))
        b = np.log1p(P(f"{ab}_2nd", stim))
        rho, _ = stats.spearmanr(a, b)
        m = np.log1p(R(sym, stim))
        y = b

        def r2(X):
            beta = np.linalg.lstsq(X, y, rcond=None)[0]
            return 1 - (y - X @ beta).var() / y.var()

        one = r2(np.column_stack([a, np.ones_like(a)]))
        two = r2(np.column_stack([a, m, np.ones_like(a)]))
        out["persistence"][sym] = dict(rho=float(rho), r2_prev=float(one),
                                       r2_prev_plus_rna=float(two))
        print(f"   {ab:7s} earlier-vs-later rho {rho:+.3f} | R2 from earlier "
              f"{one:.3f} -> +own mRNA {two:.3f} (gain {two-one:+.3f})")

    print("\n[Tracing] later-window secretion by mRNA level (median, 10hr cells)")
    out["bins"] = {}
    for sym, ab in PAIRS_TR:
        x = R(sym, stim)
        y = P(f"{ab}_2nd", stim)
        zero = x == 0
        cells = [("mRNA=0", float(np.median(y[zero])), int(zero.sum()))]
        pos = x[~zero]
        q = np.quantile(pos, [0.25, 0.5, 0.75])
        masks = [(~zero) & (x <= q[0]), (~zero) & (x > q[0]) & (x <= q[1]),
                 (~zero) & (x > q[1]) & (x <= q[2]), (~zero) & (x > q[2])]
        for i, mm in enumerate(masks, 1):
            cells.append((f"Q{i}", float(np.median(y[mm])), int(mm.sum())))
        out["bins"][sym] = cells
        fold = cells[-1][1] / max(cells[0][1], 1e-9)
        print(f"   {ab:7s} " + "  ".join(f"{n}:{v:7.1f}(n={k})" for n, v, k in cells)
              + f"   top/zero = {fold:.1f}x")
    return out


def report() -> dict:
    print("=" * 78)
    print("M24  does cytokine mRNA predict what the same cell secretes?")
    print("     TRAPS-seq, GEO GSE200690 (Wu et al., Nat Methods 2023)")
    print("=" * 78)
    ensure_data()
    A_end = _align(END, PAIRS_END)
    print(f"\nEndTimePoint: {A_end['n']} cells aligned across RNA / secretion / "
          f"surface / hashtag")
    res_end = end_timepoint(A_end)
    A_tr = _align(TRACE, PAIRS_TR)
    print(f"\nTracing: {A_tr['n']} cells aligned")
    res_tr = tracing(A_tr)

    print("\n" + "-" * 78)
    print("WHAT THIS MEANS FOR A peptide -> pathway -> RNA -> cytokine MODEL")
    print("  The RNA -> secretion arrow is real and monotone, so the chain can")
    print("  RANK conditions. Per cell it is weak (R2 0.02-0.25) and secretion is")
    print("  better predicted by the same cell's earlier secretion than by its")
    print("  mRNA, so the chain cannot output a per-cell quantity.")
    print("-" * 78)
    return {"end": res_end, "tracing": res_tr, "n_end": A_end["n"], "n_tracing": A_tr["n"]}


main = report

if __name__ == "__main__":
    report()
