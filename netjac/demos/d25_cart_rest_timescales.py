"""
D25  CAR-T exhaustion: does rest restore everything at one speed, or is there
a slow mode (the "false recovery" window)?

The hypothesis (scenario 6 of the indication catalogue). Function relaxes fast
and the epigenetic state relaxes slowly, so after a short rest a CAR-T product
can look functional while still being locked in, and collapse on re-challenge.
In the two-timescale form

    df/dt = -(f - F(e)) / tau_f        tau_f ~ hours to days   (function)
    de/dt = -(e - e_inf) / tau_e       tau_e ~ a week          (chromatin)

the decoupling window is where f has returned but e has not.

Data: Weber et al. 2021, "Transient rest induces functional reinvigoration and
epigenetic remodeling in exhausted CAR-T cells" (GSE164950):
  GSE164949  RNA-seq, 3 donors x 11 arms: HA-28z CAR-T at days 7 / 11 / 15,
             always ON (tonic signalling, exhausting), always OFF, transient
             rest started at day 7 or day 11, and always ON + nivolumab.
             Rest durations available: 4 days (twice, from different start
             days) and 8 days.
  GSE164947  ATAC-seq and H3K27me3 ChIP-seq.

What this deposit does and does not allow (checked, not assumed):
  - The RNA arm is a proper rest time course with matched always-ON and
    always-OFF controls on the same day and in the same donors. Everything
    below is measured there.
  - The chromatin arm CANNOT answer the cross-modality question. Its deposited
    ATAC matrix has only always-on / CD19-control / CD19-stimulated groups -
    no rested arm. Its H3K27me3 matrix has 18 samples labelled by construct
    (D11 / HA / SI) which cannot be mapped onto the GEO sample titles
    (rested / always-off / always-on) from the deposit alone, and guessing
    that mapping would invert the result. So the chromatin side is reported as
    a design gap, not as a number.

Questions on the RNA arm:
  Q1  Is there one exhaustion axis, and is it stable across days and donors?
  Q2  How much does rest restore, as a function of rest duration (4 d, 8 d)?
  Q3  Positive/negative controls: nivolumab (PD-1 blockade) should NOT restore
      the axis the way rest does; the housekeeping module should be undefined.
  Q4  Do modules recover at different speeds - i.e. is there a slow mode inside
      the transcriptome at all? Per-module time constants from the two rest
      durations, with the two independent 4-day arms as a reproducibility check.

Run:  python3 -m netjac.demos.d25_cart_rest_timescales
"""
import os
import pickle
import re

import numpy as np

from netjac import DATA_CACHE, FIG

RES = os.path.join(DATA_CACHE, "d25_results.pkl")
DIR = os.path.join(DATA_CACHE, "weber_rest")
RNA = os.path.join(DIR, "rna.csv.gz")
RNA_URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/GSE164949/suppl/"
           "GSE164949_Processed_log2TPM_Tuning-CAR-dataset.csv.gz")
ATAC = os.path.join(DIR, "fig_s6b_atacseq_peak_matrix.csv.gz")

# Pre-specified modules (function, not chosen by how they recover).
MOD = {
    "effector function": ["IFNG", "IL2", "TNF", "GZMB", "GZMA", "PRF1", "CSF2",
                          "IL2RA"],
    "exhaustion receptors": ["PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "CD244",
                             "ENTPD1", "CD38"],
    "exhaustion TF programme": ["TOX", "TOX2", "NR4A1", "NR4A2", "NR4A3", "IRF4",
                                "BATF", "EOMES"],
    "proliferation": ["MKI67", "TOP2A", "CCNB1", "CDK1", "BUB1", "TYMS", "RRM2",
                      "CCNA2"],
    "memory / stemness (down in exhaustion)": ["TCF7", "SELL", "CCR7", "IL7R",
                                               "LEF1", "BACH2", "KLF2"],
    "epigenetic machinery": ["EZH2", "SUZ12", "EED", "DNMT1", "DNMT3A", "DNMT3B",
                             "KDM6B", "TET2", "HDAC1", "UHRF1"],
    "housekeeping (negative control)": ["ACTB", "GAPDH", "PGK1", "TBP", "HPRT1",
                                        "PPIA", "RPLP0", "GUSB"],
}
DONORS = ["DN0022", "DN3047", "DN614"]
# arm -> (day measured, rest duration in days); ON/OFF/Nivo have no rest
ARMS = {"D7_A_ON": (7, None), "D7_A_OFF": (7, None),
        "D11_A_ON": (11, None), "D11_A_OFF": (11, None),
        "D11_A_ON_Nivo": (11, None), "D11_OFF_D7": (11, 4),
        "D15_A_ON": (15, None), "D15_A_OFF": (15, None),
        "D15_A_ON_Nivo": (15, None), "D15_OFF_D11": (15, 4), "D15_OFF_D7": (15, 8)}


# ------------------------------------------------------------------ loading
def load():
    import pandas as pd
    import subprocess
    from netjac.data.keratinocyte import _ens2sym
    if not os.path.exists(RNA):
        os.makedirs(DIR, exist_ok=True)
        subprocess.run(["curl", "-sSL", "--max-time", "900", "-o", RNA, RNA_URL],
                       check=True)
    X = pd.read_csv(RNA, index_col=0)
    sym = _ens2sym()
    ens = [i.split(".")[0] for i in X.index]
    keep = [i for i, e in enumerate(ens) if e in sym]
    X = X.iloc[keep]
    X.index = [sym[ens[i]] for i in keep]
    X = X.groupby(level=0).max()
    cols = {}
    for c in X.columns:
        m = re.match(r"(.+)_(DN\d+)$", c)
        if m and m.group(1) in ARMS:
            cols[(m.group(1), m.group(2))] = c
    return {"X": X, "cols": cols}


# ------------------------------------------------------------------ helpers
def _axis(D, day, genes_idx=None):
    """Exhaustion displacement at one day: mean over donors of ON - OFF."""
    X, C = D["X"], D["cols"]
    on = np.array([X[C[(f"D{day}_A_ON", d)]].values for d in DONORS])
    off = np.array([X[C[(f"D{day}_A_OFF", d)]].values for d in DONORS])
    return (on - off).mean(0)


def _phi(D, arm, day, idx, u):
    """Per-donor position on the day's exhaustion axis: 1 = always ON, 0 = OFF."""
    X, C = D["X"], D["cols"]
    uu = u[idx] / np.linalg.norm(u[idx])
    out = []
    for d in DONORS:
        off = X[C[(f"D{day}_A_OFF", d)]].values
        on = X[C[(f"D{day}_A_ON", d)]].values
        x = X[C[(arm, d)]].values
        den = float((on - off)[idx] @ uu)
        if abs(den) < 1e-9:
            out.append(np.nan)
        else:
            out.append(float((x - off)[idx] @ uu) / den)
    return np.array(out)


def _cos(a, b):
    a, b = a - a.mean(), b - b.mean()
    return float(a @ b / max(np.linalg.norm(a) * np.linalg.norm(b), 1e-12))


# ------------------------------------------------------------------ analysis
def run():
    from scipy import stats
    D = load()
    X = D["X"]
    genes = np.array(X.index)
    gi = {g: i for i, g in enumerate(genes)}
    out = {"n_genes": len(genes), "donors": DONORS}

    # ---- Q1 one axis?
    ax = {d: _axis(D, d) for d in (7, 11, 15)}
    out["q1_axis_cos"] = {f"day {a} vs day {b}": _cos(ax[a], ax[b])
                          for a, b in ((7, 11), (7, 15), (11, 15))}
    # per-donor axes at day 11 for a within-day scale
    per = []
    for d in DONORS:
        on = X[D["cols"][("D11_A_ON", d)]].values
        off = X[D["cols"][("D11_A_OFF", d)]].values
        per.append(on - off)
    out["q1_donor_cos_day11"] = float(np.median(
        [_cos(per[i], per[j]) for i in range(3) for j in range(i + 1, 3)]))
    out["q1_n_disrupted"] = {d: int((np.abs(ax[d]) > 0.585).sum()) for d in ax}

    # ---- module coordinates per arm
    res = {}
    for mn, gl in MOD.items():
        row = {"present": sum(g in gi for g in gl), "total": len(gl)}
        idx = np.array([gi[g] for g in gl if g in gi and abs(ax[11][gi[g]]) > 0.585],
                       dtype=int)
        row["n_disrupted"] = len(idx)
        if len(idx) < 3:
            row["defined"] = False
            res[mn] = row
            continue
        row["defined"] = True
        row["genes"] = [g for g in gl if g in gi and abs(ax[11][gi[g]]) > 0.585]
        for arm, (day, rest) in ARMS.items():
            if arm.endswith("_A_ON") or arm.endswith("_A_OFF"):
                continue
            row[arm] = {"phi": _phi(D, arm, day, idx, ax[day]).tolist(),
                        "rest_days": rest}
        # sanity: the ON arm must sit at 1 and OFF at 0 by construction
        row["check_on_day11"] = _phi(D, "D11_A_ON", 11, idx, ax[11]).tolist()
        res[mn] = row
    out["modules"] = res

    # ---- Q2 restoration vs rest duration (whole-axis, all disrupted genes)
    idx_all = np.where(np.abs(ax[11]) > 0.585)[0]
    whole = {}
    for arm, (day, rest) in ARMS.items():
        if rest is None and "Nivo" not in arm:
            continue
        p = _phi(D, arm, day, idx_all, ax[day])
        whole[arm] = {"phi": p.tolist(), "median": float(np.nanmedian(p)),
                      "restoration": float(1 - np.nanmedian(p)), "rest_days": rest}
    out["q2_whole_axis"] = whole

    # ---- Q3 nivolumab vs rest at the same day, paired by donor
    q3 = {}
    for day in (11, 15):
        rest_arm = "D11_OFF_D7" if day == 11 else "D15_OFF_D7"
        nv = np.array(whole[f"D{day}_A_ON_Nivo"]["phi"])
        rs = np.array(whole[rest_arm]["phi"])
        q3[f"day {day}"] = {"nivolumab_phi": nv.tolist(), "rest_phi": rs.tolist(),
                            "paired_diff": float(np.mean(nv - rs)),
                            "p": float(stats.ttest_rel(nv, rs).pvalue)}
    out["q3_nivolumab"] = q3

    # ---- Q4 per-module time constants from the two rest durations
    def tau(phi, t):
        phi = min(max(phi, 1e-3), 0.999)
        return float(-t / np.log(phi))
    q4 = {}
    for mn, row in res.items():
        if not row.get("defined"):
            continue
        e = {}
        for arm in ("D11_OFF_D7", "D15_OFF_D11", "D15_OFF_D7"):
            p = float(np.nanmedian(row[arm]["phi"]))
            e[arm] = {"phi": p, "restoration": 1 - p,
                      "tau_days": tau(p, row[arm]["rest_days"])}
        e["tau_4d_reproducibility"] = (e["D11_OFF_D7"]["tau_days"],
                                       e["D15_OFF_D11"]["tau_days"])
        q4[mn] = e
    out["q4_modules"] = q4
    taus = [v["D15_OFF_D7"]["tau_days"] for v in q4.values()]
    if taus:
        out["q4_tau_spread"] = {"min": float(np.min(taus)), "max": float(np.max(taus)),
                                "ratio": float(np.max(taus) / max(np.min(taus), 1e-9))}

    # ---- chromatin: what the deposit contains
    if os.path.exists(ATAC):
        import pandas as pd
        A = pd.read_csv(ATAC, index_col=0)
        g = {"always_on": [0, 1, 2], "cd19_ctr": [3, 4, 5], "cd19_stim": [6, 7, 8]}
        L = np.log1p(A.values / A.values.sum(0, keepdims=True) * 1e6)
        d_exh = L[:, g["always_on"]].mean(1) - L[:, g["cd19_ctr"]].mean(1)
        d_stim = L[:, g["cd19_stim"]].mean(1) - L[:, g["cd19_ctr"]].mean(1)
        out["chromatin"] = {
            "n_peaks": int(A.shape[0]), "groups": {k: len(v) for k, v in g.items()},
            "sd_exhaustion_axis": float(np.std(d_exh)),
            "n_peaks_moved_1log2": int((np.abs(d_exh) > 1).sum()),
            "cos_exhaustion_vs_stim": _cos(d_exh, d_stim),
            "note": ("ATAC deposit has no rested arm; the H3K27me3 matrix is "
                     "labelled by construct (D11/HA/SI) and cannot be mapped to "
                     "rested/always-off/always-on from the deposit alone")}
    return out


# ------------------------------------------------------------------ report
def report(R):
    print("=" * 78)
    print("D25  CAR-T rest: one timescale or two? (Weber 2021, GSE164949/47)")
    print("=" * 78)
    print(f"  {R['n_genes']} genes with symbols, donors {', '.join(R['donors'])}")

    print("\n  Q1  is there one exhaustion axis?")
    for k, v in R["q1_axis_cos"].items():
        print(f"      {k:22s} cos {v:+.3f}")
    print(f"      between donors within day 11: cos {R['q1_donor_cos_day11']:+.3f}"
          "   (the scale for 'the same axis')")
    print("      disrupted genes (|log2 ON-OFF| > 0.585): " +
          ", ".join(f"day {d}: {n}" for d, n in R["q1_n_disrupted"].items()))

    print("\n  Q2  restoration of the whole axis (1 = back to always-OFF)")
    for arm, v in sorted(R["q2_whole_axis"].items(),
                         key=lambda kv: (kv[1]["rest_days"] or 0)):
        tag = f"{v['rest_days']} d rest" if v["rest_days"] else "no rest"
        print(f"      {arm:18s} {tag:10s} phi {v['median']:+.3f}  "
              f"restoration {v['restoration']:+.3f}")

    print("\n  Q3  rest vs PD-1 blockade, paired by donor")
    for k, v in R["q3_nivolumab"].items():
        print(f"      {k}: nivolumab phi {np.median(v['nivolumab_phi']):+.3f} vs "
              f"rest phi {np.median(v['rest_phi']):+.3f}   "
              f"difference {v['paired_diff']:+.3f} (p = {v['p']:.4f})")

    print("\n  Q4  per-module restoration and implied time constants")
    print(f"      {'module':40s} {'4d(a)':>7s} {'4d(b)':>7s} {'8d':>7s} "
          f"{'tau 8d':>8s} {'tau 4d a/b':>14s}")
    for mn, v in sorted(R["q4_modules"].items(),
                        key=lambda kv: -kv[1]["D15_OFF_D7"]["tau_days"]):
        a, b = v["tau_4d_reproducibility"]
        print(f"      {mn:40s} {v['D11_OFF_D7']['restoration']:+7.2f} "
              f"{v['D15_OFF_D11']['restoration']:+7.2f} "
              f"{v['D15_OFF_D7']['restoration']:+7.2f} "
              f"{v['D15_OFF_D7']['tau_days']:8.1f} {a:6.1f} /{b:6.1f}")
    for mn, row in R["modules"].items():
        if not row.get("defined"):
            print(f"      {mn:40s} undefined ({row['n_disrupted']} disrupted of "
                  f"{row['present']}/{row['total']} present)")
    if "q4_tau_spread" in R:
        s = R["q4_tau_spread"]
        print(f"      spread of module time constants at 8 d: {s['min']:.1f} to "
              f"{s['max']:.1f} days (ratio {s['ratio']:.2f})")

    c = R.get("chromatin")
    if c:
        print(f"\n  chromatin deposit: {c['n_peaks']} peaks, groups "
              + ", ".join(f"{k}={v}" for k, v in c["groups"].items()))
        print(f"      exhaustion axis exists in chromatin: "
              f"{c['n_peaks_moved_1log2']} peaks move > 1 log2, "
              f"cos(exhaustion, stimulation) {c['cos_exhaustion_vs_stim']:+.3f}")
        print(f"      LIMIT: {c['note']}")


def figures(R):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    fig, ax = plt.subplots(1, 3, figsize=(17, 5.2))

    a = ax[0]
    for k, (mn, v) in enumerate(R["q4_modules"].items()):
        p4a, p4b = v["D11_OFF_D7"]["phi"], v["D15_OFF_D11"]["phi"]
        p8 = v["D15_OFF_D7"]["phi"]
        c = f"C{k}"
        a.plot([0, 4, 8], [1.0, np.mean([p4a, p4b]), p8], "-", lw=1.8, color=c,
               label=mn.split(" (")[0][:26], alpha=0.9)
        a.plot([4, 4], [p4a, p4b], "o", ms=5, color=c, alpha=0.9)
        a.plot([0, 8], [1.0, p8], "o", ms=5, color=c, alpha=0.9)
    nv = np.median(R["q3_nivolumab"]["day 15"]["nivolumab_phi"])
    a.plot([8], [nv], "D", ms=8, color="k", label="nivolumab, 8 d (no rest)")
    a.axhline(1, color="k", ls=":", lw=1)
    a.axhline(0, color="k", lw=0.8)
    a.set_xticks([0, 4, 8])
    a.set_xlabel("days of rest  (0 = always ON)")
    a.set_ylabel("position on the exhaustion axis (1 = always ON)")
    a.set_title("Q2/Q4  every module collapses onto one fast\ntimescale; the two "
                "4-day arms agree", fontsize=10)
    a.legend(fontsize=7, loc="upper right")

    a = ax[1]
    ks = list(R["q4_modules"])
    tau = [R["q4_modules"][k]["D15_OFF_D7"]["tau_days"] for k in ks]
    o = np.argsort(tau)
    a.barh(range(len(ks)), [tau[i] for i in o],
           color=["0.6" if ks[i].startswith("housekeeping") else "C0" for i in o])
    a.set_yticks(range(len(ks)), [ks[i].split(" (")[0][:28] for i in o], fontsize=8)
    a.set_xlabel("implied time constant from the 8-day arm (days)")
    a.set_title("Q4  no slow mode inside the transcriptome", fontsize=10)

    a = ax[2]
    lbl, val = [], []
    for day in (11, 15):
        v = R["q3_nivolumab"][f"day {day}"]
        lbl += [f"rest\nday {day}", f"nivolumab\nday {day}"]
        val += [np.median(v["rest_phi"]), np.median(v["nivolumab_phi"])]
    a.bar(range(len(val)), val, 0.6,
          color=["C2", "C3", "C2", "C3"])
    a.set_xticks(range(len(lbl)), lbl, fontsize=8.5)
    a.axhline(1, color="k", ls=":", lw=1)
    a.axhline(0, color="k", lw=0.8)
    a.set_ylabel("position on the exhaustion axis")
    a.set_title("Q3  specificity control:\nPD-1 blockade does not substitute for rest",
                fontsize=10)

    fig.tight_layout()
    fn = os.path.join(FIG, "d25_cart_rest_timescales.png")
    fig.savefig(fn, dpi=140)
    print(f"\nfigure: {fn}")


def main(force=False):
    if os.path.exists(RES) and not force:
        R = pickle.load(open(RES, "rb"))
    else:
        R = run()
        pickle.dump(R, open(RES, "wb"))
    report(R)
    figures(R)
    return R


if __name__ == "__main__":
    main(force=bool(os.environ.get("D25_FORCE")))
