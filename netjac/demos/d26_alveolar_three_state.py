"""
D26  Alveolar regeneration stalls: can a composition snapshot prove that the
EXIT is blocked rather than the entry accelerated?

The hypothesis (scenario 10 of the indication catalogue). Repair runs
AT2 -> transitional (Krt8+ ADI / "DATP") -> AT1. In fibrosis the transitional
state piles up. Two stories fit the pile-up: the entry was accelerated
(more cells pushed in) or the exit was blocked (cells cannot leave). The
claim is that a composition snapshot can rule out the first.

    dA/dt = s  - k1 A        A = AT2
    dD/dt = k1 A - k2 D      D = transitional
    dT/dt = k2 D - d  T      T = AT1

PART A  the algebra, then a numerical check of it (no data needed):
  A1  at steady state A : D : T = 1/k1 : 1/k2 : 1/d, so a snapshot identifies
      only RATIOS of residence times - never absolute days.
  A2  D/T = d/k2 exactly, independent of k1 and of s. No acceleration of entry
      can move the transitional-to-AT1 ratio at steady state. This is the sharp
      version of the scenario's claim.
  A3  the fraction ceiling: with k2 and d held at their healthy values, driving
      k1 -> infinity gives f_D -> f_D/(f_D + f_T) evaluated on the healthy
      composition. (The catalogue's "ceiling 0.200" has no source; the ceiling
      is whatever the healthy composition says it is, and it is computed here.)

PART B  human IPF (Habermann 2020, GSE135893): AT2 / Transitional AT2 / AT1
  counts per donor, bootstrapped over donors.
PART C  mouse bleomycin time course (Strunz 2020, GSE141259, epithelium-
  enriched): the same ratio day by day from day 2 to day 54.
PART D  the measurement check that decides whether any of this is usable:
  AT1 cells are large and fragile and are lost in dissociation. If their
  recovery differs between conditions, the ratio test is not interpretable.
  Quantified here, with a sensitivity sweep over AT1 recovery.

Run:  python3 -m netjac.demos.d26_alveolar_three_state
"""
import os
import pickle

import numpy as np

from netjac import DATA_CACHE, FIG

RES = os.path.join(DATA_CACHE, "d26_results.pkl")
DIR = os.path.join(DATA_CACHE, "alveolar")
HUMAN = os.path.join(DIR, "habermann_meta.csv.gz")
MOUSE = os.path.join(DIR, "strunz_hires_cellinfo.csv.gz")
MOUSE_WHOLE = os.path.join(DIR, "strunz_cellinfo.csv.gz")
URLS = {
    HUMAN: ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135893/suppl/"
            "GSE135893_IPF_metadata.csv.gz"),
    MOUSE: ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141259/suppl/"
            "GSE141259_HighResolution_cellinfo.csv.gz"),
    MOUSE_WHOLE: ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141259/suppl/"
                  "GSE141259_WholeLung_cellinfo.csv.gz"),
}
# compartment -> annotation labels in each dataset
HUMAN_MAP = {"AT2": ["AT2", "Transitional AT2"], "D": ["Transitional AT2"],
             "T": ["AT1"]}
MOUSE_MAP = {"A": ["AT2", "AT2 activated"], "D": ["Krt8+ ADI"], "T": ["AT1"]}


def _fetch(path):
    if not os.path.exists(path):
        os.makedirs(DIR, exist_ok=True)
        import subprocess
        subprocess.run(["curl", "-sSL", "--max-time", "900", "-o", path, URLS[path]],
                       check=True)
    return path


# ============================================================ part A: algebra
def part_a(seed=0):
    """Verify the three statements numerically; they are what the data is then
    tested against."""
    from scipy.linalg import expm
    rng = np.random.default_rng(seed)
    out = {}

    def steady(s, k1, k2, d):
        """Steady state solved from the linear system M y + b = 0, i.e. without
        using the hand-derived formula that is being checked."""
        M = np.array([[-k1, 0.0, 0.0], [k1, -k2, 0.0], [0.0, k2, -d]])
        b = np.array([s, 0.0, 0.0])
        return np.linalg.solve(M, -b)

    def trajectory_reaches(s, k1, k2, d, t):
        """Independent check that the solved point is actually approached:
        y(t) = expm(Mt) y0 + (I - expm(Mt)) y*, started from zero."""
        M = np.array([[-k1, 0.0, 0.0], [k1, -k2, 0.0], [0.0, k2, -d]])
        y_ss = steady(s, k1, k2, d)
        E = expm(M * t)
        return (np.eye(3) - E) @ y_ss

    # A1 fractions are proportional to residence times
    err = []
    for _ in range(200):
        s, k1, k2, d = rng.uniform(0.2, 5, 4)
        y = steady(s, k1, k2, d)
        pred = np.array([1 / k1, 1 / k2, 1 / d]) * s
        err.append(float(np.max(np.abs(y - pred) / pred)))
    out["a1_max_rel_error"] = float(np.max(err))
    # and the steady state is actually reached from zero
    s0, k1, k2, d = 1.0, 0.5, 0.7, 0.4
    y_ss = steady(s0, k1, k2, d)
    y_t = trajectory_reaches(s0, k1, k2, d, 200.0)
    out["a1_trajectory_rel_gap"] = float(np.max(np.abs(y_t - y_ss) / y_ss))

    # A2 D/T is independent of k1 (and of s)
    k2, d = 0.7, 0.4
    ratios = []
    for k1 in np.logspace(-2, 3, 12):
        for s in (0.3, 3.0, 30.0):
            y = steady(s, k1, k2, d)
            ratios.append(y[1] / y[2])
    out["a2_ratio_mean"] = float(np.mean(ratios))
    out["a2_ratio_spread"] = float(np.max(ratios) - np.min(ratios))
    out["a2_predicted_d_over_k2"] = float(d / k2)

    # A3 the fraction ceiling as k1 -> infinity, from a healthy composition
    def ceiling_from_composition(fA, fD, fT):
        return float(fD / (fD + fT))
    s, k1h, k2h, dh = 1.0, 0.5, 0.7, 0.4
    yh = steady(s, k1h, k2h, dh)
    fh = yh / yh.sum()
    ceil_pred = ceiling_from_composition(*fh)
    big = steady(s, 1e9, k2h, dh)
    out["a3"] = {"healthy_fractions": fh.tolist(), "ceiling_predicted": ceil_pred,
                 "ceiling_simulated": float(big[1] / big.sum()),
                 "note": "k2 and d held at their healthy values"}
    return out


# ============================================================ part B: human
def part_b(n_boot=4000, seed=0):
    import pandas as pd
    H = pd.read_csv(_fetch(HUMAN), index_col=0)
    H = H[H.population == "Epithelial"]
    rows = []
    for (donor, status), g in H.groupby(["Sample_Name", "Status"], observed=True):
        c = g.celltype.value_counts()
        a = int(c.get("AT2", 0))
        dd = int(c.get("Transitional AT2", 0))
        t = int(c.get("AT1", 0))
        rows.append({"donor": donor, "status": status, "A": a, "D": dd, "T": t,
                     "KRT5_KRT17": int(c.get("KRT5-/KRT17+", 0))})
    R = pd.DataFrame(rows)
    out = {"per_donor": R.to_dict("records")}
    rng = np.random.default_rng(seed)

    def stats(sub):
        A, D, T = sub.A.sum(), sub.D.sum(), sub["T"].sum()
        tot = A + D + T
        return {"A": int(A), "D": int(D), "T": int(T),
                "f_D": float(D / tot) if tot else np.nan,
                "f_T": float(T / tot) if tot else np.nan,
                "D_over_T": float(D / T) if T else np.nan,
                "n_donors": len(sub)}

    for st in ("Control", "Disease"):
        sub = R[R.status == st]
        s = stats(sub)
        # bootstrap over DONORS, not cells
        bs = [stats(sub.iloc[rng.integers(0, len(sub), len(sub))])["D_over_T"]
              for _ in range(n_boot)]
        bs = np.array([b for b in bs if b == b and np.isfinite(b)])
        s["D_over_T_ci"] = [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]
        out[st] = s
    # the ceiling implied by the controls, and whether disease exceeds it
    c = out["Control"]
    out["ceiling_from_controls"] = float(c["f_D"] / (c["f_D"] + c["f_T"]))
    out["disease_f_D"] = out["Disease"]["f_D"]
    out["fraction_ceiling_exceeded"] = bool(out["disease_f_D"] >
                                            out["ceiling_from_controls"])
    # ratio fold change with a donor bootstrap
    ctrl, dis = R[R.status == "Control"], R[R.status == "Disease"]
    fc = []
    for _ in range(n_boot):
        a = stats(ctrl.iloc[rng.integers(0, len(ctrl), len(ctrl))])["D_over_T"]
        b = stats(dis.iloc[rng.integers(0, len(dis), len(dis))])["D_over_T"]
        if a and b and np.isfinite(a) and np.isfinite(b):
            fc.append(b / a)
    out["ratio_fold_change"] = {"obs": out["Disease"]["D_over_T"] / out["Control"]["D_over_T"],
                                "ci": [float(np.percentile(fc, 2.5)),
                                       float(np.percentile(fc, 97.5))]}
    return out


# ============================================================ part C: mouse
def part_c():
    import pandas as pd
    M = pd.read_csv(_fetch(MOUSE), sep="\t", index_col=0)
    rows = []
    for tp, g in M.groupby("time_point", observed=True):
        c = g.cell_type.value_counts()
        a = int(sum(c.get(x, 0) for x in MOUSE_MAP["A"]))
        d = int(sum(c.get(x, 0) for x in MOUSE_MAP["D"]))
        t = int(sum(c.get(x, 0) for x in MOUSE_MAP["T"]))
        day = 0 if "PBS" in str(tp) else int(str(tp).split()[-1])
        rows.append({"time_point": str(tp), "day": day, "A": a, "D": d, "T": t,
                     "n_cells": len(g),
                     "f_D": d / max(a + d + t, 1), "f_T": t / max(a + d + t, 1),
                     "D_over_T": (d / t) if t else np.nan})
    R = pd.DataFrame(rows).sort_values("day")
    out = {"time_course": R.to_dict("records")}
    pbs = R[R.day == 0].iloc[0]
    out["pbs"] = pbs.to_dict()
    out["ceiling_from_pbs"] = float(pbs.f_D / (pbs.f_D + pbs.f_T)) if (pbs.f_D + pbs.f_T) else np.nan
    late = R[(R.day >= 21)]
    out["late"] = {"days": late.day.tolist(), "f_D": late.f_D.tolist(),
                   "D_over_T": late.D_over_T.tolist()}

    # whole-lung version, to show the compartment is unmeasurable there
    W = pd.read_csv(_fetch(MOUSE_WHOLE), index_col=0)
    c = W["cell.type"].value_counts()
    out["whole_lung"] = {"n_cells": int(len(W)), "AT1": int(c.get("AT1 cells", 0)),
                         "Krt8_ADI": int(c.get("Krt8 ADI", 0)),
                         "AT2": int(c.get("AT2 cells", 0) + c.get("Activated AT2 cells", 0))}
    return out


# ============================================================ part D: can we?
def part_d(B, C):
    """Is the AT1 compartment recovered comparably across conditions? If not,
    the ratio test is not interpretable, whatever it says."""
    import pandas as pd
    R = pd.DataFrame(C["time_course"])
    pbs = R[R.day == 0].iloc[0]
    inj = R[(R.day >= 2) & (R.day <= 21)]
    out = {"mouse_AT1_fraction_pbs": float(pbs["T"] / pbs.n_cells),
           "mouse_AT1_fraction_injury_median": float(np.median(inj["T"] / inj.n_cells)),
           "mouse_AT1_counts": {"pbs": int(pbs["T"]),
                                "injury_range": [int(inj["T"].min()), int(inj["T"].max())]},
           "human_AT1_counts": {"control": B["Control"]["T"], "disease": B["Disease"]["T"]},
           "human_AT1_fraction": {
               "control": B["Control"]["T"] / (B["Control"]["A"] + B["Control"]["D"] + B["Control"]["T"]),
               "disease": B["Disease"]["T"] / (B["Disease"]["A"] + B["Disease"]["D"] + B["Disease"]["T"])}}
    out["mouse_AT1_recovery_changes"] = bool(
        out["mouse_AT1_fraction_injury_median"] > 2 * out["mouse_AT1_fraction_pbs"])
    # sensitivity: if AT1 is recovered with efficiency r in disease relative to
    # control, what does the observed ratio fold change become?
    obs = B["ratio_fold_change"]["obs"]
    out["sensitivity"] = {str(r): float(obs * r) for r in (0.25, 0.5, 1.0, 2.0, 4.0)}
    out["r_that_cancels"] = float(1.0 / obs)
    return out


# ============================================================ report
def report(R):
    A, B, C, D = R["A"], R["B"], R["C"], R["D"]
    print("=" * 78)
    print("D26  alveolar three-state: what a composition snapshot can and cannot prove")
    print("=" * 78)

    print("\n  PART A  the algebra, checked numerically")
    print(f"    A1  fractions proportional to residence times (1/k1 : 1/k2 : 1/d):"
          f" max relative error {A['a1_max_rel_error']:.2e} over 200 random parameter sets")
    print(f"        => a snapshot identifies only RATIOS of residence times, never days")
    print(f"    A2  D/T across k1 spanning 1e-2 to 1e3 and s spanning 100x: "
          f"mean {A['a2_ratio_mean']:.6f}, spread {A['a2_ratio_spread']:.2e}; "
          f"predicted d/k2 = {A['a2_predicted_d_over_k2']:.6f}")
    print("        => accelerating the entry cannot change the transitional-to-AT1 ratio")
    a3 = A["a3"]
    print(f"    A3  ceiling as k1 -> infinity: predicted {a3['ceiling_predicted']:.4f}, "
          f"simulated {a3['ceiling_simulated']:.4f} ({a3['note']})")

    print("\n  PART B  human IPF (Habermann 2020)")
    for st in ("Control", "Disease"):
        s = B[st]
        print(f"    {st:8s} AT2 {s['A']:5d}  transitional {s['D']:4d}  AT1 {s['T']:4d}   "
              f"f_D {s['f_D']:.3f}  D/T {s['D_over_T']:.2f} "
              f"[{s['D_over_T_ci'][0]:.2f}, {s['D_over_T_ci'][1]:.2f}]  "
              f"({s['n_donors']} donors)")
    print(f"    fraction ceiling implied by the controls: "
          f"{B['ceiling_from_controls']:.3f}; disease f_D {B['disease_f_D']:.3f} -> "
          f"ceiling {'EXCEEDED' if B['fraction_ceiling_exceeded'] else 'not exceeded'}")
    fc = B["ratio_fold_change"]
    print(f"    D/T fold change disease vs control: {fc['obs']:.2f}x "
          f"[{fc['ci'][0]:.2f}, {fc['ci'][1]:.2f}] (donor bootstrap)")

    print("\n  PART C  mouse bleomycin time course (Strunz 2020, epithelium-enriched)")
    print(f"    {'day':>5s} {'AT2':>6s} {'Krt8+ADI':>9s} {'AT1':>5s} {'f_D':>7s} {'D/T':>7s}")
    for r in C["time_course"]:
        flag = "  <- too few cells to use" if (r["n_cells"] < 200 or r["T"] < 10) else ""
        print(f"    {r['day']:5d} {r['A']:6d} {r['D']:9d} {r['T']:5d} "
              f"{r['f_D']:7.3f} {r['D_over_T']:7.2f}{flag}")
    print(f"    ceiling implied by PBS: {C['ceiling_from_pbs']:.3f}")
    w = C["whole_lung"]
    print(f"    the same question in the WHOLE-LUNG deposit: AT1 = {w['AT1']} cells of "
          f"{w['n_cells']} (Krt8 ADI {w['Krt8_ADI']}, AT2 {w['AT2']}) -> unmeasurable")

    print("\n  PART D  is the AT1 compartment recovered comparably across conditions?")
    print(f"    mouse AT1 as a fraction of recovered cells: PBS "
          f"{D['mouse_AT1_fraction_pbs']:.4f} vs injury median "
          f"{D['mouse_AT1_fraction_injury_median']:.4f}  "
          f"(counts {D['mouse_AT1_counts']['pbs']} vs "
          f"{D['mouse_AT1_counts']['injury_range'][0]}-{D['mouse_AT1_counts']['injury_range'][1]})")
    print(f"    human AT1 counts: control {D['human_AT1_counts']['control']} vs disease "
          f"{D['human_AT1_counts']['disease']}  (fractions "
          f"{D['human_AT1_fraction']['control']:.3f} vs {D['human_AT1_fraction']['disease']:.3f})")
    print(f"    recovery appears condition-dependent in mouse: "
          f"{D['mouse_AT1_recovery_changes']}")
    print("    sensitivity - observed D/T fold change x relative AT1 recovery r:")
    print("      " + ", ".join(f"r={k}: {v:.2f}x" for k, v in D["sensitivity"].items()))
    print(f"      a relative recovery of r = {D['r_that_cancels']:.2f} would erase the "
          "human result entirely")


def figures(R):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
    os.makedirs(FIG, exist_ok=True)
    A, B, C, D = R["A"], R["B"], R["C"], R["D"]
    fig, ax = plt.subplots(1, 3, figsize=(17, 5.2))

    a = ax[0]
    T = pd.DataFrame(C["time_course"]).sort_values("day")
    inj = T[(T.day > 0) & (T.n_cells >= 200)]     # day 28 has 48 cells, 0 AT1
    a.plot(inj.day, inj.f_D, "o-", color="C0", label="transitional fraction f_D")
    a.plot(inj.day, inj.f_T, "s-", color="C2", label="AT1 fraction f_T")
    a.axhline(T[T.day == 0].f_D.iloc[0], color="C0", ls=":", lw=1.2, label="PBS f_D")
    a.axhline(C["ceiling_from_pbs"], color="C3", ls="--", lw=1.4,
              label="ceiling from PBS composition")
    a.set_xlabel("days after bleomycin")
    a.set_ylabel("fraction of the AT2 / transitional / AT1 pool")
    a.set_title("C  mouse: the transitional fraction does exceed\nthe PBS-implied ceiling "
                "- but PBS rests on 3 cells", fontsize=10)
    a.legend(fontsize=7.5)

    a = ax[1]
    labs, vals, cis = [], [], []
    for st in ("Control", "Disease"):
        labs.append(st)
        vals.append(B[st]["D_over_T"])
        cis.append(B[st]["D_over_T_ci"])
    a.bar(range(2), vals, 0.5, color=["C2", "C3"],
          yerr=[[v - c[0] for v, c in zip(vals, cis)],
                [c[1] - v for v, c in zip(vals, cis)]], capsize=5)
    a.set_xticks(range(2), labs)
    a.set_ylabel("transitional / AT1 ratio  (D/T)")
    fc = B["ratio_fold_change"]
    a.set_title(f"B  human IPF: D/T rises {fc['obs']:.1f}x\n"
                f"[{fc['ci'][0]:.1f}, {fc['ci'][1]:.1f}] — but see panel 3", fontsize=10)

    a = ax[2]
    xs = [0, 1, 2, 3]
    vals = [D["mouse_AT1_fraction_pbs"], D["mouse_AT1_fraction_injury_median"],
            D["human_AT1_fraction"]["control"], D["human_AT1_fraction"]["disease"]]
    a.bar(xs, vals, 0.55, color=["C0", "C1", "C0", "C1"])
    a.set_xticks(xs, ["mouse\nPBS", "mouse\ninjury", "human\ncontrol", "human\ndisease"],
                 fontsize=8.5)
    a.set_ylabel("AT1 share of the recovered pool")
    a.set_title("D  AT1 recovery moves with condition,\nand in opposite directions",
                fontsize=10)

    fig.tight_layout()
    fn = os.path.join(FIG, "d26_alveolar_three_state.png")
    fig.savefig(fn, dpi=140)
    print(f"\nfigure: {fn}")


def main(force=False):
    if os.path.exists(RES) and not force:
        R = pickle.load(open(RES, "rb"))
    else:
        A = part_a()
        B = part_b()
        C = part_c()
        R = {"A": A, "B": B, "C": C, "D": part_d(B, C)}
        pickle.dump(R, open(RES, "wb"))
    report(R)
    figures(R)
    return R


if __name__ == "__main__":
    main(force=bool(os.environ.get("D26_FORCE")))
