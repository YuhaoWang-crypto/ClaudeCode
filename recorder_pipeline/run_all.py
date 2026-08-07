"""
Run the whole recorder pipeline and write figures to figures/.

    python3 -m recorder_pipeline.run_all
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from . import r1_kernel as R1
from . import r2_pairing as R2
from . import r3_individuality as R3
from . import r4_catalog as R4
from . import r5_datacases as R5
from . import r6_validation as R6
from . import r7_crossdisease as R7
from . import r8_audit as R8
from . import r9_discovery as R9

FIGDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      "figures")


def fig_kernels():
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))

    # (a) kernel shapes
    x = np.linspace(0, 130, 4000)
    ax[0].plot(x, R1.kernel_deterministic(x, 120), lw=2.2,
               label="RBC, fixed lifespan 120 d\n(HbA1c) — linear ramp, HARD horizon")
    ax[0].plot(x, R1.kernel_exponential(x, 19 / R1.LN2), lw=2.2, ls="--",
               label="albumin, random clearance\n(glycated albumin) — exponential tail")
    ax[0].plot(x, R1.kernel_exponential(x, 10 / R1.LN2), lw=1.8, ls=":",
               label="transferrin (%CDT)")
    ax[0].set_xlabel("lag before sampling (days)")
    ax[0].set_ylabel("K(x) = P(carrier age > x)")
    ax[0].set_title("(a) The kernel is the carrier AGE survivor\nnot its decay curve")
    ax[0].legend(fontsize=7.5, loc="upper right")
    ax[0].grid(alpha=0.3)

    # (b) HbA1c weight by month, derived
    L = 120.0
    xx = np.linspace(0, L, 20001)
    K = R1.kernel_deterministic(xx, L)
    fr = R1.window_fractions(K, xx, [0, 30, 60, 90, 120])
    bars = ax[1].bar(["0-30 d", "30-60 d", "60-90 d", "90-120 d"],
                     [f * 100 for f in fr], color="#4c72b0")
    for b, f in zip(bars, fr):
        ax[1].text(b.get_x() + b.get_width() / 2, f * 100 + 1,
                   f"{f*100:.1f}%", ha="center", fontsize=9)
    ax[1].axhline(50, color="crimson", ls="--", lw=1.2,
                  label="clinically quoted ~50% for recent month")
    ax[1].set_ylabel("% of HbA1c signal")
    ax[1].set_title("(b) HbA1c time weighting derived from\none number (RBC lifespan)")
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=0.3, axis="y")

    # (c) step response
    t = np.linspace(0, 180, 1801)
    for name, kind, param, tm in R1.CARRIERS:
        if "AGE" in name or "TAT" in name:
            continue
        xx2, KK = R1.build(name, kind, param, tm)
        traj = R1.step_response(KK, xx2, t, drop=0.5)
        ax[2].plot(t, traj, lw=2, label=name.split(" (")[0])
    ax[2].axvline(28, color="k", ls=":", lw=1.2)
    ax[2].text(30, 0.95, "4-week\nendpoint", fontsize=8)
    ax[2].axhline(0.5, color="gray", ls="--", lw=0.8)
    ax[2].set_xlabel("days after disease activity halves")
    ax[2].set_ylabel("recorder / baseline")
    ax[2].set_title("(c) Which recorder can serve a\n4-week treatment endpoint")
    ax[2].legend(fontsize=7.5)
    ax[2].grid(alpha=0.3)

    fig.tight_layout()
    p = os.path.join(FIGDIR, "recorder_r1_kernels.png")
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p


def fig_pairing():
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))

    # (a) admissibility map: gain of the fixed ratio over (rho, kappa)
    rho = np.linspace(0, 0.99, 300)
    kappa = np.linspace(0.05, 3.0, 300)
    RR, KK = np.meshgrid(rho, kappa)
    gain = 1.0 / np.sqrt(1.0 + KK**2 - 2 * RR * KK)
    m = ax[0].pcolormesh(RR, KK, gain, cmap="RdBu", vmin=0.3, vmax=1.7,
                         shading="auto")
    ax[0].plot(rho, 2 * rho, "k-", lw=2.5, label=r"admissibility  $\rho=\kappa/2$")
    ax[0].set_xlabel(r"$\rho$  (correlation of non-disease variation)")
    ax[0].set_ylabel(r"$\kappa=\sigma_B/\sigma_A$  (denominator noise)")
    ax[0].set_title("(a) Fixed ratio A/B: gain vs harm\nRED = the ratio destroys signal")
    ax[0].set_ylim(0.05, 3.0)
    ax[0].legend(fontsize=8, loc="upper left")
    fig.colorbar(m, ax=ax[0], label="separation gain")

    # (b) residual gain bound
    r = np.linspace(0, 0.97, 400)
    ax[1].plot(r, R2.residual_gain(r), lw=2.4, color="#c44e52")
    ax[1].axhline(1.4, color="k", ls="--", lw=1)
    ax[1].axvline(0.7, color="k", ls="--", lw=1)
    ax[1].text(0.72, 1.05, "below ρ≈0.7 a second\nassay buys <40%", fontsize=8.5)
    ax[1].set_xlabel(r"$\rho$")
    ax[1].set_ylabel("separation gain of fitted residual")
    ax[1].set_title(r"(b) Residual bound $1/\sqrt{1-\rho^2}$" "\nnever below 1 — never harmful")
    ax[1].set_ylim(1, 4)
    ax[1].grid(alpha=0.3)

    # (c) the four denominator levels
    labels = ["L0\nalbumin\n(unrelated)", "L3\nsame\nsource cell", "L1\nsame\nparent",
              "L4\nopposing\naxis"]
    setups = [(0.0, 0.5, 0.15), (0.0, 1.0, 0.60), (0.0, 1.0, 0.90), (-1.0, 1.0, 0.60)]
    gains = []
    for d_b, kap, rh in setups:
        d0 = R2.dprime(0.0, 1.0, d_b, 1.0, kap, rh)
        d1 = R2.dprime(1.0, 1.0, d_b, 1.0, kap, rh)
        gains.append(d1 / d0)
    cols = ["#c44e52" if g < 1 else "#55a868" for g in gains]
    bars = ax[2].bar(labels, gains, color=cols)
    for b, g in zip(bars, gains):
        ax[2].text(b.get_x() + b.get_width() / 2, g + 0.03, f"{g:.2f}x",
                   ha="center", fontsize=9)
    ax[2].axhline(1.0, color="k", lw=1.2)
    ax[2].set_ylabel("separation vs the marker alone")
    ax[2].set_title("(c) The denominator hierarchy,\nderived not asserted")
    ax[2].grid(alpha=0.3, axis="y")

    fig.tight_layout()
    p = os.path.join(FIGDIR, "recorder_r2_pairing.png")
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p


def fig_individuality():
    rows = []
    for name, cvi, cvg, cva, note in R3.BV_TABLE:
        rows.append((name, R3.index_of_individuality(cvi, cvg),
                     R3.personal_dividend(cvi, cvg, cva)))
    rows.sort(key=lambda r: r[1])
    names = [r[0] for r in rows]
    iis = [r[1] for r in rows]
    divs = [r[2] for r in rows]

    fig, ax = plt.subplots(1, 2, figsize=(12.5, 5.2))
    cols = ["#c44e52" if v < 0.6 else "#8c8c8c" for v in iis]
    ax[0].barh(names, iis, color=cols)
    ax[0].axvline(0.6, color="k", ls="--", lw=1.4)
    ax[0].text(0.62, 0.3, "II < 0.6:\npopulation reference\ninterval uninformative",
               fontsize=8.5)
    ax[0].set_xlabel("index of individuality  CV_I / CV_G")
    ax[0].set_title("(a) Which analytes a population cutoff wastes")
    ax[0].grid(alpha=0.3, axis="x")

    ax[1].barh(names, divs, color=cols)
    ax[1].axvline(1.0, color="k", lw=1.2)
    ax[1].set_xlabel("separation gain from a personal baseline")
    ax[1].set_title("(b) The personal-baseline dividend\n(single prior sample)")
    ax[1].grid(alpha=0.3, axis="x")

    fig.tight_layout()
    p = os.path.join(FIGDIR, "recorder_r3_individuality.png")
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p


def fig_dataspace():
    """The actual clinical data space: where normal and disease sit, and
    whether the assay is good enough to tell them apart."""
    from scipy.stats import norm as _norm
    fig, ax = plt.subplots(1, 2, figsize=(14.5, 5.6))

    # (a) normal vs disease distributions on a log axis
    rows = []
    for c in R5.CASES:
        r = R5.analyse(c)
        rows.append((c, r))
    rows = sorted(rows, key=lambda t: t[1]["delta"])

    for i, (c, r) in enumerate(rows):
        mu_n = np.log(c["normal_median"])
        mu_d = np.log(c["disease_median"])
        for mu, s, col, lab in [(mu_n, r["sigma_norm"], "#4c72b0", "normal"),
                                (mu_d, r["sigma_dis"], "#c44e52", "disease")]:
            lo, hi = mu - 1.96 * s, mu + 1.96 * s
            ax[0].plot([lo, hi], [i, i], color=col, lw=7, alpha=0.55,
                       solid_capstyle="butt",
                       label=lab if i == 0 else None)
            ax[0].plot([mu], [i], "o", color=col, ms=7, zorder=3)
        # annotate fold change
        ax[0].text(max(mu_n, mu_d) + 0.35, i, f"{r['fold']:.1f}x",
                   va="center", fontsize=8.5)
    def _wrap(name):
        if " (" in name:
            head, tail = name.split(" (", 1)
            return f"{head}\n({tail}"
        return name
    ax[0].set_yticks(range(len(rows)))
    ax[0].set_yticklabels([_wrap(c["name"]) for c, _ in rows], fontsize=7.5)
    ax[0].set_xlabel("ln(concentration)  — units differ per analyte")
    ax[0].set_title("(a) The data space: published normal vs disease\n"
                    "bars = central 95%, dot = median")
    ax[0].legend(fontsize=8, loc="lower right")
    ax[0].grid(alpha=0.3, axis="x")

    # (b) feasibility plane: biological effect vs assay noise
    eff = np.linspace(0.15, 4.0, 300)
    cvs = np.linspace(0.01, 0.35, 300)
    E, C = np.meshgrid(eff, cvs)
    sig_bio = 0.45                      # typical pooled biological ln-sd
    sig_tot = np.sqrt(sig_bio**2 + np.log(1 + C**2))
    AUC = _norm.cdf((E / sig_tot) / np.sqrt(2))
    m = ax[1].contourf(E, C * 100, AUC, levels=np.linspace(0.5, 1.0, 21),
                       cmap="viridis")
    ax[1].contour(E, C * 100, AUC, levels=[0.70, 0.80, 0.90],
                  colors="w", linewidths=1.6)
    SHORT = {"GFAP (TBI, CT-/MRI+ vs healthy)": "GFAP",
             "Serum uromodulin (CKD G3 vs non-CKD)": "uromodulin",
             "Carbamylated albumin (ESRD vs non-uremic)": "C-Alb (vs normal)",
             "C-Alb for 1-yr mortality WITHIN ESRD": "C-Alb (mortality)",
             "PRO-C3 (F4 vs F0/F1, NAFLD)": "PRO-C3",
             "TAT (DIC vs reference)": "TAT"}
    offsets = {"C-Alb (mortality)": (6, -16), "C-Alb (vs normal)": (6, 8),
               "uromodulin": (6, -16), "TAT": (8, 4), "GFAP": (-46, 8),
               "PRO-C3": (8, 6)}
    for c, r in rows:
        lab = SHORT.get(c["name"], c["name"])
        ax[1].plot(r["delta"], c["cv_analytical"] * 100, "o",
                   ms=9, mec="w", mew=1.6, color="#d62728")
        ax[1].annotate(lab, (r["delta"], c["cv_analytical"] * 100),
                       textcoords="offset points",
                       xytext=offsets.get(lab, (7, 6)), fontsize=8,
                       color="w", fontweight="bold")
    ax[1].set_xlabel("biological effect  |Δ ln(median)|")
    ax[1].set_ylabel("assay analytical CV (%)")
    ax[1].set_title("(b) Feasibility plane (pooled biological sd fixed at 0.45)\n"
                    "white contours = AUC 0.70 / 0.80 / 0.90")
    fig.colorbar(m, ax=ax[1], label="achievable AUC")

    fig.tight_layout()
    p = os.path.join(FIGDIR, "recorder_r5_dataspace.png")
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p


def fig_validation():
    """Primary-data validation: source specificity and within-person behaviour."""
    hpa = R6.load_hpa()
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10))
    ax = axes.ravel()

    # (a) how many cell types express each parent
    labels, counts, clean = [], [], []
    for lab in R6.CLAIMS:
        rec = hpa[lab]
        _, blood_clean, n_ct = R6.source_verdict(rec)
        labels.append(lab.split(" (")[0])
        counts.append(n_ct)
        clean.append(blood_clean)
    order = np.argsort(counts)
    labels = [labels[i] for i in order]
    counts = [counts[i] for i in order]
    clean = [clean[i] for i in order]
    cols = ["#55a868" if c else "#c44e52" for c in clean]
    ax[0].barh(labels, counts, color=cols)
    ax[0].set_xlabel("# cell types with enriched/enhanced expression (HPA)")
    ax[0].set_title("(a) Source specificity of each parent\n"
                    "green = no leukocyte expression, red = leukocyte source")
    ax[0].grid(alpha=0.3, axis="x")
    ax[0].tick_params(labelsize=8)

    # (b) AD panel within-person CV, ratio highlighted
    names = [r[0] for r in R6.AD_BV]
    cvi = [r[1] for r in R6.AD_BV]
    cols = ["#dd8452"] * len(names)
    cols[names.index("Abeta42/Abeta40")] = "#4c72b0"
    bars = ax[1].bar(range(len(names)), cvi, color=cols)
    ax[1].set_xticks(range(len(names)))
    ax[1].set_xticklabels(names, rotation=40, ha="right", fontsize=8)
    ax[1].set_ylabel("within-person CV_I (%)")
    ax[1].set_title("(b) 20 healthy adults, weekly x 10 wk\n"
                    "the RATIO halves its own components' noise")
    for i in (0, 1, 2):
        ax[1].text(i, cvi[i] + 0.4, f"{cvi[i]}", ha="center", fontsize=8.5)
    ax[1].annotate("", xy=(2, 3.8), xytext=(0.5, 6.9),
                   arrowprops=dict(arrowstyle="->", lw=1.8, color="#4c72b0"))
    ax[1].text(0.9, 8.4, r"measured $\rho$=0.87" "\n1.97x gain", fontsize=8.5,
               color="#4c72b0")
    ax[1].grid(alpha=0.3, axis="y")

    # (c) creatinine normalisation: works for one, not the other
    names = [r[0] for r in R6.URINE_BV]
    raw = [r[1] for r in R6.URINE_BV]
    nrm = [r[2] for r in R6.URINE_BV]
    x = np.arange(len(names))
    ax[2].bar(x - 0.19, raw, 0.38, label="raw concentration", color="#8c8c8c")
    ax[2].bar(x + 0.19, nrm, 0.38, label="/ urine creatinine", color="#4c72b0")
    for i, (r_, n_) in enumerate(zip(raw, nrm)):
        ax[2].text(i + 0.19, n_ + 0.8, f"−{100*(r_-n_)/r_:.0f}%", ha="center",
                   fontsize=8.5)
    ax[2].set_xticks(x)
    ax[2].set_xticklabels(names, fontsize=9)
    ax[2].set_ylabel("within-person CV_I (%)")
    ax[2].set_title("(c) Same denominator, opposite verdicts\n"
                    "23 stable CKD patients, daily x 10 d")
    ax[2].legend(fontsize=8)
    ax[2].grid(alpha=0.3, axis="y")

    # (d) CELLxGENE: healthy vs CKD vs AKI, within dataset & cell type
    rows = R6._read_tsv(R6.CENSUS_PATH)
    order = ["normal", "chronic kidney disease", "acute kidney failure"]
    short = ["normal", "CKD", "AKI"]
    panels = [("HAVCR1", "proximal tubule", "KIM-1 (reproducible)", "-"),
              ("LCN2", "proximal tubule", "NGAL (NOT reproducible)", "--"),
              ("UMOD", "thick ascending limb", "UMOD (flat)", ":")]
    cols = {"a12ccb9b": "#4c72b0", "dea717d4": "#c44e52"}
    xs = np.arange(3)
    for gene, ctk, lab, ls in panels:
        for ds in ["a12ccb9b", "dea717d4"]:
            sel = {r["disease"]: r for r in rows if r["gene"] == gene
                   and r["dataset"] == ds and ctk in r["cell_type"]}
            if len(sel) < 3:
                continue
            base = float(sel["normal"]["mean_expr"])
            y = [float(sel[d]["mean_expr"]) / base for d in order]
            ax[3].plot(xs, y, ls, marker="o", color=cols[ds], lw=2,
                       label=f"{gene} · {ds}")
    ax[3].set_yscale("log")
    ax[3].axhline(1.0, color="k", lw=1)
    ax[3].set_xticks(xs)
    ax[3].set_xticklabels(short)
    ax[3].set_ylabel("fold change vs normal (log scale)")
    ax[3].set_title("(d) Healthy vs indication, single cell\n"
                    "within dataset & cell type (CELLxGENE, 1.05M cells)")
    ax[3].legend(fontsize=7.5, ncol=2)
    ax[3].grid(alpha=0.3)

    fig.tight_layout()
    p = os.path.join(FIGDIR, "recorder_r6_validation.png")
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p


def fig_crossdisease():
    """SLE, brain, and the individual-level creatinine result."""
    import pandas as pd
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))

    # (a) SLE: IFN control up, CR1 down -> CR1 is not a valid denominator
    s = R7._load("census_sle.tsv")

    def piv(gene):
        t = s.pivot(index="cell_type", columns="disease", values=f"{gene}_pct").dropna()
        t.columns = ["SLE" if "lupus" in c else "normal" for c in t.columns]
        return t
    cts = ["B cell", "plasmablast", "classical monocyte",
           "CD4-positive, alpha-beta T cell"]
    ifn, cr1 = piv("IFI44L"), piv("CR1")
    lab = [c.split(",")[0][:16] for c in cts if c in ifn.index and c in cr1.index]
    fi = [ifn.loc[c, "SLE"] / ifn.loc[c, "normal"] for c in cts
          if c in ifn.index and c in cr1.index]
    fc = [cr1.loc[c, "SLE"] / cr1.loc[c, "normal"] for c in cts
          if c in ifn.index and c in cr1.index]
    x = np.arange(len(lab))
    ax = axes[0]
    ax.bar(x - 0.19, fi, 0.38, color="#dd8452", label="IFI44L (disease control)")
    ax.bar(x + 0.19, fc, 0.38, color="#4c72b0", label="CR1 (proposed denominator)")
    ax.axhline(1.0, color="k", lw=1.2)
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(lab, rotation=25, ha="right", fontsize=8)
    ax.set_ylabel("SLE / normal (% positive)")
    ax.set_title("(a) SLE, 1.26M cells, one dataset\n"
                 "CR1 falls in SLE — a denominator must not")
    ax.legend(fontsize=7.5)
    ax.grid(alpha=0.3, axis="y")

    # (b) brain: GFAP flat per cell, VIM up
    b = R7._load("census_brain_astro.tsv")
    ax = axes[1]
    for gene, col in [("GFAP", "#c44e52"), ("VIM", "#55a868"),
                      ("CAPN1", "#8c8c8c")]:
        t = b.pivot(index="ds", columns="disease", values=f"{gene}_mean").dropna()
        ad = [c for c in t.columns if "Alzheimer" in c][0]
        f = (t[ad] / t["normal"]).values
        ax.scatter(np.full(len(f), gene), f, s=55, color=col, alpha=0.8,
                   edgecolor="w", zorder=3)
        ax.scatter([gene], [np.median(f)], marker="_", s=900, color=col, zorder=4)
    ax.axhline(1.0, color="k", lw=1.2)
    ax.set_yscale("log")
    ax.set_ylabel("Alzheimer / control, per astrocyte")
    ax.set_title("(b) Astrocytes, 6 paired datasets\n"
                 "reactivity (VIM) rises, GFAP synthesis does not")
    ax.grid(alpha=0.3, axis="y")

    # (c) individual-level creatinine trajectories
    cr = pd.read_csv(os.path.join(R7._D, "mimic_creatinine.tsv.gz"), sep="\t",
                     parse_dates=["charttime"])
    ax = axes[2]
    shown = 0
    for sid, g in cr.groupby("subject_id"):
        g = g.sort_values("charttime")
        if len(g) < 12 or shown >= 25:
            continue
        h = (g.charttime - g.charttime.iloc[0]).dt.total_seconds() / 86400
        if h.max() > 30:
            continue
        ax.plot(h, g.valuenum, lw=1.0, alpha=0.65)
        shown += 1
    ax.axhspan(0.5, 1.2, color="#55a868", alpha=0.18)
    ax.text(0.4, 1.25, "population reference interval 0.5–1.2 mg/dL",
            fontsize=8, color="#2d6a4f")
    ax.set_xlabel("days from first result")
    ax.set_ylabel("serum creatinine (mg/dL)")
    ax.set_title(f"(c) {shown} real patient trajectories (MIMIC-IV demo)\n"
                 "19.2% of KDIGO AKI events stay inside the green band")
    ax.grid(alpha=0.3)

    fig.tight_layout()
    p = os.path.join(FIGDIR, "recorder_r7_crossdisease.png")
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    for mod in (R1, R2, R3, R4, R5, R6, R7, R8, R9):
        mod.main()
        print()
    paths = [fig_kernels(), fig_pairing(), fig_individuality(), fig_dataspace(),
             fig_validation(), fig_crossdisease()]
    tsv = R4.export_tsv(os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "recorder_catalog.tsv"))
    print("=" * 78)
    print("figures written:")
    for p in paths:
        print("   ", os.path.relpath(p))
    print("   ", os.path.relpath(tsv))


if __name__ == "__main__":
    main()
