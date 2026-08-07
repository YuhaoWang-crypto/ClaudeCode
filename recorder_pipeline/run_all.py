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


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    for mod in (R1, R2, R3, R4):
        mod.main()
        print()
    paths = [fig_kernels(), fig_pairing(), fig_individuality()]
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
