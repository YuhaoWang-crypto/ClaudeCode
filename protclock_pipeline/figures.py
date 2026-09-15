"""Figures for the proteomic-aging-clock reproduction."""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = os.path.join("figures", "protclock")
ARM_ORDER = ["placebo", "rento_30mg_QD", "rento_30mg_BID", "rento_60mg_QD"]
ARM_LABEL = {"placebo": "placebo", "rento_30mg_QD": "30 QD",
             "rento_30mg_BID": "30 BID", "rento_60mg_QD": "60 QD"}


def _ensure():
    os.makedirs(OUT, exist_ok=True)


def fig_trajectories(scored, fname="f1_accel_trajectories.png"):
    """Mean age acceleration over time, per arm, one panel per clock."""
    _ensure()
    clocks = sorted(scored["clock"].unique())
    fig, axes = plt.subplots(2, 3, figsize=(14, 7), sharex=True)
    for ax, clk in zip(axes.ravel(), clocks):
        d = scored[scored["clock"] == clk]
        for arm in ARM_ORDER:
            g = d[d["arm"] == arm].groupby("week")["accel_z"]
            m, se = g.mean(), g.std() / np.sqrt(g.count())
            ax.errorbar(m.index, m.values, yerr=se.values, marker="o",
                        capsize=3, label=ARM_LABEL[arm],
                        lw=2.2 if arm != "placebo" else 1.6,
                        ls="-" if arm != "placebo" else "--")
        ax.axhline(0, color="0.6", lw=0.8)
        ax.set_title(clk, fontsize=10)
        ax.set_xlabel("week")
        ax.set_ylabel("age accel (ref SD)")
    axes.ravel()[0].legend(fontsize=7, ncol=2)
    fig.suptitle("Age acceleration over the trial  (SIMULATED DATA)",
                 fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, fname), dpi=140)
    plt.close(fig)
    return fname


def fig_forest(vs_placebo, fname="f2_forest.png"):
    """Week-12 change vs placebo, per clock and arm, with 95% CIs."""
    _ensure()
    arms = [a for a in ARM_ORDER if a != "placebo"]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.6), sharex=True, sharey=True)
    for ax, arm in zip(axes, arms):
        d = vs_placebo[vs_placebo["arm"] == arm].sort_values("clock")
        y = np.arange(len(d))
        ax.errorbar(d["diff_z"], y,
                    xerr=[d["diff_z"] - d["ci_low"], d["ci_high"] - d["diff_z"]],
                    fmt="o", capsize=3, color="#22577a")
        ax.axvline(0, color="0.4", lw=1)
        ax.set_yticks(y)
        ax.set_yticklabels(d["clock"], fontsize=8)
        ax.set_title(ARM_LABEL[arm], fontsize=10)
        ax.set_xlabel("change vs placebo (ref SD)")
    fig.suptitle("Week-12 age acceleration vs placebo; left of 0 = younger"
                 "  (SIMULATED)", fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, fname), dpi=140)
    plt.close(fig)
    return fname


def fig_power(sweep, fname="f3_power.png"):
    """Detection and concordance against planted effect size."""
    _ensure()
    fig, ax = plt.subplots(figsize=(7, 4.6))
    ax.plot(sweep["planted_years"], sweep["detect_rate"], marker="o",
            label="fraction of tests p<0.05", color="#c1121f")
    ax.plot(sweep["planted_years"], sweep["concordant_frac"], marker="s",
            label="replicates with >=5/6 clocks agreeing", color="#22577a")
    ax.plot(sweep["planted_years"], sweep["best_arm_correct"], marker="^",
            label="strongest arm identified correctly", color="#588157")
    ax.axhline(0.05, color="0.5", ls=":", label="nominal alpha = 0.05")
    ax.set_xlabel("planted effect (proteomic-age-equivalent years, undiluted)")
    ax.set_ylabel("rate")
    ax.set_ylim(-0.03, 1.05)
    ax.legend(fontsize=8)
    ax.set_title("What 42 patients can detect  (SIMULATED)",
                 fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, fname), dpi=140)
    plt.close(fig)
    return fname


def fig_volcano(differential, aging, fname="f4_volcano.png"):
    """Effect vs significance, aging proteins highlighted."""
    _ensure()
    eff = differential["effect_npx"].to_numpy()
    q = differential["q"].to_numpy()
    is_age = aging["is_aging"].to_numpy()
    y = -np.log10(np.maximum(q, 1e-12))

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(eff[~is_age], y[~is_age], s=6, alpha=0.35, color="0.6",
               label="not age-associated")
    ax.scatter(eff[is_age], y[is_age], s=9, alpha=0.75, color="#c1121f",
               label="age-associated (reference cohort)")
    ax.axhline(-np.log10(0.05), color="0.3", ls="--", lw=1,
               label="BH q = 0.05")
    ax.set_xlabel("treated change from baseline (NPX)")
    ax.set_ylabel("-log10 BH q")
    ax.legend(fontsize=8)
    ax.set_title("Differential abundance at week 12  (SIMULATED)",
                 fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, fname), dpi=140)
    plt.close(fig)
    return fname


def fig_pathways(pathways, fname="f5_pathways.png"):
    """Aging-aligned shift per pathway."""
    _ensure()
    d = pathways.sort_values("aligned_npx")
    colors = ["#588157" if v < 0 else "#c1121f" for v in d["aligned_npx"]]
    fig, ax = plt.subplots(figsize=(7.5, 3.8))
    ax.barh(d["pathway"], d["aligned_npx"], color=colors)
    for i, (v, qv) in enumerate(zip(d["aligned_npx"], d["aligned_q"])):
        ax.text(v, i, f"  q={qv:.3g}", va="center", fontsize=8,
                ha="right" if v < 0 else "left")
    ax.axvline(0, color="0.3", lw=1)
    ax.set_xlabel("aging-aligned shift (NPX); negative = toward younger")
    ax.set_title("Pathway response  (SIMULATED)", fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, fname), dpi=140)
    plt.close(fig)
    return fname


def fig_clock_corr(corr, eff, fname="f6_clock_correlation.png"):
    """REAL cross-clock correlation from the paper's supplementary Table S2."""
    _ensure()
    m = corr.to_numpy()
    labels = list(corr.columns)
    fig, ax = plt.subplots(figsize=(6.6, 5.6))
    im = ax.imshow(m, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=8)
    for i in range(len(labels)):
        for j in range(len(labels)):
            ax.text(j, i, f"{m[i, j]:.2f}", ha="center", va="center",
                    fontsize=7.5,
                    color="white" if abs(m[i, j]) > 0.6 else "black")
    fig.colorbar(im, ax=ax, shrink=0.8, label="Pearson r")
    ax.set_title(f"Six clocks on the REAL trial samples\n"
                 f"effective independent clocks: {eff['li_ji']:.1f} of 6",
                 fontweight="bold", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, fname), dpi=140)
    plt.close(fig)
    return fname


def write_all(pre, scored, r5, r6, r7, sweep=None, supp=None):
    made = [fig_trajectories(scored), fig_forest(r5["vs_placebo"])]
    made.append(fig_volcano(r6["paired_treated"]["differential"], r6["aging"]))
    made.append(fig_pathways(r7))
    if sweep is not None:
        made.append(fig_power(sweep))
    if supp is not None:
        made.append(fig_clock_corr(supp["corr"], supp["effective"]))
    print(f"\n  figures written to {OUT}/:")
    for f in made:
        print(f"    {f}")
    return made
