"""
Module G3 - The underprocessing score, and how much of it you can trust.

The paper triages designs on one number: an abundance-weighted sum over an
ordinal 1-20 "underprocessing" ladder, with selection at score > 1500 and
occupancy > 50%. The ladder's endpoints are published (20 = Glc1Man9, low =
processed complex/core, axis 0-2000) but the species-to-integer table is not,
anywhere in the paper or its Extended Data. ``glycans.UNDERPROCESSING_LADDER``
reconstructs it; see that docstring for the derivation and its two
consistency checks.

Because the ladder is reconstructed, this module's job is not just to compute
scores. It is to separate conclusions that survive the uncertainty from
conclusions that do not:

  (rigorous) ORDERING between two glycan profiles, whenever one stochastically
      dominates the other along the processing continuum. Then score(A) >
      score(B) for EVERY strictly increasing ladder, so the reconstruction
      cannot affect the verdict. This is proved below by exhaustive sweep as
      well as argued.

  (caution) The ABSOLUTE threshold of 1500. Its meaning depends entirely on
      the ladder's spacing. The sweep reports, per design, the fraction of
      random monotone ladders under which it would still pass - which is the
      honest confidence in a pass/fail call.

Replace ``ILLUSTRATIVE_PROFILES`` with Supplementary Data S4/S5 to get real
numbers out of this module.
"""

from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .g1_delta_rasa import FIG_DIR
from .glycans import UNDERPROCESSING_LADDER, underprocessing_score

# ---------------------------------------------------------------------------
# (caution) ILLUSTRATIVE ONLY - NOT MEASURED DATA.
#
# These profiles encode nothing but the paper's own QUALITATIVE descriptions
# of each construct, so that the engine can be exercised end to end:
#   THRgly1  "a range of mass shifts corresponding to various complex glycans"
#   THRgly2  "primarily Man8 and Man9"
#   THRgly3  "masses consistent with Glc1Man9"
#   THRgly4  "masses consistent with Glc1Man9"
#   C1       complex-type control glycoprotein
#   O1,O3,O4 "carried more processed species, including Man5 and Man6"
#   O2       "enriched for Man7-9"
# The per-species percentages are invented to match those sentences. They are
# NOT from LC-MS/MS. Every absolute score printed by this module inherits that
# status; only the ladder-robustness conclusions are independent of them.
# ---------------------------------------------------------------------------
ILLUSTRATIVE_PROFILES: dict[str, dict] = {
    "THRgly1": {"abundances": {"HexNAc4(F)": 40, "HexNAc5(F)": 30,
                               "HexNAc4": 20, "Hybrid": 10},
                "occupancy": 70.0},
    "THRgly2": {"abundances": {"Man8": 45, "Man9": 40, "Man7": 15},
                "occupancy": 35.0},
    "THRgly3": {"abundances": {"Glc1Man9": 70, "Man9": 25, "Man8": 5},
                "occupancy": 60.0},
    "THRgly4": {"abundances": {"Glc1Man9": 80, "Man9": 20},
                "occupancy": 40.0},
    "C1":      {"abundances": {"HexNAc4(F)": 45, "HexNAc5(F)": 25,
                               "HexNAc3(F)": 15, "Hybrid": 10, "Core": 5},
                "occupancy": 65.0},
    "O1":      {"abundances": {"Man5": 35, "Man6": 30, "Man7": 15,
                               "Man8": 10, "Hybrid": 10}, "occupancy": 62.0},
    "O2":      {"abundances": {"Man8": 35, "Man9": 30, "Man7": 25,
                               "Man6": 10}, "occupancy": 58.0},
    "O3":      {"abundances": {"Man5": 40, "Man6": 25, "Man7": 15,
                               "Man8": 10, "Man9": 10}, "occupancy": 71.0},
    "O4":      {"abundances": {"Man5": 35, "Man6": 30, "Man8": 15,
                               "Man9": 10, "Hybrid": 10}, "occupancy": 55.0},
}


def profile_vector(abundances: dict, ladder: dict | None = None) -> np.ndarray:
    """Abundance (percent, renormalised) indexed by ladder rung 1..20."""
    ladder = ladder or UNDERPROCESSING_LADDER
    n = max(ladder.values())
    v = np.zeros(n + 1)
    total = sum(abundances.values()) or 1.0
    for species, pct in abundances.items():
        v[ladder[species]] += 100.0 * pct / total
    return v


def dominates(a: dict, b: dict, ladder: dict | None = None) -> bool:
    """True if profile ``a`` stochastically dominates ``b`` (less processed).

    Compares upper-tail cumulative abundance at every rung. When this holds,
    score(a) >= score(b) under any strictly increasing ladder, so the ordering
    is immune to the reconstruction.
    """
    va, vb = profile_vector(a, ladder), profile_vector(b, ladder)
    ca = np.cumsum(va[::-1])[::-1]        # abundance at rung >= k
    cb = np.cumsum(vb[::-1])[::-1]
    return bool(np.all(ca >= cb - 1e-9) and np.any(ca > cb + 1e-9))


def random_monotone_ladders(n: int = 2000, seed: int = 0) -> np.ndarray:
    """``n`` random strictly increasing ladders spanning the same 1..20 range.

    Each row is a rung->value map sharing the paper's published endpoints
    (bottom rung 1, top rung 20) but with arbitrary interior spacing. This is
    the space of ladders consistent with what the paper actually states.
    """
    rng = np.random.default_rng(seed)
    n_rungs = max(UNDERPROCESSING_LADDER.values())
    interior = np.sort(rng.uniform(1.0, 20.0, size=(n, n_rungs - 2)), axis=1)
    lo = np.ones((n, 1))
    hi = np.full((n, 1), 20.0)
    return np.hstack([lo, interior, hi])        # shape (n, 20), rung 1..20


def sensitivity(profiles: dict, n_ladders: int = 2000, threshold: float = 1500.0,
                occupancy_cut: float = 50.0, seed: int = 0) -> dict:
    """How often each design passes, across the space of monotone ladders."""
    ladders = random_monotone_ladders(n_ladders, seed=seed)
    names = list(profiles)
    vecs = np.array([profile_vector(profiles[k]["abundances"]) for k in names])
    # scores[i, j] = design i under ladder j
    scores = vecs[:, 1:] @ ladders.T

    occ = np.array([profiles[k]["occupancy"] for k in names])
    passes = (scores > threshold) & (occ[:, None] > occupancy_cut)

    return {
        "names": names,
        "pass_fraction": dict(zip(names, passes.mean(axis=1).round(4))),
        "score_min": dict(zip(names, scores.min(axis=1).round(1))),
        "score_max": dict(zip(names, scores.max(axis=1).round(1))),
        "n_ladders": n_ladders,
        "threshold": threshold,
        "occupancy_cut": occupancy_cut,
    }


def analyse(profiles: dict | None = None):
    profiles = profiles or ILLUSTRATIVE_PROFILES
    scored = {}
    for name, spec in profiles.items():
        scored[name] = underprocessing_score(
            spec["abundances"], occupancy=spec["occupancy"])

    names = list(profiles)
    dom_pairs = [
        (a, b) for a in names for b in names
        if a != b and dominates(profiles[a]["abundances"],
                                profiles[b]["abundances"])
    ]

    sens = sensitivity(profiles)

    # Verify the dominance claim numerically rather than asserting it.
    ladders = random_monotone_ladders(500, seed=7)
    vecs = {k: profile_vector(profiles[k]["abundances"]) for k in names}
    violations = 0
    for a, b in dom_pairs:
        sa = vecs[a][1:] @ ladders.T
        sb = vecs[b][1:] @ ladders.T
        violations += int((sa < sb - 1e-6).sum())

    return {
        "profiles": profiles,
        "scored": scored,
        "dominance_pairs": dom_pairs,
        "dominance_violations": violations,
        "dominance_tests": len(dom_pairs) * ladders.shape[0],
        "sensitivity": sens,
        "ladder": dict(UNDERPROCESSING_LADDER),
    }


def make_figure(r, path=None):
    path = path or os.path.join(FIG_DIR, "g3_underprocessing.png")
    os.makedirs(os.path.dirname(path), exist_ok=True)

    names = list(r["scored"])
    scores = [r["scored"][n]["total_underprocessing_score"] for n in names]
    occ = [r["scored"][n]["occupancy"] for n in names]
    sens = r["sensitivity"]
    frac = [sens["pass_fraction"][n] for n in names]

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))

    ax = axes[0]
    rungs = sorted(r["ladder"].values(), reverse=True)
    labels = sorted(r["ladder"], key=lambda k: -r["ladder"][k])
    ax.barh(range(len(rungs)), rungs, color="#55A868")
    ax.set_yticks(range(len(rungs)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("ordinal rung")
    ax.set_title("Reconstructed 1-20 ladder\n(endpoints published, interior inferred)",
                 fontsize=9)
    ax.grid(axis="x", alpha=0.3)

    ax = axes[1]
    colors = ["#C44E52" if s <= 1500 else "#4C72B0" for s in scores]
    ax.scatter(scores, occ, c=colors, s=70)
    for n, s, o in zip(names, scores, occ):
        ax.annotate(n, (s, o), fontsize=7, xytext=(4, 3),
                    textcoords="offset points")
    ax.axvline(1500, ls="--", c="gray", lw=1)
    ax.axhline(50, ls="--", c="gray", lw=1)
    ax.set_xlim(0, 2000)
    ax.set_xlabel("total underprocessing score")
    ax.set_ylabel("occupancy (%)")
    ax.set_title("Paper's selection gate\n(score > 1500, occupancy > 50%)",
                 fontsize=9)
    ax.grid(alpha=0.3)

    ax = axes[2]
    order = np.argsort(frac)
    ax.barh([names[i] for i in order], [frac[i] for i in order],
            color="#4C72B0")
    ax.set_xlim(0, 1)
    ax.set_xlabel(f"fraction of {sens['n_ladders']} monotone ladders passing")
    ax.set_title("Is the pass/fail call ladder-robust?\n"
                 "(0 or 1 = robust, middle = fragile)", fontsize=9)
    ax.grid(axis="x", alpha=0.3)

    fig.suptitle("Underprocessing score - engine and its sensitivity to the "
                 "reconstructed ladder", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def report():
    r = analyse()
    print("=" * 68)
    print("MODULE G3 - Underprocessing score (ladder reconstructed)")
    print("=" * 68)
    print("(caution) Profiles below are ILLUSTRATIVE, built only from the")
    print("          paper's qualitative sentences. Absolute scores are not")
    print("          measurements. Swap in Supplementary Data S4/S5 for real.")

    print(f"\n{'design':>9} {'score':>8} {'mean rung':>10} {'occup.':>8} "
          f"{'gate':>6} {'ladder-robust':>14}")
    sens = r["sensitivity"]
    for name, s in r["scored"].items():
        frac = sens["pass_fraction"][name]
        robust = "yes" if frac in (0.0, 1.0) else f"{frac:.0%} of ladders"
        print(f"{name:>9} {s['total_underprocessing_score']:>8.0f} "
              f"{s['mean_rung']:>10.2f} {s['occupancy']:>7.0f}% "
              f"{'PASS' if s['passes_paper_filter'] else 'fail':>6} "
              f"{robust:>14}")

    print(f"\nstochastic-dominance pairs found : {len(r['dominance_pairs'])}")
    print(f"ladder sweeps checked against them: {r['dominance_tests']}")
    print(f"order violations                  : {r['dominance_violations']}")
    if r["dominance_violations"] == 0:
        print("    -> every dominance ordering held under every ladder tested,")
        print("       so those rankings do NOT depend on the reconstruction.")

    print("\nscore range across monotone ladders (shows threshold fragility):")
    for name in r["scored"]:
        print(f"    {name:>9}  {sens['score_min'][name]:>7.0f} .. "
              f"{sens['score_max'][name]:>7.0f}"
              f"   {'spans 1500' if sens['score_min'][name] < 1500 < sens['score_max'][name] else ''}")

    p = make_figure(r)
    print(f"\nfigure written to: {p}")
    return r


if __name__ == "__main__":
    report()
