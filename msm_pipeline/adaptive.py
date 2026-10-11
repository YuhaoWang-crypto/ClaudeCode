"""
Seed-selection strategies for adaptive sampling, shared by the controlled
benchmark (m9) and the production OpenMM driver (adaptive_md.py).

The adaptive sampling loop is:

    round r:  run N short walkers  ->  append to the dataset
              re-discretise, re-estimate the MSM
              pick N seed STATES by some strategy
              pick an actual stored CONFIGURATION from each state
    round r+1: restart the walkers from those configurations

Everything hinges on which states you pick.  The strategies here:

  min_counts   pick the least-visited states.  Pure exploration; the classic
               choice when you do not yet know where the interesting states
               are.
  bottleneck   pick states whose forward committor q+ is near 0.5, i.e. the
               transition state.  Pure exploitation; the theoretically
               motivated choice for reducing the variance of a RATE estimate,
               because reactive trajectories are born there.
  hybrid       min_counts for the first half of the rounds, bottleneck after.
  random       uniform over visited states (control).
  population   inverse-population weighting (soft min_counts).

WHAT THE BENCHMARK ACTUALLY SHOWS (m9, 24 paired repeats, known ground truth):
adaptive reseeding of any kind beats a few long trajectories decisively, but
the four strategies are NOT statistically distinguishable from each other.
Do not over-engineer the seeding rule; do run rounds.
"""
from __future__ import annotations

import numpy as np

STRATEGIES = ("min_counts", "bottleneck", "hybrid", "random", "population")


def _committor_weights(msm, active, state_a, state_b, width=0.15):
    """Gaussian weight peaked at q+ = 0.5; None if A or B not yet sampled."""
    from deeptime.markov import reactive_flux
    act = list(active)
    if state_a not in act or state_b not in act:
        return None
    ia = act.index(state_a)
    ib = act.index(state_b)
    try:
        q = reactive_flux(msm.transition_matrix, [ia], [ib]).forward_committor
    except Exception:
        return None
    w = np.exp(-((q - 0.5) ** 2) / (2.0 * width ** 2))
    s = w.sum()
    return None if not np.isfinite(s) or s <= 0 else w / s


def select_seed_states(n_seeds, dtraj_all, rng, strategy="min_counts",
                       msm=None, active=None, state_a=None, state_b=None,
                       round_index=0, n_rounds=1, switch=0.5):
    """
    -> array of n_seeds state indices to restart from.

    dtraj_all : 1-D array of discrete states, all data so far
    msm/active: current MSM and its active-set state labels (for 'bottleneck')
    state_a/b : the source / sink states defining the committor

    Falls back to min_counts whenever the committor is unavailable (e.g. the
    sink has not been discovered yet), which is the behaviour you want: you
    cannot exploit a barrier you have not found.
    """
    if strategy not in STRATEGIES:
        raise ValueError(f"unknown strategy {strategy!r}; pick from {STRATEGIES}")
    uniq, cnt = np.unique(dtraj_all, return_counts=True)

    eff = strategy
    if strategy == "hybrid":
        eff = "min_counts" if round_index < switch * n_rounds else "bottleneck"

    if eff == "bottleneck" and msm is not None and active is not None:
        w = _committor_weights(msm, active, state_a, state_b)
        if w is not None:
            return rng.choice(np.asarray(active), size=n_seeds, p=w, replace=True)
        eff = "min_counts"          # barrier not locatable yet -> explore

    if eff == "random":
        return rng.choice(uniq, size=n_seeds, replace=True)
    if eff == "population":
        w = 1.0 / cnt
        return rng.choice(uniq, size=n_seeds, p=w / w.sum(), replace=True)

    # min_counts: least-visited first, cycling if fewer states than seeds
    order = uniq[np.argsort(cnt)]
    return np.array([order[i % len(order)] for i in range(n_seeds)])


def pick_configurations(seed_states, dtraj_all, rng):
    """-> frame indices into the pooled dataset, one per seed state."""
    out = []
    for s in seed_states:
        idx = np.flatnonzero(dtraj_all == s)
        out.append(int(rng.choice(idx)) if len(idx)
                   else int(rng.integers(len(dtraj_all))))
    return np.array(out, dtype=np.int64)


def its_converged(series, tol=1.5):
    """
    series: list of slowest implied timescales, ordered by increasing lag.
    Converged if the value stops growing over the upper half of the lag range.
    """
    s = np.asarray(series, float)
    s = s[np.isfinite(s)]
    if len(s) < 3:
        return False, np.nan
    mid = s[len(s) // 2]
    drift = s[-1] / mid if mid > 0 else np.inf
    return bool(drift < tol), float(drift)
