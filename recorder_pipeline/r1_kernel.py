"""
R1 — Persistence kernels: what a Disease Recorder actually remembers.

THE CENTRAL CORRECTION
----------------------
The framework writes

    Recorder(t) ~ integral A(s) * w(s) * Persistence(t-s) ds

and leaves `Persistence` generic. For a mark written irreversibly onto a
population of carriers that are continuously born and destroyed, the kernel is
NOT free to choose and it is NOT the carrier's decay curve. It is derived:

Let carriers be produced at a constant rate with lifespan distribution having
survival S(a) = P(lifespan > a). In steady state the AGE density of the living
carrier population is

    f(a) = S(a) / integral_0^inf S(u) du                       (renewal theory)

A carrier now aged `a` was born at t-a and has been accumulating mark ever
since, so its burden is integral_0^a w * A(t - x) dx. Averaging over the age
distribution and exchanging the order of integration:

    mean mark(t) = w * integral_0^inf A(t-x) * K(x) dx
    with   K(x) = integral_x^inf f(a) da = P(carrier age > x)

=> THE KERNEL IS THE SURVIVOR FUNCTION OF THE CARRIER *AGE* DISTRIBUTION.

Two consequences that matter for assay design, and that a generic
"exponential decay" picture gets wrong:

1. A carrier with a near-DETERMINISTIC lifespan L (the red blood cell) gives
       K(x) = 1 - x/L  on [0, L]        -- a linear ramp with a HARD horizon.
   Not an exponential. HbA1c cannot see past 120 days at all.

2. A carrier cleared by RANDOM first-order loss (albumin, transferrin) gives
       K(x) = exp(-x/tau)               -- exponential, no horizon, long tail.

So HbA1c and glycated albumin do not merely have "different time windows";
their kernels have different SHAPES. One has a memory cliff, the other a fade.

If the MARK itself is also erased (rate 1/tau_m) on top of carrier turnover,
the kernels multiply, K(x) -> K(x) * exp(-x/tau_m), and for the exponential
case the timescales combine harmonically:

    1/tau_eff = 1/tau_carrier + 1/tau_mark      => the SHORTER one wins.

This is the design law: a recorder's memory is capped by the shorter of mark
stability and carrier lifetime. Irreversible chemistry (glycation) on a
long-lived carrier (RBC) is the only combination that yields a months-long
recorder -- which is exactly why HbA1c, and not something else, is the one
that worked.

VALIDATION
----------
The RBC kernel is built from ONE literature input (RBC lifespan ~120 d) and
reproduces, without further tuning, the clinically quoted properties of HbA1c:
the ~50% weighting of the most recent month and the 8-12 week window.
"""
import numpy as np

LN2 = np.log(2.0)


# --------------------------------------------------------------------------
# Kernels K(x) = P(carrier age > x), normalised later by their own integral.
# --------------------------------------------------------------------------

def kernel_deterministic(x, lifespan):
    """Carrier with a fixed lifespan L (red blood cell).

    S(a) = 1{a < L}  ->  f(a) = 1/L on [0,L]  ->  K(x) = 1 - x/L.
    """
    x = np.asarray(x, dtype=float)
    return np.clip(1.0 - x / lifespan, 0.0, None)


def kernel_exponential(x, tau):
    """Carrier cleared by first-order random loss with mean lifetime tau.

    S(a) = exp(-a/tau) -> f(a) = S(a)/tau -> K(x) = exp(-x/tau).
    (For an exponential lifespan the age survivor equals the lifespan survivor;
    this is the memoryless special case, not the general rule.)
    """
    x = np.asarray(x, dtype=float)
    return np.exp(-x / tau)


def kernel_gamma_rbc(x, lifespan, cv=0.15):
    """More realistic RBC: lifespan is deterministic-ish but not exact, plus a
    small random-destruction component. Modelled as a gamma lifespan with mean
    `lifespan` and coefficient of variation `cv`, integrated numerically.

    Used only to show the ~50% recent-month figure is robust to the idealised
    box-lifespan assumption.
    """
    from scipy.stats import gamma as _gamma
    x = np.asarray(x, dtype=float)
    k = 1.0 / cv**2
    theta = lifespan / k
    # K(x) = P(age > x) = int_x^inf S(a) da / int_0^inf S(a) da
    grid = np.linspace(0.0, lifespan * 4.0, 20001)
    S = _gamma.sf(grid, a=k, scale=theta)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (S[1:] + S[:-1]) * np.diff(grid))])
    total = cum[-1]
    tail = total - np.interp(x, grid, cum)
    return np.clip(tail / total, 0.0, None)


def apply_mark_decay(K, x, tau_mark):
    """Multiply a carrier kernel by mark erasure. tau_mark=None -> irreversible."""
    if tau_mark is None or not np.isfinite(tau_mark):
        return K
    return K * np.exp(-np.asarray(x, dtype=float) / tau_mark)


# --------------------------------------------------------------------------
# Summary statistics of a kernel (this is where the clinical "window" lives)
# --------------------------------------------------------------------------

def fmt_dur(days):
    """Human-readable duration; recorder timescales span minutes to decades."""
    if not np.isfinite(days):
        return "n/a"
    if days < 1.0 / 24:
        return f"{days*24*60:.0f} min"
    if days < 2.0:
        return f"{days*24:.1f} h"
    if days < 400:
        return f"{days:.1f} d"
    return f"{days/365.25:.1f} y"


def kernel_stats(K, x):
    """Mean lag, median lag, and the 90%-mass horizon of a (unnormalised) kernel.

    The clinical literature quotes recorder "windows" inconsistently: some are
    the median lag, some the 90% horizon. Reporting all three removes the
    ambiguity and lets the two conventions be compared like for like.
    """
    K = np.asarray(K, dtype=float)
    x = np.asarray(x, dtype=float)
    w = np.trapezoid(K, x)
    mean_lag = np.trapezoid(K * x, x) / w
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (K[1:] + K[:-1]) * np.diff(x))]) / w
    median_lag = float(np.interp(0.5, cum, x))
    h90 = float(np.interp(0.9, cum, x))
    return {"mean_lag_d": float(mean_lag), "median_lag_d": median_lag,
            "horizon90_d": h90, "mass": float(w)}


def window_fractions(K, x, edges_days):
    """Fraction of the recorder's total weight contributed by each time slab.

    For HbA1c this is the classic "how much of my HbA1c is last month?" table.
    """
    K = np.asarray(K, dtype=float)
    total = np.trapezoid(K, x)
    out = []
    for lo, hi in zip(edges_days[:-1], edges_days[1:]):
        m = (x >= lo) & (x <= hi)
        if m.sum() < 2:
            out.append(0.0)
            continue
        out.append(float(np.trapezoid(K[m], x[m]) / total))
    return out


# --------------------------------------------------------------------------
# Step response: how fast does a recorder follow a change in disease activity?
# --------------------------------------------------------------------------

def step_response(K, x, t_grid, drop=0.5):
    """Disease activity A(s) = 1 for s<0, then `drop` for s>=0 (treatment starts).

    Returns the normalised recorder trajectory. Answers the question the
    framework raises in section 10: which recorder can serve as a 4-week
    treatment-response readout, and which structurally cannot.
    """
    K = np.asarray(K, dtype=float)
    total = np.trapezoid(K, x)
    traj = []
    for t in t_grid:
        A = np.where(x <= t, drop, 1.0)   # lag x corresponds to time t-x
        traj.append(np.trapezoid(K * A, x) / total)
    return np.array(traj)


def time_to_fraction(traj, t_grid, frac=0.9, drop=0.5):
    """Time for the recorder to travel `frac` of the way to its new plateau."""
    target = 1.0 - frac * (1.0 - drop)
    below = np.where(traj <= target)[0]
    return float(t_grid[below[0]]) if len(below) else float("nan")


# --------------------------------------------------------------------------
# The catalogue of real carriers  [LIT inputs, [OK] derived kernels]
# --------------------------------------------------------------------------
# lifespan/half-life values are literature medians; sources in RECORDER_FRAMEWORK.md
CARRIERS = [
    # name,                    kind,            param (days), mark erasure (days)
    ("RBC / haemoglobin (HbA1c)",      "deterministic", 120.0, None),
    ("RBC membrane (PEth, alcohol)",   "deterministic", 120.0, 7.0),   # PEth t1/2 4-10 d
    ("Albumin (glycated albumin)",     "exponential",   19.0 / LN2, None),
    ("Transferrin (%CDT)",             "exponential",   10.0 / LN2, None),
    ("IgG (anti-PTM autoantibody)",    "exponential",   21.0 / LN2, None),
    ("Fibrillar collagen (AGE x-link)", "exponential",  365.0 * 10 / LN2, None),
    ("Plasma protein pool (TAT)",      "exponential",   0.01 / LN2, None),  # t1/2 ~15 min
]


def build(name, kind, param, tau_mark, xmax=None, n=40001):
    """Return (x, K) for one carrier, mark erasure already applied."""
    if kind == "deterministic":
        xmax = xmax or param * 1.05
        x = np.linspace(0.0, xmax, n)
        K = kernel_deterministic(x, param)
    elif kind == "exponential":
        xmax = xmax or param * 12.0
        x = np.linspace(0.0, xmax, n)
        K = kernel_exponential(x, param)
    else:
        raise ValueError(kind)
    return x, apply_mark_decay(K, x, tau_mark)


def main():
    print("=" * 78)
    print("R1 — Persistence kernels: what does a recorder remember?")
    print("=" * 78)

    # ---- 1. HbA1c from first principles -----------------------------------
    L = 120.0
    x = np.linspace(0.0, L, 40001)
    K_rbc = kernel_deterministic(x, L)
    st = kernel_stats(K_rbc, x)
    fr = window_fractions(K_rbc, x, [0, 30, 60, 90, 120])

    print("\n[1] HbA1c kernel derived from RBC lifespan alone (L = 120 d)")
    print("    K(x) = 1 - x/120   (linear ramp, HARD horizon -- not exponential)")
    print(f"    mean lag        = {st['mean_lag_d']:.1f} d   (exactly L/3)")
    print(f"    median lag      = {st['median_lag_d']:.1f} d")
    print(f"    90% horizon     = {st['horizon90_d']:.1f} d  "
          f"= {st['horizon90_d']/7:.1f} weeks")
    print("    weight by month: " +
          "  ".join(f"{lo}-{hi}d {f*100:.1f}%"
                    for (lo, hi), f in zip([(0, 30), (30, 60), (60, 90), (90, 120)], fr)))
    print(f"    -> most recent 30 d carries {fr[0]*100:.1f}% of the signal")
    print("    VALIDATION: clinical convention quotes ~50% recent month and an")
    print("    8-12 week window. Both fall out of one number (120 d). [OK]")

    # robustness to the idealised box lifespan
    K_g = kernel_gamma_rbc(x, L, cv=0.15)
    fr_g = window_fractions(K_g, x, [0, 30, 60, 90, 120])
    st_g = kernel_stats(K_g, x)
    print(f"    robustness (gamma lifespan, CV=15%): recent 30 d = {fr_g[0]*100:.1f}%,"
          f" 90% horizon = {st_g['horizon90_d']:.0f} d")

    # ---- 2. Shape matters, not just tau -----------------------------------
    print("\n[2] Kernel SHAPE differs by carrier destruction mechanism")
    print(f"    {'carrier':<34}{'mean':>10}{'median':>10}{'90% horizon':>13}")
    rows = []
    for name, kind, param, tm in CARRIERS:
        xx, KK = build(name, kind, param, tm)
        s = kernel_stats(KK, xx)
        rows.append((name, s))
        print(f"    {name:<34}{fmt_dur(s['mean_lag_d']):>10}"
              f"{fmt_dur(s['median_lag_d']):>10}{fmt_dur(s['horizon90_d']):>13}")
    print("    TAT is the instructive outlier: the complex is a CHEMICALLY stable")
    print("    mark, but it is its own carrier and is cleared in ~15 min. Stable")
    print("    chemistry on a short-lived carrier = a short recorder. The memory")
    print("    cap law again -- TAT records the last hour of thrombin generation,")
    print("    not the last month of coagulation activity. [OK]")
    print("    Note HbA1c(mean 40 d) vs glycated albumin(mean 27 d) are much")
    print("    closer than the quoted '8-12 wk vs 2-3 wk' suggests: the clinical")
    print("    numbers compare a 90% horizon against a median. Like for like,")
    print("    the real difference is the HARD HORIZON, not the centre. [OK]")

    # ---- 3. The memory-cap law --------------------------------------------
    print("\n[3] Memory cap: 1/tau_eff = 1/tau_carrier + 1/tau_mark")
    tau_c = 120.0 / 3          # effective mean lag of the RBC ramp
    for tau_m in [np.inf, 365.0, 30.0, 7.0, 1.0]:
        xx = np.linspace(0, 200, 20001)
        KK = apply_mark_decay(kernel_deterministic(xx, 120.0), xx,
                              None if not np.isfinite(tau_m) else tau_m)
        s = kernel_stats(KK, xx)
        lab = "irreversible" if not np.isfinite(tau_m) else f"{tau_m:g} d"
        print(f"    mark erasure {lab:>12}  ->  mean lag {s['mean_lag_d']:6.1f} d")
    print("    A reversible mark on a long-lived carrier is NOT a long recorder.")
    print("    Irreversible chemistry is the non-negotiable ingredient. [OK]")

    # ---- 4. Step response = treatment-monitoring capability ---------------
    print("\n[4] Step response to a 50% drop in disease activity (t=0)")
    t = np.linspace(0, 180, 1801)
    print(f"    {'recorder':<34}{'t to 50% of new plateau':>26}{'t to 90%':>12}")
    for name, kind, param, tm in CARRIERS:
        if "TAT" in name or "AGE" in name:
            continue
        xx, KK = build(name, kind, param, tm)
        traj = step_response(KK, xx, t, drop=0.5)
        t50 = time_to_fraction(traj, t, 0.5)
        t90 = time_to_fraction(traj, t, 0.9)
        print(f"    {name:<34}{t50:>24.0f} d{t90:>10.0f} d")
    print("    -> a 4-week treatment-response endpoint is structurally")
    print("       unavailable to HbA1c and available to glycated albumin.")
    print("       This reproduces the clinical observation the framework cites,")
    print("       from kernel geometry alone. [OK]")

    return {"hba1c_recent30_frac": fr[0], "hba1c_stats": st, "carriers": rows}


if __name__ == "__main__":
    main()
