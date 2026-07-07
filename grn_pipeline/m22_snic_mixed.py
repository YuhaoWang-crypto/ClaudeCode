"""
Module 22 — Mixed bifurcation: Saddle-Node on an Invariant Circle (SNIC).

M19 covered saddle-node (bistable switch, no oscillation); M21 covered Hopf
(oscillation, finite frequency, amplitude -> 0). Some pathways -- notably cell-
cycle models (Tyson-Novak) and type-I excitable systems -- pass through a SNIC
(a.k.a. SNIPER): a saddle-node that sits ON a limit cycle. It combines both:
below SNIC there is a stable rest state (excitable / arrested); above it a
FINITE-amplitude oscillation whose period diverges (frequency -> 0) at the
bifurcation.

Distinctive early-warning signature (this is the point):
  * saddle-node : variance up, autocorrelation up, NO oscillation.
  * Hopf        : variance up, spectral peak at a FIXED frequency (amp -> 0).
  * SNIC        : FINITE-amplitude spikes whose PERIOD -> infinity
                  (frequency -> 0); the spectral peak slides toward 0 Hz and
                  the inter-spike interval and its variability blow up.

Canonical SNIC normal form (theta / Ermentrout-Kopell model):
    dtheta/dt = (1 - cos theta) + (1 + cos theta) * I
  I < 0 : two fixed points on the circle (stable node + saddle) -> rest.
  I = 0 : SNIC.
  I > 0 : rotation (spiking) with period T ~ pi / sqrt(I) -> infinity.
Observable v = (1 - cos theta)/2 spikes once per rotation.
Biology: near the cycling/arrest boundary the cell cycle slows continuously
(long, variable periods) -- the SNIC hallmark, unlike a Hopf clock.
"""
import numpy as np


def dtheta(theta, I):
    return (1 - np.cos(theta)) + (1 + np.cos(theta)) * I


def integrate(I, T=4000.0, dt=0.01, sigma=0.0, seed=0):
    n = int(T / dt); th = 0.1; rng = np.random.default_rng(seed)
    sq = np.sqrt(dt); rec = np.empty(n)
    for i in range(n):
        th = th + dt * dtheta(th, I) + sigma * sq * rng.normal()
        rec[i] = th
    return rec


def period(I, T=6000.0, dt=0.005):
    """Deterministic spike period from theta crossing multiples of 2*pi."""
    if I <= 0:
        return np.inf                       # rest state: no oscillation
    th = integrate(I, T=T, dt=dt)
    crossings = np.where(np.diff(np.floor(th / (2 * np.pi))) >= 1)[0]
    if len(crossings) < 3:
        return np.inf
    times = crossings * dt
    return float(np.median(np.diff(times)))


def spikes_isi(I, sigma, T=8000.0, dt=0.01, seed=1):
    """Noisy inter-spike intervals near SNIC (on the oscillating side)."""
    th = integrate(I, T=T, dt=dt, sigma=sigma, seed=seed)
    v = (1 - np.cos(th)) / 2
    # spike = crossing 2*pi boundary
    idx = np.where(np.diff(np.floor(th / (2 * np.pi))) >= 1)[0]
    isi = np.diff(idx * dt)
    return v, isi


def report(fig_path="figures/m22_snic.png"):
    print("=" * 70)
    print("MODULE 22 — SNIC mixed bifurcation (saddle-node ON a limit cycle)")
    print("=" * 70)
    print("theta model dθ/dt=(1-cosθ)+(1+cosθ)I ; SNIC at I=0\n")

    print("period vs drive I (approaching SNIC from the oscillating side):")
    print(f"  {'I':>8} {'period T':>10} {'frequency':>10}  regime")
    Is = [0.6, 0.3, 0.15, 0.08, 0.04, 0.02, 0.01, -0.05]
    per = []
    for I in Is:
        T = period(I)
        f = 1.0 / T if np.isfinite(T) else 0.0
        per.append((I, T, f))
        reg = "rest (excitable/arrest)" if I <= 0 else "spiking"
        print(f"  {I:8.3f} {T:10.1f} {f:10.4f}  {reg}")
    print("  -> as I -> 0+, period -> infinity and frequency -> 0 while spike")
    print("     amplitude stays finite: the SNIC hallmark (T ~ pi/sqrt(I)).")

    # power law check T ~ 1/sqrt(I)
    osc = [(I, T) for I, T, f in per if I > 0 and np.isfinite(T)]
    Iv = np.array([x[0] for x in osc]); Tv = np.array([x[1] for x in osc])
    slope = np.polyfit(np.log(Iv), np.log(Tv), 1)[0]
    print(f"\n  log-log slope of T vs I = {slope:.2f}  (SNIC predicts -0.5)")

    print("\nnoisy inter-spike intervals approaching SNIC:")
    print(f"  {'I':>6} {'mean ISI':>9} {'ISI CV':>7}")
    isi_rows = []
    for I in [0.3, 0.1, 0.03, 0.012]:
        v, isi = spikes_isi(I, sigma=0.05)
        if len(isi) > 3:
            m, cv = np.mean(isi), np.std(isi) / np.mean(isi)
            isi_rows.append((I, m, cv, v))
            print(f"  {I:6.3f} {m:9.1f} {cv:7.3f}")
    print("  -> mean ISI and its coefficient of variation both grow (period")
    print("     lengthening + rising variability) = the SNIC early warning.")

    print("\nBIFURCATION-TYPE CONTRAST (unified biomarker framework):")
    print("  saddle-node (switch) : variance^, autocorr^, no oscillation")
    print("  Hopf (clock)         : variance^, spectral peak at FIXED freq")
    print("  SNIC (excitable/cycle): finite-amp spikes, period->inf (freq->0)")
    print("  -> read variance+autocorr for switches, the SPECTRUM for clocks,")
    print("     and the PERIOD/ISI trend for SNIC (e.g. slow-cycling cells near")
    print("     the cell-cycle arrest boundary).")
    _figure(per, isi_rows, slope, fig_path)
    print(f"\nfigure written to: {fig_path}")
    return {"snic": 0.0, "slope": slope}


def _figure(per, isi_rows, slope, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    osc = [(I, T) for I, T, f in per if I > 0 and np.isfinite(T)]
    Iv = [x[0] for x in osc]; Tv = [x[1] for x in osc]
    ax[0].loglog(Iv, Tv, "o-", color="#e67e22")
    ax[0].set_xlabel("drive I (-> SNIC at 0)"); ax[0].set_ylabel("period T")
    ax[0].set_title(f"Period diverges at SNIC\nlog-log slope {slope:.2f} "
                    f"(theory -0.5)")

    # example spike trains at three I, vertically offset (amplitude fixed,
    # spacing grows as I -> SNIC)
    dt = 0.01
    for k, (I, col) in enumerate([(0.3, "#95a5a6"), (0.08, "#2980b9"),
                                  (0.02, "#8e44ad")]):
        th = integrate(I, T=350, dt=dt)
        v = (1 - np.cos(th)) / 2
        ax[1].plot(np.arange(len(v)) * dt, v + 1.25 * k, color=col, lw=0.9,
                   label=f"I={I}")
    ax[1].set_xlabel("time"); ax[1].set_yticks([])
    ax[1].set_title("Spikes: same amplitude,\nspacing grows toward SNIC")
    ax[1].legend(fontsize=7, loc="upper right")

    if isi_rows:
        Ii = [r[0] for r in isi_rows]
        ax[2].plot(Ii, [r[1] for r in isi_rows], "o-", color="#16a085",
                   label="mean ISI")
        axb = ax[2].twinx()
        axb.plot(Ii, [r[2] for r in isi_rows], "s--", color="#c0392b",
                 label="ISI CV")
        axb.set_ylabel("ISI CV", color="#c0392b")
        ax[2].set_xlabel("drive I"); ax[2].set_ylabel("mean ISI",
                                                       color="#16a085")
        ax[2].set_title("Noisy ISI: period + variability\nboth grow toward SNIC")
    fig.suptitle("SNIC: a mixed saddle-node + oscillation — finite-amplitude "
                 "spikes with diverging period (frequency -> 0)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(path, dpi=130); plt.close(fig)


if __name__ == "__main__":
    report()
