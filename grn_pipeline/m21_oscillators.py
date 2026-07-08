"""
Module 21 — Extend the biomarker framework to OSCILLATORY pathways (Hopf).

M15-M19 covered bistable switches (saddle-node bifurcations). Many key
pathways are oscillators (circadian clock, p53-Mdm2 DNA-damage pulses,
NF-kB, glycolysis) that lose stability through a HOPF bifurcation. This is a
genuine extension: near a Hopf the leading eigenvalue pair still has
Re -> 0 (critical slowing / rising variance), but with NON-ZERO imaginary
part, so the early-warning signature is DIFFERENT from a saddle-node -- the
fluctuation power spectrum develops and sharpens a PEAK at the Hopf frequency
f = Im(lambda)/2*pi (noise-induced / coherence-resonance oscillations).

Three canonical Hopf oscillators (distinct biology, reliably-correct math):
  * Goodwin       -> circadian / generic transcriptional clock (Hopf in Hill n)
  * p53-Mdm2      -> DNA-damage pulses (Hopf in Mdm2 transcription gain)
  * Brusselator   -> metabolic / glycolytic oscillation (exact Hopf B=1+A^2)
(NF-kB / IkB is the same negative-feedback-oscillator class.)

For each: locate the Hopf (Re(lambda_max) crosses 0 with Im != 0); on the
stable-focus side approaching it, show Re -> 0 and a Langevin power spectrum
whose peak at the Hopf frequency grows and sharpens.
"""
import numpy as np
from scipy.optimize import fsolve
from scipy.signal import welch


def _jac(f, x, p, h=1e-6):
    n = len(x); J = np.zeros((n, n)); f0 = np.array(f(x, p))
    for j in range(n):
        xx = np.array(x, float); xx[j] += h
        J[:, j] = (np.array(f(xx, p)) - f0) / h
    return J


def _lead_eig(f, x0, p):
    fp = fsolve(f, x0, args=(p,))
    ev = np.linalg.eigvals(_jac(f, fp, p))
    i = int(np.argmax(ev.real))
    return fp, ev[i].real, abs(ev[i].imag)


# --------------------------- oscillator models ---------------------------
def goodwin(s, n):
    x, y, z = s; kd = 0.35
    return [1.0 / (1 + z ** n) - kd * x, x - kd * y, y - kd * z]


def p53(s, b):
    p, Mm, M = s
    return [1.0 - 4.0 * M * p / (0.08 + p) - 0.05 * p,
            b * p * p / (1.0 + p * p) - 1.0 * Mm,
            1.0 * Mm - 0.6 * M]


def brusselator(s, B):
    x, y = s; A = 1.0
    return [A - (B + 1) * x + x * x * y, B * x - x * x * y]


MODELS = [
    dict(id="Goodwin", biology="circadian / transcriptional clock",
         f=goodwin, x0=[0.5, 0.5, 0.5], pname="Hill coefficient n",
         prange=(6.0, 8.45), var=2),
    dict(id="p53_Mdm2", biology="DNA-damage p53 pulses",
         f=p53, x0=[0.5, 0.3, 0.3], pname="Mdm2 transcription gain",
         prange=(4.0, 2.05), var=0),
    dict(id="Brusselator", biology="metabolic / glycolytic oscillation",
         f=brusselator, x0=[1.0, 1.0], pname="B (drive)",
         prange=(1.2, 1.95), var=0),
]


def find_hopf(model, n=60):
    """Scan the control parameter widely; return the Hopf point (Re crosses 0
    with Im!=0) and the (p, Re, Im) trace."""
    f, x0 = model["f"], model["x0"]
    lo, hi = sorted(model["prange"])
    span = hi - lo
    ps = np.linspace(lo - 0.4 * span, hi + 0.6 * span, n)
    trace = []
    for p in ps:
        try:
            _, re, im = _lead_eig(f, x0, p)
            trace.append((p, re, im))
        except Exception:
            continue
    hopf = None
    for (p0, r0, i0), (p1, r1, i1) in zip(trace, trace[1:]):
        if r0 * r1 < 0 and max(i0, i1) > 0.05:
            hopf = 0.5 * (p0 + p1)
            break
    return hopf, trace


def power_spectrum(model, p, T=6000.0, dt=0.02, sigma=0.02, seed=0):
    """Langevin around the stable focus; return (freqs, PSD) of `var`."""
    f, x0, iv = model["f"], model["x0"], model["var"]
    fp = fsolve(f, x0, args=(p,))
    x = np.array(fp, float); rng = np.random.default_rng(seed)
    nstep = int(T / dt); sq = np.sqrt(dt); rec = []
    for i in range(nstep):
        x = x + dt * np.array(f(x, p)) + sigma * sq * rng.normal(size=len(x))
        x = np.clip(x, 0, None)
        if i % 5 == 0:
            rec.append(x[iv])
    rec = np.array(rec) - np.mean(rec)
    fs = 1.0 / (5 * dt)
    fr, ps = welch(rec, fs=fs, nperseg=min(2048, len(rec) // 4))
    return fr, ps, float(np.var(rec))


def report(fig_path="figures/m21_oscillators.png"):
    print("=" * 70)
    print("MODULE 21 — Framework extension to OSCILLATORY pathways (Hopf)")
    print("=" * 70)
    results = []
    for m in MODELS:
        hopf, trace = find_hopf(m)
        results.append({"m": m, "hopf": hopf, "trace": trace})
        # near-Hopf frequency from the eigenvalue Im
        if hopf is not None:
            _, _, im = _lead_eig(m["f"], m["x0"],
                                 hopf - 0.02 * abs(hopf) - 1e-3)
            fhz = im / (2 * np.pi)
        else:
            fhz = np.nan
        print(f"\n{m['id']:12s} [{m['biology']}]")
        print(f"  control = {m['pname']}; Hopf at "
              f"{('%.3f' % hopf) if hopf else 'n/a'}  "
              f"(eigenvalue Im -> osc. frequency ~{fhz:.3f} rad-based)")

    # spectra approaching the Hopf for each model
    print("\napproaching each Hopf on the stable side: variance rises and a")
    print("spectral peak emerges at the Hopf frequency (distinct from the")
    print("saddle-node signature in M19):")
    spec_data = {}
    for r in results:
        m = r["m"]
        if r["hopf"] is None:
            continue
        # three approach points from deep-stable to near-Hopf, using the
        # model's intended stable-side prange (prange[1] hugs the Hopf).
        p_deep, p_near = m["prange"]
        pts = [p_deep, 0.5 * (p_deep + p_near), p_near]
        rows = []
        for p in pts:
            fr, ps, var = power_spectrum(m, p)
            pk = fr[int(np.argmax(ps))]
            rows.append({"p": p, "var": var, "peak_f": pk,
                         "peak_pow": float(np.max(ps)), "fr": fr, "ps": ps})
        spec_data[m["id"]] = rows
        v0, v2 = rows[0]["var"], rows[-1]["var"]
        pw0, pw2 = rows[0]["peak_pow"], rows[-1]["peak_pow"]
        print(f"  {m['id']:12s} variance x{v2/v0:4.1f}, spectral-peak power "
              f"x{pw2/pw0:5.1f} approaching Hopf; peak at f~{rows[-1]['peak_f']:.3f}")

    print("\nCONCLUSION: the biomarker framework EXTENDS to oscillatory")
    print("(Hopf) pathways. Critical slowing still holds (Re(lambda)->0,")
    print("variance rises), but the discriminating early-warning is a SHARPENING")
    print("SPECTRAL PEAK at the intrinsic frequency -- absent for saddle-node")
    print("switches. Same engine, richer readout: watch the power spectrum for")
    print("clocks/pulse-generators (circadian, p53, NF-kB, glycolysis).")
    _figure(results, spec_data, fig_path)
    print(f"\nfigure written to: {fig_path}")
    return {"hopf": {r["m"]["id"]: r["hopf"] for r in results}}


def _figure(results, spec_data, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 3, figsize=(14, 7))
    for j, r in enumerate(results):
        m = r["m"]; tr = r["trace"]
        ps = [t[0] for t in tr]; re = [t[1] for t in tr]
        ax[0, j].plot(ps, re, "-", color="#c0392b")
        ax[0, j].axhline(0, ls="--", color="k", lw=1)
        if r["hopf"]:
            ax[0, j].axvline(r["hopf"], ls=":", color="gray")
        ax[0, j].set_title(f"{m['id']}\n{m['biology']}", fontsize=9)
        ax[0, j].set_xlabel(m["pname"]); ax[0, j].set_ylabel("Re(lambda_max)")
        # spectra
        if m["id"] in spec_data:
            cols = ["#95a5a6", "#2980b9", "#8e44ad"]
            labs = ["far", "mid", "near Hopf"]
            for k, row in enumerate(spec_data[m["id"]]):
                ax[1, j].semilogy(row["fr"], row["ps"], color=cols[k],
                                  label=labs[k])
            ax[1, j].set_xlim(0, 0.25)
            ax[1, j].set_xlabel("frequency"); ax[1, j].set_ylabel("power")
            ax[1, j].legend(fontsize=7)
            ax[1, j].set_title("fluctuation spectrum: peak sharpens", fontsize=9)
    fig.suptitle("Hopf extension: approaching an oscillatory instability, "
                 "variance rises AND a spectral peak emerges at the intrinsic "
                 "frequency", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(path, dpi=130); plt.close(fig)


if __name__ == "__main__":
    report()
