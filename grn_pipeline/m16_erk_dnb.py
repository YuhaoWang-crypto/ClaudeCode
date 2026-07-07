"""
Module 16 — Critical slowing down & DNB on the REAL ERK switch.

Connects M15 (exact Markevich MM ERK bistable switch) back to the M4/M9
early-warning machinery, but now at the genuine saddle-node boundaries
[MAPKK] = 39.25 / 57.38 nM instead of a toy fold.

For each stable branch we compute:
  * leading Jacobian eigenvalue lambda_max  (= largest Lyapunov exponent at a
    fixed point) -> 0^- as MAPKK approaches its saddle-node boundary;
  * recovery time  tau = -1/lambda_max  -> diverges (critical slowing down);
  * from a Langevin (Euler-Maruyama) simulation around the stable state:
    fluctuation SD, lag-1 autocorrelation, and a DNB index -> all rise near
    the boundary.

Reproduces the report's MAPK critical-slowing/DNB table (Table 9): near a
boundary lambda_max ~ -2e-4, tau ~ 4700 s, DNB proxy peaks.
"""
import numpy as np
from grn_pipeline import m15_markevich_mm as mk

MTOT = mk.MTOT


def stable_branches(K):
    ss = mk.steady_states(K, n=700)
    stab = [r for r in ss if r["stable"]]
    off = min(stab, key=lambda r: r["Mpp"]) if stab else None
    on = max(stab, key=lambda r: r["Mpp"]) if stab else None
    return off, on


def slowing_table():
    """lambda_max and tau along each branch across the bistable window."""
    rows = []
    for K in [35, 39.3, 42, 46, 50, 54, 57.3, 60]:
        off, on = stable_branches(K)
        for label, st in (("OFF/low", off), ("ON/high", on)):
            if st is None:
                continue
            lam = st["lmax"]
            tau = -1.0 / lam if lam < 0 else np.inf
            rows.append({"K": K, "branch": label, "Mpp": st["Mpp"],
                         "lam": lam, "tau": tau})
    return rows


def langevin_dnb(K, branch="OFF/low", sigma=2.0, T=4000.0, dt=0.02, seed=1):
    """Euler-Maruyama around a stable ERK state; return SD, lag-1 autocorr,
    cross-correlation and a DNB index."""
    off, on = stable_branches(K)
    st = off if branch.startswith("OFF") else on
    if st is None:
        return None
    x = np.array([st["M"], st["Mpp"]], float)
    rng = np.random.default_rng(seed)
    n = int(T / dt); sq = np.sqrt(dt)
    rec_every = int(5.0 / dt)               # sample every 5 s so lag-1
    Ms, Ps = [], []                          # autocorrelation is meaningful
    for i in range(n):
        drift = mk.rhs(x, K)
        x = x + dt * drift + sigma * sq * rng.normal(size=2)
        x[0] = min(max(x[0], 0.0), MTOT)
        x[1] = min(max(x[1], 0.0), MTOT)
        if i % rec_every == 0:
            Ms.append(x[0]); Ps.append(x[1])
    Ms, Ps = np.array(Ms), np.array(Ps)
    burn = len(Ms) // 4
    Ms, Ps = Ms[burn:], Ps[burn:]
    sd = 0.5 * (Ms.std() + Ps.std())
    ac = np.corrcoef(Ps[:-1], Ps[1:])[0, 1]
    corr = abs(np.corrcoef(Ms, Ps)[0, 1])
    dnb = sd * corr
    return {"sd": sd, "autocorr": ac, "corr": corr, "dnb": dnb}


def report(fig_path="figures/m16_erk_dnb.png"):
    print("=" * 68)
    print("MODULE 16 — Critical slowing down & DNB on the real ERK switch")
    print("=" * 68)
    print("saddle-node boundaries (M15): MAPKK = 39.25 and 57.38 nM\n")

    rows = slowing_table()
    print(f"{'MAPKK':>6} {'branch':>8} {'Mpp':>8} {'lambda_max':>11} "
          f"{'tau (s)':>10}")
    for r in rows:
        print(f"{r['K']:6.1f} {r['branch']:>8} {r['Mpp']:8.1f} "
              f"{r['lam']:11.5f} {r['tau']:10.1f}")
    print("\n-> near 39.3 nM the ON/high branch has lambda_max -> 0^- "
          "(tau huge);")
    print("   near 57.3 nM the OFF/low branch has lambda_max -> 0^- "
          "(critical slowing).")
    print("   (report Table 9: lambda_max ~ -2.1e-4, tau ~ 4720 s at the edges)")

    # DNB along OFF branch approaching the upper boundary 57.38
    print(f"\nDNB early-warning along OFF/low branch -> upper boundary "
          f"57.38 nM:")
    print(f"  {'MAPKK':>6} {'Mpp':>7} {'SD':>7} {'autocorr':>9} {'DNB':>8}")
    dnb_track = []
    for K in [50, 53, 55, 56, 57]:
        d = langevin_dnb(K, "OFF/low")
        if d:
            off, _ = stable_branches(K)
            print(f"  {K:6.1f} {off['Mpp']:7.1f} {d['sd']:7.2f} "
                  f"{d['autocorr']:9.3f} {d['dnb']:8.2f}")
            dnb_track.append((K, d))
    print("  -> SD, lag-1 autocorrelation and DNB all climb as the OFF state")
    print("     nears its saddle-node -> model-free early warning of the")
    print("     ERK off->on tipping (biomarker: pERK-DUSP time-series stats).")

    _figure(rows, dnb_track, fig_path)
    print(f"\nfigure written to: {fig_path}")
    return {"slowing": rows, "dnb_track": dnb_track}


def _figure(rows, dnb_track, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))

    for branch, col in (("OFF/low", "#2980b9"), ("ON/high", "#c0392b")):
        pts = [(r["K"], r["lam"]) for r in rows if r["branch"] == branch]
        if pts:
            xs, ys = zip(*sorted(pts))
            ax[0].plot(xs, ys, "o-", color=col, label=branch)
    ax[0].axhline(0, ls="--", color="k", lw=1)
    for b in (39.25, 57.38):
        ax[0].axvline(b, ls=":", color="gray")
    ax[0].set_xlabel("[MAPKK] (nM)"); ax[0].set_ylabel("lambda_max (1/s)")
    ax[0].set_title("Leading eigenvalue (=LLE) -> 0\nat real saddle-nodes")
    ax[0].legend(fontsize=8)

    for branch, col in (("OFF/low", "#2980b9"), ("ON/high", "#c0392b")):
        pts = [(r["K"], min(r["tau"], 1e4)) for r in rows
               if r["branch"] == branch and np.isfinite(r["tau"])]
        if pts:
            xs, ys = zip(*sorted(pts))
            ax[1].semilogy(xs, ys, "o-", color=col, label=branch)
    for b in (39.25, 57.38):
        ax[1].axvline(b, ls=":", color="gray")
    ax[1].set_xlabel("[MAPKK] (nM)"); ax[1].set_ylabel("recovery time tau (s)")
    ax[1].set_title("Critical slowing down\ntau = -1/lambda_max diverges")
    ax[1].legend(fontsize=8)

    if dnb_track:
        ks = [k for k, _ in dnb_track]
        ax[2].plot(ks, [d["dnb"] for _, d in dnb_track], "o-",
                   color="#8e44ad", label="DNB index")
        ax[2].plot(ks, [d["sd"] for _, d in dnb_track], "s--",
                   color="#16a085", label="SD")
        ax[2].axvline(57.38, ls=":", color="gray", label="boundary 57.38")
        ax[2].set_xlabel("[MAPKK] (nM)  (OFF branch -> boundary)")
        ax[2].set_title("DNB early-warning rises\ntoward the tipping point")
        ax[2].legend(fontsize=8)
    fig.suptitle("Real ERK switch (Markevich): DNB / critical slowing / "
                 "Lyapunov at genuine saddle-node boundaries", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(path, dpi=130); plt.close(fig)


if __name__ == "__main__":
    report()
