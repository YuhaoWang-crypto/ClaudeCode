"""
Module 23 — Chromatin memory as a continuous-time Markov chain (CTMC):
the A -> R -> I "three-state" model of Bintu et al., Science 2016.

Paper (exact anchor):
  Bintu, Yong, Antebi, McCue, Kazuki, Uno, Oshimura, Elowitz,
  "Dynamics of epigenetic regulation at the single-cell level",
  Science 351:720-724 (2016). doi:10.1126/science.aab2956 (PMC5108652).

What the paper established (and what we encode here):
  Recruiting a chromatin regulator (CR) to a reporter drives single cells,
  in all-or-none stochastic events, through three phenomenological states:

      A  (active / open)  --k_AR-->  R  (reversibly silent)  --k_RI-->  I
                             ^                                   (irreversibly
                             |   (after release of the CR)        silent)
                             +------------- k_RA ----------------- R

  * During recruitment: forward silencing A->R->I is on; its RATE is set by
    the CR identity and by recruitment STRENGTH (dox concentration).
  * After release: forward rates become negligible; R reverts to A at a rate
    k_RA that depends ONLY on the CR identity (NOT on dox) — a key, clean
    experimental result. I never reverts: it is an ABSORBING state.

  The four CRs map onto distinct memory classes:
    HDAC4  : short-term  -> no stable I (k_RI ~ 0); A<->R only; all reactivate.
    DNMT3B : long-term   -> commits to I (no reactivating R); A->I absorbing.
    EED/KRAB: HYBRID      -> full A->R->I; a fraction reverts (from R), the
                             rest is locked permanently (reached I).

This module turns that verbal model into THE core data structure of the whole
chromatin extension: a generator matrix Q (a flow / diffusion of probability
on the state graph, dp/dt = p Q), from which everything else is computed:

  (1) irreversible-trajectory classification  = recurrent (closed) communicating
      classes of the OFF-phase chain, found by graph SCC decomposition;
  (2) memory retention                         = permanent silent fraction after
      release  (= probability absorbed in I)   and its growth with recruitment
      duration (the paper's "duration-dependent fractional control");
  (3) extinction / decay speed (消退速度)       = the slowest OFF-phase relaxation
      rate = spectral gap of Q_off on the transient states (= reactivation rate).

RIGOR:
  ✅ The model FORM (states, allowed transitions, dox-modulates-silencing-only,
     I absorbing), the generator algebra, the master-equation integration, the
     SCC state classification, the absorption probabilities and the spectral
     gap are all exact, deterministic computations.
  ⚠️ The specific rate CONSTANTS are representative values calibrated to the
     timescales the paper REPORTS in text (50%-silencing at ~15 h KRAB/HDAC4,
     ~35 h EED, ~62 h DNMT3B; reactivation over days vs weeks; DNMT3B stable
     for 30 days). They are NOT the exact fitted constants from the paper's
     supplement (those live in figure panels we cannot transcribe to the
     decimal). The QUALITATIVE class structure they produce is the paper's.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

STATES = ("A", "R", "I")   # active, reversibly silent, irreversibly silent
A, R, I = 0, 1, 2


# --------------------------------------------------------------------------
# representative rate table (per hour), calibrated to the paper's REPORTED
# timescales. dox in [0,1] = relative recruitment strength (1 = saturating).
# --------------------------------------------------------------------------
# k_AR0 : A->R forward silencing at saturating dox (sets 50%-silencing time)
# k_RI  : R->I commitment to the irreversible state (sets permanent fraction)
# k_RA  : R->A reactivation after release (dox-INDEPENDENT; sets decay speed)
# A direct A->I for pure DNA-methylation (DNMT3B) is modelled as going via a
# very short-lived R with huge k_RI, i.e. essentially A->I.
CR_PARAMS = {
    # CR      k_AR0     k_RI     k_RA     note
    "HDAC4":  dict(k_AR0=0.046, k_RI=0.000, k_RA=0.020,
                   note="histone deacetylation -> short-term, fully reversible"),
    "KRAB":   dict(k_AR0=0.046, k_RI=0.010, k_RA=0.0040,
                   note="H3K9me3 -> hybrid memory (fast silence, weeks reactivation)"),
    "EED":    dict(k_AR0=0.020, k_RI=0.012, k_RA=0.0035,
                   note="PRC2/H3K27me3 -> hybrid memory (slower silence)"),
    "DNMT3B": dict(k_AR0=0.011, k_RI=0.300, k_RA=0.0000,
                   note="DNA methylation -> long-term, irreversible (I absorbing)"),
}


def generator(k_AR, k_RI, k_RA):
    """Build a 3x3 CTMC generator Q (rows sum to 0).  dp/dt = p Q.
    A->R (k_AR), R->I (k_RI), R->A (k_RA).  I is absorbing (no out-rates)."""
    Q = np.zeros((3, 3))
    Q[A, R] = k_AR
    Q[R, I] = k_RI
    Q[R, A] = k_RA
    for s in range(3):
        Q[s, s] = -Q[s].sum()
    return Q


def Q_recruit(cr, dox=1.0):
    """Recruitment phase (dox ON): forward silencing active, strength-scaled.
    Reactivation is negligible WHILE the CR is bound (k_RA ~ 0 during dox)."""
    p = CR_PARAMS[cr]
    return generator(k_AR=p["k_AR0"] * dox, k_RI=p["k_RI"], k_RA=0.0)


def Q_release(cr):
    """Release phase (dox OFF): forward silencing off; R reverts at k_RA;
    I absorbing.  k_RA depends only on CR identity (paper's clean result)."""
    p = CR_PARAMS[cr]
    return generator(k_AR=0.0, k_RI=0.0, k_RA=p["k_RA"])


# --------------------------------------------------------------------------
# master-equation flow:  p(t) = p0 · expm(Q t)   (probability diffusing on the
# state graph).  We use a dense matrix exponential via eigen-decomposition-free
# scaling-and-squaring (scipy) to stay exact.
# --------------------------------------------------------------------------
def evolve(p0, Q, t):
    from scipy.linalg import expm
    return p0 @ expm(Q * t)


def silencing_curve(cr, dox=1.0, t_max=120.0, n=240):
    """Fraction of cells silenced (not in A) vs time during recruitment."""
    Q = Q_recruit(cr, dox)
    ts = np.linspace(0, t_max, n)
    p0 = np.array([1.0, 0.0, 0.0])
    frac_silent = np.array([1.0 - evolve(p0, Q, t)[A] for t in ts])
    return ts, frac_silent


def t_half_silencing(cr, dox=1.0):
    """Time to 50% silenced during recruitment (interp)."""
    ts, fs = silencing_curve(cr, dox, t_max=400.0, n=2000)
    if fs[-1] < 0.5:
        return np.inf
    k = np.argmax(fs >= 0.5)
    t0, t1, f0, f1 = ts[k - 1], ts[k], fs[k - 1], fs[k]
    return t0 + (0.5 - f0) * (t1 - t0) / (f1 - f0)


# --------------------------------------------------------------------------
# (1) irreversible-trajectory classification via GRAPH THEORY.
# The OFF-phase chain's directed transition graph decomposes into strongly
# connected components (SCCs). A CLOSED SCC (no edges leaving it) is a
# recurrent class = a state you can never leave = an irreversible attractor.
# --------------------------------------------------------------------------
def classify_states(Q, tol=1e-12):
    import networkx as nx
    G = nx.DiGraph()
    G.add_nodes_from(range(3))
    for i in range(3):
        for j in range(3):
            if i != j and Q[i, j] > tol:
                G.add_edge(i, j)
    sccs = list(nx.strongly_connected_components(G))
    classes = []
    for comp in sccs:
        # closed (recurrent) iff no edge leaves the component
        leaves = any(v not in comp for u in comp for v in G.successors(u))
        classes.append({"states": sorted(STATES[i] for i in comp),
                        "recurrent": not leaves})
    return classes


def absorption_fraction(p_release):
    """After release, with forward rates off, the permanent (irreversibly
    silent) fraction = the probability mass already in I at release; the R
    mass all returns to A.  Returns (permanent_silent, reactivatable)."""
    return p_release[I], p_release[R]


# --------------------------------------------------------------------------
# (3) extinction / decay speed = slowest OFF-phase relaxation rate.
# Restrict Q_off to the TRANSIENT states {A,R} (I is absorbing). Its dominant
# (least-negative) eigenvalue lambda_slow sets reactivation: tau = -1/lambda.
# --------------------------------------------------------------------------
def decay_rate(cr):
    Q = Q_release(cr)
    sub = Q[np.ix_([A, R], [A, R])]      # transient block
    ev = np.linalg.eigvals(sub).real
    ev = ev[ev < -1e-12]
    if len(ev) == 0:                      # no decay (fully locked)
        return 0.0, np.inf
    lam = ev.max()                        # least negative = slowest
    return lam, -1.0 / lam                # (rate, timescale hours)


# --------------------------------------------------------------------------
# (2) duration-dependent fractional control: permanent silent fraction after
# release vs recruitment DURATION (the paper's central prediction, tested in
# their Fig. 3).
# --------------------------------------------------------------------------
def duration_scan(cr, dox=1.0, durations=None):
    if durations is None:
        durations = np.array([24, 48, 72, 96, 120, 168, 240])  # hours
    p0 = np.array([1.0, 0.0, 0.0])
    Qon = Q_recruit(cr, dox)
    out = []
    for T in durations:
        p_rec = evolve(p0, Qon, T)
        perm, react = absorption_fraction(p_rec)
        out.append({"T_rec_h": float(T), "silent_at_release": float(1 - p_rec[A]),
                    "permanent": float(perm), "reactivatable": float(react)})
    return out


def memory_class(cr):
    """Operational class label from the computed structure."""
    perm120 = duration_scan(cr, durations=[120])[0]["permanent"]
    _, tau = decay_rate(cr)
    if perm120 < 0.02:
        return "short-term (reversible)"
    if tau == np.inf:
        return "long-term (irreversible / absorbing)"
    return "hybrid (partial permanent + slow reactivation)"


# --------------------------------------------------------------------------
# figure
# --------------------------------------------------------------------------
def make_figure(path):
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    colors = {"HDAC4": "#16a085", "KRAB": "#8e44ad",
              "EED": "#2980b9", "DNMT3B": "#c0392b"}

    # panel 1: silencing curves during recruitment
    for cr, c in colors.items():
        ts, fs = silencing_curve(cr, t_max=120)
        ax[0].plot(ts, 100 * fs, color=c, lw=2, label=cr)
    ax[0].axhline(50, ls=":", color="k", lw=1)
    ax[0].set_title("Silencing during recruitment\n(all-or-none, single-rate)")
    ax[0].set_xlabel("time (h)"); ax[0].set_ylabel("% cells silenced")
    ax[0].legend(fontsize=8)

    # panel 2: permanent (memory) fraction vs recruitment duration
    for cr, c in colors.items():
        ds = duration_scan(cr)
        ax[1].plot([d["T_rec_h"] for d in ds],
                   [100 * d["permanent"] for d in ds],
                   "o-", color=c, lw=2, label=cr)
    ax[1].set_title("Epigenetic memory vs recruitment duration\n"
                    "(duration-dependent fractional control)")
    ax[1].set_xlabel("recruitment duration (h)")
    ax[1].set_ylabel("% permanently (irreversibly) silent")
    ax[1].legend(fontsize=8)

    # panel 3: reactivation (extinction) after a 120 h pulse
    for cr, c in colors.items():
        Qoff = Q_release(cr)
        p_rec = evolve(np.array([1.0, 0, 0]), Q_recruit(cr), 120.0)
        ts = np.linspace(0, 30 * 24, 200)        # 30 days in hours
        silent = np.array([1 - evolve(p_rec, Qoff, t)[A] for t in ts])
        ax[2].plot(ts / 24, 100 * silent, color=c, lw=2, label=cr)
    ax[2].set_title("Reactivation after release of a 120 h pulse\n"
                    "(decay speed set by CR identity only)")
    ax[2].set_xlabel("days after release"); ax[2].set_ylabel("% still silent")
    ax[2].legend(fontsize=8)

    fig.suptitle("M23 — chromatin memory as an A->R->I Markov chain "
                 "(Bintu et al., Science 2016)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(path, dpi=130)
    plt.close(fig)


def report(fig_path="figures/m23_chromatin_memory.png"):
    print("=" * 68)
    print("MODULE 23 — chromatin memory CTMC (A->R->I, Bintu 2016)")
    print("=" * 68)

    print("\ncalibration: 50%-silencing time during recruitment (dox=1):")
    print(f"  {'CR':>7} {'t_half(h)':>10}  paper (text)")
    paper_th = {"KRAB": "~15 (<=1 cell cycle)", "HDAC4": "~15 (<=1 cell cycle)",
                "EED": "~35", "DNMT3B": "~62"}
    for cr in CR_PARAMS:
        th = t_half_silencing(cr)
        print(f"  {cr:>7} {th:10.1f}  {paper_th[cr]}")

    print("\n(1) irreversible-trajectory classification (OFF-phase SCCs):")
    for cr in CR_PARAMS:
        cls = classify_states(Q_release(cr))
        rec = [c["states"] for c in cls if c["recurrent"]]
        tra = [c["states"] for c in cls if not c["recurrent"]]
        print(f"  {cr:>7}: recurrent(absorbing)={rec}  transient={tra}")
    print("  -> {I} is a closed (recurrent) class for every hybrid/long CR:")
    print("     the irreversible branch is a graph-theoretic absorbing class,")
    print("     NOT a tuned number — reach it and the trajectory never leaves.")

    print("\n(2) duration-dependent memory (permanent silent fraction):")
    print(f"  {'CR':>7} " + " ".join(f"{int(T)}h" for T in
          [24, 72, 120, 240]))
    for cr in CR_PARAMS:
        ds = {int(d["T_rec_h"]): d["permanent"]
              for d in duration_scan(cr, durations=[24, 72, 120, 240])}
        print(f"  {cr:>7} " + " ".join(f"{ds[T]:.2f}" for T in
              [24, 72, 120, 240]))
    print("  -> longer recruitment -> larger permanent fraction (hybrid CRs);")
    print("     HDAC4 ~0 (no memory), DNMT3B saturates high (full memory).")

    print("\n(3) extinction / decay speed after release (消退速度):")
    print(f"  {'CR':>7} {'lambda_slow(1/h)':>16} {'tau(days)':>10} "
          f"{'class':>34}")
    results = {}
    for cr in CR_PARAMS:
        lam, tau = decay_rate(cr)
        taud = tau / 24 if tau != np.inf else np.inf
        mc = memory_class(cr)
        results[cr] = {"lambda_slow": lam, "tau_h": tau, "class": mc}
        tau_s = "inf" if taud == np.inf else f"{taud:.1f}"
        print(f"  {cr:>7} {lam:16.4f} {tau_s:>10} {mc:>34}")
    print("  -> decay speed = OFF-phase spectral gap (slowest relaxation);")
    print("     DNMT3B: no transient decay (locked) -> tau=inf (long memory).")

    make_figure(fig_path)
    print(f"\nfigure written to: {fig_path}")
    return {"t_half": {cr: t_half_silencing(cr) for cr in CR_PARAMS},
            "classes": {cr: classify_states(Q_release(cr)) for cr in CR_PARAMS},
            "duration": {cr: duration_scan(cr) for cr in CR_PARAMS},
            "decay": results, "params": CR_PARAMS}


if __name__ == "__main__":
    report()
