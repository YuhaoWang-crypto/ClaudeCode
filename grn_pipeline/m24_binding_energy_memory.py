"""
Module 24 — Binding energy / efficiency / concentration -> transition rate ->
memory & extinction speed.

This answers the quantitative question directly: does a transcription factor's
BINDING EFFICIENCY / BINDING ENERGY / CONCENTRATION change (a) how much
epigenetic memory forms, and (b) how fast it decays (消退速度)?

Two layers, both with REAL data where possible:

  Layer 1 — thermodynamic occupancy from a REAL motif.
    A JASPAR position-frequency matrix (PFM) is read as a Berg-von Hippel
    ENERGY matrix: eps(b,i) = -ln[(f(b,i)+p) / (f_consensus(i)+p)]  (units kT).
    A site's mismatch energy ddG(seq) = sum_i eps(seq_i, i) >= 0 (consensus=0).
    Dissociation constant  Kd(seq) = Kd_consensus * exp(ddG).
    Occupancy  theta = [TF] / ([TF] + Kd)  (Hill-1 / Langmuir).
    We use the REAL REST motif (MA0138.2 — a repressor that recruits
    CoREST/HDAC, i.e. directly on the silencing axis) and the REAL KRAB-ZNF
    motif ZNF263 (MA0528.2). Energies and the occupancy curve are REAL
    consequences of the measured PFM.

  Layer 2 — occupancy -> rate -> memory, through the M23 chain.
    The paper's clean result: recruitment STRENGTH modulates the SILENCING
    rate, but NOT the reactivation rate. We encode k_AR = k_AR_max * theta and
    keep k_RA fixed, then ask how the permanent memory fraction and the
    population decay behave as a function of occupancy (hence of [TF] and ddG).

RIGOR:
  ✅ The energy matrix, ddG distribution, Kd, occupancy curve are exact
     computations from the REAL measured JASPAR PFM.
  ✅ The occupancy->rate->memory propagation through the M23 generator is exact.
  ✅ The qualitative split (concentration/energy move the SILENCING rate and the
     amount of memory; the intrinsic DECAY rate is occupancy-invariant) is the
     paper's experimental finding, reproduced here structurally.
  ⚠️ Kd_consensus, [TF] range and k_AR_max are representative scales (labelled);
     the shapes and the invariance result are what matter, not absolute nM.
"""
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from grn_pipeline import m23_chromatin_memory as m23

CACHE = "figures/_jaspar_cache/{id}.json"
MIRROR = "https://jaspar.elixir.no/api/v1/matrix/{id}/?format=json"


def load_pfm(matrix_id):
    """Load a JASPAR PFM, from local cache or (fallback) the JASPAR API."""
    path = CACHE.format(id=matrix_id)
    if not os.path.exists(path):
        import urllib.request
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with urllib.request.urlopen(MIRROR.format(id=matrix_id), timeout=30) as r:
            data = r.read().decode()
        with open(path, "w") as f:
            f.write(data)
    d = json.loads(open(path).read())
    pfm = np.array([d["pfm"][b] for b in "ACGT"], dtype=float)  # 4 x L
    return d.get("name", matrix_id), pfm


def energy_matrix(pfm, pseudo=1.0):
    """Berg-von Hippel energy matrix (kT). Consensus base per column = 0;
    mismatches cost positive energy.  eps = -ln[(f+p)/(f_max+p)]."""
    f = pfm + pseudo
    fmax = f.max(axis=0, keepdims=True)
    eps = -np.log(f / fmax)                   # 4 x L, >= 0, consensus col = 0
    return eps


def ddG_random_sites(eps, n=200_000, seed=0):
    """Distribution of mismatch binding energies ddG over random genomic
    sequence (i.i.d. uniform bases) — the real spread the motif imposes.
    By construction the CONSENSUS site has ddG=0 (the global minimum); random
    genomic sites are tens of kT worse, which is exactly why a motif confers
    specificity."""
    rng = np.random.default_rng(seed)
    L = eps.shape[1]
    idx = rng.integers(0, 4, size=(n, L))
    ddg = eps[idx, np.arange(L)].sum(axis=1)   # n,
    return ddg


def ddG_consensus(eps):
    """True consensus ddG = 0 (best base at every position)."""
    return float(eps.min(axis=0).sum())        # == 0 by construction


def ddG_worst(eps):
    """Worst possible site = max mismatch energy at every position."""
    return float(eps.max(axis=0).sum())


def ddG_k_mismatch(eps, k):
    """Representative ddG of a near-consensus site with the k costliest
    single-base mismatches applied to the consensus (a realistic 'weak site')."""
    per_pos_cost = np.sort(eps.max(axis=0) - eps.min(axis=0))[::-1]
    return float(per_pos_cost[:k].sum())


def occupancy(ddG, tf_conc, Kd_consensus=1.0):
    """Langmuir occupancy theta = C/(C + Kd), Kd = Kd_cons * exp(ddG)."""
    Kd = Kd_consensus * np.exp(ddG)
    return tf_conc / (tf_conc + Kd)


# --------------------------------------------------------------------------
# occupancy -> silencing rate -> memory, via the M23 chain
# --------------------------------------------------------------------------
def memory_vs_occupancy(cr="KRAB", thetas=None, T_rec=120.0):
    """For a hybrid CR, scale k_AR by occupancy theta (strength) and read the
    permanent memory fraction after a fixed recruitment window.  k_RA fixed."""
    if thetas is None:
        thetas = np.linspace(0.02, 1.0, 40)
    p = m23.CR_PARAMS[cr]
    p0 = np.array([1.0, 0.0, 0.0])
    out = []
    for th in thetas:
        Q = m23.generator(k_AR=p["k_AR0"] * th, k_RI=p["k_RI"], k_RA=0.0)
        p_rec = m23.evolve(p0, Q, T_rec)
        out.append({"theta": float(th),
                    "silent_at_release": float(1 - p_rec[m23.A]),
                    "permanent": float(p_rec[m23.I])})
    return out


def decay_is_strength_invariant(cr="KRAB", thetas=(0.2, 0.5, 1.0), T_rec=120.0):
    """Show the OFF-phase decay RATE is the same for every recruitment strength
    (only the amplitude/permanent fraction differs) — the paper's result."""
    p = m23.CR_PARAMS[cr]
    p0 = np.array([1.0, 0.0, 0.0])
    Qoff = m23.Q_release(cr)
    ts = np.linspace(0, 30 * 24, 240)
    curves = {}
    rates = {}
    for th in thetas:
        Qon = m23.generator(k_AR=p["k_AR0"] * th, k_RI=p["k_RI"], k_RA=0.0)
        p_rec = m23.evolve(p0, Qon, T_rec)
        silent = np.array([1 - m23.evolve(p_rec, Qoff, t)[m23.A] for t in ts])
        curves[th] = (ts, silent)
        # fit the decaying part (toward the permanent floor) to get a rate
        floor = p_rec[m23.I]
        y = silent - floor
        mask = y > 1e-4
        if mask.sum() > 5:
            slope = np.polyfit(ts[mask], np.log(y[mask]), 1)[0]
            rates[th] = -slope        # 1/h
        else:
            rates[th] = np.nan
    return curves, rates


def make_figure(name_rest, eps_rest, ddg_rest, name_znf, eps_znf, ddg_znf, path):
    fig, ax = plt.subplots(2, 3, figsize=(15, 8.4))

    # row 1: REAL motif -> energy -> occupancy
    im = ax[0, 0].imshow(eps_rest, aspect="auto", cmap="viridis")
    ax[0, 0].set_yticks(range(4)); ax[0, 0].set_yticklabels(list("ACGT"))
    ax[0, 0].set_title(f"{name_rest} energy matrix (kT)\n(real JASPAR MA0138.2)")
    ax[0, 0].set_xlabel("motif position")
    fig.colorbar(im, ax=ax[0, 0], fraction=0.046)

    ax[0, 1].hist(ddg_rest, bins=60, color="#2980b9", alpha=0.85)
    ax[0, 1].axvline(ddg_rest.min(), color="k", ls=":",
                     label=f"consensus ddG={ddg_rest.min():.1f}")
    ax[0, 1].set_title(f"{name_rest}: binding-energy spread\nover random sites")
    ax[0, 1].set_xlabel("ddG (kT)  [lower = tighter binding]")
    ax[0, 1].set_ylabel("# sites"); ax[0, 1].legend(fontsize=8)

    for ddg, nm, c in [(0.0, "consensus", "#c0392b"),
                       (ddG_k_mismatch(eps_rest, 1), "1 mismatch", "#27ae60"),
                       (ddG_k_mismatch(eps_rest, 2), "2 mismatch", "#8e44ad")]:
        C = np.logspace(-1, 5, 200)
        ax[0, 2].semilogx(C, occupancy(ddg, C), color=c, lw=2,
                          label=f"{nm} (ddG={ddg:.1f})")
    ax[0, 2].set_title("Occupancy vs [TF] concentration\n"
                       "theta = C/(C+Kd), Kd=exp(ddG)")
    ax[0, 2].set_xlabel("[TF] (a.u.)"); ax[0, 2].set_ylabel("occupancy theta")
    ax[0, 2].legend(fontsize=8)

    # row 2: occupancy -> memory, and decay invariance
    mv = memory_vs_occupancy("KRAB")
    ax[1, 0].plot([m["theta"] for m in mv], [100 * m["silent_at_release"] for m in mv],
                  "o-", color="#e67e22", label="silent at release")
    ax[1, 0].plot([m["theta"] for m in mv], [100 * m["permanent"] for m in mv],
                  "s-", color="#c0392b", label="permanent (memory)")
    ax[1, 0].set_title("Occupancy -> memory (KRAB, 120 h window)\n"
                       "higher occupancy -> more permanent silencing")
    ax[1, 0].set_xlabel("occupancy theta (set by [TF] & ddG)")
    ax[1, 0].set_ylabel("% cells"); ax[1, 0].legend(fontsize=8)

    curves, rates = decay_is_strength_invariant("KRAB")
    cmap = {0.2: "#95a5a6", 0.5: "#3498db", 1.0: "#2c3e50"}
    for th, (ts, silent) in curves.items():
        ax[1, 1].plot(ts / 24, 100 * silent, color=cmap[th], lw=2,
                      label=f"theta={th} (rate={rates[th]:.4f}/h)")
    ax[1, 1].set_title("Decay after release at 3 strengths\n"
                       "amplitude differs, RATE is invariant")
    ax[1, 1].set_xlabel("days after release"); ax[1, 1].set_ylabel("% silent")
    ax[1, 1].legend(fontsize=8)

    ax[1, 2].imshow(eps_znf, aspect="auto", cmap="viridis")
    ax[1, 2].set_yticks(range(4)); ax[1, 2].set_yticklabels(list("ACGT"))
    ax[1, 2].set_title(f"{name_znf} energy matrix (kT)\n"
                       f"(real KRAB-ZNF, JASPAR MA0528.2)")
    ax[1, 2].set_xlabel("motif position")

    fig.suptitle("M24 — binding energy / concentration -> silencing rate -> "
                 "memory  (decay rate is strength-invariant)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(path, dpi=130)
    plt.close(fig)


def report(fig_path="figures/m24_binding_energy.png"):
    print("=" * 68)
    print("MODULE 24 — binding energy / concentration -> rate -> memory")
    print("=" * 68)

    name_rest, pfm_rest = load_pfm("MA0138.2")
    name_znf, pfm_znf = load_pfm("MA0528.2")
    eps_rest = energy_matrix(pfm_rest)
    eps_znf = energy_matrix(pfm_znf)
    ddg_rest = ddG_random_sites(eps_rest)
    ddg_znf = ddG_random_sites(eps_znf)

    print(f"\nreal motifs (JASPAR):")
    print(f"  {name_rest:>8} (MA0138.2) L={eps_rest.shape[1]}  "
          f"consensus ddG={ddG_consensus(eps_rest):.1f}, "
          f"worst site ddG={ddG_worst(eps_rest):.1f} kT, "
          f"median random site={np.median(ddg_rest):.1f} kT")
    print(f"  {name_znf:>8} (MA0528.2) L={eps_znf.shape[1]}  "
          f"consensus ddG={ddG_consensus(eps_znf):.1f}, "
          f"worst site ddG={ddG_worst(eps_znf):.1f} kT, "
          f"median random site={np.median(ddg_znf):.1f} kT")

    # occupancy at a [TF] scaled so the consensus site is ~half-occupied
    # (Kd_consensus=1), i.e. [TF]=100 a.u. to also fill near-consensus sites.
    tf = 100.0
    print(f"\noccupancy at [TF]={tf:.0f} (a.u., Kd_consensus=1) for REST sites "
          f"of varying quality:")
    for q, ddg in [("consensus (0 mm)", 0.0),
                   ("1 mismatch", ddG_k_mismatch(eps_rest, 1)),
                   ("2 mismatch", ddG_k_mismatch(eps_rest, 2)),
                   ("median random", float(np.median(ddg_rest)))]:
        print(f"  {q:>16}: ddG={ddg:5.1f} kT -> Kd={np.exp(ddg):.2e} -> "
              f"theta={occupancy(ddg, tf):.3f}")

    print("\noccupancy -> permanent memory (KRAB, 120 h recruitment window):")
    mv = memory_vs_occupancy("KRAB", thetas=[0.1, 0.3, 0.5, 0.8, 1.0])
    print(f"  {'theta':>6} {'silent@release':>15} {'permanent':>10}")
    for m in mv:
        print(f"  {m['theta']:6.2f} {m['silent_at_release']:15.3f} "
              f"{m['permanent']:10.3f}")

    print("\nextinction/decay RATE vs recruitment strength (KRAB):")
    curves, rates = decay_is_strength_invariant("KRAB")
    for th in sorted(rates):
        print(f"  theta={th:.1f}: population decay rate = {rates[th]:.5f} /h")
    rate_vals = [r for r in rates.values() if not np.isnan(r)]
    spread = (max(rate_vals) - min(rate_vals)) / np.mean(rate_vals)
    print(f"  -> decay rate spread across strengths = {spread*100:.1f}% "
          f"(≈invariant): CONCENTRATION sets how MUCH memory forms, NOT how")
    print(f"     fast it decays. Decay speed is a property of the CR/eraser")
    print(f"     identity (k_RA), not of recruiter concentration — the paper's")
    print(f"     result, reproduced structurally.")

    make_figure(name_rest, eps_rest, ddg_rest, name_znf, eps_znf, ddg_znf,
                fig_path)
    print(f"\nfigure written to: {fig_path}")
    return {"rest": {"name": name_rest, "max_ddG": float(eps_rest.sum(0).max()),
                     "median_ddG": float(np.median(ddg_rest))},
            "znf": {"name": name_znf, "max_ddG": float(eps_znf.sum(0).max())},
            "memory_vs_occupancy": mv,
            "decay_rates": {float(k): float(v) for k, v in rates.items()},
            "decay_rate_spread": float(spread)}


if __name__ == "__main__":
    report()
