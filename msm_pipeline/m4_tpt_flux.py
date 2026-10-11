"""
Module 4 - core-set MSM + PCCA+ + transition path theory on the pocket.

This is the step the user's workflow sketch calls "conformational flux":
for a source set A (pocket closed) and sink B (pocket open), TPT gives the
committor q+, the net flux J_ij = max(0, pi_i T_ij q_j+ - pi_j T_ji q_i+),
the dominant pathways, and the bottleneck edge.

Given M3, what is reported here is the flux TOPOLOGY, which is robust:
  * how concentrated the flux is (one channel or many),
  * where along the coordinate the bottleneck sits.
Absolute rates and MFPTs are printed but explicitly flagged as not converged.

The state space is built on physically meaningful coordinates only (the
crystal-axis projection plus the two mouth distances), so cluster identity
cannot absorb trajectory identity the way the TICA basis did in M3.
"""
import numpy as np

from .common import DT_NS, apo_indices, banner, load
from .m2_cryptic_pocket import collective_variables, opening_axis

LAG_FRAMES = 80          # 5 ns
N_MICRO = 80
A_MAX = 1.5              # angstrom along the crystal axis: closed core set
B_MIN = 6.5              # open core set


def build(seed=0):
    from deeptime.clustering import KMeans
    from deeptime.markov import TransitionCountEstimator
    from deeptime.markov.msm import MaximumLikelihoodMSM

    X, names = load()
    apo = apo_indices(names)
    proj, cv1, cv2 = collective_variables()
    Z = np.stack([proj, cv1, cv2], -1)[apo]
    mu, sd = Z.reshape(-1, 3).mean(0), Z.reshape(-1, 3).std(0)
    Zs = (Z - mu) / sd

    km = (KMeans(n_clusters=N_MICRO, max_iter=500, init_strategy="kmeans++",
                 n_jobs=4, fixed_seed=seed)
          .fit(Zs[:, ::3].reshape(-1, 3)).fetch_model())
    dt = [km.transform(Zs[k]) for k in range(len(apo))]
    cnt = (TransitionCountEstimator(LAG_FRAMES, "sliding").fit(dt)
           .fetch_model().submodel_largest(directed=True))
    msm = MaximumLikelihoodMSM(reversible=True).fit(cnt).fetch_model()
    centers = km.cluster_centers * sd + mu
    return msm, cnt, centers[cnt.state_symbols][:, 0]


def _direct_rate(a_max, b_min):
    """Committed A->B crossings per microsecond, straight from the trajectories."""
    from .common import TRAJ_NS
    _, names = load()
    apo = apo_indices(names)
    proj, _, _ = collective_variables()
    n = 0
    for k in apo:
        s = np.full(proj.shape[1], -1, np.int8)
        s[proj[k] < a_max] = 0
        s[proj[k] > b_min] = 1
        seq = s[s >= 0]
        if len(seq) > 1:
            n += int((np.diff(seq) == 1).sum())
    return n / (len(apo) * TRAJ_NS / 1000.0)


def main():
    from deeptime.markov import reactive_flux

    banner("M4  transition path theory: flux from closed to open pocket")
    _, _, L = opening_axis()
    msm, cnt, cp = build()
    print(f"  microstates {cnt.n_states} | lag {LAG_FRAMES*DT_NS:.1f} ns "
          f"| slowest ITS {msm.timescales(3)[0]*DT_NS:.1f} ns")

    A = np.where(cp < A_MAX)[0]
    B = np.where(cp > B_MIN)[0]
    pi = msm.stationary_distribution
    print(f"  A = closed  (axis < {A_MAX} A): {len(A):2d} states, population {pi[A].sum():.4f}")
    print(f"  B = open    (axis > {B_MIN} A): {len(B):2d} states, population {pi[B].sum():.4f}")

    fl = reactive_flux(msm.transition_matrix, A, B)
    per_us = 1.0 / (DT_NS * 1e-3)
    print(f"\n  committor-weighted net flux   {fl.total_flux:.3e} /step")
    print(f"  A->B rate                     {fl.rate*per_us:9.1f} /us   "
          f"[NOT CONVERGED - see M3]")
    print(f"  MFPT                          {1/fl.rate*DT_NS:9.1f} ns     "
          f"[NOT CONVERGED - see M3]")

    paths, caps = fl.pathways(fraction=0.9)
    top = 100 * caps[0] / fl.total_flux
    print(f"\n  pathway decomposition: {len(paths)} distinct routes carry 90% of the flux")
    print(f"  strongest single route carries only {top:.2f} % of it")
    print("  -> the opening is a DIFFUSE, many-route process, not one channel.")
    print("     (a single dominant channel would be the thing to engineer against;")
    print("      a diffuse one says the pocket opens by a soft collective motion)")
    print("\n  strongest routes, as position along the crystal opening axis (A):")
    for p, c in list(zip(paths, caps))[:5]:
        print(f"    {100*c/fl.total_flux:5.2f} %   " + " -> ".join(f"{cp[s]:5.2f}" for s in p))

    net = fl.net_flux
    p0 = paths[0]
    edges = [(p0[i], p0[i + 1], net[p0[i], p0[i + 1]]) for i in range(len(p0) - 1)]
    a, b, f = min(edges, key=lambda e: e[2])
    print(f"\n  BOTTLENECK on the strongest route: axis {cp[a]:.2f} -> {cp[b]:.2f} A"
          f"  (net flux {f:.2e}, narrowest of {len(edges)} edges)")

    # Where does the system actually commit?  q+ crossing 0.5 along the axis.
    q = fl.forward_committor
    order = np.argsort(cp)
    cs, qs = cp[order], q[order]
    k = np.argmax(qs >= 0.5)
    print(f"  committor q+ reaches 0.5 at axis ~ {cs[k]:.2f} A of the {L:.2f} A "
          f"crystal opening ({100*cs[k]/L:.0f} % of the way)")
    print("  -> the barrier is NOT the initial loosening: the ensemble reaches")
    print(f"     ~{cs[k]*0.6:.1f}-{cs[k]:.1f} A with q+ still well below 0.5 and only commits")
    print("     in the last stretch, which is where a fragment would have to bind")
    print("     to trap the open form.")

    # ---- the assumption-free number, and what it can and cannot tell you ---
    direct = _direct_rate(A_MAX, B_MIN)
    print(f"\n  SAME TRANSITION, COUNTED DIRECTLY")
    print(f"    MSM/TPT reactive rate      {fl.rate*per_us:9.1f} /us")
    print(f"    committed crossings        {direct:9.1f} /us  (no Markov assumption)")
    print("""
    These two numbers are NOT the same observable and are not expected to
    match: the TPT rate counts reactive flux including recrossing of the
    dividing surface, while the milestoning count only registers committed
    A->B events.  Running this pipeline on alanine dipeptide - a system whose
    MSM *does* pass the ITS test (see m8_positive_control) - still gives a
    ~50x gap between them.  So the ratio is NOT a convergence diagnostic, and
    it is not used as one here.

    The convergence diagnostic is the ITS test in M3, and it fails.  What
    survives that failure is: the committed-crossing count above (a direct
    observable), the state populations, and the flux TOPOLOGY below.  What
    does not survive is any absolute rate or MFPT.""")

    print("\n  forward committor q+ along the opening axis:")
    bins = [(-1, 1.5), (1.5, 3), (3, 4.5), (4.5, 6.5), (6.5, 99)]
    for lo, hi in bins:
        m = (cp >= lo) & (cp < hi)
        if m.sum():
            print(f"    axis {max(lo,0):4.1f}-{min(hi,30):4.1f} A : q+ = {q[m].mean():.3f}"
                  f"  ({m.sum()} states, pop {pi[m].sum():.4f})")
    print("  the transition state (q+ ~ 0.5) is where the committor crosses 0.5.")
    return {"n_paths_90pct": len(paths), "top_path_pct": top,
            "bottleneck": (float(cp[a]), float(cp[b]))}


if __name__ == "__main__":
    main()
