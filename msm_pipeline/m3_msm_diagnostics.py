"""
Module 3 - build the textbook MSM (TICA -> k-means -> transition matrix) and
test whether it converges.  It does not, and the reason is instructive.

The implied-timescale (ITS) test is the honest gate on every MSM number: if
t_i(tau) still climbs with the lag time tau, the discretisation is not
Markovian and the eigenvalues are not rates.  Reporting a flux network without
showing this plot is the single most common way MSM analyses mislead.

Diagnosis here: with 32 *independent* 500 ns runs that each start from a
different crystal structure, the slowest apparent "process" is the failure of
the runs to interconvert at all.  TICA happily picks that up -- 60-85 % of the
variance along the leading components is BETWEEN trajectories rather than
within them -- and the MSM then reports the inverse of an almost-zero rate,
i.e. a timescale of 10^5 ns from 500 ns of data.

Takeaway for the user's own runs: aggregate length is not the binding
constraint, *connectivity* is.  Many short trajectories seeded from the
barrier region beat a few long ones.
"""
import numpy as np

from .common import DT_NS, apo_indices, banner, load

LAGS_FRAMES = [10, 20, 40, 80, 160, 320]


def _featurize(X, residues, cutoff, stride_ref=50):
    """Pairwise C-alpha distances among `residues`, kept if mean dist < cutoff."""
    from .common import amb_to_idx
    R = np.array(sorted({amb_to_idx(a) for a in residues}))
    I, J = np.triu_indices(X.shape[2], 3)
    keep = np.isin(I, R) & np.isin(J, R)
    I, J = I[keep], J[keep]
    ref = np.asarray(X[:32, ::stride_ref]).reshape(-1, X.shape[2], 3)
    d = np.linalg.norm(ref[:, I] - ref[:, J], axis=-1).mean(0)
    I, J = I[d < cutoff], J[d < cutoff]
    return I, J


def _its(dtrajs, lags, n=3):
    from deeptime.markov import TransitionCountEstimator
    from deeptime.markov.msm import MaximumLikelihoodMSM
    out = []
    for lag in lags:
        c = (TransitionCountEstimator(lag, "sliding").fit(dtrajs)
             .fetch_model().submodel_largest(directed=True))
        m = MaximumLikelihoodMSM(reversible=True).fit(c).fetch_model()
        out.append((lag, c.n_states, m.timescales(n) * DT_NS))
    return out


def main():
    from deeptime.clustering import KMeans
    from deeptime.decomposition import TICA

    banner("M3  MSM diagnostics - the implied-timescale test, and why it fails")
    X, names = load()
    apo = apo_indices(names)

    for tag, residues, cutoff, dim in [
        ("global (whole protein, <10 A contacts)",
         [a for a in range(26, 291) if a not in (239, 253)], 10.0, 4),
        ("pocket-local (H10-H11/H11 + H12 only)",
         list(range(214, 232)) + list(range(269, 287)), 18.0, 2),
    ]:
        I, J = _featurize(X, residues, cutoff)
        F = [np.linalg.norm(np.asarray(X[k])[:, I] - np.asarray(X[k])[:, J],
                            axis=-1).astype(np.float32) for k in apo]
        tica = TICA(lagtime=100, dim=max(dim, 3), scaling="kinetic_map")
        for f in F:
            tica.partial_fit(f)
        m = tica.fetch_model()
        Y = np.stack([m.transform(f) for f in F])

        btw = [Y[:, :, i].mean(1).var() / Y[:, :, i].var() for i in range(3)]
        print(f"\n  {tag}")
        print(f"    features {len(I)} distance pairs")
        print(f"    fraction of TIC variance that is BETWEEN trajectories: "
              f"{np.round(btw, 2)}")

        km = (KMeans(n_clusters=100, max_iter=300, init_strategy="kmeans++",
                     n_jobs=4, fixed_seed=1)
              .fit(Y[:, ::3, :dim].reshape(-1, dim)).fetch_model())
        dt = [km.transform(Y[k, :, :dim]) for k in range(len(apo))]
        print(f"    {'lag (ns)':>9s} {'states':>7s}   slowest 3 implied timescales (ns)")
        for lag, ns, ts in _its(dt, LAGS_FRAMES):
            print(f"    {lag*DT_NS:9.2f} {ns:7d}   " + ", ".join(f"{v:10.1f}" for v in ts))

    print("""
  READING THIS TABLE
    A converged MSM shows timescales that FLATTEN as the lag grows.  Neither
    discretisation does.  The global one is catastrophic (10^2 -> 10^5 ns): its
    slowest 'process' is trajectory identity, not physics.  The pocket-local one
    is far better behaved but still climbing, so its eigenvalues are lower
    bounds, not rates.

  CONSEQUENCE, APPLIED HONESTLY
    * state POPULATIONS and COMMITTED CROSSING COUNTS are direct observables -
      no Markov assumption - and are used for the quantitative claims (M5).
    * the TPT flux TOPOLOGY (where the bottleneck sits, how many routes carry
      the flux) is robust to this failure and is reported in M4.
    * absolute rates / MFPTs from this dataset are NOT trustworthy and are
      labelled as such.""")


if __name__ == "__main__":
    main()
