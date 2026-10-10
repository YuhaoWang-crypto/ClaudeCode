"""
Module 23 — Conformational flux: MSM + Transition Path Theory on real MD.

RMSD measures how far a frame is from a reference structure; it cannot tell
a molecule rattling inside one well from one crossing a barrier into a new
state. Conformational flux can: discretise the trajectory into states, estimate
the transition matrix T(tau), and use Transition Path Theory (TPT) to compute
the net probability current from a source set A to a target set B.

Pipeline (deeptime):
  1. featurize   : internal coordinates (backbone phi/psi/chi1 as sin/cos,
                   optional C-alpha contact distances) -- no superposition,
                   no reference frame
  2. TICA        : slowest linear combinations (kinetic-map scaling)
  3. k-means     : micro-states in TICA space -> discrete trajectories
  4. MSM         : reversible max-likelihood T(tau); lag chosen from implied
                   timescales; validated with a Chapman-Kolmogorov test
  5. PCCA+       : metastable macro-states
  6. TPT         : committors q+, q-, net flux f+_ij, total flux, rate k_AB,
                   pathway decomposition, bottleneck edge
  7. validation  : MSM mean first-passage time vs. MFPT counted directly in
                   the trajectories (independent of the Markov assumption)

Net flux (Metzner, Schuette & Vanden-Eijnden 2009):
    f_ij  = pi_i q-_i T_ij q+_j          (i != j)
    f+_ij = max(0, f_ij - f_ji)
For a reversible MSM q- = 1 - q+, and this reduces exactly to
    f+_ij = max(0, pi_i T_ij (q+_j - q+_i)),
i.e. the "committor-difference" form. The module checks both numerically.

Demo system: alanine dipeptide (ACE-ALA-NME), amber99sb-ildn + OBC implicit
solvent, 300 K, simulated here with OpenMM. Runs on any GROMACS/AMBER
trajectory via the CLI:

    python3 -m grn_pipeline.m23_msm_tpt --top conf.gro --traj md.xtc \
        --dt-ps 10 --features backbone --n-macro 4

(.tpr is not readable by mdtraj: export a .gro/.pdb with `gmx editconf`, and
make molecules whole first with `gmx trjconv -pbc mol -center`.)
"""
import argparse
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "..", "figures", "_md_cache")


# ---------------------------------------------------------------------------
# 0. demo trajectory: build + simulate alanine dipeptide with OpenMM
# ---------------------------------------------------------------------------
def _place(a, b, c, bond, angle, torsion):
    """NeRF: position of D from A,B,C given |CD| (A), angle BCD and dihedral
    ABCD (degrees)."""
    angle, torsion = np.radians(angle), np.radians(torsion)
    bc = (c - b) / np.linalg.norm(c - b)
    n = np.cross(b - a, bc); n /= np.linalg.norm(n)
    m = np.cross(n, bc)
    d = np.array([-bond * np.cos(angle),
                  bond * np.sin(angle) * np.cos(torsion),
                  bond * np.sin(angle) * np.sin(torsion)])
    return c + d[0] * bc + d[1] * m + d[2] * n


def build_alanine_dipeptide(path):
    """Write an all-atom ACE-ALA-NME PDB built from standard internal
    coordinates (extended phi=-150, psi=150; L-chirality via the
    N-C-CA-CB dihedral of +122.6 deg), hydrogens added by OpenMM."""
    import openmm as mm
    import openmm.app as app
    import openmm.unit as u
    X = {"ACE_CH3": np.zeros(3), "ACE_C": np.array([1.52, 0.0, 0.0])}
    X["ALA_N"] = X["ACE_C"] + 1.33 * np.array(
        [np.cos(np.radians(64)), np.sin(np.radians(64)), 0.0])
    X["ALA_CA"] = _place(X["ACE_CH3"], X["ACE_C"], X["ALA_N"], 1.46, 122, 180)
    X["ALA_C"] = _place(X["ACE_C"], X["ALA_N"], X["ALA_CA"], 1.52, 111, -150)
    X["NME_N"] = _place(X["ALA_N"], X["ALA_CA"], X["ALA_C"], 1.33, 116, 150)
    X["NME_C"] = _place(X["ALA_CA"], X["ALA_C"], X["NME_N"], 1.46, 122, 180)
    X["ACE_O"] = _place(X["ALA_N"], X["ACE_CH3"], X["ACE_C"], 1.23, 121, 180)
    X["ALA_O"] = _place(X["NME_N"], X["ALA_CA"], X["ALA_C"], 1.23, 121, 180)
    X["ALA_CB"] = _place(X["ALA_N"], X["ALA_C"], X["ALA_CA"], 1.53, 110.1,
                         122.6)
    top = app.Topology(); chain = top.addChain(); pos = []
    for rn, names in [("ACE", ["CH3", "C", "O"]),
                      ("ALA", ["N", "CA", "CB", "C", "O"]),
                      ("NME", ["N", "C"])]:
        res = top.addResidue(rn, chain)
        for nm in names:
            top.addAtom(nm, app.element.Element.getBySymbol(nm[0]), res)
            pos.append(X[f"{rn}_{nm}"] / 10.0)
    top.createStandardBonds()
    ff = app.ForceField("amber99sbildn.xml", "amber99_obc.xml")
    mod = app.Modeller(top, u.Quantity([mm.Vec3(*p) for p in pos],
                                       u.nanometer))
    mod.addHydrogens(ff)
    system = ff.createSystem(mod.topology, nonbondedMethod=app.NoCutoff,
                             constraints=app.HBonds)
    sim = app.Simulation(mod.topology, system,
                         mm.VerletIntegrator(1 * u.femtosecond),
                         mm.Platform.getPlatformByName("Reference"))
    sim.context.setPositions(mod.positions)
    sim.minimizeEnergy()
    with open(path, "w") as fh:
        app.PDBFile.writeFile(mod.topology, sim.context.getState(
            getPositions=True).getPositions(), fh)
    return path


def _run_one(args):
    pdb_path, out_xtc, ns, seed, save_ps = args
    import openmm as mm
    import openmm.app as app
    import openmm.unit as u
    pdb = app.PDBFile(pdb_path)
    ff = app.ForceField("amber99sbildn.xml", "amber99_obc.xml")
    # HMR (H mass 4 amu) allows a 4 fs step; it rescales absolute kinetics by
    # a modest factor but leaves the thermodynamics unchanged.
    system = ff.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff,
                             constraints=app.HBonds, hydrogenMass=4 * u.amu)
    dt_fs = 4.0
    integ = mm.LangevinMiddleIntegrator(300 * u.kelvin, 1 / u.picosecond,
                                        dt_fs * u.femtosecond)
    integ.setRandomNumberSeed(seed)
    sim = app.Simulation(pdb.topology, system, integ,
                         mm.Platform.getPlatformByName("Reference"))
    sim.context.setPositions(pdb.positions)
    sim.context.setVelocitiesToTemperature(300 * u.kelvin, seed)
    sim.step(25000)                                  # 100 ps equilibration
    every = int(round(save_ps * 1000 / dt_fs))
    sim.reporters.append(app.XTCReporter(out_xtc, every))
    sim.step(int(round(ns * 1e6 / dt_fs)))
    return out_xtc


def simulate_alanine_dipeptide(n_traj=4, ns=40.0, save_ps=1.0, out=CACHE):
    """n_traj independent Langevin runs in parallel processes. Returns
    (topology_pdb, [xtc paths])."""
    from multiprocessing import Pool
    os.makedirs(out, exist_ok=True)
    pdb_path = os.path.join(out, "ala2.pdb")
    if not os.path.exists(pdb_path):
        build_alanine_dipeptide(pdb_path)
    jobs = [(pdb_path, os.path.join(out, f"ala2_run{i}.xtc"), ns, 1000 + i,
             save_ps) for i in range(n_traj)]
    todo = [j for j in jobs if not os.path.exists(j[1])]
    if todo:
        with Pool(min(len(todo), os.cpu_count() or 1)) as pool:
            pool.map(_run_one, todo)
    return pdb_path, [j[1] for j in jobs]


# ---------------------------------------------------------------------------
# 1. featurization: internal coordinates only (no superposition)
# ---------------------------------------------------------------------------
def featurize(top, trajs, features="backbone", stride=1, ca_min_sep=3):
    """Returns (feature arrays per trajectory, feature names, phi/psi of the
    first residue per trajectory in degrees or None)."""
    import mdtraj as md
    out, rama, names = [], [], None
    for path in trajs:
        t = md.load(path, top=top, stride=stride)
        cols, nm = [], []
        if features in ("backbone", "both"):
            for fn, lab in [(md.compute_phi, "phi"), (md.compute_psi, "psi"),
                            (md.compute_chi1, "chi1")]:
                idx, ang = fn(t)
                if ang.shape[1]:
                    cols += [np.sin(ang), np.cos(ang)]
                    nm += [f"sin_{lab}{k}" for k in range(ang.shape[1])]
                    nm += [f"cos_{lab}{k}" for k in range(ang.shape[1])]
        if features in ("contacts", "both"):
            ca = t.topology.select("name CA")
            pairs = [(ca[a], ca[b]) for a in range(len(ca))
                     for b in range(a + ca_min_sep, len(ca))]
            if pairs:
                cols.append(md.compute_distances(t, pairs))
                nm += [f"d_CA{i}-CA{j}" for i, j in pairs]
        if not cols:
            raise ValueError(f"no features of kind '{features}' for {path}")
        out.append(np.hstack(cols).astype(np.float64))
        names = nm
        _, phi = md.compute_phi(t); _, psi = md.compute_psi(t)
        rama.append(np.degrees(np.c_[phi[:, 0], psi[:, 0]])
                    if phi.shape[1] and psi.shape[1] else None)
    if any(r is None for r in rama):
        rama = None
    return out, names, rama


# ---------------------------------------------------------------------------
# 2-4. TICA -> k-means -> MSM
# ---------------------------------------------------------------------------
def tica_cluster(feats, tica_lag, n_clusters=100, var_cutoff=0.95, seed=7,
                 method="regspace"):
    """TICA projection, then micro-state clustering.
    regspace: centres at least dmin apart (dmin bisected to give about
    n_clusters), so sparse barrier regions and rare states get their own
    micro-states -- this is where TPT needs resolution. kmeans: centres
    follow density (fine basins, coarse barriers)."""
    from deeptime.clustering import KMeans, RegularSpace
    from deeptime.decomposition import TICA
    tica = TICA(lagtime=tica_lag, var_cutoff=var_cutoff,
                scaling="kinetic_map").fit(feats).fetch_model()
    Y = [tica.transform(x) for x in feats]
    cat = np.concatenate(Y)
    fit = cat[::max(1, len(cat) // 50000)]
    if method == "kmeans":
        km = KMeans(n_clusters, init_strategy="kmeans++", fixed_seed=seed,
                    max_iter=300).fit(fit).fetch_model()
    else:
        lo, hi = 1e-4, float(np.ptp(fit, 0).max())
        for _ in range(30):
            dmin = 0.5 * (lo + hi)
            km = RegularSpace(dmin=dmin, max_centers=10 * n_clusters) \
                .fit(fit).fetch_model()
            if km.n_clusters > n_clusters:
                lo = dmin
            else:
                hi = dmin
            if abs(km.n_clusters - n_clusters) <= 0.05 * n_clusters:
                break
    dtrajs = [km.transform(y).astype(np.int32) for y in Y]
    return tica, Y, km, dtrajs


def estimate_msm(dtrajs, lag, reversible=True):
    """Reversible MLE MSM on the largest connected set. Returns (msm,
    dtrajs re-indexed to the active set with -1 for inactive frames)."""
    from deeptime.markov import TransitionCountEstimator
    from deeptime.markov.msm import MaximumLikelihoodMSM
    counts = TransitionCountEstimator(lag, "sliding").fit(dtrajs) \
        .fetch_model().submodel_largest()
    msm = MaximumLikelihoodMSM(reversible=reversible).fit(counts) \
        .fetch_model()
    active = counts.transform_discrete_trajectories_to_submodel(dtrajs)
    return msm, active


def implied_timescales(dtrajs, lags, n_its=5, n_samples=0):
    """ITS (in frames) per lag; optional Bayesian 95% intervals."""
    rows, lo, hi = [], [], []
    for lag in lags:
        msm, _ = estimate_msm(dtrajs, lag)
        ts = msm.timescales(n_its)
        rows.append(np.pad(ts, (0, n_its - len(ts)), constant_values=np.nan))
        if n_samples:
            from deeptime.markov import TransitionCountEstimator
            from deeptime.markov.msm import BayesianMSM
            c = TransitionCountEstimator(lag, "effective").fit(dtrajs) \
                .fetch_model().submodel_largest()
            post = BayesianMSM(n_samples=n_samples).fit(c).fetch_model()
            st = post.gather_stats("timescales", k=n_its)
            lo.append(st.L); hi.append(st.R)
    return (np.array(rows), np.array(lo) if lo else None,
            np.array(hi) if hi else None)


def choose_lag(lags, its, tol=0.10):
    """Smallest lag after which the slowest ITS changes by < tol."""
    t1 = its[:, 0]
    for k in range(len(lags) - 1):
        if abs(t1[k + 1] - t1[k]) / t1[k] < tol:
            return lags[k]
    return lags[-1]


def choose_n_macro(timescales, max_sets=6):
    """Largest spectral gap t_k / t_(k+1) among the leading timescales."""
    t = np.asarray(timescales[:max_sets])
    ratio = t[:-1] / t[1:]
    return int(np.argmax(ratio)) + 2


def ck_test(dtrajs, msm, lag, memberships, mults=(1, 2, 3, 4, 5)):
    """Chapman-Kolmogorov test on PCCA+ memberships: predicted
    P(k*tau) = coarse-grained T(tau)^k vs. MSMs re-estimated at k*tau.
    Both use the active set of the tau model (states missing at k*tau
    contribute zero weight)."""
    from deeptime.markov import TransitionCountEstimator
    from deeptime.markov.msm import MaximumLikelihoodMSM
    pi, T, M = msm.stationary_distribution, msm.transition_matrix, memberships
    w = (pi[:, None] * M)
    norm = w.sum(0)
    syms = msm.count_model.state_symbols
    pred, est = [], []
    for k in mults:
        Tk = np.linalg.matrix_power(T, k)
        pred.append((w.T @ Tk @ M) / norm[:, None])
        c = TransitionCountEstimator(lag * k, "sliding").fit(dtrajs) \
            .fetch_model().submodel_largest()
        m = MaximumLikelihoodMSM().fit(c).fetch_model()
        # map memberships onto the k*tau active set via state symbols
        pos = {s: i for i, s in enumerate(syms)}
        idx = np.array([pos.get(s, -1) for s in c.state_symbols])
        Mk = np.zeros((len(idx), M.shape[1]))
        Mk[idx >= 0] = M[idx[idx >= 0]]
        wk = m.stationary_distribution[:, None] * Mk
        est.append((wk.T @ m.transition_matrix @ Mk) / wk.sum(0)[:, None])
    return np.array(mults), np.array(pred), np.array(est)


# ---------------------------------------------------------------------------
# 5-6. PCCA+ and Transition Path Theory
# ---------------------------------------------------------------------------
def manual_tpt(T, pi, A, B):
    """Committors and net flux from first principles (general, not assuming
    reversibility), plus the 'committor-difference' shortcut that is exact
    only under detailed balance."""
    n = len(pi)
    A, B = np.asarray(A), np.asarray(B)
    I = np.setdiff1d(np.arange(n), np.r_[A, B])
    qp = np.zeros(n); qp[B] = 1.0
    L = T - np.eye(n)
    qp[I] = np.linalg.solve(L[np.ix_(I, I)], -L[np.ix_(I, B)].sum(1))
    Tb = (T.T * pi[None, :]) / pi[:, None]          # time-reversed chain
    qm = np.zeros(n); qm[A] = 1.0
    Lb = Tb - np.eye(n)
    qm[I] = np.linalg.solve(Lb[np.ix_(I, I)], -Lb[np.ix_(I, A)].sum(1))
    f = pi[:, None] * qm[:, None] * T * qp[None, :]
    np.fill_diagonal(f, 0.0)
    fnet = np.maximum(0.0, f - f.T)
    g = pi[:, None] * T * qp[None, :]
    np.fill_diagonal(g, 0.0)
    short = np.maximum(0.0, g - g.T)
    F = fnet[A][:, np.setdiff1d(np.arange(n), A)].sum()
    return dict(qp=qp, qm=qm, fnet=fnet, short=short, F=F,
                rate=F / (pi @ qm))


def direct_rate(labels, a, b, dt):
    """Trajectory-level A->B rate with no Markov assumption: number of
    A->B transitions divided by the time spent with A as the last-visited
    set (frames outside A and B keep the previous label -- milestoning)."""
    n_ab, t_a = 0, 0
    for d in labels:
        last = None
        for s in d:
            if s == a:
                last = "A"
            elif s == b:
                if last == "A":
                    n_ab += 1
                last = "B"
            if last == "A":
                t_a += 1
    return (n_ab / (t_a * dt) if n_ab else float("nan")), n_ab


def analyze(feats, dt_ps, rama=None, tica_lag=None, n_clusters=100,
            lags=None, msm_lag=None, n_macro=None, source=None, target=None,
            its_samples=50, seed=7, cluster="regspace", run_tpt=True,
            verbose=True):
    """Full MSM + TPT analysis. Returns a dict of results (times in ps)."""
    say = print if verbose else (lambda *a, **k: None)
    n_frames = sum(len(x) for x in feats)
    if lags is None:
        lags = [1, 2, 5, 10, 20, 30, 50, 75, 100]
    tica_lag = tica_lag or lags[3]
    tica, Y, km, dtrajs = tica_cluster(feats, tica_lag, n_clusters, seed=seed,
                                       method=cluster)
    say(f"data: {len(feats)} trajectories, {n_frames} frames, dt={dt_ps} ps "
        f"({n_frames * dt_ps / 1000:.1f} ns), {feats[0].shape[1]} features")
    say(f"TICA (lag {tica_lag * dt_ps:g} ps, kinetic map): "
        f"{tica.output_dimension} dims keep 95% kinetic variance")
    say(f"{cluster} clustering: {km.n_clusters} micro-states")

    its, lo, hi = implied_timescales(dtrajs, lags, n_samples=its_samples)
    say("\nimplied timescales (ps) vs lag:")
    say("  lag(ps) " + " ".join(f"{'t' + str(i + 1):>8}" for i in range(3)))
    for lag, row in zip(lags, its):
        say(f"  {lag * dt_ps:7g} " + " ".join(f"{t * dt_ps:8.1f}"
                                             for t in row[:3]))
    lag = msm_lag or choose_lag(lags, its)
    msm, active = estimate_msm(dtrajs, lag)
    ts = msm.timescales(6) * dt_ps
    gap = choose_n_macro(ts)
    say(f"\nMSM lag = {lag * dt_ps:g} ps "
        f"({'user-set' if msm_lag else 'slowest ITS converged'}); "
        f"{msm.n_states} active states; leading timescales "
        + ", ".join(f"{t:.0f}" for t in ts[:4]) + " ps")
    say(f"largest spectral gap suggests {gap} macro-states; using "
        f"{n_macro or gap}")
    n_macro = n_macro or gap

    pcca = msm.pcca(n_macro)
    M, assign = pcca.memberships, pcca.assignments
    pis = pcca.coarse_grained_stationary_probability
    mults, pred, est = ck_test(dtrajs, msm, lag, M)
    ck_err = float(np.max(np.abs(pred - est)))
    say(f"Chapman-Kolmogorov test (tau..{mults[-1]}tau): max |pred-est| = "
        f"{ck_err:.3f}")

    # macro labels for plotting / direct counting: crisp PCCA+ assignment
    macro_dtrajs = [np.where(a >= 0, assign[np.maximum(a, 0)], -1)
                    for a in active]
    centers = None
    if rama is not None:
        cat_r = np.concatenate(rama); cat_a = np.concatenate(active)
        centers = np.array([_circ_mean(cat_r[cat_a == s])
                            for s in range(msm.n_states)])
    macro_center = None
    if centers is not None:
        macro_center = np.array([_circ_mean(centers[assign == k],
                                            weights=msm.stationary_distribution
                                            [assign == k])
                                 for k in range(n_macro)])

    # source/target from the slowest eigenvector psi2 (pi-orthogonal to 1, so
    # its sign splits the slowest process): B = the minority side's extreme,
    # A = the most populated macro-state on the opposite side.
    if source is None or target is None:
        psi2 = msm.eigenvectors_right()[:, 1]
        w = msm.stationary_distribution[:, None] * M
        score = (w * psi2[:, None]).sum(0) / w.sum(0)
        lo_, hi_ = int(np.argmin(score)), int(np.argmax(score))
        side = np.sign(score)
        b_ = hi_ if pis[side > 0].sum() < pis[side < 0].sum() else lo_
        opp = np.where(side == -side[b_])[0]
        a_ = int(opp[np.argmax(pis[opp])])
        source = a_ if source is None else source
        target = b_ if target is None else target

    res = dict(lags=np.array(lags), its=its, its_lo=lo, its_hi=hi,
               dt_ps=dt_ps, lag=lag, msm=msm, n_macro=n_macro, pcca=pcca,
               assign=assign, pis=pis, ck=(mults, pred, est), ck_err=ck_err,
               active=active, macro_dtrajs=macro_dtrajs, centers=centers,
               macro_center=macro_center, Y=Y, rama=rama, timescales_ps=ts)
    kT = 0.0083144626 * 300.0
    say(f"\nmacro-states (PCCA+):")
    for k in range(n_macro):
        loc = (f" at (phi,psi)=({macro_center[k][0]:.0f},"
               f"{macro_center[k][1]:.0f})" if macro_center is not None
               else "")
        say(f"  M{k}: pi={pis[k]:.4f}  dG={-kT * np.log(pis[k] / pis.max()):5.2f} "
            f"kJ/mol{loc}")
    if run_tpt:
        res.update(tpt_analysis(res, source, target, verbose=verbose))
    return res


def tpt_analysis(r, source, target, core=0.9, verbose=True):
    """TPT between macro-states `source` (A) and `target` (B) of an
    analyze() result: committors, net flux, rates both ways (MSM and direct
    trajectory counting), q+=1/2 dividing-surface edges, macro pathways.

    A and B are CORE sets: micro-states whose PCCA+ membership in the source
    / target macro-state is >= `core`. The barrier region in between is left
    to the committor, so q+ varies smoothly and the q+=1/2 surface is
    resolved. (Crisp PCCA+ partitions would make A and B adjacent: q+ would
    jump 0->1 on every edge.) The direct rate is counted on the raw frames
    with transition-based (milestoning) assignment to the same cores, which
    is independent of the MSM lag -- a genuine test of the Markov model."""
    say = print if verbose else (lambda *a, **k: None)
    msm, assign, lag, dt_ps = r["msm"], r["assign"], r["lag"], r["dt_ps"]
    macro_dtrajs, centers, n_macro = (r["macro_dtrajs"], r["centers"],
                                      r["n_macro"])
    M = r["pcca"].memberships
    A = np.where(M[:, source] >= core)[0]; B = np.where(M[:, target] >= core)[0]
    if not len(A):
        A = np.where(assign == source)[0]
    if not len(B):
        B = np.where(assign == target)[0]
    lab = np.full(msm.n_states, -1); lab[A] = source; lab[B] = target
    core_dtrajs = [np.where(a >= 0, lab[np.maximum(a, 0)], -1)
                   for a in r["active"]]
    flux = msm.reactive_flux(A, B)
    back = msm.reactive_flux(B, A)
    man = manual_tpt(msm.transition_matrix, msm.stationary_distribution, A, B)
    rate_ps = flux.rate / (lag * dt_ps)
    rate_back_ps = back.rate / (lag * dt_ps)
    d_lib = float(np.max(np.abs(man["fnet"] - flux.net_flux)))
    d_short = float(np.max(np.abs(man["short"] - man["fnet"])))
    k_direct, n_events = direct_rate(core_dtrajs, source, target, dt_ps)
    k_back_direct, n_back = direct_rate(core_dtrajs, target, source, dt_ps)

    # rate-limiting location: every A->B reactive trajectory crosses the
    # q+ = 1/2 isocommittor surface, so the net flux through the edges that
    # cut it equals F exactly (flux conservation). The edges carrying most
    # of it are the transition-state region.
    q = flux.forward_committor; f = flux.net_flux
    cut = [(i, j) for i, j in zip(*np.nonzero(f))
           if q[i] < 0.5 <= q[j]]
    cut.sort(key=lambda e: -f[e])
    cut_total = float(sum(f[e] for e in cut))
    ts_states = np.unique([k for e in cut[:5] for k in e])

    sets = [list(np.where(assign == k)[0]) for k in range(n_macro)]
    cg_sets, cg = flux.coarse_grain(sets)
    cg_order = [int(assign[next(iter(s_))]) for s_ in cg_sets]
    Aset, Bset = set(A.tolist()), set(B.tolist())
    cg_lab = ["A" if set(s_) <= Aset else "B" if set(s_) <= Bset
              else f"M{k}~" for s_, k in zip(cg_sets, cg_order)]
    cg_net = np.zeros((n_macro, n_macro))
    for i, si in enumerate(cg_order):
        for j, sj in enumerate(cg_order):
            if si != sj:
                cg_net[si, sj] += cg.net_flux[i, j]
    cpaths, ccaps = cg.pathways(fraction=0.999)
    cpaths = [[cg_lab[k] for k in pth] for pth in cpaths]
    ccaps = np.asarray(ccaps) / cg.total_flux

    say(f"\nTPT  A=core(M{source}) -> B=core(M{target}): |A|={len(A)}, "
        f"|B|={len(B)}, intermediate={msm.n_states - len(A) - len(B)} "
        f"micro-states (core = PCCA+ membership >= {core})")
    say(f"  total reactive flux F = {flux.total_flux:.3e} per lag step")
    say(f"  MSM-TPT : MFPT A->B = {_fmt_t(1 / rate_ps):>9}   "
        f"B->A = {_fmt_t(1 / rate_back_ps):>9}")
    say(f"  direct  : MFPT A->B = {_fmt_t(1 / k_direct):>9}   "
        f"B->A = {_fmt_t(1 / k_back_direct):>9}   ({n_events:.0f} A->B "
        f"core-to-core events in the raw {dt_ps:g} ps frames; no Markov "
        f"assumption, no lag)")
    say(f"  manual TPT vs deeptime net flux: max|diff| = {d_lib:.1e}")
    say(f"  committor-difference form pi_i T_ij (q+_j - q+_i) vs general "
        f"form pi_i q-_i T_ij q+_j: max|diff| = {d_short:.1e} "
        f"(equal: reversible MSM)")
    say(f"  flux through the q+=1/2 surface = {cut_total / flux.total_flux:.4f}"
        f" x F over {len(cut)} edges (conservation check)")
    for e in cut[:3]:
        loc = (f" (phi,psi) ({centers[e[0]][0]:.0f},{centers[e[0]][1]:.0f})"
               f" -> ({centers[e[1]][0]:.0f},{centers[e[1]][1]:.0f})"
               if centers is not None else "")
        say(f"    edge {e[0]:>3}->{e[1]:<3} q+ {q[e[0]]:.2f}->{q[e[1]]:.2f} "
            f"carries {f[e] / flux.total_flux * 100:4.1f}% of F{loc}")
    say("  coarse-grained pathways (A/B = cores, Mk~ = non-core part of "
        "macro-state k; fraction of F):")
    for pth, c in zip(cpaths, ccaps):
        if c > 0.005:
            say("    " + " -> ".join(pth) + f" : {c:.2f}")

    if n_events < 10:
        say(f"  WARNING: only {n_events:.0f} A->B events were sampled; this "
            f"rate is order-of-magnitude only (an MSM cannot create "
            f"transitions the data never visited)")
    toy = shortcut_counterexample()
    say(f"  driven chains without detailed balance (20 random): the "
        f"committor-difference form keeps the total rate but misplaces edge "
        f"fluxes by {toy[0] * 100:.0f}% median / {toy[1] * 100:.0f}% max "
        f"-> pathways need q- in general")

    return dict(source=source, target=target, A=A, B=B, flux=flux,
                rate_ps=rate_ps,
                rate_back_ps=rate_back_ps, k_direct=k_direct,
                k_back_direct=k_back_direct, n_events=n_events, d_lib=d_lib,
                d_short=d_short, d_toy=toy, cut=cut, cut_total=cut_total,
                ts_states=ts_states, cg_net=cg_net, cpaths=cpaths,
                ccaps=ccaps)


def shortcut_counterexample(n_chains=20, n=6, seed=0):
    """Generic driven chains (random row-stochastic T, no detailed balance),
    A={0}, B={n-1}. The committor-difference form always reproduces the total
    flux out of A (there q+=0, q-=1), but not the edge-by-edge flux map.
    Returns (median, max) relative edge error; general TPT is checked
    against deeptime on every chain."""
    from deeptime.markov.msm import MarkovStateModel
    rng = np.random.default_rng(seed)
    errs = []
    for _ in range(n_chains):
        T = rng.dirichlet(np.full(n, 0.5), size=n)
        w, v = np.linalg.eig(T.T)
        pi = np.real(v[:, np.argmax(np.real(w))]); pi /= pi.sum()
        m = manual_tpt(T, pi, [0], [n - 1])
        ref = MarkovStateModel(T).reactive_flux([0], [n - 1])
        assert np.allclose(ref.net_flux, m["fnet"], atol=1e-12)
        assert np.isclose(m["short"][0].sum(), m["fnet"][0].sum())
        errs.append(np.abs(m["short"] - m["fnet"]).max() / m["fnet"].max())
    return float(np.median(errs)), float(np.max(errs))


def _circ_mean(x, weights=None):
    """Circular mean of angles in degrees, column-wise."""
    x = np.radians(np.atleast_2d(x))
    w = np.ones(len(x)) if weights is None else np.asarray(weights)
    s = (w[:, None] * np.sin(x)).sum(0); c = (w[:, None] * np.cos(x)).sum(0)
    return np.degrees(np.arctan2(s, c))


# ---------------------------------------------------------------------------
# 7. figure, demo report, CLI
# ---------------------------------------------------------------------------
_COLORS = ["#2a6fdb", "#e4572e", "#29a36a", "#9b59b6", "#f2a541", "#17becf",
           "#7f8c8d", "#c0392b"]


def figure(r, path, rmsd=None, title="MSM + TPT conformational flux"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch
    dt, lag, n = r["dt_ps"], r["lag"], r["n_macro"]
    msm, flux, assign = r["msm"], r["flux"], r["assign"]
    rama = r["rama"] is not None
    fig, ax = plt.subplots(2, 3, figsize=(16, 9.6))
    ax = ax.ravel()

    # (a) free energy landscape with macro-state centres
    if rama:
        R = np.concatenate(r["rama"])
        H, xe, ye = np.histogram2d(R[:, 0], R[:, 1], bins=90,
                                   range=[[-180, 180], [-180, 180]])
        xlab, ylab = "phi (deg)", "psi (deg)"
        pts = r["centers"]
    else:
        R = np.concatenate(r["Y"])[:, :2]
        H, xe, ye = np.histogram2d(R[:, 0], R[:, 1], bins=90)
        xlab, ylab = "TIC 1", "TIC 2"
        cat_a = np.concatenate(r["active"])
        pts = np.array([R[cat_a == s].mean(0) for s in range(msm.n_states)])
    with np.errstate(divide="ignore"):
        F = -np.log(H.T / H.sum())
    F -= F[np.isfinite(F)].min()
    F *= 0.0083144626 * 300
    cs = ax[0].contourf(0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1]),
                        F, levels=np.linspace(0, 20, 21), cmap="Greys_r",
                        extend="max")
    for k in range(n):
        sel = assign == k
        ax[0].scatter(pts[sel, 0], pts[sel, 1], s=14, c=_COLORS[k % 8],
                      label=f"M{k} (pi={r['pis'][k]:.2f})", edgecolor="none")
    fig.colorbar(cs, ax=ax[0], label="free energy (kJ/mol)")
    ax[0].set_xlabel(xlab); ax[0].set_ylabel(ylab)
    ax[0].legend(fontsize=7, loc="lower right")
    ax[0].set_title("(a) free energy + PCCA+ macro-states\n"
                    "(dots = micro-state centres)")

    # (b) implied timescales
    lags = r["lags"] * dt
    for i in range(min(4, r["its"].shape[1])):
        ax[1].semilogy(lags, r["its"][:, i] * dt, "o-", ms=3,
                       label=f"t{i + 1}")
        if r["its_lo"] is not None:
            ax[1].fill_between(lags, r["its_lo"][:, i] * dt,
                               r["its_hi"][:, i] * dt, alpha=0.2)
    ax[1].fill_between(lags, 1e-3, lags, color="k", alpha=0.15)
    ax[1].axvline(lag * dt, color="r", ls="--", lw=1)
    ax[1].set_ylim(bottom=max(0.5 * lags.min(), 1e-2))
    ax[1].set_xlabel("lag time (ps)"); ax[1].set_ylabel("implied timescale (ps)")
    ax[1].legend(fontsize=7)
    ax[1].set_title(f"(b) implied timescales (95% Bayesian band)\n"
                    f"chosen lag {lag * dt:g} ps (red); grey = unresolvable")

    # (c) Chapman-Kolmogorov test
    mults, pred, est = r["ck"]
    for k in range(n):
        ax[2].plot(mults * lag * dt, pred[:, k, k], "-", c=_COLORS[k % 8],
                   label=f"M{k}->M{k} MSM prediction")
        ax[2].plot(mults * lag * dt, est[:, k, k], "o", c=_COLORS[k % 8],
                   mfc="none")
    ax[2].set_xlabel("time (ps)"); ax[2].set_ylabel("P(stay in macro-state)")
    ax[2].set_ylim(0, 1.05); ax[2].legend(fontsize=7)
    ax[2].set_title(f"(c) Chapman-Kolmogorov test\nlines = T(tau)^k, circles "
                    f"= re-estimated; max err {r['ck_err']:.3f}")

    # (d) committor
    q = flux.forward_committor
    sc = ax[3].scatter(pts[:, 0], pts[:, 1], c=q, cmap="coolwarm", s=22,
                       vmin=0, vmax=1, edgecolor="none")
    ts = r["ts_states"]
    ax[3].scatter(pts[ts, 0], pts[ts, 1], s=70, facecolor="none",
                  edgecolor="k", lw=1, label="ends of top q+=1/2 crossing edges")
    for e in r["cut"][:3]:
        ax[3].annotate("", xy=pts[e[1]], xytext=pts[e[0]],
                       arrowprops=dict(arrowstyle="->", color="k", lw=1.5))
    fig.colorbar(sc, ax=ax[3], label="forward committor q+")
    if rama:
        ax[3].set_xlim(-180, 180); ax[3].set_ylim(-180, 180)
    ax[3].set_xlabel(xlab); ax[3].set_ylabel(ylab)
    ax[3].legend(fontsize=7, loc="lower right")
    ax[3].set_title(f"(d) committor q+ for M{r['source']} -> M{r['target']}\n"
                    "arrows = edges carrying most flux across q+ = 1/2")

    # (e) coarse-grained flux network
    a = ax[4]; cg = r["cg_net"] / flux.total_flux
    if r["macro_center"] is not None:
        P = r["macro_center"].copy()
    else:
        ang = np.linspace(0, 2 * np.pi, n, endpoint=False)
        P = np.c_[np.cos(ang), np.sin(ang)]
    for k in range(n):
        a.scatter(*P[k], s=6000 * r["pis"][k] + 200, c=_COLORS[k % 8],
                  alpha=0.75, edgecolor="k", zorder=3)
        a.text(*P[k], f"M{k}\n{r['pis'][k]:.2f}", ha="center", va="center",
               fontsize=8, zorder=4)
    for i in range(n):
        for j in range(n):
            if cg[i, j] > 0.01:
                a.add_patch(FancyArrowPatch(
                    P[i], P[j], arrowstyle="-|>", mutation_scale=18,
                    lw=1 + 10 * cg[i, j], color="#444", alpha=0.8,
                    connectionstyle="arc3,rad=0.15", shrinkA=22, shrinkB=22,
                    zorder=2))
                mid = 0.5 * (P[i] + P[j])
                a.text(*mid, f"{cg[i, j] * 100:.0f}%", fontsize=8,
                       color="#a00", ha="center")
    pad = 0.35 * (np.ptp(P, 0).max() + 1e-9)
    a.set_xlim(P[:, 0].min() - pad, P[:, 0].max() + pad)
    a.set_ylim(P[:, 1].min() - pad, P[:, 1].max() + pad)
    a.set_xlabel(xlab if rama else ""); a.set_ylabel(ylab if rama else "")
    a.set_title(f"(e) net reactive flux network M{r['source']} -> "
                f"M{r['target']}\nnode = population, arrow = % of total flux; "
                f"MFPT A->B {_fmt_t(1 / r['rate_ps'])}")

    # (f) RMSD is degenerate: distinct macro-states share RMSD values
    a = ax[5]
    if rmsd is not None:
        Rm = np.concatenate(rmsd) * 10
        lab = np.concatenate(r["macro_dtrajs"])
        bins = np.linspace(0, Rm.max(), 60)
        for k in range(n):
            if (lab == k).sum():
                a.hist(Rm[lab == k], bins=bins, density=True, alpha=0.5,
                       color=_COLORS[k % 8], label=f"M{k}")
        ov = rmsd_overlap(rmsd, r["macro_dtrajs"], n)
        i, j = np.unravel_index(np.argmax(ov), ov.shape)
        a.set_xlabel("heavy-atom RMSD to the starting structure (A)")
        a.set_ylabel("density"); a.legend(fontsize=7)
        a.set_title(f"(f) RMSD cannot tell states apart:\nM{i} and M{j} "
                    f"overlap {ov[i, j] * 100:.0f}% in RMSD")
    else:
        a.axis("off")
    fig.suptitle(title, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(path, dpi=120); plt.close(fig)


def _fmt_t(ps):
    return f"{ps / 1e3:.2f} ns" if ps >= 1e3 else f"{ps:.0f} ps"


def _rmsd(top, xtcs):
    """Heavy-atom RMSD (nm) of every frame to the topology structure."""
    import mdtraj as md
    ref = md.load(top)
    heavy = ref.topology.select("not element H")
    return [md.rmsd(md.load(x, top=top), ref, atom_indices=heavy)
            for x in xtcs]


def rmsd_overlap(rmsd, macro_dtrajs, n, bins=60):
    """Overlap coefficient sum_x min(p_a(x), p_b(x)) of the RMSD histograms
    of each pair of macro-states (1 = indistinguishable by RMSD)."""
    R = np.concatenate(rmsd); lab = np.concatenate(macro_dtrajs)
    edges = np.linspace(R.min(), R.max(), bins + 1)
    H = [np.histogram(R[lab == k], edges)[0] / max(1, (lab == k).sum())
         for k in range(n)]
    ov = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            ov[i, j] = ov[j, i] = np.minimum(H[i], H[j]).sum()
    return ov


def _ala2_basins(macro_center):
    """Label alanine-dipeptide macro-states by their (phi, psi) centre."""
    lab = {}
    for k, (phi, psi) in enumerate(macro_center):
        if phi > 0:
            lab[k] = "alphaL"
        elif -120 < psi < 50:
            lab[k] = "alphaR"
        else:
            lab[k] = "beta/PII"
    return lab


def _ala2_channels(r):
    """Split the q+=1/2 crossing flux of beta->alphaR by where psi crosses:
    descending through psi ~ +60..120 or wrapping through psi = +-180."""
    f, c = r["flux"].net_flux, r["centers"]
    out = {"psi descends through +60..+120 (direct)": 0.0,
           "psi wraps through +-180 (via -120)": 0.0}
    for i, j in r["cut"]:
        psi_mid = _circ_mean(np.array([c[i], c[j]]))[1]
        key = ("psi descends through +60..+120 (direct)" if psi_mid > 0
               else "psi wraps through +-180 (via -120)")
        out[key] += f[i, j]
    tot = sum(out.values())
    return {k: v / tot for k, v in out.items()}


def report(fig_path=os.path.join(HERE, "..", "figures", "m23_msm_tpt.png"),
           n_traj=4, ns=50.0):
    print("=" * 70)
    print("MODULE 23 — Conformational flux: MSM + TPT on real MD")
    print("=" * 70)
    print("system: alanine dipeptide, amber99sb-ildn + OBC implicit solvent, "
          "300 K, OpenMM\n")
    top, xtcs = simulate_alanine_dipeptide(n_traj=n_traj, ns=ns)
    feats, names, rama = featurize(top, xtcs, "backbone")
    print("features: " + ", ".join(names))
    # 3 macro-states (beta/PII, alphaR, alphaL): the 2-state split chosen by
    # the largest spectral gap would lump beta and alphaR together.
    r = analyze(feats, dt_ps=1.0, rama=rama, n_clusters=150, n_macro=3,
                run_tpt=False)
    lab = _ala2_basins(r["macro_center"])
    print("  basins: " + ", ".join(f"M{k}={v}" for k, v in lab.items()))
    idx = {v: k for k, v in lab.items()}
    beta, aR, aL = idx.get("beta/PII"), idx.get("alphaR"), idx.get("alphaL")

    print("\n--- primary: beta/PII -> alphaR (well sampled) ---")
    r.update(tpt_analysis(r, beta, aR))
    ch = _ala2_channels(r)
    print("  reaction channels (flux across q+=1/2, by psi at the crossing):")
    for k, v in ch.items():
        print(f"    {k:<38}: {v * 100:5.1f}% of F")
    rare = None
    if aL is not None:
        print("\n--- secondary: beta/PII -> alphaL (rare event) ---")
        rare = tpt_analysis(r, beta, aL)
    rmsd = _rmsd(top, xtcs)
    ov = rmsd_overlap(rmsd, r["macro_dtrajs"], r["n_macro"])
    print("\nRMSD degeneracy (overlap of per-state RMSD histograms, "
          "1 = indistinguishable):")
    for i in range(r["n_macro"]):
        for j in range(i + 1, r["n_macro"]):
            print(f"  {lab[i]:>8} vs {lab[j]:<8}: {ov[i, j]:.2f}")
    figure(r, fig_path, rmsd=rmsd,
           title="Alanine dipeptide: MSM + Transition Path Theory "
                 "(conformational flux, not RMSD)")
    print(f"\nfigure written to: {fig_path}")
    return dict(lag_ps=r["lag"] * r["dt_ps"], labels=lab,
                timescales_ps=r["timescales_ps"], pis=r["pis"],
                ck_err=r["ck_err"], d_short=r["d_short"],
                mfpt_ps=1 / r["rate_ps"], mfpt_direct_ps=1 / r["k_direct"],
                n_events=r["n_events"], channels=ch, rmsd_overlap=ov,
                cut_total=r["cut_total"] / r["flux"].total_flux,
                rare_mfpt_ns=(1 / rare["rate_ps"] / 1e3 if rare else None),
                rare_events=(rare["n_events"] if rare else None))


def main(argv=None):
    p = argparse.ArgumentParser(
        description="MSM + TPT conformational flux from MD trajectories "
                    "(GROMACS .xtc/.trr, AMBER .nc/.dcd, ...).")
    p.add_argument("--top", required=True,
                   help="topology readable by mdtraj (.gro/.pdb/.prmtop)")
    p.add_argument("--traj", nargs="+", required=True,
                   help="one or more trajectory files")
    p.add_argument("--dt-ps", type=float, required=True,
                   help="time between saved frames (ps), before --stride")
    p.add_argument("--stride", type=int, default=1)
    p.add_argument("--features", default="backbone",
                   choices=["backbone", "contacts", "both"])
    p.add_argument("--n-clusters", type=int, default=200)
    p.add_argument("--cluster", default="regspace",
                   choices=["regspace", "kmeans"])
    p.add_argument("--tica-lag", type=int, default=None,
                   help="TICA lag in frames (default: 4th entry of --lags)")
    p.add_argument("--lags", type=int, nargs="+", default=None,
                   help="lag times (frames) scanned for implied timescales")
    p.add_argument("--msm-lag", type=int, default=None,
                   help="MSM lag (frames); default: ITS convergence")
    p.add_argument("--n-macro", type=int, default=None,
                   help="PCCA+ macro-states; default: largest spectral gap")
    p.add_argument("--source", type=int, default=None,
                   help="source macro-state index A")
    p.add_argument("--target", type=int, default=None,
                   help="target macro-state index B")
    p.add_argument("--out", default="msm_tpt.png")
    a = p.parse_args(argv)
    feats, names, rama = featurize(a.top, a.traj, a.features, a.stride)
    r = analyze(feats, dt_ps=a.dt_ps * a.stride, rama=rama,
                tica_lag=a.tica_lag, n_clusters=a.n_clusters, lags=a.lags,
                msm_lag=a.msm_lag, n_macro=a.n_macro, source=a.source,
                target=a.target, cluster=a.cluster)
    figure(r, a.out, rmsd=_rmsd(a.top, a.traj) if a.stride == 1 else None)
    print(f"\nfigure written to: {a.out}")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        main()
    else:
        report()
