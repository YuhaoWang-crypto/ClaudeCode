"""
Module 25 — Retrospective test: does the MSM perturbation scan (M24) find
KNOWN allosteric drug sites on real proteins?

Systems (two directions of allosteric control):
  SHP2 phosphatase  : inactive autoinhibited 4DGP (apo) <-> active open 6CRF.
                      Known allosteric INHIBITOR SHP099 (PDB 5EHR) -- used
                      ONLY to score the prediction, never to build the model.
  Glucokinase (GK)  : inactive super-open 1V4T <-> active closed 1V4S.
                      Known allosteric ACTIVATOR (MRK, in 1V4S). The only
                      closed structure carries the activator, so the active
                      basin is built from a ligand-bound geometry (mild
                      leakage; protein coordinates only, ligand removed).

Ground truth = residues with any heavy atom within 4.5 A of the drug in the
deposited co-crystal structure.

Model: all-atom explicit-solvent MSMs of ~450-500 residue proteins need
microseconds on GPUs; this module therefore uses a dual-basin C-alpha
structure-based (Go) model built from the two experimental states
(Okazaki-Takada exponential mixing). It tests whether flux / sensitivity
analysis locates the known sites GIVEN the conformational landscape; it does
not test force-field accuracy.

Pipeline: simulate many transitions -> M23 MSM on contact features ->
2 macro-states (inactive / active) -> for every residue r, a virtual
ligand that stabilises the local contact network of the state the drug is
expected to prefer (inactive for an inhibitor, active for an activator) ->
predicted d ln K(active/inactive) and d ln k(inactive->active) by M24
transition-state reweighting -> rank residues -> ROC-AUC / enrichment
against the drug-contact residues, vs. structural and RMSF baselines and a
random-site null.
"""
import os
import urllib.request

import numpy as np

from grn_pipeline import m23_msm_tpt as m23
from grn_pipeline import m24_msm_validation as m24

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "..", "figures")
CACHE = os.path.join(FIG, "_retro_cache")

SYSTEMS = {
    "SHP2": dict(inactive=("4DGP", "A"), active=("6CRF", "A"),
                 holo=("5EHR", "A", "5OD"), kind="inhibitor"),
    "GK": dict(inactive=("1V4T", "A"), active=("1V4S", "A"),
               holo=("1V4S", "A", "MRK"), kind="activator"),
}


# ---------------------------------------------------------------------------
# structures
# ---------------------------------------------------------------------------
def fetch(pdb_id):
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, f"{pdb_id}.pdb")
    if not os.path.exists(path):
        urllib.request.urlretrieve(
            f"https://files.rcsb.org/download/{pdb_id}.pdb", path)
    return path


def read_chain(pdb_id, chain):
    """Residue number -> dict(ca=xyz, heavy=[xyz...]) for the first model,
    altloc A/blank, standard residues only (ATOM records)."""
    res = {}
    for line in open(fetch(pdb_id)):
        if line.startswith("ENDMDL"):
            break
        if not line.startswith("ATOM") or line[21] != chain:
            continue
        if line[16] not in " A" or line[76:78].strip() == "H":
            continue
        rn = int(line[22:26]); name = line[12:16].strip()
        xyz = np.array([float(line[30:38]), float(line[38:46]),
                        float(line[46:54])]) / 10.0          # nm
        d = res.setdefault(rn, dict(ca=None, heavy=[]))
        d["heavy"].append(xyz)
        if name == "CA":
            d["ca"] = xyz
    return {k: v for k, v in res.items() if v["ca"] is not None}


def ligand_site(pdb_id, chain, resname, cutoff=0.45):
    """Residue numbers with a heavy atom within `cutoff` nm of the ligand."""
    lig = []
    for line in open(fetch(pdb_id)):
        if line.startswith("HETATM") and line[17:20] == resname \
                and line[21] == chain:
            lig.append([float(line[30:38]), float(line[38:46]),
                        float(line[46:54])])
    lig = np.array(lig) / 10.0
    prot = read_chain(pdb_id, chain)
    return sorted(r for r, d in prot.items()
                  if np.min(np.linalg.norm(np.array(d["heavy"])[:, None]
                                           - lig[None], axis=2)) < cutoff)


# ---------------------------------------------------------------------------
# dual-basin C-alpha structure-based model
# ---------------------------------------------------------------------------
def _dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n; w = b2 - np.dot(b2, b1n) * b1n
    return np.arctan2(np.dot(np.cross(b1n, v), w), np.dot(v, w))


def _angle(p0, p1, p2):
    u, v = p0 - p1, p2 - p1
    return np.arccos(np.clip(np.dot(u, v) / np.linalg.norm(u)
                             / np.linalg.norm(v), -1, 1))


def _contacts(res, ids, cutoff=0.45, min_sep=3):
    """Native contacts: any heavy-atom pair within `cutoff` nm."""
    ca = np.array([res[r]["ca"] for r in ids])
    heavy = [np.array(res[r]["heavy"]) for r in ids]
    out = set()
    for i in range(len(ids)):
        near = np.where(np.linalg.norm(ca - ca[i], axis=1) < 1.4)[0]
        for j in near[near > i]:
            if ids[j] - ids[i] < min_sep:
                continue
            d = np.linalg.norm(heavy[i][:, None] - heavy[j][None], axis=2)
            if d.min() < cutoff:
                out.add((i, j))
    return out


class DualBasin:
    """Structure data + contact sets for a two-state Go model."""

    def __init__(self, name):
        s = SYSTEMS[name]
        a, b = read_chain(*s["inactive"]), read_chain(*s["active"])
        self.name, self.kind = name, s["kind"]
        self.ids = sorted(set(a) & set(b))
        self.X1 = np.array([a[r]["ca"] for r in self.ids])     # inactive
        self.X2 = np.array([b[r]["ca"] for r in self.ids])     # active
        c1, c2 = _contacts(a, self.ids), _contacts(b, self.ids)
        d1 = lambda i, j: np.linalg.norm(self.X1[i] - self.X1[j])
        d2 = lambda i, j: np.linalg.norm(self.X2[i] - self.X2[j])
        self.common = [(i, j, 0.5 * (d1(i, j) + d2(i, j)))
                       for (i, j) in sorted(c1 & c2)
                       if abs(d1(i, j) - d2(i, j)) < 0.1]
        cset = {(i, j) for i, j, _ in self.common}
        self.spec1 = [(i, j, d1(i, j)) for (i, j) in sorted(c1 - cset)]
        self.spec2 = [(i, j, d2(i, j)) for (i, j) in sorted(c2 - cset)]
        self.site = ligand_site(*s["holo"])
        self.site_idx = [self.ids.index(r) for r in self.site
                         if r in self.ids]

    def summary(self):
        return (f"{self.name}: {len(self.ids)} beads, common contacts "
                f"{len(self.common)}, inactive-specific {len(self.spec1)}, "
                f"active-specific {len(self.spec2)}, drug site "
                f"{len(self.site_idx)} residues")


def build_system(m, eps=1.0, kT=2.494, bmix=None, dV=0.0, kd=1.0):
    """OpenMM System: bonds/angles (mean of both states), common contacts,
    and exponential mixing of the two state-specific potentials
        V_mix = -1/b ln(exp(-b V1) + exp(-b (V2 + dV)))
    written in a numerically stable min/log form. eps, kd in units of kT."""
    import openmm as mm
    e = eps * kT
    n = len(m.ids)
    sys_ = mm.System()
    for _ in range(n):
        sys_.addParticle(110.0)
    ids = m.ids
    bond = mm.HarmonicBondForce()
    ang = mm.HarmonicAngleForce()
    excl = set()
    for k in range(n - 1):
        gap = ids[k + 1] - ids[k]
        if gap > 12:
            continue
        r0 = 0.5 * (np.linalg.norm(m.X1[k + 1] - m.X1[k])
                    + np.linalg.norm(m.X2[k + 1] - m.X2[k]))
        bond.addBond(k, k + 1, r0, 2e4 if gap == 1 else 2e3)
        excl.add((k, k + 1))
    for k in range(n - 2):
        if ids[k + 2] - ids[k] == 2:
            th = 0.5 * (_angle(*m.X1[k:k + 3]) + _angle(*m.X2[k:k + 3]))
            ang.addAngle(k, k + 1, k + 2, th, 40.0 * kT)
            excl.add((k, k + 2))
    sys_.addForce(bond); sys_.addForce(ang)

    lj = "e*(5*(s/r)^12-6*(s/r)^10)"

    def contact_force(pairs):
        f = mm.CustomBondForce(lj)
        f.addGlobalParameter("e", e); f.addPerBondParameter("s")
        for i, j, s_ in pairs:
            f.addBond(i, j, [s_])
        return f

    sys_.addForce(contact_force(m.common))

    def state_force(pairs, X):
        f1 = contact_force(pairs)
        dih = mm.CustomTorsionForce(
            "kd*((1-cos(theta-t0))+0.5*(1-cos(3*(theta-t0))))")
        dih.addGlobalParameter("kd", kd * kT); dih.addPerTorsionParameter("t0")
        for k in range(n - 3):
            if ids[k + 3] - ids[k] == 3:
                dih.addTorsion(k, k + 1, k + 2, k + 3,
                               [_dihedral(*X[k:k + 4])])
        return f1, dih

    bmix = bmix if bmix is not None else 1.0 / (20 * kT)
    cv = mm.CustomCVForce("min(V1,V2p)-log(1+exp(-b*abs(V1-V2p)))/b;"
                          "V2p=V2+dV; V1=c1+t1; V2=c2+t2")
    c1, t1 = state_force(m.spec1, m.X1)
    c2, t2 = state_force(m.spec2, m.X2)
    for nm, f in [("c1", c1), ("t1", t1), ("c2", c2), ("t2", t2)]:
        cv.addCollectiveVariable(nm, f)
    cv.addGlobalParameter("b", bmix); cv.addGlobalParameter("dV", dV * kT)
    sys_.addForce(cv)

    rep = mm.CustomNonbondedForce("er*(0.4/r)^12")
    rep.addGlobalParameter("er", e)
    rep.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic)
    rep.setCutoffDistance(0.8)
    for _ in range(n):
        rep.addParticle([])
    for i, j, _ in m.common + m.spec1 + m.spec2:
        excl.add((i, j))
    for i, j in excl:
        rep.addExclusion(i, j)
    sys_.addForce(rep)
    return sys_


def build_gauss(m, n_copies=1, eps=1.0, lam=1.0, kT=2.494, w=0.05,
                kd=0.5, spacing=60.0, alpha=1.0, ligand=None):
    """Faster dual-basin form without CustomCVForce (Lammert, Schug &
    Onuchic 2009 Gaussian contacts). State-specific contacts are scaled by
    `alpha` (< 1 lowers the switching barrier while the common contacts keep
    the fold):
      common / single-state contact :  e[(1+(r_rep/r)^12)(1-G(r1)) - 1]
      contact formed in BOTH states at different distances:
                                       e[(1+(r_rep/r)^12)(1-G(r1))(1-G(r2))-1]
    G(rk) = exp(-(r-rk)^2/(2 w^2)). Active-state-specific contacts are
    scaled by `lam` to set the inactive/active balance. Torsions get a
    double-well -kd ln(exp(-U1)+exp(-U2)). `n_copies` non-interacting copies
    share one OpenMM context (amortises per-step overhead)."""
    import openmm as mm
    e = eps * kT
    n = len(m.ids); ids = m.ids
    sys_ = mm.System()
    for _ in range(n * n_copies):
        sys_.addParticle(110.0)
    bond = mm.HarmonicBondForce(); ang = mm.HarmonicAngleForce()
    dih = mm.CustomTorsionForce(
        "-kd*log(exp(-(1-cos(theta-t1))-0.5*(1-cos(3*(theta-t1))))"
        "+exp(-(1-cos(theta-t2))-0.5*(1-cos(3*(theta-t2)))))")
    dih.addGlobalParameter("kd", kd * kT)
    dih.addPerTorsionParameter("t1"); dih.addPerTorsionParameter("t2")
    con = mm.CustomBondForce(
        "e*s*((1+(0.4/r)^12)*(1-exp(-(r-r1)^2/(2*w^2)))"
        "*(1-two*exp(-(r-r2)^2/(2*w^2)))-1)")
    con.addGlobalParameter("e", e); con.addGlobalParameter("w", w)
    for p in ("r1", "r2", "two", "s"):
        con.addPerBondParameter(p)
    rep = mm.CustomNonbondedForce("er*(0.4/r)^12")
    rep.addGlobalParameter("er", e)
    rep.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic)
    rep.setCutoffDistance(0.6)
    for _ in range(n * n_copies):
        rep.addParticle([])
    # template terms for one copy
    bonds, angles, tors, excl = [], [], [], set()
    for k in range(n - 1):
        gap = ids[k + 1] - ids[k]
        if gap <= 12:
            r0 = 0.5 * (np.linalg.norm(m.X1[k + 1] - m.X1[k])
                        + np.linalg.norm(m.X2[k + 1] - m.X2[k]))
            bonds.append((k, k + 1, r0, 2e4 if gap == 1 else 2e3))
            excl.add((k, k + 1))
    for k in range(n - 2):
        if ids[k + 2] - ids[k] == 2:
            th = 0.5 * (_angle(*m.X1[k:k + 3]) + _angle(*m.X2[k:k + 3]))
            angles.append((k, k + 1, k + 2, th)); excl.add((k, k + 2))
    for k in range(n - 3):
        if ids[k + 3] - ids[k] == 3:
            tors.append((k, k + 1, k + 2, k + 3, _dihedral(*m.X1[k:k + 4]),
                         _dihedral(*m.X2[k:k + 4])))
    d1 = {(i, j): d for i, j, d in m.spec1}
    d2 = {(i, j): d for i, j, d in m.spec2}
    contacts = [(i, j, r, r, 0.0, 1.0) for i, j, r in m.common]
    for key in set(d1) | set(d2):
        if key in d1 and key in d2:      # formed in both, different lengths
            contacts.append((*key, d1[key], d2[key], 1.0, alpha))
        elif key in d1:
            contacts.append((*key, d1[key], d1[key], 0.0, alpha))
        else:
            contacts.append((*key, d2[key], d2[key], 0.0, alpha * lam))
    for i, j, *_ in contacts:
        excl.add((i, j))
    for c in range(n_copies):
        o = c * n
        for i, j, r0, k in bonds:
            bond.addBond(o + i, o + j, r0, k)
        for i, j, l, th in angles:
            ang.addAngle(o + i, o + j, o + l, th, 40.0 * kT)
        for i, j, l, q, t1, t2 in tors:
            dih.addTorsion(o + i, o + j, o + l, o + q, [t1, t2])
        for i, j, r1, r2, two, s_ in contacts:
            con.addBond(o + i, o + j, [r1, r2, two, s_])
        for i, j in excl:
            rep.addExclusion(o + i, o + j)
    forces = [bond, ang, dih, con, rep]
    if ligand is not None:               # explicit virtual ligand (check)
        P, R0, eps_l, wl = ligand
        lig = mm.CustomBondForce("-el*exp(-(r-r0)^2/(2*wl^2))")
        lig.addGlobalParameter("el", eps_l / max(1, len(P)))
        lig.addGlobalParameter("wl", wl); lig.addPerBondParameter("r0")
        for c in range(n_copies):
            for (i, j), r0 in zip(P, R0):
                lig.addBond(c * n + int(i), c * n + int(j), [float(r0)])
        forces.append(lig)
    for f in forces:
        sys_.addForce(f)
    m.contacts = contacts
    return sys_


def copies_positions(m, n_copies, which=None, spacing=60.0):
    X = []
    for c in range(n_copies):
        base = m.X1 if (which or ("1" if c % 2 == 0 else "2")) == "1" \
            else m.X2
        X.append(base - base.mean(0) + np.array([c * spacing, 0, 0]))
    return np.concatenate(X)


def order_params(m, X):
    """Q_common, Q_inactive-specific, Q_active-specific for coordinates X
    (..., n, 3): fraction of contacts within 1.2x their native length."""
    def q(pairs):
        i = np.array([p[0] for p in pairs]); j = np.array([p[1] for p in pairs])
        r0 = np.array([p[2] for p in pairs])
        r = np.linalg.norm(X[..., i, :] - X[..., j, :], axis=-1)
        return (r < 1.2 * r0).mean(-1)
    return q(m.common), q(m.spec1), q(m.spec2)


def state_coordinate(m, X):
    """dRMS to the inactive minus dRMS to the active structure, over the
    pairs whose C-alpha distance differs by > 0.3 nm between the states
    (< 0: inactive-like, > 0: active-like). Shape (..., )."""
    if not hasattr(m, "_disc"):
        i, j = np.triu_indices(len(m.ids), 4)
        d1 = np.linalg.norm(m.X1[i] - m.X1[j], axis=1)
        d2 = np.linalg.norm(m.X2[i] - m.X2[j], axis=1)
        k = (np.abs(d1 - d2) > 0.3) & (np.minimum(d1, d2) < 1.5)
        m._disc = (i[k], j[k], d1[k], d2[k])
    i, j, d1, d2 = m._disc
    r = np.linalg.norm(X[..., i, :] - X[..., j, :], axis=-1)
    return (np.sqrt(((r - d1) ** 2).mean(-1))
            - np.sqrt(((r - d2) ** 2).mean(-1)))


def _run_copies(args):
    """Worker: simulate n_copies in one context; save C-alpha frames.
    args = dict(name, eps, lam, w, n_copies, n_frames, stride, seed, out,
    start)."""
    a = dict(lam=1.0, w=0.05, start=None, eps=1.0, alpha=1.0, ligand=None)
    a.update(args)
    name, eps, lam, n_copies = a["name"], a["eps"], a["lam"], a["n_copies"]
    n_frames, stride, seed, out, start = (a["n_frames"], a["stride"],
                                          a["seed"], a["out"], a["start"])
    import openmm as mm
    m = DualBasin(name)
    s = build_gauss(m, n_copies=n_copies, eps=eps, lam=lam, w=a["w"],
                    alpha=a["alpha"], ligand=a["ligand"])
    integ = mm.LangevinMiddleIntegrator(300, 1.0, 0.02)
    integ.setRandomNumberSeed(seed)
    ctx = mm.Context(s, integ, mm.Platform.getPlatformByName("CPU"),
                     {"Threads": "1"})
    ctx.setPositions(copies_positions(m, n_copies, which=start))
    mm.LocalEnergyMinimizer.minimize(ctx)
    ctx.setVelocitiesToTemperature(300, seed)
    n = len(m.ids)
    traj = np.empty((n_frames, n_copies, n, 3), dtype=np.float32)
    for f in range(n_frames):
        integ.step(stride)
        X = ctx.getState(getPositions=True).getPositions(asNumpy=True)._value
        traj[f] = X.reshape(n_copies, n, 3)
    np.save(out, traj)
    return out


PARAMS = {"GK": dict(eps=1.2, alpha=0.5, lam=1.5, w=0.05),
          "SHP2": dict(eps=1.2, alpha=0.3, lam=0.6, w=0.05)}


def produce(name, n_proc=4, n_copies=8, n_frames=200, stride=2500,
            seed0=100):
    """Production: n_proc processes x n_copies copies, started alternately
    from the inactive and active structures. Frames every `stride` steps
    (dt = 20 fs CG time). Cached as .npy in figures/_retro_cache/."""
    from multiprocessing import Pool
    os.makedirs(CACHE, exist_ok=True)
    jobs = [dict(PARAMS[name], name=name, n_copies=n_copies,
                 n_frames=n_frames, stride=stride, seed=seed0 + k,
                 out=os.path.join(CACHE, f"{name}_prod{k}.npy"))
            for k in range(n_proc)]
    todo = [j for j in jobs if not os.path.exists(j["out"])]
    if todo:
        with Pool(min(len(todo), os.cpu_count() or 1)) as p:
            p.map(_run_copies, todo)
    return [j["out"] for j in jobs]


# ---------------------------------------------------------------------------
# analysis: MSM, virtual-ligand site scan, evaluation
# ---------------------------------------------------------------------------
FRAME_NS = 2500 * 0.02 / 1000           # 50 ps of CG time per frame


def load_copies(paths, burn=0.1):
    """List of per-copy trajectories (n_frames, n, 3), burn-in removed."""
    out = []
    for p in paths:
        T = np.load(p)
        b = int(len(T) * burn)
        out += [T[b:, c] for c in range(T.shape[1])]
    return out


def msm_features(m, trajs, n_pairs=300, seed=0):
    """Distances of a fixed random subset of the state-discriminating
    pairs (internal coordinates, no superposition)."""
    state_coordinate(m, m.X1[None])            # builds m._disc
    i, j, _, _ = m._disc
    rng = np.random.default_rng(seed)
    k = rng.choice(len(i), size=min(n_pairs, len(i)), replace=False)
    return [np.linalg.norm(t[:, i[k]] - t[:, j[k]], axis=-1).astype(float)
            for t in trajs]


def build_msm(m, trajs, lags=(1, 2, 3, 5, 8, 12, 20), n_clusters=100):
    feats = msm_features(m, trajs)
    r = m23.analyze(feats, dt_ps=FRAME_NS, n_clusters=n_clusters, n_macro=2,
                    its_samples=0, run_tpt=False, verbose=False,
                    lags=list(lags))
    z = np.concatenate([state_coordinate(m, t) for t in trajs])
    a = np.concatenate(r["active"])
    zmac = [z[np.isin(a, np.where(r["assign"] == k)[0])].mean()
            for k in range(2)]
    inact, act = int(np.argmin(zmac)), int(np.argmax(zmac))
    r.update(m23.tpt_analysis(r, inact, act, verbose=False))
    r["inactive"], r["active_macro"] = inact, act
    r["z"] = z
    return r


def site_contacts(m, target):
    """For every residue: the contacts (i, j, r_target) among residues
    within 0.8 nm of it in the target-state structure."""
    X = m.X1 if target == "inactive" else m.X2
    D = np.linalg.norm(X[:, None] - X[None], axis=2)
    pairs = [(i, j, D[i, j]) for i, j, _ in
             m.common + (m.spec1 if target == "inactive" else m.spec2)]
    P = np.array([(i, j) for i, j, _ in pairs])
    R0 = np.array([d for *_, d in pairs])
    sites = []
    for r_ in range(len(m.ids)):
        near = D[r_] < 0.8
        k = np.where(near[P[:, 0]] & near[P[:, 1]])[0]
        sites.append((P[k], R0[k]))
    return sites


def ligand_energy(trajs, site, eps_l=2.0 * 2.494, width=0.1):
    """Virtual ligand: V = -eps_l * mean_k exp(-(r_k - r0_k)^2 / 2 w^2)
    over the site's target-state contacts (kJ/mol, per frame)."""
    P, R0 = site
    if len(P) == 0:
        return [np.zeros(len(t)) for t in trajs]
    out = []
    for t in trajs:
        r = np.linalg.norm(t[:, P[:, 0]] - t[:, P[:, 1]], axis=-1)
        out.append(-eps_l * np.exp(-(r - R0) ** 2 / (2 * width ** 2))
                   .mean(-1))
    return out


def scan_sites(m, trajs, r, target, beta=1 / 2.494):
    """Predicted d ln K(active/inactive) and d ln k(inactive -> active)
    for a virtual ligand at every residue (reweighting one MSM)."""
    sites = site_contacts(m, target)
    msm = r["msm"]; n = msm.n_states
    pi0 = msm.stationary_distribution
    act = r["assign"] == r["active_macro"]
    lnK0 = np.log(pi0[act].sum() / pi0[~act].sum())
    k0 = r["flux"].rate
    dK, dk = np.zeros(len(sites)), np.zeros(len(sites))
    for s_idx, site in enumerate(sites):
        V = ligand_energy(trajs, site)
        d = m24.state_shifts(V, r["active"], n, beta)
        Vb = m24.crossing_shifts(V, r["active"], r["lag"], n, beta)
        kp, mp = m24.predicted_rate(msm, r["A"], r["B"], d, beta, Vb)
        pp = mp.stationary_distribution
        dK[s_idx] = np.log(pp[act].sum() / pp[~act].sum()) - lnK0
        dk[s_idx] = np.log(kp / k0)
    return dK, dk, sites


def rmsf(trajs):
    """Per-residue C-alpha RMSF after Kabsch superposition on the mean."""
    X = np.concatenate(trajs)[::5]
    ref = X.mean(0)
    for _ in range(2):
        al = []
        for x in X:
            P, Q = x - x.mean(0), ref - ref.mean(0)
            U, _, Vt = np.linalg.svd(P.T @ Q)
            dd = np.sign(np.linalg.det(U @ Vt))
            al.append(P @ (U @ np.diag([1, 1, dd]) @ Vt))
        X = np.array(al); ref = X.mean(0)
    return np.sqrt(((X - ref) ** 2).sum(-1).mean(0))


def structural_baseline(m, target):
    """Static two-structure comparison: mean |d_inactive - d_active| over
    pairs within 0.8 nm of each residue (local structural change)."""
    X = m.X1 if target == "inactive" else m.X2
    D = np.linalg.norm(X[:, None] - X[None], axis=2)
    D1 = np.linalg.norm(m.X1[:, None] - m.X1[None], axis=2)
    D2 = np.linalg.norm(m.X2[:, None] - m.X2[None], axis=2)
    return np.array([np.abs(D1 - D2)[np.ix_(D[r] < 0.8, D[r] < 0.8)].mean()
                     for r in range(len(m.ids))])


def evaluate(score, m, target, n_null=2000, seed=0):
    """ROC-AUC (site residues vs rest), top-10% enrichment, best rank, and
    the percentile of the true site's mean score among random spatial
    patches of the same size (patch = k nearest residues to a random
    centre in the target-state structure)."""
    from scipy.stats import rankdata
    pos = np.zeros(len(score), bool); pos[m.site_idx] = True
    ranks = rankdata(score)
    auc = (ranks[pos].sum() - pos.sum() * (pos.sum() + 1) / 2) / \
        (pos.sum() * (~pos).sum())
    top = np.argsort(-score)[:max(1, len(score) // 10)]
    enrich = pos[top].mean() / pos.mean()
    best = int(np.min(len(score) - ranks[pos] + 1))
    X = m.X1 if target == "inactive" else m.X2
    D = np.linalg.norm(X[:, None] - X[None], axis=2)
    rng = np.random.default_rng(seed); k = pos.sum()
    null = np.array([score[np.argsort(D[c])[:k]].mean()
                     for c in rng.integers(0, len(score), n_null)])
    pct = float((null < score[pos].mean()).mean())
    return dict(auc=float(auc), enrich_top10=float(enrich), best_rank=best,
                patch_percentile=pct)


def ligand_checks(m, r, trajs, sites, tags, n_copies=8, n_frames=100,
                  stride=2500, eps_l=2.0 * 2.494, beta=1 / 2.494):
    """Simulate WITH the virtual ligand applied explicitly (2 processes per
    site, all sites concurrently) and compare the measured d ln K with the
    prediction from reweighting the unbiased frames. Populations are split
    on the order parameter z (independent of the MSM)."""
    from multiprocessing import Pool
    jobs = []
    for (P, R0), tag in zip(sites, tags):
        jobs += [dict(PARAMS[m.name], name=m.name, n_copies=n_copies,
                      n_frames=n_frames, stride=stride, seed=500 + k,
                      ligand=(P, R0, eps_l, 0.1),
                      out=os.path.join(CACHE, f"{m.name}_lig_{tag}_{k}.npy"))
                 for k in range(2)]
    todo = [j for j in jobs if not os.path.exists(j["out"])]
    if todo:
        with Pool(min(len(todo), os.cpu_count() or 1)) as p:
            p.map(_run_copies, todo)
    z0 = r["z"]
    act = r["assign"] == r["active_macro"]
    a = np.concatenate(r["active"])
    thr = 0.5 * (np.median(z0[np.isin(a, np.where(act)[0])])
                 + np.median(z0[np.isin(a, np.where(~act)[0])]))

    def lnK(z, w=None):
        w = np.ones_like(z) if w is None else w
        return np.log((w * (z > thr)).sum() / (w * (z <= thr)).sum())

    out = []
    rng = np.random.default_rng(0)
    for site, tag in zip(sites, tags):
        bt = load_copies([j["out"] for j in jobs if f"_lig_{tag}_" in
                          j["out"]])
        zc = [state_coordinate(m, t) for t in bt]
        zb = np.concatenate(zc)
        bs = [lnK(np.concatenate([zc[i] for i in rng.integers(0, len(zc),
                                                             len(zc))]))
              for _ in range(300)]
        w = np.exp(-beta * np.concatenate(ligand_energy(trajs, site, eps_l)))
        out.append(dict(tag=tag, dlnK_sim=float(lnK(zb) - lnK(z0)),
                        ci=[float(x - lnK(z0))
                            for x in np.percentile(bs, [2.5, 97.5])],
                        dlnK_pred=float(lnK(z0, w) - lnK(z0))))
    return out


def run_system(name, checks=True):
    m = DualBasin(name)
    print(m.summary())
    trajs = load_copies(produce(name))
    r = build_msm(m, trajs)
    z = r["z"]
    flips = sum(int(np.abs(np.diff(np.sign(state_coordinate(m, t))))
                    .sum() // 2) for t in trajs)
    target = "inactive" if m.kind == "inhibitor" else "active"
    sign = -1.0 if m.kind == "inhibitor" else 1.0
    dK, dk, sites = scan_sites(m, trajs, r, target)
    scores = {"MSM: d ln K": sign * dK, "MSM: d ln k (rate)": sign * dk,
              "static structure difference": structural_baseline(m, target),
              "RMSF (conventional MD)": rmsf(trajs)}
    ev = {k: evaluate(v, m, target) for k, v in scores.items()}
    act = r["assign"] == r["active_macro"]
    pi = r["msm"].stationary_distribution
    res = dict(name=name, kind=m.kind, n_beads=len(m.ids),
               n_site=len(m.site_idx), n_frames=int(sum(len(t)
                                                        for t in trajs)),
               transitions=flips, lag_ns=r["lag"] * FRAME_NS,
               timescales_ns=r["timescales_ps"][:3].tolist(),
               ck_err=r["ck_err"], p_active=float(pi[act].sum()),
               mfpt_activation_ns=float(1 / r["rate_ps"]),
               eval=ev, ids=m.ids, site=[m.ids[i] for i in m.site_idx],
               scores={k: v.tolist() for k, v in scores.items()})
    top = np.argsort(-scores["MSM: d ln K"])[:10]
    res["top10"] = [int(m.ids[i]) for i in top]
    print(f"  MSM: lag {res['lag_ns']:.2f} ns, CK err {r['ck_err']:.3f}, "
          f"p(active) {res['p_active']:.2f}, {flips} transitions")
    for k, v in ev.items():
        print(f"  {k:30s} AUC {v['auc']:.2f}  top10% enrichment "
              f"{v['enrich_top10']:.1f}x  best rank {v['best_rank']:3d}  "
              f"site vs random patches {v['patch_percentile'] * 100:.0f}th pct")
    print(f"  top-10 MSM residues: {res['top10']}; drug site: {res['site']}")
    if checks:
        i_top = int(top[0])
        X = m.X1 if target == "inactive" else m.X2
        far = np.linalg.norm(X - X[i_top], axis=1) > 2.0
        cand = np.where(far)[0]
        i_ctl = int(cand[np.argmin(np.abs(scores["MSM: d ln K"][cand]))])
        res["checks"] = ligand_checks(m, r, trajs,
                                      [sites[i_top], sites[i_ctl]],
                                      [f"top{m.ids[i_top]}",
                                       f"ctl{m.ids[i_ctl]}"])
        for c in res["checks"]:
            print(f"  explicit-ligand check {c['tag']}: predicted d ln K "
                  f"{c['dlnK_pred']:+.2f}, simulated {c['dlnK_sim']:+.2f} "
                  f"[{c['ci'][0]:+.2f}, {c['ci'][1]:+.2f}]")
    return res




# ---------------------------------------------------------------------------
# figure + report
# ---------------------------------------------------------------------------
def figure(results, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cols = {"MSM: d ln K": "#2a6fdb", "MSM: d ln k (rate)": "#9b59b6",
            "static structure difference": "#e4572e",
            "RMSF (conventional MD)": "#7f8c8d"}
    nsys = len(results)
    fig, ax = plt.subplots(nsys, 3, figsize=(16, 4.6 * nsys), squeeze=False)
    for row, res in enumerate(results):
        ids = np.array(res["ids"]); site = set(res["site"])
        pos = np.array([i in site for i in ids])
        a = ax[row, 0]
        for k in ("MSM: d ln K", "static structure difference"):
            v = np.array(res["scores"][k])
            v = (v - v.mean()) / (v.std() + 1e-12)
            a.plot(ids, v, lw=0.9, c=cols[k], label=k)
        a.scatter(ids[pos], np.full(pos.sum(), a.get_ylim()[1] * 0.92),
                  marker="v", c="k", s=18, label="drug-contact residue")
        a.set_xlabel("residue"); a.set_ylabel("score (z-scored)")
        a.legend(fontsize=7, loc="lower right")
        a.set_title(f"({'abc'[0]}) {res['name']} ({res['kind']}): site "
                    f"scores along the sequence")
        a = ax[row, 1]
        for k, v in res["scores"].items():
            v = np.array(v); o = np.argsort(-v)
            tpr = np.cumsum(pos[o]) / pos.sum()
            fpr = np.cumsum(~pos[o]) / (~pos).sum()
            a.plot(fpr, tpr, c=cols[k],
                   label=f"{k}  AUC {res['eval'][k]['auc']:.2f}")
        a.plot([0, 1], [0, 1], "k:", lw=0.8)
        a.set_xlabel("false-positive rate"); a.set_ylabel("true-positive rate")
        a.legend(fontsize=7, loc="lower right")
        a.set_title(f"(b) ROC against the {res['n_site']} drug-contact "
                    f"residues")
        a = ax[row, 2]
        if res.get("checks"):
            for i, c in enumerate(res["checks"]):
                a.errorbar([i - 0.1], [c["dlnK_sim"]],
                           yerr=[[c["dlnK_sim"] - c["ci"][0]],
                                 [c["ci"][1] - c["dlnK_sim"]]],
                           fmt="o", c="#29a36a", capsize=5,
                           label="simulated with ligand" if i == 0 else None)
                a.plot([i + 0.1], [c["dlnK_pred"]], "s", c="#2a6fdb",
                       label="predicted (reweighting)" if i == 0 else None)
            a.axhline(0, c="k", lw=0.8)
            a.set_xticks(range(len(res["checks"])),
                         [c["tag"] for c in res["checks"]])
            a.set_ylabel("d ln K(active/inactive)")
            a.legend(fontsize=8)
            a.set_title("(c) explicit virtual-ligand simulations\n"
                        "top-ranked site vs distant control site")
        else:
            a.axis("off")
    fig.suptitle("M25 — retrospective test on proteins with known "
                 "allosteric drugs (dual-basin C-alpha model)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96]); fig.savefig(path, dpi=120)
    plt.close(fig)


def report(systems=("GK", "SHP2")):
    import json
    print("=" * 70)
    print("MODULE 25 — retrospective allosteric-site test")
    print("=" * 70)
    results = [run_system(s) for s in systems]
    figure(results, os.path.join(FIG, "m25_allosteric_retro.png"))
    with open(os.path.join(FIG, "m25_summary.json"), "w") as fh:
        json.dump(results, fh, indent=1, default=lambda o: (
            o.tolist() if hasattr(o, "tolist") else float(o)))
    return results


if __name__ == "__main__":
    import sys
    report(tuple(sys.argv[1:]) or ("GK", "SHP2"))
