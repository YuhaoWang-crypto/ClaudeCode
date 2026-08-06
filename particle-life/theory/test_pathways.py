"""
Self-checks for plife.bio.sites and plife.bio.pathways.  Run with
`python3 test_pathways.py`.

Every closed-form condition stated in the pathways module docstring is checked
against the full Wertheim solver, which knows nothing about those formulas.
"""

import sys

import numpy as np

from plife.bio import pathways as PW, phases as P, sites as SI, thermo as T, units as U

fails = []
run = []


def check(name, cond, detail=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"   {detail}" if detail else ""))
    run.append(name)
    if not cond:
        fails.append(name)


def cross_mix(kd_uM, nA, nB, d=5.0):
    """Two species, one site class each, binding only across."""
    kd = np.array([[np.inf, kd_uM], [kd_uM, np.inf]]) * 1e-6
    return SI.SiteMixture.from_valence(["A", "B"], [d, d], [nA, nB], kd)


def cmax_uM(d_nm, phi=0.55):
    return float(U.density_to_conc(phi / ((4 / 3) * np.pi * (d_nm / 2) ** 3)))


# --------------------------------------------------------------------------- #
print("site-resolved model reduces to the single-class model")
kd = np.array([[3e-4, 3e-5], [3e-5, np.inf]])
mix = T.Mixture(names=["A", "B"], diameters=[5.2, 4.3], valence=[5.0, 3.0], kd=kd)
sm = SI.SiteMixture.from_valence(["A", "B"], [5.2, 4.3], [5.0, 3.0], kd)
rho = U.conc_to_density(np.array([[3.0, 7.0], [100.0, 250.0], [1e-3, 1e-3]]))
check("unbonded fractions identical",
      np.abs(T.unbonded_fractions(rho, mix) - SI.unbonded_fractions(rho, sm)).max() == 0)
check("free energy identical",
      np.abs(T.free_energy(rho, mix) - T.free_energy(rho, sm)).max() == 0)
check("stability eigenvalues identical",
      np.abs(T.spinodal_eig(rho, mix) - T.spinodal_eig(rho, sm)).max() == 0)

print("\ncondition 1+2: Flory-Stockmayer gel point vs the Wertheim solver")
worst = 0.0
for nA, nB in [(3, 3), (4, 4), (6, 3), (4, 2), (5, 5), (8, 4)]:
    for kdv in (1.0, 30.0, 350.0):
        m = cross_mix(kdv, nA, nB)
        S = np.geomspace(1e-2, 1e5, 800)
        lam = np.array([SI.percolates([s / nA, s / nB], m)[1] for s in S])
        got = S[lam >= 1.0][0] if (lam >= 1.0).any() else np.inf
        want = float(PW.gel_site_conc(kdv, nA, nB))
        worst = max(worst, abs(got / want - 1.0))
check("closed-form S* matches the solver", worst < 0.08, f"max rel err {worst:.1%}")

print("\ncondition 1: the valence gate is a gate, not a tendency")
S = np.geomspace(1e-3, 1e5, 400)
lam22 = max(SI.percolates([s / 2, s / 2], cross_mix(1e-4, 2, 2))[1] for s in S)
check("2x2 approaches but never reaches percolation, even at Kd = 0.1 nM",
      lam22 < 1.0, f"max branching eigenvalue {lam22:.6f}")
check("...and the formula says so too", not np.isfinite(PW.percolation_p(2, 2)))
check("3x3 does percolate", SI.percolates([40.0, 40.0], cross_mix(30.0, 3, 3))[0])

print("\ncondition 3: the stoichiometry window")
worst_w = 0.0
for nA, nB in [(3, 3), (4, 4), (6, 3), (5, 3)]:
    G = PW.stoichiometry_window(nA, nB)
    m = cross_mix(10.0, nA, nB)
    cA, cB = 300.0, np.geomspace(1e-2, 3e5, 800)
    lam = np.array([SI.percolates([cA, c], m)[1] for c in cB])
    r = nA * cA / (nB * cB[lam >= 1.0])
    worst_w = max(worst_w, abs(r.max() / G - 1), abs((1 / r.min()) / G - 1))
check("window edges reach (nA-1)(nB-1)", worst_w < 0.10, f"max rel err {worst_w:.1%}")

m44 = cross_mix(10.0, 4, 4)
lam_mid = SI.percolates([60.0, 60.0], m44)[1]
lam_hi = SI.percolates([60.0, 60.0 * 40], m44)[1]
check("excess partner dissolves the network (re-entrance)", lam_mid >= 1 > lam_hi,
      f"branching {lam_mid:.2f} -> {lam_hi:.2f} on 40x partner")

print("\ncondition 4: miscibility, chi vs the common-tangent construction")


def dense_phases(kd_aa, kd_bb, kd_ab, n=4.0, d=5.0, cq=100.0):
    mx = T.Mixture(names=["A", "B"], diameters=[d, d], valence=[n, n],
                   kd=np.array([[kd_aa, kd_ab], [kd_ab, kd_bb]]) * 1e-6)
    c = np.geomspace(0.3, cmax_uM(d), 40)
    rho2, shape = P.grid2d(c, c)
    res = T.coexistence(mx, rho2)
    i = int(np.argmin(np.abs(c - cq)))
    sp = P._facet_split(res, i * shape[1] + i, rho2)
    if sp is None:
        return None
    v = sp[0][sp[1] > 1e-3]
    return U.density_to_conc(v[np.argsort(v.sum(axis=1))])


ph = dense_phases(30.0, 30.0, 3.0)          # cross 10x tighter than self
check("cross-binding beats self -> ONE shared body",
      len(ph) == 2 and abs(ph[-1][0] - ph[-1][1]) / max(ph[-1]) < 0.05,
      f"{len(ph)} phases, dense = ({ph[-1][0]:.0f}, {ph[-1][1]:.0f}) uM")
check("...and chi agrees", PW.miscible(30.0, 30.0, 3.0))

ph = dense_phases(30.0, 30.0, 300.0)        # cross 10x weaker than self
check("self-binding beats cross -> TWO pure bodies",
      len(ph) == 3 and min(ph[-1][1], ph[-2][0]) / max(ph[-1][0], ph[-2][1]) < 0.05,
      f"{len(ph)} phases, dense = " + " and ".join(f"({p[0]:.0f},{p[1]:.0f})" for p in ph[1:]))
check("...and chi agrees", not PW.miscible(30.0, 30.0, 300.0))

print("\nclient recruitment and steric exclusion")
sc = dense_phases(30.0, np.inf, 10.0)
check("a monovalent client is still recruited", sc is not None and sc[-1][1] / sc[0][1] > 10,
      f"K_client = {sc[-1][1] / sc[0][1]:.0f}")


def exclusion(d_client):
    mx = T.Mixture(names=["A", "B"], diameters=[5.0, d_client], valence=[4.0, 1.0],
                   kd=np.array([[30.0, np.inf], [np.inf, np.inf]]) * 1e-6)
    cA = np.geomspace(0.3, cmax_uM(5.0), 40)
    cB = np.geomspace(0.3, cmax_uM(d_client), 40)
    rho2, shape = P.grid2d(cA, cB)
    res = T.coexistence(mx, rho2)
    i = int(np.argmin(np.abs(cA - 100.0))), int(np.argmin(np.abs(cB - 100.0)))
    sp = P._facet_split(res, i[0] * shape[1] + i[1], rho2)
    v = sp[0][sp[1] > 1e-3]
    v = v[np.argsort(v.sum(axis=1))]
    return v[-1][1] / v[0][1]


k_small, k_big = exclusion(3.0), exclusion(8.0)
check("big inert clients are excluded more than small ones", k_big < k_small,
      f"K(3nm) = {k_small:.3f}  vs  K(8nm) = {k_big:.4f}")

# --------------------------------------------------------------------------- #
print("\nsite classes matter: one SH2 is not two")
d = [2 * U.radius_idp(210), 2 * U.radius_from_mw(25), 2 * U.radius_from_mw(150)]
c0 = np.array([1.58, 2.24, 0.158])
collapsed = SI.SiteMixture.from_valence(
    ["pLAT", "Grb2", "SOS1"], d, [4.0, 3.0, 4.0],
    np.array([[np.inf, 1.0, np.inf], [1.0, np.inf, 10.0],
              [np.inf, 10.0, np.inf]]) * 1e-6)
resolved = PW.to_mixture(PW.get("lat"))


def threshold(m, lo=1.0, hi=1e6):
    if SI.percolates(c0 * hi, m)[1] < 1.0:
        return np.inf
    while hi / lo > 1.001:
        mid = np.sqrt(lo * hi)
        lo, hi = (lo, mid) if SI.percolates(c0 * mid, m)[0] else (mid, hi)
    return np.sqrt(lo * hi)


t_c, t_r = threshold(collapsed), threshold(resolved)
check("collapsing Grb2 to 3 equivalent sites lowers the threshold", t_c < t_r / 2,
      f"collapsed switches on at {t_c:.1f}x, site-resolved needs {t_r:.1f}x")
between = np.sqrt(t_c * t_r)
check("...so the two models disagree over a whole band of concentrations",
      SI.percolates(c0 * between, collapsed)[0]
      and not SI.percolates(c0 * between, resolved)[0],
      f"at {between:.1f}x: collapsed says network, site-resolved says soup")

print("\nthe atlas")
for s in PW.ATLAS:
    a = PW.assess(s)
    if not (a["conc_uM"].min() > 0 and np.isfinite(a["branching"])):
        check(f"{s.key} well-formed", False)
check("every entry assessable", not any("well-formed" in f for f in fails))

pka = PW.assess(PW.get("pka"))
check("PKA-AKAP: 1 nM affinity but valence 1 -> no network",
      not pka["gate_open"] and not np.isfinite(PW.boost_to_percolate(PW.get("pka"))),
      f"Kd = {pka['kd_uM'] * 1e3:.0f} nM, branching = {pka['branching']:.3f}")
check("MAPK core is a chain, not a network",
      not np.isfinite(PW.boost_to_percolate(PW.get("mapk"))))
wnt = PW.assess(PW.get("wnt"))
check("Wnt DIX head-to-tail: no 3-D network at any concentration",
      not np.isfinite(PW.boost_to_percolate(PW.get("wnt"))),
      f"branching = {wnt['branching']:.3f}")
check("MAVS filament likewise 1-D",
      not np.isfinite(PW.boost_to_percolate(PW.get("mavs"))))
cgas = PW.assess(PW.get("cgas"))
check("cGAS + long dsDNA: gate wide open", cgas["percolates"] and cgas["G"] > 20,
      f"G = {cgas['G']:.0f}, p* = {cgas['p_star']:.3f}, branching = {cgas['branching']:.2f}")

print("\nlocalisation is the switch")
boost = float(PW.membrane_boost())
check("2-D recruitment concentrates by ~10^2", 50 < boost < 500, f"{boost:.0f}x")
for key in ("lat", "nephrin"):
    a = PW.assess(PW.get(key))
    need = PW.boost_to_percolate(PW.get(key))
    am = PW.assess(PW.get(key), local_boost=boost)
    check(f"{key}: soup in solution, network once localised",
          not a["percolates"] and need < boost and am["percolates"],
          f"needs {need:.0f}x, membrane gives {boost:.0f}x "
          f"(branching {a['branching']:.2f} -> {am['branching']:.2f})")

print("\nthe pool")
pool = PW.pool_analysis()
live = set(pool["live"])
check("most modules are OFF at bulk cytosolic concentration", len(live) <= 6,
      f"{len(live)}/{len(PW.ATLAS)} percolate: {sorted(live)}")
check("orthogonal modules sort into one body each",
      pool["n_bodies"] == len(live), f"{pool['n_bodies']} bodies from {len(live)} modules")
merged = PW.pool_analysis(cross={("psd", "spop"): 0.05})
check("a promiscuous cross-interaction merges two bodies",
      merged["n_bodies"] == pool["n_bodies"] - 1,
      f"{pool['n_bodies']} -> {merged['n_bodies']} bodies")

print()
if fails:
    print(f"{len(fails)}/{len(run)} FAILED: {fails}")
    sys.exit(1)
print(f"all {len(run)} checks passed")
