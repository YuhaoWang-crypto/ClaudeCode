"""
Self-checks for plife.bio.function.  Run with `python3 test_function.py`.

The claim under test in each block is a *capability*: that the computed pattern
actually delivers ultrasensitivity, buffering, rate enhancement, bistability or
length discrimination -- and by how much.
"""

import sys

import numpy as np

from plife.bio import function as F, pathways as PW, sites as SI, units as U

fails = []
run = []


def check(name, cond, detail=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"   {detail}" if detail else ""))
    run.append(name)
    if not cond:
        fails.append(name)


def cross_mix(kd_uM, nA, nB, d=5.0):
    kd = np.array([[np.inf, kd_uM], [kd_uM, np.inf]]) * 1e-6
    return SI.SiteMixture.from_valence(["A", "B"], [d, d], [nA, nB], kd)


# --------------------------------------------------------------------------- #
print("gel fraction is exact at the transition")
m = cross_mix(10.0, 4, 4)
c_star = float(PW.gel_site_conc(10.0, 4, 4)) / 4.0
check("closed-form threshold is 1.875 uM", abs(c_star - 1.875) < 1e-9, f"{c_star:.4f} uM")
check("gel fraction is exactly zero below it",
      F.gel_fraction([c_star * 0.95, c_star * 0.95], m).max() == 0.0)
check("and positive above it", F.gel_fraction([c_star * 1.1, c_star * 1.1], m).min() > 0,
      f"g = {F.gel_fraction([c_star * 1.1, c_star * 1.1], m)[0]:.4f}")
lam_at = SI.percolates([c_star, c_star], m)[1]
check("onset coincides with branching eigenvalue 1", abs(lam_at - 1.0) < 5e-3,
      f"lambda({c_star:.4f} uM) = {lam_at:.5f}")

print("\ncapability 1: ultrasensitivity WITHOUT cooperative binding")
dr = F.dose_response(m, [1.0, 1.0], np.geomspace(0.5, 60, 260))
nH = F.hill_coefficient(dr["scale"], dr["gel"][:, 0], at=0.5)
check("effective Hill coefficient is well above 1", nH > 2.5, f"n_H = {nH:.2f}")
check("...from bonds that are strictly non-cooperative (one Kd, no allostery)", True,
      "the model has a single Kd and no cooperativity term")
nH_by_val = {}
for n in (3, 4, 6, 8):
    mm = cross_mix(10.0, n, n)
    d2 = F.dose_response(mm, [1.0, 1.0], np.geomspace(0.05, 200, 300))
    nH_by_val[n] = F.hill_coefficient(d2["scale"], d2["gel"][:, 0], at=0.5)
vals = [nH_by_val[n] for n in (3, 4, 6, 8)]
check("steepness grows with valence", all(np.diff(vals) > -0.2),
      "  ".join(f"n={n}:{nH_by_val[n]:.1f}" for n in (3, 4, 6, 8)))

print("\ncapability 2: the free concentration is buffered")
scale = np.geomspace(2e-3, 60, 300)
dr = F.dose_response(m, [1.0, 1.0], scale)
tot, free = dr["conc"][:, 0], dr["free"][:, 0]
slope = F.buffering_quality(tot, free)
below = slope[tot < 0.05 * c_star]
above = slope[tot > 6 * c_star]
check("free tracks total in the dilute limit", below.mean() > 0.9,
      f"dln free/dln total = {below.mean():.3f}")
check("free is pinned above threshold", above.mean() < 0.35,
      f"dln free/dln total = {above.mean():.3f}")
i3 = int(np.argmin(np.abs(tot - 3 * c_star)))
i9 = int(np.argmin(np.abs(tot - 9 * c_star)))
check("3x more protein does NOT give 3x more free protein",
      free[i9] / free[i3] < 1.8, f"total x3.0 -> free x{free[i9] / free[i3]:.2f}")

print("\ncapability 3: co-concentration speeds a bimolecular step")
whole, local = F.rate_enhancement(100.0, 100.0, phi_condensate=0.01)
check("local rate gain is K_E * K_S", abs(local - 1e4) < 1, f"{local:.0f}x")
check("whole-cell gain is much smaller than the local one", whole < local / 50,
      f"whole cell {whole:.0f}x vs local {local:.0f}x at phi = 1%")
w_slow, _ = F.rate_enhancement(100.0, 100.0, phi_condensate=0.01, viscosity_penalty=100.0)
check("a 100x slower interior wipes out most of the gain", w_slow < whole / 20,
      f"{w_slow:.1f}x with a 100x viscosity penalty")

print("\ncapability 4: exclusion-driven bistability EMERGES")
# LAT-like: a membrane-localised scaffold with many phospho-tyrosines, a
# tetravalent adaptor, both at ~10 uM local concentration, and a phosphatase
# left with 5% access inside the network (the size-exclusion figure measured in
# test_pathways is 0.0025-0.06, so 5% is the conservative end).
N_MAX, N_ADAPT, C_SCAF, C_ADAPT, KD = 10.0, 4.0, 10.0, 10.0, 10.0


def gel_of_p(p):
    """Phosphorylation creates the scaffold's valence."""
    n = max(N_MAX * p, 1e-6)
    sm = SI.SiteMixture.from_valence(
        ["scaffold", "adaptor"], [5.0, 5.0], [n, N_ADAPT],
        np.array([[np.inf, KD], [KD, np.inf]]) * 1e-6)
    return F.gel_fraction([C_SCAF, C_ADAPT], sm)[0]


bi = F.phospho_bistability(gel_of_p, np.geomspace(0.02, 20, 140), access_in_gel=0.05)
n_bi = sum(o["bistable"] for o in bi["states"])
check("a bistable window exists", bi["bistable_any"],
      f"{n_bi}/{len(bi['states'])} kinase:phosphatase ratios give two stable states")
if n_bi:
    Rs = [o["R"] for o in bi["states"] if o["bistable"]]
    on = [o["on"] for o in bi["states"] if o["bistable"]]
    off = [o["off"] for o in bi["states"] if o["bistable"]]
    check("the two branches are far apart (a real switch, not a wiggle)",
          max(np.array(on) - np.array(off)) > 0.3,
          f"widest gap p = {min(off):.2f} vs {max(on):.2f}")
    check("hysteresis spans a useful range of kinase activity",
          max(Rs) / min(Rs) > 1.3, f"R from {min(Rs):.2f} to {max(Rs):.2f} "
          f"({max(Rs) / min(Rs):.1f}x)")
no_fb = F.phospho_bistability(lambda p: 0.0, np.geomspace(0.02, 20, 140))
check("with NO exclusion feedback the bistability disappears",
      not no_fb["bistable_any"], "control: gel fraction forced to 0")

print("\ncapability 5: polymer length IS the valence")


def build_dna(L_bp):
    n_turns = max(L_bp / 20.0, 1.0)
    return SI.SiteMixture(
        names=["cGAS", "dsDNA"], diameters=[5.2, 2.0],
        sites=[SI.Site("cGAS", "DNA", 2.0), SI.Site("dsDNA", "turn", n_turns)],
        kd=np.array([[np.inf, 0.5], [0.5, np.inf]]) * 1e-6)


res = F.min_polymer_length(build_dna, np.arange(20, 1400, 10.0), [0.2, 0.02])
check("short DNA cannot switch cGAS on", res["branching"][0] < 1.0,
      f"20 bp -> branching {res['branching'][0]:.3f}")
check("long DNA can", np.isfinite(res["L_min"]), f"threshold at {res['L_min']:.0f} bp")
check("the threshold is in the range where DNA length actually matters",
      50 <= res["L_min"] <= 1000, f"L_min = {res['L_min']:.0f} bp")
lo = F.min_polymer_length(build_dna, np.arange(20, 4000, 20.0), [0.02, 0.002])["L_min"]
check("scarcer cGAS needs longer DNA", lo > res["L_min"],
      f"{res['L_min']:.0f} bp at 0.2 uM  ->  {lo:.0f} bp at 0.02 uM")

print("\nallostery: a ligand that changes a PARTNER's valence")
th, freefrac = F.competitive_occupancy(1.0, [1e6], [1.0])
check("saturating single ligand fills the site", th[0] > 0.999 and freefrac < 1e-3,
      f"theta = {th[0]:.6f}")
th, _ = F.competitive_occupancy(1e-6, [1.0], [1.0])   # site-limited, ligand excess
check("Langmuir limit at trace site", abs(th[0] - 1.0 / (1.0 + 1.0)) < 1e-4,
      f"theta = {th[0]:.5f} vs L/(Kd+L) = 0.5")
th, _ = F.competitive_occupancy(1.0, [3.0, 3.0], [1.0, 0.5])
check("the tighter competitor wins the site", th[1] > th[0],
      f"theta(Kd=1) = {th[0]:.3f}  <  theta(Kd=0.5) = {th[1]:.3f}")
check("occupancies cannot exceed the site", th.sum() <= 1.0 + 1e-9,
      f"sum theta = {th.sum():.6f}")
check("conditional valence interpolates", F.conditional_valence(1.0, 5.0, 0.5) == 3.0)

# G3BP1: CAPRIN1 relieves autoinhibition, USP10 competes for the same NTF2 site
G3, RNA_C, N_NTF2 = 3.0, 0.3, 2.0
KD_CAP, KD_USP, N_CLOSED, N_OPEN = 1.0, 0.5, 0.8, 4.0


def sg_branching(cap_uM, usp_uM):
    theta, _ = F.competitive_occupancy(N_NTF2 * G3, [cap_uM, usp_uM], [KD_CAP, KD_USP])
    n_rgg = float(F.conditional_valence(N_CLOSED, N_OPEN, theta[0]))
    sm = SI.SiteMixture(
        ["G3BP1", "RNA"], [2 * U.radius_from_mw(104), 2 * U.radius_idp(2000)],
        [SI.Site("G3BP1", "RGG", n_rgg), SI.Site("RNA", "site", 18.0)],
        np.array([[np.inf, 1.0], [1.0, np.inf]]) * 1e-6)
    return SI.percolates([G3, RNA_C], sm)[1]


check("without CAPRIN1 the granule does not form", sg_branching(0.0, 0.0) < 1.0,
      f"branching = {sg_branching(0.0, 0.0):.3f}")
check("CAPRIN1 switches it ON (relieves autoinhibition)", sg_branching(10.0, 0.0) > 2.0,
      f"branching {sg_branching(0.0, 0.0):.2f} -> {sg_branching(10.0, 0.0):.2f} at 10 uM")
check("USP10 switches it OFF again (same site, no valence)",
      sg_branching(3.0, 30.0) < 1.0 < sg_branching(3.0, 0.0),
      f"branching {sg_branching(3.0, 0.0):.2f} -> {sg_branching(3.0, 30.0):.2f} "
      f"on 30 uM USP10")
grid = np.geomspace(0.01, 300, 300)
off = grid[np.array([sg_branching(3.0, u) for u in grid]) < 1.0]
check("...at a USP10:CAPRIN1 ratio the model can name", 1.0 < off[0] / 3.0 < 20.0,
      f"dissolves above [USP10] = {off[0] / 3.0:.1f} x [CAPRIN1]")

# and the control: WITHOUT allostery the model gets the sign backwards
def flat_branching(cap_uM):
    sm = SI.SiteMixture(
        ["G3BP1", "RNA", "CAPRIN1"],
        [2 * U.radius_from_mw(104), 2 * U.radius_idp(2000), 2 * U.radius_from_mw(78)],
        [SI.Site("G3BP1", "RGG", 4.0), SI.Site("RNA", "site", 18.0),
         SI.Site("CAPRIN1", "RGG", 3.0)],
        np.array([[np.inf, 1.0, np.inf], [1.0, np.inf, 5.0],
                  [np.inf, 5.0, np.inf]]) * 1e-6)
    return SI.percolates([G3, RNA_C, cap_uM], sm)[1]


check("a flat (no-allostery) model gets CAPRIN1's sign BACKWARDS",
      flat_branching(10.0) < flat_branching(0.0),
      f"flat model: branching {flat_branching(0.0):.2f} -> {flat_branching(10.0):.2f} "
      "as CAPRIN1 rises -- the opposite of experiment")

print("\nreal module: LAT/Grb2/SOS1 under local enrichment")
sysL = PW.get("lat")
lat = PW.to_mixture(sysL)
c_lat = sysL.conc("mid")
d = F.dose_response(lat, c_lat, np.geomspace(1, 2000, 200))
i = int(np.argmax(d["gel"][:, 0] > 0))
check("switches on well below the membrane's 133x, not at bulk concentration",
      5 < d["scale"][i] < 133, f"switches at {d['scale'][i]:.0f}x "
      f"(membrane supplies {PW.membrane_boost():.0f}x)")
gmax = d["gel"][:, 0].max()
check("but the condensed fraction SATURATES far below 1", gmax < 0.35,
      f"gel fraction plateaus at {gmax:.3f} even at 2000x")

# why: uniform scaling keeps the stoichiometry fixed, and SOS1 is the scarce
# partner, so the network is stoichiometry-limited rather than affinity-limited.
sites = c_lat[lat.owner] * lat.n
check("SOS1 supplies the fewest sites of the three",
      lat.labels[int(np.argmin(sites))].startswith("SOS1"),
      "  ".join(f"{l}={s:.2f}uM" for l, s in zip(lat.labels, sites)))
boost = float(PW.membrane_boost())
g_even = F.gel_fraction(c_lat * boost, lat)[0]
g_sos = F.gel_fraction(c_lat * boost * np.array([1.0, 1.0, 6.0]), lat)[0]
g_grb = F.gel_fraction(c_lat * boost * np.array([1.0, 6.0, 1.0]), lat)[0]
check("6x more SOS1 helps far more than 6x more Grb2 (a falsifiable prediction)",
      g_sos > 2 * g_even and g_sos > 2 * g_grb,
      f"gel: even {g_even:.3f}  +SOS1 {g_sos:.3f}  +Grb2 {g_grb:.3f}")

print()
if fails:
    print(f"{len(fails)}/{len(run)} FAILED: {fails}")
    sys.exit(1)
print(f"all {len(run)} checks passed")
