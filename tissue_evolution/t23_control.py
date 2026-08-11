"""
T23 - flux control coefficients: which of these enzymes actually decides anything.

Every rate comparison so far has treated the enzymes of a pathway as
interchangeable units of "a protein in a pathway". They are not. In metabolic
control analysis the influence of enzyme i on pathway flux J is

    C_i^J  =  d ln J / d ln a_i

i.e. change the enzyme's activity by 1%, see what fraction of a percent the
flux moves. The summation theorem says these sum to 1 across the whole system,
so control is a budget: if one enzyme has 0.6, everything else shares 0.4.
An enzyme with C_i^J ~ 0 is not a control point no matter how tissue-specific,
how highly expressed, or how fast-evolving it is.

This matters for the evolutionary question because the usual story ("this
enzyme is tissue-enriched, therefore it is what the tissue tuned") assumes
expression implies leverage. C^J is the test of that assumption.

MODEL. Eight TCA reactions plus Complex I and a lumped ubiquinol oxidase
(Complex III+IV), with convenience kinetics (Liebermeister & Klipp 2006) so
that every rate law is thermodynamically consistent: v = 0 exactly when the
mass-action ratio reaches Keq, and no reaction can be driven backwards through
its own equilibrium. Conserved moieties (NAD/NADH, Q/QH2, CoA, and the carbon
in the cycle intermediates) are imposed as constraints rather than left to the
solver, which is what makes the steady state unique.

The equilibrium constants are the interesting part and they are not free
parameters - they are thermodynamics:

    CS      ~2e5    irreversible, pulls the cycle
    ACO2    0.067   sits against citrate
    IDH3    ~1      driven by NAD and CO2 removal
    OGDH    ~1e5    irreversible
    SUCL    3.8     near equilibrium
    SDH     ~1      NEAR EQUILIBRIUM - this is why complex II can reverse
    FH      4.4
    MDH2    2.8e-5  strongly uphill; only runs because CS drains OAA

Two of those lines already answer questions from the candidate list. SDH at
Keq ~ 1 is the structural reason complex II is the one that switches direction
when the quinone pool reduces. MDH2 running uphill is why the cycle is a cycle
and not a chain.

⚠️ HONESTY. Km and Vmax are literature-typical mitochondrial values, not fitted
to any one tissue. A single control coefficient computed from one parameter set
is not a measurement. So nothing here is reported as a point estimate: every
coefficient is recomputed over parameter sets sampled from log-normal
distributions around those values, and what gets reported is the distribution
and the rank. A conclusion that survives 2-fold parameter scatter is about the
network topology and thermodynamics, which are known; a conclusion that does
not survive it is about my guesses, and is labelled as such.

The summation theorem is used as a unit test, not as a result: if the numerics
are right the coefficients must sum to 1.0, and that catches sign errors,
bad finite differences and unconverged steady states in one number.
"""
import os

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

from .common import DERIVED, FIGDIR, apply_figure_style, PALETTE

# ---------------------------------------------------------------- model

SPECIES_ORDER = ["Cit", "Iso", "aKG", "SucCoA", "Suc", "Fum", "Mal", "OAA",
                 "NADH", "QH2"]

# reaction: (name, gene, {substrate: Km}, {product: Km}, Keq, Vmax)
# concentrations in mM, rates in mM/s.
#
# Two things here are physics, not fitting, and both were wrong in my first
# pass in ways that produced a nonsense steady state:
#
#   * Km for Q and QH2 must be comparable to the quinone POOL (3 mM), not to
#     a free-solution ligand concentration. Setting them to 0.3 uM saturates
#     every membrane complex permanently, and a saturated pool cannot poise -
#     the redox state of the quinone pool then carries no information and
#     complex II can never reverse.
#   * Complex I's equilibrium constant is not the solution-chemistry value.
#     CI pumps 4 H+ against the protonmotive force, and at dp ~ 180 mV that
#     back-pressure (~69 kJ/mol) very nearly cancels the NADH -> Q drop. The
#     EFFECTIVE Keq is therefore near unity, which is precisely why complex I
#     sits close to equilibrium in vivo and can run backwards (reverse
#     electron transport). Using the uncoupled Keq ~ 1e3 makes CI an
#     irreversible sink that pins NADH/NAD at 0.005 instead of ~0.25.
REACTIONS = [
    ("CS",   "CS",     {"OAA": 0.005, "AcCoA": 0.016}, {"Cit": 1.6},            2.0e5, 1.6),
    ("ACO2", "ACO2",   {"Cit": 0.48},                  {"Iso": 0.12},           0.067, 4.0),
    ("IDH3", "IDH3A",  {"Iso": 0.13, "NAD": 0.08},     {"aKG": 0.25, "NADH": 0.05}, 1.0, 1.3),
    ("OGDH", "OGDH",   {"aKG": 0.40, "NAD": 0.07, "CoA": 0.007},
                       {"SucCoA": 0.05, "NADH": 0.05},                          1.0e5, 1.1),
    ("SUCL", "SUCLA2", {"SucCoA": 0.04},               {"Suc": 0.60},           3.8,   2.0),
    ("SDH",  "SDHA",   {"Suc": 0.13, "Q": 1.0},        {"Fum": 0.50, "QH2": 1.0},   1.0, 1.4),
    ("FH",   "FH",     {"Fum": 0.05},                  {"Mal": 0.25},           4.4,   5.0),
    ("MDH2", "MDH2",   {"Mal": 0.50, "NAD": 0.06},     {"OAA": 0.04, "NADH": 0.02}, 2.8e-5, 8.0),
    ("CI",   "NDUFS1", {"NADH": 0.02, "Q": 1.0},       {"NAD": 0.5, "QH2": 1.0},    3.0, 3.0),
    ("OX",   "MT-CO1", {"QH2": 1.0, "O2": 0.001},      {"Q": 1.0},              1.0e6, 4.0),
]

# Absolute enzyme concentrations in vivo are not known to me, so the three
# free capacity scales (complex I, the oxidase, and a common multiplier on the
# eight TCA enzymes, whose RATIOS come from literature Vmax/Km) are calibrated
# to reproduce three measured steady-state observables rather than guessed:
#     NADH/(NAD+NADH) = 0.25,  QH2/(Q+QH2) = 0.40,  TCA turnover = 0.05 mM/s.
# calibrate() below re-derives these; they are frozen here so the module is
# deterministic and does not re-fit on every import.
CALIBRATION = dict(vmax_CI=0.79267, vmax_OX=0.38001, tca_scale=2.499)
TARGETS = dict(nadh_frac=0.25, qh2_frac=0.40, J=0.05)

POOLS = dict(NAD_tot=3.0, Q_tot=3.0, CoA_tot=1.3, AcCoA=0.05, O2=0.05)

# Physiological mitochondrial concentrations (mM). OAA is ~1e-4, four orders
# below malate, and that is not a typo - it is what Keq(MDH2) = 2.8e-5 forces.
# Starting the solver anywhere near "0.01 mM OAA" puts MDH2 five orders of
# magnitude out of balance and no root finder recovers from that.
Y0 = np.array([0.35, 0.02, 0.15, 0.02, 0.50, 0.08, 0.60, 7e-5, 0.60, 1.00])


def _conc(y, pools):
    """State vector + conservation relations -> full concentration dict."""
    c = dict(zip(SPECIES_ORDER[:-2], y[:8]))
    c["NADH"], c["QH2"] = y[8], y[9]
    c["NAD"] = pools["NAD_tot"] - c["NADH"]
    c["Q"] = pools["Q_tot"] - c["QH2"]
    c["AcCoA"] = pools["AcCoA"]
    c["CoA"] = pools["CoA_tot"] - c["AcCoA"] - c["SucCoA"]
    c["O2"] = pools["O2"]
    return c


def rate(rxn, c, vmax_scale=1.0):
    """Convenience kinetics, with the Haldane relation applied correctly.

    The subtlety that bit me: writing the numerator as

        prod(s_i/Ks_i) - prod(p_j/Kp_j) / Keq

    does NOT make the reaction stop at Keq. It stops where

        prod(p_j) / prod(s_i)  =  Keq * prod(Kp_j) / prod(Ks_i)

    so the stated equilibrium constant is silently multiplied by the ratio of
    the Michaelis constants. For MDH2 that factor is (0.5*0.06)/(0.04*0.02) =
    37.5, i.e. the model was running malate dehydrogenase at an equilibrium
    constant 37-fold away from the thermodynamic one, and complex II was
    sitting at an effective Keq of 3.85 rather than 1.0 - which is precisely
    the number that decides whether it can reverse.

    The Haldane relation fixes it: divide the stated Keq by the Km ratio
    before using it, so that v = 0 lands exactly on the thermodynamic
    equilibrium in concentration units.
    """
    _, _, subs, prods, keq, vmax = rxn
    num_f = 1.0
    num_r = 1.0
    den_f = 1.0
    den_r = 1.0
    km_s = 1.0
    km_p = 1.0
    for s, km in subs.items():
        x = max(c[s], 1e-12) / km
        num_f *= x
        den_f *= (1.0 + x)
        km_s *= km
    for p, km in prods.items():
        x = max(c[p], 1e-12) / km
        num_r *= x
        den_r *= (1.0 + x)
        km_p *= km
    keq_eff = keq * km_s / km_p
    return vmax * vmax_scale * (num_f - num_r / keq_eff) / (den_f + den_r - 1.0)


# stoichiometry: species -> {reaction: coefficient}
STOICH = {
    "Cit":    {"CS": +1, "ACO2": -1},
    "Iso":    {"ACO2": +1, "IDH3": -1},
    "aKG":    {"IDH3": +1, "OGDH": -1},
    "SucCoA": {"OGDH": +1, "SUCL": -1},
    "Suc":    {"SUCL": +1, "SDH": -1},
    "Fum":    {"SDH": +1, "FH": -1},
    "Mal":    {"FH": +1, "MDH2": -1},
    "OAA":    {"MDH2": +1, "CS": -1},
    "NADH":   {"IDH3": +1, "OGDH": +1, "MDH2": +1, "CI": -1},
    "QH2":    {"SDH": +1, "CI": +1, "OX": -1},
}


def derivatives(y, pools, scales, reactions=REACTIONS):
    c = _conc(y, pools)
    v = {r[0]: rate(r, c, scales.get(r[0], 1.0)) for r in reactions}
    d = np.zeros(len(SPECIES_ORDER))
    for i, sp in enumerate(SPECIES_ORDER):
        d[i] = sum(k * v[r] for r, k in STOICH[sp].items())
    return d, v


def steady_state(pools=None, scales=None, y0=None, reactions=None,
                 t_end=2.0e4):
    """Integrate to steady state, then polish the root.

    Solving dy/dt = 0 algebraically fails here: the Jacobian spans five orders
    of magnitude in concentration (OAA ~ 1e-4 mM against citrate ~ 0.35) and
    the rate laws are stiff near equilibrium, so a root finder started even a
    little off wanders into negative concentrations. Integrating the ODE
    instead is unconditionally stable and has a second advantage - the sum of
    the eight cycle intermediates is EXACTLY conserved by the stoichiometry, so
    forward integration carries the carbon constraint for free rather than
    having it imposed as an extra equation. The total carbon of the run is
    therefore whatever y0 sums to.

    least_squares afterwards is only a polish; the tolerance check is what
    decides whether the point is accepted.
    """
    pools = dict(POOLS, **(pools or {}))
    scales = scales or {}
    reactions = calibrated() if reactions is None else reactions
    y_init = Y0.copy() if y0 is None else np.asarray(y0, float)

    def rhs(_, y):
        d, _v = derivatives(np.maximum(y, 1e-14), pools, scales, reactions)
        return d

    sol = solve_ivp(rhs, (0.0, t_end), y_init, method="BDF",
                    rtol=1e-10, atol=1e-14)
    y = np.maximum(sol.y[:, -1], 1e-14)

    def resid(z):
        d, _ = derivatives(np.exp(z), pools, scales, reactions)
        out = list(d)
        out[7] = float(np.sum(np.exp(z[:8])) - np.sum(y_init[:8]))
        return np.array(out)

    try:
        pol = least_squares(resid, np.log(y), method="lm",
                            xtol=1e-15, ftol=1e-15, max_nfev=5000)
        if np.max(np.abs(resid(pol.x))) < np.max(np.abs(resid(np.log(y)))):
            y = np.exp(pol.x)
    except Exception:                             # noqa: BLE001
        pass

    d, v = derivatives(y, pools, scales, reactions)
    # scale-free convergence test: every net rate must be negligible against
    # the fluxes actually running through the system
    scale = max(abs(v["OX"]), 1e-12)
    ok = bool(np.max(np.abs(d)) < 1e-6 * scale) and bool(np.all(y > 0))
    return y, v, ok


def control_coefficients(pools=None, reactions=None, eps=0.01, flux="OX",
                         only=None, warm=None):
    """C_i^J = dlnJ/dln(Vmax_i) by central difference in log space."""
    reactions = calibrated() if reactions is None else reactions
    y0, v0, ok = steady_state(pools, reactions=reactions, y0=warm)
    if not ok or v0[flux] <= 0:
        return None, None, None
    out = {}
    for r in reactions:
        name = r[0]
        if only is not None and name != only:
            continue
        _, vp, okp = steady_state(pools, {name: 1 + eps}, y0=y0, reactions=reactions)
        _, vm, okm = steady_state(pools, {name: 1 - eps}, y0=y0, reactions=reactions)
        if not (okp and okm) or vp[flux] <= 0 or vm[flux] <= 0:
            out[name] = np.nan
            continue
        out[name] = ((np.log(vp[flux]) - np.log(vm[flux]))
                     / (np.log(1 + eps) - np.log(1 - eps)))
    return out, v0, y0


def calibrated(cal=None, reactions=REACTIONS):
    """Apply the capacity calibration to a reaction list."""
    c = dict(CALIBRATION, **(cal or {}))
    out = []
    for name, gene, subs, prods, keq, vmax in reactions:
        if name == "CI":
            vmax = c["vmax_CI"]
        elif name == "OX":
            vmax = c["vmax_OX"]
        else:
            vmax = vmax * c["tca_scale"]
        out.append((name, gene, subs, prods, keq, vmax))
    return out


def calibrate(x0=(0.25, 0.35, 50.0)):
    """Re-derive CALIBRATION from the three steady-state targets."""
    def resid(z):
        vci, vox, tca = np.exp(z)
        rxn = calibrated(dict(vmax_CI=vci, vmax_OX=vox, tca_scale=tca))
        y, v, ok = steady_state(reactions=rxn)
        if not ok:
            return np.array([10.0, 10.0, 10.0])
        return np.array([y[8] / POOLS["NAD_tot"] - TARGETS["nadh_frac"],
                         y[9] / POOLS["Q_tot"] - TARGETS["qh2_frac"],
                         np.log(max(v["CS"], 1e-9) / TARGETS["J"])])

    sol = least_squares(resid, np.log(x0), diff_step=0.05,
                        xtol=1e-13, ftol=1e-13)
    vci, vox, tca = np.exp(sol.x)
    return dict(vmax_CI=vci, vmax_OX=vox, tca_scale=tca), \
        float(np.max(np.abs(resid(sol.x))))


def sample_parameters(rng, spread=2.0):
    """Log-normal jitter on every Km and Vmax; Keq is thermodynamics, untouched."""
    s = np.log(spread) / 1.96                     # 95% within [1/spread, spread]
    out = []
    for name, gene, subs, prods, keq, vmax in calibrated():
        j = lambda x: float(x * np.exp(rng.normal(0, s)))   # noqa: E731
        out.append((name, gene,
                    {k: j(x) for k, x in subs.items()},
                    {k: j(x) for k, x in prods.items()},
                    keq, j(vmax)))
    return out


def robustness(n=300, seed=0, pools=None):
    """Distribution of every C^J over sampled parameter sets."""
    rng = np.random.default_rng(seed)
    rows, sums = [], []
    for _ in range(n):
        rxn = sample_parameters(rng)
        cc, v0, _ = control_coefficients(pools, reactions=rxn)
        if cc is None or any(np.isnan(list(cc.values()))):
            continue
        rows.append(cc)
        sums.append(sum(cc.values()))
    return pd.DataFrame(rows), np.array(sums)


def reserve_scan(fracs=None, thresh=0.1, flux="OX"):
    """How far must each enzyme be knocked down before it controls anything?

    A single C_i^J at one operating point is not a property of the enzyme, it
    is a property of the enzyme AT THAT CAPACITY. Reporting "SUCL has zero
    control" from the calibrated model would be close to circular: the
    calibration put the eight TCA enzymes in ~57-fold excess in order to
    reproduce a reduced NAD pool at physiological flux, and anything in
    57-fold excess has no control by definition.

    The non-circular question is the one an evolving lineage actually faces:
    how much of this enzyme can you LOSE before flux notices? So sweep each
    enzyme's own capacity down and record two thresholds - where C_i^J crosses
    `thresh`, and where flux itself has fallen 10%. An enzyme that tolerates a
    20-fold knockdown before either happens has enormous reserve, and
    tissue-to-tissue expression differences smaller than that cannot be about
    controlling this flux.
    """
    fracs = fracs or [1.0, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05, 0.03, 0.02, 0.01,
                      0.005, 0.002]
    y_base, v0, _ok = steady_state()
    j0 = v0[flux]
    rows = []
    for name, gene, *_ in REACTIONS:
        c_cross = j_cross = np.nan
        warm = y_base
        for f in fracs:
            rxn = [(n, g, sb, pr, k, vm * (f if n == name else 1.0))
                   for n, g, sb, pr, k, vm in calibrated()]
            cc, v, y = control_coefficients(reactions=rxn, flux=flux,
                                            only=name, warm=warm)
            if cc is None:
                continue
            warm = y                      # next frac starts from this state
            ci, jr = cc[name], v[flux] / j0
            rows.append(dict(enzyme=name, gene=gene, frac=f, C=ci, J_rel=jr))
            if np.isnan(c_cross) and ci >= thresh:
                c_cross = f
            if np.isnan(j_cross) and jr <= 0.9:
                j_cross = f
        rows.append(dict(enzyme=name, gene=gene, frac=np.nan,
                         C=c_cross, J_rel=j_cross))
    df = pd.DataFrame(rows)
    curves = df[df["frac"].notna()]
    thr = df[df["frac"].isna()].rename(
        columns={"C": "frac_at_C_thresh", "J_rel": "frac_at_J_minus10pct"})
    return curves, thr.drop(columns=["frac"]).set_index("enzyme")


def allocation_scan(scales=(0.3, 0.6, 1, 2.499, 6, 20, 60)):
    """Sweep TCA capacity against respiratory capacity.

    The calibrated model sets the eight TCA enzymes at ~2.5x the literature
    scale, and a reserve measured at that point partly reflects that choice.
    But the 2.5 is itself an inference, not a free knob: it is what the three measured
    observables (NAD poise, quinone poise, turnover) force given literature
    Km values. This sweep reports the whole curve anyway, so the reader can
    see at which allocation each enzyme would acquire control - and therefore
    how far from that allocation a real tissue would have to sit before a
    cycle enzyme became a decision point.
    """
    rows = []
    for sc in scales:
        cal = dict(CALIBRATION, tca_scale=sc)
        rxn = calibrated(cal)
        cc, v, _y = control_coefficients(reactions=rxn)
        if cc is None:
            continue
        rows.append(dict(tca_scale=sc, J=v["OX"], **cc))
    return pd.DataFrame(rows).set_index("tca_scale")


def hypoxia_scan(levels=(1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01)):
    """Drop O2 and watch complex II approach, then cross, zero."""
    rows = []
    for f in levels:
        pools = dict(POOLS, O2=POOLS["O2"] * f)
        y, v, ok = steady_state(pools)
        if not ok:
            continue
        c = _conc(y, pools)
        cc, _, _ = control_coefficients(pools)
        # Complex II runs backwards when [Fum][QH2]/([Suc][Q]) exceeds Keq,
        # i.e. when Fum/Suc rises above Keq * Q/QH2. Reporting the observed
        # ratio against that moving threshold says exactly how far from
        # reversal the system is, and what it would take to get there.
        keq_sdh = dict((r[0], r[4]) for r in REACTIONS)["SDH"]
        need = keq_sdh * c["Q"] / max(c["QH2"], 1e-12)
        have = c["Fum"] / max(c["Suc"], 1e-12)
        rows.append(dict(o2_frac=f, J=v["OX"], v_SDH=v["SDH"], v_CI=v["CI"],
                         QH2_frac=c["QH2"] / (c["QH2"] + c["Q"]),
                         NADH_frac=c["NADH"] / (c["NADH"] + c["NAD"]),
                         Suc=c["Suc"], Fum=c["Fum"],
                         fum_suc=have, fum_suc_needed=need,
                         fold_short=need / max(have, 1e-12),
                         C_SDH=(cc or {}).get("SDH", np.nan),
                         C_CS=(cc or {}).get("CS", np.nan),
                         C_OX=(cc or {}).get("OX", np.nan)))
    return pd.DataFrame(rows)


def make_figure(dist, hyp, curves=None):
    plt = apply_figure_style()
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.4))

    order = dist.median().sort_values(ascending=False).index
    ax = axes[0]
    if curves is not None:
        for i, (name, g) in enumerate(curves.groupby("enzyme")):
            g = g.sort_values("frac")
            ax.plot(g["frac"], g["C"], "o-", ms=3, lw=1.2,
                    color=PALETTE[i % len(PALETTE)], label=name)
        ax.set_xscale("log")
        ax.axhline(0.1, color="k", ls=":", lw=0.9)
        ax.set_xlabel("enzyme capacity (fraction of calibrated)")
        ax.set_ylabel(r"$C_i^{J}$")
        ax.set_title("reserve: how far each enzyme falls\nbefore it controls anything")
        ax.legend(fontsize=6, ncol=2)
    else:
        ax.boxplot([dist[c].dropna() for c in order], vert=False,
                   tick_labels=list(order), showfliers=False,
                   medianprops=dict(color=PALETTE[1]))
        ax.set_xlabel(r"$C_i^{J}$")

    ax = axes[1]
    ax.plot(hyp["o2_frac"], hyp["v_SDH"], "o-", color=PALETTE[0], label="complex II flux")
    ax.axhline(0, color="k", lw=1.0)
    ax.set_xscale("log")
    ax.invert_xaxis()
    ax.set_xlabel(r"O$_2$ (fraction of normoxic)")
    ax.set_ylabel("flux (mM/s)")
    ax2 = ax.twinx()
    ax2.plot(hyp["o2_frac"], hyp["QH2_frac"], "s--", color=PALETTE[2], ms=3.5,
             label="quinone pool reduced")
    ax2.set_ylabel("QH$_2$ / (Q + QH$_2$)")
    ax.set_title("complex II direction vs the quinone pool")
    ax.legend(fontsize=7, loc="lower left")
    ax2.legend(fontsize=7, loc="upper right")

    ax = axes[2]
    for i, col in enumerate(["C_CS", "C_SDH", "C_OX"]):
        ax.plot(hyp["o2_frac"], hyp[col], "o-", ms=3.5,
                color=PALETTE[i], label=col.replace("C_", r"$C^J$ "))
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xscale("log")
    ax.invert_xaxis()
    ax.set_xlabel(r"O$_2$ (fraction of normoxic)")
    ax.set_ylabel(r"$C_i^{J}$")
    ax.set_title("control migrates as oxygen falls")
    ax.legend(fontsize=7)

    p = os.path.join(FIGDIR, "t23_control.png")
    fig.savefig(p)
    plt.close(fig)
    return p


def report():
    print("=" * 72)
    print("T23 - flux control coefficients: which enzymes are decisive")
    print("=" * 72)

    cc, v0, _y0 = control_coefficients()
    print("\nbaseline steady state, calibrated to measured poise:")
    y, v, _ = steady_state()
    print(f"  NADH/(NAD+NADH) = {y[8] / POOLS['NAD_tot']:.3f}   "
          f"QH2/(Q+QH2) = {y[9] / POOLS['Q_tot']:.3f}   "
          f"TCA turnover = {v['CS']:.4f} mM/s   OAA = {y[7]:.2e} mM")
    tot = sum(cc.values())
    print(f"\n  SUMMATION THEOREM: sum C_i^J = {tot:.6f}  (must be 1.0)  -> "
          f"{'PASS' if abs(tot - 1) < 5e-3 else 'FAIL'}")
    print("  (this is a unit test of the numerics, not a result)")

    cap = {r[0]: r[5] for r in calibrated()}
    print(f"\n{'enzyme':<8}{'gene':<10}{'C^J':>9}{'flux':>10}{'Vmax/flux':>11}")
    for name, gene, *_ in REACTIONS:
        print(f"{name:<8}{gene:<10}{cc[name]:>9.4f}{v0[name]:>10.4f}"
              f"{cap[name] / max(abs(v0[name]), 1e-9):>11.1f}")
    print("  numerical resolution of the finite difference is ~1e-3;")
    print("  anything printed as 0.0000 means |C| below that, not exactly zero.")

    print("\nreserve: how far must each enzyme fall before it matters?")
    curves, thr = reserve_scan()
    print(f"{'enzyme':<8}{'gene':<10}{'C>0.1 at':>11}{'J-10% at':>11}")
    for k, r in thr.iterrows():
        a = ("never" if np.isnan(r["frac_at_C_thresh"])
             else f"{1 / r['frac_at_C_thresh']:.0f}x down")
        b = ("never" if np.isnan(r["frac_at_J_minus10pct"])
             else f"{1 / r['frac_at_J_minus10pct']:.0f}x down")
        print(f"{k:<8}{r['gene']:<10}{a:>11}{b:>11}")
    print("  'never' = survived a 500-fold knockdown without crossing.")

    print("\ncontrol vs capacity allocation (TCA capacity / respiratory capacity)")
    alloc = allocation_scan()
    cols = ["CS", "IDH3", "OGDH", "SUCL", "SDH", "MDH2", "CI", "OX"]
    print("  tca_scale  " + "".join(f"{c:>8}" for c in cols) + f"{'J':>9}")
    for k, r in alloc.iterrows():
        print(f"  {k:>9.0f}  " + "".join(f"{r[c]:>8.3f}" for c in cols)
              + f"{r['J']:>9.4f}")

    print("\nhypoxia: does complex II reverse?")
    hyp = hypoxia_scan()
    print(f"{'O2':>7}{'J':>10}{'v_SDH':>10}{'QH2frac':>9}{'NADHfrac':>10}"
          f"{'Fum/Suc':>10}{'needed':>9}{'short':>8}")
    for _, r in hyp.iterrows():
        print(f"{r['o2_frac']:>7.2f}{r['J']:>10.4f}{r['v_SDH']:>10.5f}"
              f"{r['QH2_frac']:>9.3f}{r['NADH_frac']:>10.3f}"
              f"{r['fum_suc']:>10.4f}{r['fum_suc_needed']:>9.3f}"
              f"{r['fold_short']:>7.0f}x")
    print("  complex II reverses only when Fum/Suc exceeds 'needed'.")

    print("\nparameter robustness (150 log-normal draws, 2-fold 95% spread)")
    dist, sums = robustness(n=150)
    print(f"  converged {len(dist)}/150;  summation theorem median "
          f"{np.median(sums):.5f}, max deviation {np.max(np.abs(sums - 1)):.1e}")
    q = dist.quantile([0.05, 0.5, 0.95]).T
    q.columns = ["p5", "median", "p95"]
    q["P(top)"] = (dist.idxmax(axis=1).value_counts()
                   / len(dist)).reindex(q.index).fillna(0.0)
    q = q.sort_values("median", ascending=False)
    print(f"\n{'enzyme':<8}{'p5':>9}{'median':>9}{'p95':>9}{'P(is top)':>11}")
    for k, r in q.iterrows():
        print(f"{k:<8}{r['p5']:>9.3f}{r['median']:>9.3f}{r['p95']:>9.3f}"
              f"{r['P(top)']:>11.2f}")

    dist.to_csv(os.path.join(DERIVED, "t23_control_samples.tsv"),
                sep="\t", index=False)
    q.to_csv(os.path.join(DERIVED, "t23_control_summary.tsv"), sep="\t")
    thr.to_csv(os.path.join(DERIVED, "t23_reserve.tsv"), sep="\t")
    curves.to_csv(os.path.join(DERIVED, "t23_reserve_curves.tsv"),
                  sep="\t", index=False)
    alloc.to_csv(os.path.join(DERIVED, "t23_allocation.tsv"), sep="\t")
    hyp.to_csv(os.path.join(DERIVED, "t23_hypoxia.tsv"), sep="\t", index=False)
    p = make_figure(dist, hyp, curves)
    print(f"\nfigure written to: {p}")
    return cc, dist, hyp, thr, alloc


if __name__ == "__main__":
    report()
