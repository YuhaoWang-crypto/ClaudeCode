"""
Diagnostics runner — what could a fit to each endpoint actually recover?

Answers, per core, how many parameter directions the assay's own readout
constrains, which parameters it cannot see at all, and whether a
time-resolved readout would help. Writes ``figures/genotox_identifiability.png``.

    python3 -m genotox.run_identifiability
"""
from __future__ import annotations

import os
from dataclasses import replace

import numpy as np

from .assay import VirtualAssay
from .core import Exposure, SOSCore
from .cytogenetic import CytogeneticCore
from .damage import TabulatedSource, demo
from .doseresponse import (COMET, GADD45A_GFP, MN_CBMN, UMU, dose_series,
                           log_doses)
from .identifiability import (check_symmetry, print_summary,
                              relative_sensitivity, summarise)
from .p53 import P53Core
from .readouts import BetaGalONPG

FIGDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      "figures")

#: measurement noise used for the practical-rank count.  5% relative is
#: optimistic for a plate assay and deliberately so: it is an upper bound on
#: what could be learned, not a typical figure.
NOISE = 0.05
FACTOR = 1.65      # "estimable" = pinned to better than this factor


def _nest(core, overrides):
    """Apply overrides, routing ``p53.<name>`` into the shared sub-model."""
    own = {k: v for k, v in overrides.items() if "." not in k}
    sub = {k.split(".", 1)[1]: v for k, v in overrides.items()
           if k.startswith("p53.")}
    core = replace(core, **own) if own else core
    if sub:
        core = replace(core, p53=replace(core.p53, **sub))
    return core


def _values(core, names):
    out = {}
    for n in names:
        out[n] = (getattr(core.p53, n.split(".", 1)[1]) if n.startswith("p53.")
                  else getattr(core, n))
    return out


# --------------------------------------------------------------------------
# Endpoint 1 — umu
# --------------------------------------------------------------------------
SOS_PARAMS = ("v_repair", "k_repair", "k_repair_lin", "K_recA", "h_recA",
              "beta_lexA", "K_lexA", "k_cleave", "kumu_ratio", "h_umu",
              "leak", "k_txn", "d_mrna", "k_tsl", "K_tox", "h_tox", "mu_max")
SOS_DOSES = log_doses(0.02, 4.0, 9)
COMP = demo("direct-acting bulky (4NQO-like)")


def sos_endpoint_observe(core):
    """What the umu kit reports: induction ratio and growth factor per dose."""
    ser = dose_series(VirtualAssay(source=TabulatedSource(), core=core),
                      COMP, SOS_DOSES, protocol=UMU)
    rows = ser["rows"][1:]
    return np.concatenate([[r["IR"] for r in rows],
                           [r["growth_factor"] for r in rows]])


def sos_timecourse_observe(core, times=(30.0, 60.0, 90.0, 120.0)):
    """A readout that samples the cascade instead of only its endpoint."""
    src = TabulatedSource()
    out = []
    ctrl = core.simulate(Exposure(flux=src.flux(COMP, 0.0),
                                  duration_min=max(times)), n_points=241)
    for d in SOS_DOSES[1:]:
        obs = core.simulate(Exposure(flux=src.flux(COMP, float(d)),
                                     duration_min=max(times)), n_points=241)
        for t in times:
            i = int(np.argmin(np.abs(obs["t"] - t)))
            out.append(obs["reporter_enzyme"][i] / ctrl["reporter_enzyme"][i])
    return np.array(out)


def umu_section() -> dict:
    print("=" * 70)
    print("ENDPOINT 1 — umu (SOS / umuDC-lacZ)")
    print("=" * 70)
    base = SOSCore()
    vals = _values(base, SOS_PARAMS)

    end = relative_sensitivity(lambda **o: _nest(base, o), sos_endpoint_observe,
                               vals, label="umu, endpoint readout (IR + gate)")
    s_end = summarise(end, NOISE, FACTOR)
    print_summary(s_end)

    tc = relative_sensitivity(lambda **o: _nest(base, o), sos_timecourse_observe,
                              vals, label="umu, time-resolved reporter (4 times)")
    s_tc = summarise(tc, NOISE, FACTOR)
    print_summary(s_tc)

    print(f"\n  time-resolved sampling moves the estimable count "
          f"{s_end['practical_rank']} -> {s_tc['practical_rank']}")

    # --- the two structural results, demonstrated rather than inferred ----
    print("\nStructural invariances")
    exp_t = Exposure(flux=TabulatedSource().flux(COMP, 1.0), duration_min=120.0)
    exp_c = Exposure(flux=TabulatedSource().flux(COMP, 0.0), duration_min=120.0)
    ot, oc = base.simulate(exp_t), base.simulate(exp_c)
    irs = [BetaGalONPG(gain=g)(ot)["signal"] / BetaGalONPG(gain=g)(oc)["signal"]
           for g in (1.0, 7.3, 1000.0)]
    gain_exact = max(abs(v - irs[0]) for v in irs) < 1e-12
    print(f"  instrument gain x1 / x7.3 / x1000 -> IR "
          f"{irs[0]:.12f} / {irs[1]:.12f} / {irs[2]:.12f}")
    print(f"    the readout is a ratio, so gain cancels exactly: {gain_exact}")

    def ir_of(core):
        return np.array([r["IR"] for r in dose_series(
            VirtualAssay(source=TabulatedSource(), core=core),
            COMP, SOS_DOSES, protocol=UMU)["rows"][1:]])

    # k_txn and k_tsl enter the reporter linearly, so a ratio readout cannot
    # see them.  Demonstrated by scaling them 1000-fold rather than by a
    # threshold on a finite-difference column, which sits at the solver's own
    # noise floor and would make this a judgement call about tolerances.
    lin = {}
    for name in ("k_txn", "k_tsl"):
        lin[name] = check_symmetry(lambda **o: _nest(base, o), ir_of,
                                   {name: (getattr(base, name), 1)},
                                   scales=(0.01, 100.0, 1000.0))
        print(f"  scaling {name} by 0.01x / 100x / 1000x -> max relative "
              f"drift in IR {lin[name]['max_rel_drift']:.1e}  "
              f"(exact: {lin[name]['exact']})")

    sym = check_symmetry(lambda **o: _nest(base, o),
                         ir_of,
                         {"beta_lexA": (base.beta_lexA, 1),
                          "K_lexA": (base.K_lexA, 1)})
    print(f"  scaling (beta_lexA, K_lexA) together by 0.5 / 1.5 / 4.0:")
    for r in sym["per_scale"]:
        print(f"    x{r['scale']:<4} max relative drift in IR "
              f"{r['max_rel_drift']:.2e}")
    print(f"    exact continuous symmetry: {sym['exact']} — the LexA "
          f"concentration scale is\n    unobservable, so one of the two can "
          f"be fixed by convention rather than fitted.")
    return {"endpoint": s_end, "timecourse": s_tc, "gain_exact": gain_exact,
            "symmetry": sym, "linear_reporter": lin}


# --------------------------------------------------------------------------
# Endpoint 2 — GADD45a-GFP
# --------------------------------------------------------------------------
P53_PARAMS = ("v_repair", "k_repair", "K_dna", "h_dna", "w_dna",
              "gfp_leak", "K_gadd", "h_gadd", "k_maturation", "d_gfp_mrna",
              "K_tox", "h_tox", "p53.k_deg", "p53.K_p53", "p53.h_p53",
              "p53.k_mdm2_txn", "p53.K_arrest")
P53_DOSES = log_doses(0.2, 4.0, 6)


def p53_observe(core):
    a = VirtualAssay(source=TabulatedSource(), core=core,
                     readouts=("reporter_gfp", "growth", "damage_probe"),
                     duration_min=48 * 60, n_points=97)
    ser = dose_series(a, COMP, P53_DOSES, protocol=GADD45A_GFP)
    rows = ser["rows"][1:]
    return np.concatenate([[r["IR"] for r in rows],
                           [r["growth_factor"] for r in rows]])


def p53_section() -> dict:
    print("\n" + "=" * 70)
    print("ENDPOINT 2 — GADD45a-GFP reporter line")
    print("=" * 70)
    base = P53Core()
    s = summarise(relative_sensitivity(
        lambda **o: _nest(base, o), p53_observe, _values(base, P53_PARAMS),
        label="GADD45a-GFP, endpoint fluorescence + density"), NOISE, FACTOR)
    print_summary(s)
    return s


# --------------------------------------------------------------------------
# Endpoint 3 — comet and micronucleus
# --------------------------------------------------------------------------
CYTO_PARAMS = ("k_ner", "k_ber", "k_ligate", "k_rejoin", "f_misrepair",
               "k_ssb_to_dsb", "ssb_spontaneous", "k_clear_acentric",
               "p_spontaneous", "p_max_clast", "K_acentric", "K_tox")
CYTO_DOSES = log_doses(0.05, 3.0, 6)


ALKYLATOR = demo("alkylating agent (MMS-like)")
ALK_DOSES = log_doses(1.0, 60.0, 6)


def _cyto_observe(core, readout, protocol, minutes, points,
                  compound=None, doses=None):
    a = VirtualAssay(source=TabulatedSource(), core=core,
                     readouts=(readout, "growth", "damage_probe"),
                     duration_min=minutes, n_points=points)
    ser = dose_series(a, compound or COMP,
                      CYTO_DOSES if doses is None else doses, protocol=protocol)
    return np.array([r["IR"] for r in ser["rows"][1:]])


def comet_observe(core):
    return _cyto_observe(core, "comet_alkaline", COMET, 4 * 60, 49)


def mn_observe(core):
    return _cyto_observe(core, "micronucleus_cbmn", MN_CBMN, 36 * 60, 73)


def both_observe(core):
    return np.concatenate([comet_observe(core), mn_observe(core)])


def two_compound_observe(core):
    """The same two endpoints, run on a bulky former AND an alkylating agent.

    Identifiability is a property of the experiment, not of the model alone:
    a parameter governing alkyl-adduct excision cannot be estimated from a
    compound that makes no alkyl adducts, however rich the readout.
    """
    return np.concatenate([
        both_observe(core),
        _cyto_observe(core, "comet_alkaline", COMET, 4 * 60, 49,
                      ALKYLATOR, ALK_DOSES),
        _cyto_observe(core, "micronucleus_cbmn", MN_CBMN, 36 * 60, 73,
                      ALKYLATOR, ALK_DOSES)])


def cyto_section() -> dict:
    print("\n" + "=" * 70)
    print("ENDPOINT 3 — comet and micronucleus")
    print("=" * 70)
    base = CytogeneticCore()
    vals = _values(base, CYTO_PARAMS)
    out = {}
    for label, fn in (("comet alone (4 h)", comet_observe),
                      ("micronucleus alone (36 h)", mn_observe),
                      ("comet + micronucleus together", both_observe)):
        s = summarise(relative_sensitivity(
            lambda **o: _nest(base, o), fn, vals,
            label=f"cytogenetic core, {label}"), NOISE, FACTOR)
        print_summary(s)
        out[label] = s
    a, b, c = (out["comet alone (4 h)"]["practical_rank"],
               out["micronucleus alone (36 h)"]["practical_rank"],
               out["comet + micronucleus together"]["practical_rank"])
    print(f"\n  combining the two endpoints: {a} + {b} -> {c} estimable "
          f"directions")
    print("  the two readouts are not redundant — they constrain different "
          "parts of\n  the same damage-processing chain, which is the "
          "quantitative version of\n  the argument for building them on one "
          "core.")

    s2 = summarise(relative_sensitivity(
        lambda **o: _nest(base, o), two_compound_observe, vals,
        label="cytogenetic core, both endpoints x (bulky + alkylating agent)"),
        NOISE, FACTOR)
    print_summary(s2)
    out["two compounds"] = s2
    one = out["comet + micronucleus together"]["column_norms"]["k_ber"]
    two = s2["column_norms"]["k_ber"]
    print(f"\n  k_ber (alkyl-adduct excision) sensitivity: {one:.2e} with the "
          f"bulky probe alone,\n  {two:.2e} once an alkylating agent joins the "
          f"design; estimable directions "
          f"{out['comet + micronucleus together']['practical_rank']} -> "
          f"{s2['practical_rank']}.")
    print("  identifiability is a property of the EXPERIMENT, not of the "
          "model: no readout\n  can estimate a rate for a lesion the probe "
          "compound never makes.")
    return out


# --------------------------------------------------------------------------
def checks(umu, p53, cyto) -> list:
    out = []

    def chk(label, ok, detail):
        out.append((label, bool(ok), detail))
        print(f"  [{'PASS' if ok else 'FAIL'}] {label}: {detail}")

    print("\nStructural checks")

    lin = umu["linear_reporter"]
    chk("reporter rate constants are invisible to a ratio readout",
        all(v["exact"] for v in lin.values()),
        "k_txn and k_tsl scaled over 1e5-fold; max relative drift in IR "
        + " / ".join(f"{n}:{v['max_rel_drift']:.0e}" for n, v in lin.items()))

    chk("instrument gain cancels exactly in the induction ratio",
        umu["gain_exact"], "IR identical to 1e-12 across a 1000x gain change")

    chk("(beta_lexA, K_lexA) is an exact continuous symmetry",
        umu["symmetry"]["exact"],
        f"max relative drift {umu['symmetry']['max_rel_drift']:.1e} over a "
        f"0.5x-4x scaling")

    for tag, s in (("umu", umu["endpoint"]), ("GADD45a-GFP", p53),
                   ("comet+MN", cyto["comet + micronucleus together"])):
        chk(f"{tag}: estimable directions are far fewer than parameters",
            s["practical_rank"] < 0.7 * s["n_parameters"],
            f"{s['practical_rank']} of {s['n_parameters']} at "
            f"{100*NOISE:.0f}% noise")

    # The result here came out against the expectation that more sampling
    # buys more parameters.  It does not, and the reason is the point: the
    # flat directions are exact symmetries, and no sampling schedule removes
    # a symmetry.  What sampling buys is precision on the directions that
    # were already visible.
    e, t = umu["endpoint"], umu["timecourse"]
    gain_ratio = float(np.median(t["spectrum"][:6] / e["spectrum"][:6]))
    c1 = cyto["comet + micronucleus together"]
    c2 = cyto["two compounds"]
    chk("a parameter is identifiable only if the probe compound exercises it",
        c1["column_norms"]["k_ber"] < 1e-6 and c2["column_norms"]["k_ber"] > 1e-2
        and c2["practical_rank"] > c1["practical_rank"],
        f"k_ber sensitivity {c1['column_norms']['k_ber']:.0e} -> "
        f"{c2['column_norms']['k_ber']:.2f} on adding an alkylating agent; "
        f"estimable {c1['practical_rank']} -> {c2['practical_rank']}")

    chk("sampling buys precision, not new directions (structural nulls "
        "survive)",
        gain_ratio > 1.2 and t["practical_rank"] == e["practical_rank"],
        f"visible singular values x{gain_ratio:.2f}, estimable count "
        f"{e['practical_rank']} -> {t['practical_rank']}")

    return out


def figure(umu, p53, cyto):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"\n(figure skipped: {type(e).__name__}: {e})")
        return None

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))

    a = ax[0]
    for s, style in ((umu["endpoint"], "-o"), (umu["timecourse"], "--s")):
        sv = np.maximum(s["spectrum"], 1e-16)
        a.semilogy(range(1, len(sv) + 1), sv, style, ms=4, label=s["label"])
    a.axhline(NOISE / np.log(FACTOR), color="k", ls=":", lw=1.2)
    a.text(0.98, NOISE / np.log(FACTOR) * 1.5, "estimable above this line",
           transform=a.get_yaxis_transform(), ha="right", fontsize=8)
    a.set_xlabel("direction"); a.set_ylabel("singular value")
    a.set_title("umu: what the readout constrains"); a.legend(fontsize=7)

    a = ax[1]
    norms = umu["endpoint"]["column_norms"]
    # The finite-difference column for a provably invariant parameter sits at
    # the solver's noise floor, not at zero.  Plotting that floor on a log
    # axis would show a small bar and imply a small influence, so the
    # demonstrated invariants are drawn at the axis edge and named.
    exact = set(umu["linear_reporter"])
    items = sorted(norms.items(), key=lambda t: t[1])
    floor = 1e-4
    a.barh([n for n, _ in items],
           [floor if n in exact else max(v, floor) for n, v in items],
           color=["#c0392b" if n in exact else "#2a5a9e" for n, _ in items])
    for i, (n, _) in enumerate(items):
        if n in exact:
            a.text(floor * 1.4, i, "exactly 0 (shown by 1e5x scaling)",
                   va="center", fontsize=6.5, color="#c0392b")
    a.set_xscale("log"); a.set_xlim(left=floor * 0.7)
    a.set_xlabel("relative sensitivity magnitude")
    a.set_title("umu: per-parameter influence"); a.tick_params(labelsize=7)

    a = ax[2]
    labels, ranks, totals = [], [], []
    for tag, s in (("umu\nendpoint", umu["endpoint"]),
                   ("umu\ntimecourse", umu["timecourse"]),
                   ("GADD45a\nGFP", p53),
                   ("comet", cyto["comet alone (4 h)"]),
                   ("MN", cyto["micronucleus alone (36 h)"]),
                   ("comet\n+MN", cyto["comet + micronucleus together"]),
                   ("+alkyl\nagent", cyto["two compounds"])):
        labels.append(tag); ranks.append(s["practical_rank"])
        totals.append(s["n_parameters"])
    x = np.arange(len(labels))
    a.bar(x, totals, color="#d9dde3", label="parameters")
    a.bar(x, ranks, color="#0e6570", label=f"estimable at {100*NOISE:.0f}% noise")
    a.set_xticks(x); a.set_xticklabels(labels, fontsize=7)
    a.set_ylabel("count"); a.set_title("how much is learnable")
    a.legend(fontsize=7)

    fig.suptitle("Identifiability of the genotox cores — what a fit to each "
                 "readout could recover", fontsize=11)
    fig.tight_layout()
    os.makedirs(FIGDIR, exist_ok=True)
    path = os.path.join(FIGDIR, "genotox_identifiability.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    print(f"\nfigure -> {path}")
    return path


def report() -> dict:
    print("GENOTOX DIAGNOSTICS — local identifiability of each endpoint")
    print(f"Relative sensitivities, {100*NOISE:.0f}% measurement noise, "
          f"'estimable' = pinned to better than a factor of {FACTOR}.")
    print("Local to one operating point; says nothing about whether the "
          "model is right.\n")
    umu = umu_section()
    p53 = p53_section()
    cyto = cyto_section()
    ck = checks(umu, p53, cyto)
    path = figure(umu, p53, cyto)
    n_ok = sum(1 for _, ok, _ in ck if ok)
    print(f"\n{n_ok}/{len(ck)} structural checks passed.")
    return {"umu": umu, "p53": p53, "cyto": cyto, "checks": ck, "figure": path}


if __name__ == "__main__":
    report()
