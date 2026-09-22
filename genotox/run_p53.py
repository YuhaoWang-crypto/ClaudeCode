"""
Step 2 runner — mammalian GADD45a-GFP reporter as a virtual assay.

Reuses the upstream chemistry layer, the readout registry and the decision
layer from step 1 unchanged; only the core is new.  Prints the p53 pulse
diagnostics, dose-response tables, an endpoint comparison against the umu
core on the *same* compounds, and five structural checks.  Writes
``figures/genotox_p53.png``.

    python3 -m genotox.run_p53
"""
from __future__ import annotations

import os

import numpy as np

from .assay import VirtualAssay
from .core import SOSCore
from .damage import DEMO_COMPOUNDS, TabulatedSource
from .doseresponse import (GADD45A_GFP, UMU, call_result, dose_series,
                           log_doses)
from .p53 import P53Core

FIGDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      "figures")

EXPOSURE_MIN = 48 * 60          # 48 h, as reporter-line protocols run
N_POINTS = 577                  # ~5 min grid; fine enough to resolve pulses

PANEL = [
    (0, False, (0.02, 20.0)),     # direct-acting bulky
    (1, True, (0.05, 50.0)),      # promutagen, +S9
    (3, False, (0.5, 500.0)),     # aneugen  <- the key contrast with umu
    (5, False, (1.0, 1000.0)),    # non-genotoxic cytotoxicant
]


def build() -> VirtualAssay:
    return VirtualAssay(source=TabulatedSource(), core=P53Core(),
                        readouts=("reporter_gfp", "growth", "damage_probe",
                                  "pulse_probe"),
                        duration_min=EXPOSURE_MIN, n_points=N_POINTS)


def build_umu() -> VirtualAssay:
    return VirtualAssay(source=TabulatedSource(), core=SOSCore())


def pulse_report(assay) -> dict:
    """Show that the undamaged state is quiet and the damaged state pulses."""
    from .readouts import REGISTRY
    print("p53 dynamics (48 h exposure)")
    print(f"  {'dose/uM':>9} {'pulses':>7} {'p53 peak':>9} {'p53 mean':>9} "
          f"{'AUC':>10} {'GFP/cell':>9}")
    rows = []
    for d in (0.0, 0.1, 0.3, 0.5, 0.7, 1.0, 2.0, 3.0):
        obs = assay.trajectory(DEMO_COMPOUNDS[0], d)
        pp = REGISTRY["pulse_probe"](obs)
        gf = REGISTRY["reporter_gfp"](obs)
        rows.append({"dose": d, **pp, **gf})
        print(f"  {d:9.2f} {pp['p53_pulses']:7d} {pp['p53_peak']:9.3f} "
              f"{pp['p53_mean']:9.3f} {pp['p53_auc']:10.1f} "
              f"{gf['fluorescence_per_cell']:9.2f}")
    return {"rows": rows}


def integrator_check(rows) -> tuple:
    """What does the fluorescence actually track?

    Three candidate predictors, tested on *induced* quantities (each minus
    its unexposed control, since a stable reporter carries a large constant
    background that would otherwise flatter every hypothesis equally):

      * the time-integral of GADD45a promoter activity,
      * the time-integral of p53,
      * the peak p53 height.

    A predictor proportional to the output gives a constant output/predictor
    ratio across doses, so the coefficient of variation of that ratio ranks
    the hypotheses.  Doses whose induction is not resolvable are excluded --
    a ratio of two near-zero numbers is noise, not evidence.
    """
    ctrl = [x for x in rows if x["dose"] == 0][0]
    r = [x for x in rows if x["dose"] > 0]

    def delta(key):
        return np.array([x[key] - ctrl[key] for x in r], float)

    gfp = delta("fluorescence_per_cell")
    cand = {"promoter integral": delta("promoter_auc"),
            "p53 integral": delta("p53_auc"),
            "p53 peak": delta("p53_peak")}

    keep = gfp > 0.02 * gfp.max()
    n_drop = int((~keep).sum())
    gfp = gfp[keep]

    print("\n  predictor of induced GFP        ratio CV   (lower = better)")
    scores = {}
    for label, x in cand.items():
        xv = x[keep]
        ratio = gfp / np.where(xv > 0, xv, np.nan)
        cv = float(np.nanstd(ratio) / np.nanmean(ratio))
        scores[label] = cv
        print(f"    {label:28s} {cv:8.3f}")
    winner = min(scores, key=scores.get)
    print(f"  -> fluorescence tracks the {winner}"
          + (f"   ({n_drop} sub-resolution dose(s) excluded)" if n_drop else ""))
    return winner, scores


def series_report(assay, comp, s9, lo, hi, protocol) -> tuple:
    ser = dose_series(assay, comp, log_doses(lo, hi, 13), s9=s9,
                      protocol=protocol)
    res = call_result(ser)
    tag = "+S9" if s9 else "-S9"
    print(f"\n{comp.name}  [{tag}]  protocol={protocol.name}")
    print(f"  {'dose/uM':>10} {'IR':>7} {'rel.dens':>9} {'lesions':>8} "
          f"{'mit.stress':>10} {'pulses':>7}  valid")
    for r in ser["rows"]:
        w = assay.well(comp, r["dose_uM"], s9=s9)
        flag = "" if r["valid"] else "  <- gated"
        print(f"  {r['dose_uM']:10.3f} {r['IR']:7.2f} "
              f"{r['growth_factor']:9.2f} {r['lesions_end']:8.2f} "
              f"{w.get('mitotic_stress_end', float('nan')):10.2f} "
              f"{r.get('p53_pulses', 0):7d}  {str(r['valid']):5s}{flag}")
    ec = res["ec_ir"]
    print(f"  -> {res['verdict']}  ({res['reason']}); "
          f"max IR(valid) = {res['max_IR_valid']:.2f}; EC-IR"
          f"{protocol.ir_threshold} = {'n/a' if ec is None else f'{ec:.3f} uM'}")
    return ser, res


def cross_endpoint(assay_p53, assay_umu) -> list:
    """The same compounds, the same upstream, two endpoints.

    This table is the argument for the whole architecture: one chemistry
    prediction, read two ways, giving two different and individually correct
    answers.
    """
    print("\nCross-endpoint comparison (identical DamageFlux inputs)")
    print(f"  {'compound':34s} {'umu':>12} {'GADD45a-GFP':>14}")
    out = {}
    for idx, s9, (lo, hi) in [(0, False, (0.02, 20.0)),
                              (3, False, (0.5, 500.0))]:
        comp = DEMO_COMPOUNDS[idx]
        u = call_result(dose_series(assay_umu, comp, log_doses(lo, hi, 13),
                                    s9=s9, protocol=UMU))
        p = call_result(dose_series(assay_p53, comp, log_doses(lo, hi, 13),
                                    s9=s9, protocol=GADD45A_GFP))
        print(f"  {comp.name:34s} {u['verdict']:>12} {p['verdict']:>14}")
        # Promoter engagement relative to solvent control, maximised over
        # the dose range rather than read at the top dose.  Both cores can
        # lose their signal at the top dose for mechanistic reasons -- the
        # mammalian mitotic-surveillance route needs cells to still be
        # cycling, and a complete mitotic block switches it off -- so "does
        # this compound ever engage the pathway" is a max, not an endpoint.
        doses = log_doses(lo, hi, 13)
        uc = assay_umu.well(comp, 0.0, s9=s9)["promoter_end"]
        pc = assay_p53.well(comp, 0.0, s9=s9)["promoter_auc"]
        ufold = [(assay_umu.well(comp, d, s9=s9)["promoter_end"] / uc, d)
                 for d in doses]
        pfold = [(assay_p53.well(comp, d, s9=s9)["promoter_auc"] / pc, d)
                 for d in doses]
        umax, udose = max(ufold)
        pmax, pdose = max(pfold)
        out[comp.name] = {
            "umu": u["verdict"], "p53": p["verdict"],
            "umu_max_any": u["max_IR_any"], "p53_max_any": p["max_IR_any"],
            "umu_promoter_fold": umax, "umu_promoter_dose": udose,
            "p53_promoter_fold": pmax, "p53_promoter_dose": pdose,
        }
    an = out["aneugen (colchicine-like)"]
    print("\n  Verdicts alone are not the contrast, because both are "
          "growth-gated.\n  Compare genuine pathway engagement instead:")
    print(f"    SOS  P_umuDC       {an['umu_promoter_fold']:6.3f} x control "
          f"(best over all doses)  -- no route, never induced")
    print(f"    GADD45a promoter  {an['p53_promoter_fold']:6.2f} x control "
          f"(peaks at {an['p53_promoter_dose']:.2f} uM)  -- mitotic")
    print("    the mammalian response peaks mid-range and falls again: a "
          "complete\n    mitotic block stops the cycling that surveillance "
          "needs.")
    print("  One scalar 'genotoxic potency' could not produce both.")
    return out


def gate_sensitivity(assay) -> dict:
    """How much of the aneugen verdict is biology and how much is protocol?

    An aneugen's p53 response and its antiproliferative action appear at
    overlapping doses, so the wells that show induction are also the wells
    that fail a strict density gate.  Re-scoring one unchanged simulation
    under different gates separates the two: if the verdict moves, the call
    is being made by the protocol constant, not by the model.
    """
    from dataclasses import replace as _replace
    comp = DEMO_COMPOUNDS[3]
    doses = log_doses(0.5, 500.0, 13)
    print("\nGate sensitivity — aneugen, one simulation scored three ways")
    out = {}
    for gate in (0.80, 0.70, 0.60):
        proto = _replace(GADD45A_GFP, growth_gate=gate)
        res = call_result(dose_series(assay, comp, doses, protocol=proto))
        out[gate] = res["verdict"]
        print(f"  density gate {gate:.2f} -> {res['verdict']:13s} "
              f"(max IR valid {res['max_IR_valid']:.2f}, "
              f"any {res['max_IR_any']:.2f})")
    print("  the underlying trajectories are identical; only the operational\n"
          "  cut-off differs.  An aneugen call from a reporter line is a\n"
          "  protocol decision as much as a biological one.")
    return out


def checks(results, pulses, integ, cross_detail, gates) -> list:
    out = []

    def chk(label, ok, detail):
        out.append((label, bool(ok), detail))
        print(f"  [{'PASS' if ok else 'FAIL'}] {label}: {detail}")

    print("\nStructural checks")

    quiet = pulses["rows"][0]
    loud = [r for r in pulses["rows"] if r["dose"] >= 1.0][0]
    chk("undamaged state is quiet, damaged state pulses",
        quiet["p53_pulses"] == 0 and loud["p53_pulses"] >= 3,
        f"0 uM -> {quiet['p53_pulses']} pulses; "
        f"{loud['dose']} uM -> {loud['p53_pulses']} pulses")

    chk("fluorescence tracks the promoter integral, not p53 itself",
        integ[0] == "promoter integral",
        "ratio CVs: " + ", ".join(f"{k}={v:.3f}" for k, v in integ[1].items()))

    d = results[("direct-acting bulky (4NQO-like)", False)]
    chk("direct-acting genotoxicant is positive",
        d["verdict"] == "POSITIVE", f"verdict={d['verdict']}")

    # Compare pathway engagement, not raw induction ratios: in a
    # growth-arrested culture the ratio moves for reasons unrelated to the
    # signal (see run_umu's artifact check), so it cannot settle this.
    an = cross_detail["aneugen (colchicine-like)"]
    chk("aneugen reaches p53 by a route the bacterial core simply lacks",
        an["umu_promoter_fold"] <= 1.001 and an["p53_promoter_fold"] > 1.0,
        f"best P_umuDC fold = {an['umu_promoter_fold']:.3f} (no route at all) "
        f"vs best GADD45a promoter fold = {an['p53_promoter_fold']:.2f} "
        f"at {an['p53_promoter_dose']:.2f} uM")

    # Recorded because it came out against expectation.  The mitotic route
    # is real but weak, and it is self-limiting: the drug that trips
    # surveillance also stops the cycling surveillance needs.  Most of what
    # the plate reader sees for an aneugen is the growth-arrest artifact, so
    # a reporter line is the wrong instrument for this class -- which is the
    # case for reading the aneugenic channel directly, at endpoint #3.
    ratio = an["p53_max_any"] / an["p53_promoter_fold"]
    chk("reporter overstates the aneugen response (artifact > surveillance)",
        ratio > 1.5,
        f"max IR {an['p53_max_any']:.2f} vs genuine promoter engagement "
        f"{an['p53_promoter_fold']:.2f}x -> {ratio:.1f}x overstatement")

    chk("the aneugen verdict is set by the density gate, not by the model",
        len(set(gates.values())) > 1,
        " / ".join(f"gate {g:.2f}: {v}" for g, v in gates.items()))

    c = results[("non-genotoxic cytotoxicant", False)]
    chk("pure cytotoxicant is not called positive",
        c["verdict"] != "POSITIVE",
        f"verdict={c['verdict']} ({c['n_valid']}/{c['n_wells']} wells "
        f"passed the gate)")

    return out


def figure(assay, series_by_key):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"\n(figure skipped: {type(e).__name__}: {e})")
        return None

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))

    a = ax[0]
    for d in (0.0, 0.3, 1.0, 3.0):
        tr = assay.trajectory(DEMO_COMPOUNDS[0], d)
        a.plot(tr["t"] / 60.0, tr["p53"], lw=1.3, label=f"{d} uM")
    a.set_xlabel("time / h"); a.set_ylabel("p53 (a.u.)")
    a.set_title("p53 pulses (Mdm2 delayed feedback)"); a.legend(fontsize=8)

    a = ax[1]
    for d in (0.0, 0.3, 1.0, 3.0):
        tr = assay.trajectory(DEMO_COMPOUNDS[0], d)
        a.plot(tr["t"] / 60.0, tr["mature_gfp"], lw=1.3, label=f"{d} uM")
    a.set_xlabel("time / h"); a.set_ylabel("mature GFP / cell")
    a.set_title("reporter integrates the pulse train"); a.legend(fontsize=8)

    a = ax[2]
    for (name, s9), ser in series_by_key.items():
        rows = ser["rows"][1:]
        d = [r["dose_uM"] for r in rows]
        ir = [r["IR"] for r in rows]
        ok = [r["valid"] for r in rows]
        line, = a.plot(d, ir, "-", lw=1.4,
                       label=f"{name.split(' (')[0]}{' +S9' if s9 else ''}")
        a.plot([x for x, v in zip(d, ok) if v],
               [y for y, v in zip(ir, ok) if v], "o", ms=4,
               color=line.get_color())
        a.plot([x for x, v in zip(d, ok) if not v],
               [y for y, v in zip(ir, ok) if not v], "x", ms=5,
               color=line.get_color())
    a.axhline(GADD45A_GFP.ir_threshold, color="k", ls="--", lw=1)
    a.set_xscale("log"); a.set_xlabel("dose / uM")
    a.set_ylabel("GFP induction ratio")
    a.set_title("dose-response (x = density-gated)"); a.legend(fontsize=7)

    fig.suptitle("Virtual GADD45a-GFP reporter assay — p53/Mdm2 core "
                 "(illustrative parameters, not fitted)", fontsize=11)
    fig.tight_layout()
    os.makedirs(FIGDIR, exist_ok=True)
    path = os.path.join(FIGDIR, "genotox_p53.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    print(f"\nfigure -> {path}")
    return path


def report() -> dict:
    print("=" * 68)
    print("GENOTOX STEP 2 — virtual GADD45a-GFP reporter (p53/Mdm2)")
    print("=" * 68)
    assay = build()
    print(assay.description)
    print("Parameters are illustrative (H): pulse period was tuned toward "
          "the reported\nfew-hour range, not fitted to a time series.\n")

    pulses = pulse_report(assay)
    integ = integrator_check(pulses["rows"])

    series_by_key, results = {}, {}
    for idx, s9, (lo, hi) in PANEL:
        comp = DEMO_COMPOUNDS[idx]
        ser, res = series_report(assay, comp, s9, lo, hi, GADD45A_GFP)
        series_by_key[(comp.name, s9)] = ser
        results[(comp.name, s9)] = res

    cross = cross_endpoint(assay, build_umu())
    gates = gate_sensitivity(assay)
    ck = checks(results, pulses, integ, cross, gates)
    path = figure(assay, series_by_key)

    n_ok = sum(1 for _, ok, _ in ck if ok)
    print(f"\n{n_ok}/{len(ck)} structural checks passed.")
    return {"pulses": pulses, "integrator": integ, "series": series_by_key,
            "results": results, "cross": cross, "gates": gates,
            "checks": ck, "figure": path}


if __name__ == "__main__":
    report()
