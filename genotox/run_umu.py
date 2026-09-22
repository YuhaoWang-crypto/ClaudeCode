"""
Step 1 runner — the umu test as a virtual assay.

Prints the core's set-point diagnostics, a time course, dose-response tables
for five compounds chosen to span the behaviours the architecture has to
distinguish, and four structural checks.  Writes ``figures/genotox_umu.png``.

    python3 -m genotox.run_umu
"""
from __future__ import annotations

import os

import numpy as np

from .assay import VirtualAssay
from .damage import DEMO_COMPOUNDS, TabulatedSource
from .doseresponse import (GROWTH_GATE, IR_THRESHOLD, call_result,
                           dose_series, log_doses)

FIGDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      "figures")

# (compound index, +/-S9, dose range) — ranges differ because potencies do
PANEL = [
    (0, False, (0.02, 20.0)),
    (1, False, (0.05, 50.0)),
    (1, True, (0.05, 50.0)),
    (2, False, (0.5, 500.0)),
    (3, False, (0.5, 500.0)),
    (4, False, (0.02, 20.0)),
    (5, False, (1.0, 1000.0)),
]


def build() -> VirtualAssay:
    return VirtualAssay(source=TabulatedSource())


def setpoint_report(assay) -> dict:
    core = assay.core
    L0 = core._L_basal
    a_basal = core.promoter(L0)
    a_full = core.promoter(0.0)
    print("SOS core set-point")
    print(f"  LexA basal                 L0      = {L0:.3f}")
    print(f"  umuDC half-repression      K_umu   = {core.K_umu:.3f} "
          f"(= {core.kumu_ratio:.2f} x L0; small ratio = tightly repressed "
          f"'late' SOS gene)")
    print(f"  promoter activity  basal / full     = {a_basal:.4f} / "
          f"{a_full:.4f}   -> max dynamic range {a_full / a_basal:.1f}x")
    print(f"  doubling time                      = "
          f"{np.log(2) / core.mu_max:.0f} min;  exposure "
          f"{assay.duration_min:.0f} min")
    return {"L0": L0, "range": a_full / a_basal}


def timecourse_report(assay, comp, dose, s9=False) -> dict:
    tr = assay.trajectory(comp, dose, s9=s9)
    print(f"\nTime course — {comp.name} @ {dose} uM "
          f"({'+' if s9 else '-'}S9)")
    print(f"  {'t/min':>6} {'lesions':>9} {'RecA*':>7} {'LexA':>7} "
          f"{'P_umu':>7} {'beta-gal':>9} {'viab':>6}")
    for i in np.linspace(0, len(tr["t"]) - 1, 7).astype(int):
        print(f"  {tr['t'][i]:6.0f} {tr['lesions'][i]:9.3f} "
              f"{tr['recA_active'][i]:7.3f} {tr['LexA'][i]:7.3f} "
              f"{tr['promoter'][i]:7.3f} {tr['reporter_enzyme'][i]:9.2f} "
              f"{tr['viability'][i]:6.3f}")
    return tr


def series_report(assay, comp, s9, lo, hi) -> tuple:
    ser = dose_series(assay, comp, log_doses(lo, hi, 13), s9=s9)
    res = call_result(ser)
    tag = "+S9" if s9 else "-S9"
    print(f"\n{comp.name}  [{tag}]")
    print(f"  channels: {comp.per_uM.nonzero() or '{} (no DNA lesion)'}")
    print(f"  {'dose/uM':>10} {'IR':>7} {'growth':>7} {'lesions':>8}  valid")
    for r in ser["rows"]:
        flag = "" if r["valid"] else "  <- gated"
        print(f"  {r['dose_uM']:10.3f} {r['IR']:7.2f} "
              f"{r['growth_factor']:7.2f} {r['lesions_end']:8.2f}"
              f"  {str(r['valid']):5s}{flag}")
    ec = res["ec_ir"]
    print(f"  -> {res['verdict']}  ({res['reason']}); "
          f"max IR(valid) = {res['max_IR_valid']:.2f}; "
          f"EC-IR{IR_THRESHOLD} = "
          f"{'n/a' if ec is None else f'{ec:.3f} uM'}")
    return ser, res


def checks(results: dict, series: dict) -> list:
    """Structural properties the architecture must have, not curve fits.

    Each one fails loudly if a layer is wired wrong, which is what makes the
    seams worth having.  None of them assert a measured value.
    """
    out = []

    def chk(label, ok, detail):
        out.append((label, bool(ok), detail))
        print(f"  [{'PASS' if ok else 'FAIL'}] {label}: {detail}")

    print("\nStructural checks")

    d = results[("direct-acting bulky (4NQO-like)", False)]
    chk("direct-acting compound is positive without S9",
        d["verdict"] == "POSITIVE", f"verdict={d['verdict']}")

    m9, p9 = (results[("promutagen (2AA-like)", False)],
              results[("promutagen (2AA-like)", True)])
    chk("promutagen call flips on metabolic activation",
        m9["verdict"] != "POSITIVE" and p9["verdict"] == "POSITIVE",
        f"-S9 {m9['verdict']} (max IR {m9['max_IR_valid']:.2f}) -> "
        f"+S9 {p9['verdict']} (max IR {p9['max_IR_valid']:.2f})")

    # Test the mechanism, not the readout.  The claim is that the aneugenic
    # channel has no route into this core; whether the raw induction ratio
    # stays low is a separate question, answered by the next check.
    assay = build()
    comp = DEMO_COMPOUNDS[3]
    ctrl = assay.well(comp, 0.0)
    hi = assay.well(comp, 500.0)
    a = results[("aneugen (colchicine-like)", False)]
    prom_fold = hi["promoter_end"] / ctrl["promoter_end"]
    # The claim is that the promoter is never INDUCED, not that it is
    # unchanged: LexA is cleared only by growth dilution in this core, so
    # arresting the culture raises LexA and deepens repression.  That the
    # promoter goes *down* while the induction ratio goes up is the whole
    # point of the following check.
    chk("aneugen never engages the SOS pathway (weight 0 on its channel)",
        hi["lesions_end"] == 0.0 and hi["recA_active_end"] == 0.0
        and prom_fold <= 1.0,
        f"at 500 uM: lesions = {hi['lesions_end']:.3f}, "
        f"RecA* = {hi['recA_active_end']:.3f}, "
        f"P_umuDC fold vs control = {prom_fold:.3f} (<= 1: repressed further, "
        f"never induced)")

    # A finding, kept as a check so it cannot quietly regress: growth arrest
    # alone manufactures an induction ratio, because beta-gal stops being
    # diluted.  Nothing in the signal path moved.  The gate is what catches
    # it -- which is the argument for reporting the gate with every number.
    gated_only = all(not r["valid"] for r in series[("aneugen (colchicine-like)", False)]["rows"]
                     if r["IR"] >= IR_THRESHOLD)
    chk("apparent aneugen 'induction' is a growth-arrest artifact, gated out",
        a["max_IR_any"] >= IR_THRESHOLD > a["max_IR_valid"] and gated_only,
        f"max IR over ALL wells = {a['max_IR_any']:.2f} while the promoter "
        f"fell to {prom_fold:.2f}x; max IR among VALID wells = "
        f"{a['max_IR_valid']:.2f}")

    c = results[("non-genotoxic cytotoxicant", False)]
    chk("pure cytotoxicant is not called positive",
        c["verdict"] != "POSITIVE",
        f"verdict={c['verdict']} ({c['n_valid']}/{c['n_wells']} wells "
        f"passed the growth gate)")

    i = results[("masked genotoxicant", False)]
    chk("threshold crossed only in gated wells -> INCONCLUSIVE, not NEGATIVE",
        i["verdict"] == "INCONCLUSIVE"
        and i["max_IR_any"] >= IR_THRESHOLD > i["max_IR_valid"],
        f"verdict={i['verdict']}; max IR valid={i['max_IR_valid']:.2f} "
        f"vs any={i['max_IR_any']:.2f}")

    return out


def figure(assay, series_by_key, traj):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"\n(figure skipped: {type(e).__name__}: {e})")
        return None

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))

    a = ax[0]
    a.plot(traj["t"], traj["lesions"] / max(traj["lesions"].max(), 1e-9),
           label="lesions (norm.)")
    a.plot(traj["t"], traj["LexA"] / traj["LexA"][0], label="LexA / LexA$_0$")
    a.plot(traj["t"], traj["promoter"] / traj["promoter"][0],
           label="P$_{umuDC}$ (fold)")
    a.plot(traj["t"], traj["reporter_enzyme"] / traj["reporter_enzyme"][0],
           label=r"$\beta$-gal (fold)")
    a.set_xlabel("time / min"); a.set_ylabel("relative")
    a.set_title("SOS cascade, 2 h exposure"); a.legend(fontsize=8)

    a = ax[1]
    for (name, s9), ser in series_by_key.items():
        rows = ser["rows"][1:]
        d = [r["dose_uM"] for r in rows]
        ir = [r["IR"] for r in rows]
        ok = [r["valid"] for r in rows]
        lbl = f"{name.split(' (')[0]}{' +S9' if s9 else ''}"
        line, = a.plot(d, ir, "-", lw=1.4, label=lbl)
        a.plot([x for x, v in zip(d, ok) if v],
               [y for y, v in zip(ir, ok) if v], "o", ms=4,
               color=line.get_color())
        a.plot([x for x, v in zip(d, ok) if not v],
               [y for y, v in zip(ir, ok) if not v], "x", ms=5,
               color=line.get_color())
    a.axhline(IR_THRESHOLD, color="k", ls="--", lw=1)
    a.set_xscale("log"); a.set_xlabel("dose / uM"); a.set_ylabel("induction ratio")
    a.set_title(f"dose-response (x = growth-gated)"); a.legend(fontsize=7)

    a = ax[2]
    for (name, s9), ser in series_by_key.items():
        rows = ser["rows"][1:]
        a.plot([r["dose_uM"] for r in rows],
               [r["growth_factor"] for r in rows], "-o", ms=3,
               label=f"{name.split(' (')[0]}{' +S9' if s9 else ''}")
    a.axhline(GROWTH_GATE, color="k", ls="--", lw=1)
    a.set_xscale("log"); a.set_xlabel("dose / uM"); a.set_ylabel("growth factor")
    a.set_title("validity gate"); a.legend(fontsize=7)

    fig.suptitle("Virtual umu test — SOS/umuDC-lacZ core "
                 "(illustrative parameters, not fitted)", fontsize=11)
    fig.tight_layout()
    os.makedirs(FIGDIR, exist_ok=True)
    path = os.path.join(FIGDIR, "genotox_umu.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    print(f"\nfigure -> {path}")
    return path


def report() -> dict:
    print("=" * 68)
    print("GENOTOX STEP 1 — virtual umu test (SOS/umuDC-lacZ)")
    print("=" * 68)
    assay = build()
    print(assay.description)
    print("Parameters are illustrative (H): curve SHAPE and relative "
          "behaviour are meaningful, absolute EC values are not.\n")

    sp = setpoint_report(assay)

    series_by_key, results = {}, {}
    for idx, s9, (lo, hi) in PANEL:
        comp = DEMO_COMPOUNDS[idx]
        ser, res = series_report(assay, comp, s9, lo, hi)
        series_by_key[(comp.name, s9)] = ser
        results[(comp.name, s9)] = res

    traj = timecourse_report(assay, DEMO_COMPOUNDS[0], 1.0)
    ck = checks(results, series_by_key)
    path = figure(assay, series_by_key, traj)

    n_ok = sum(1 for _, ok, _ in ck if ok)
    print(f"\n{n_ok}/{len(ck)} structural checks passed.")
    return {"setpoint": sp, "series": series_by_key, "results": results,
            "checks": ck, "figure": path}


if __name__ == "__main__":
    report()
