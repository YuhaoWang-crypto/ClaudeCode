"""
Step 3 runner — comet and micronucleus from ONE core.

Two assays over the same :class:`~genotox.cytogenetic.CytogeneticCore`,
differing only in exposure length and readout: comet at 4 h with no division,
CBMN at 36 h scored per binucleated cell.  Prints the damage-processing time
course, both dose-response tables, the mechanism-classification table that
endpoint #2 could not produce, two experiments on the model, and seven
structural checks.  Writes ``figures/genotox_comet_mn.png``.

    python3 -m genotox.run_comet_mn
"""
from __future__ import annotations

import os
from dataclasses import replace as dc_replace

import numpy as np

from .assay import VirtualAssay
from .cytogenetic import CytogeneticCore
from .damage import Compound, DamageFlux, TabulatedSource, demo
from .doseresponse import (COMET, MN_CBMN, call_result, dose_series,
                           log_doses)

FIGDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      "figures")

COMET_MIN = 4 * 60        # short: repair intermediates still open
MN_MIN = 36 * 60          # long enough for ~1 division under cytochalasin B

#: A clastogen whose damage is frank strand breakage rather than adducts,
#: so the comet/MN contrast with the aneugen is not confounded by repair
#: kinetics.  Illustrative (H).
DIRECT_CLASTOGEN = Compound(
    name="direct clastogen (bleomycin-like)",
    per_uM=DamageFlux(dsb=0.020, ssb=0.120),
    direct_fraction=1.0,
    note="frank breaks; comet-positive and MN-positive, centromere-negative",
)

PANEL = [demo("direct-acting bulky (4NQO-like)"),
         demo("alkylating agent (MMS-like)"),
         DIRECT_CLASTOGEN,
         demo("aneugen (colchicine-like)"),
         demo("mixed clastogen/aneugen"),
         demo("non-genotoxic cytotoxicant")]

#: centromere-positive fraction above which an aneugenic component is
#: called.  Operational, like every other threshold in this package: it lives
#: here rather than in the biology so it can be moved without anyone
#: suspecting the model changed.
CPOS_CALL = 20.0

RANGES = {
    "direct-acting bulky (4NQO-like)": (0.02, 20.0),
    "alkylating agent (MMS-like)": (0.5, 500.0),
    "direct clastogen (bleomycin-like)": (0.5, 500.0),
    "aneugen (colchicine-like)": (0.05, 50.0),
    "mixed clastogen/aneugen": (0.05, 50.0),
    "non-genotoxic cytotoxicant": (1.0, 1000.0),
}


def build_comet(core=None) -> VirtualAssay:
    return VirtualAssay(source=TabulatedSource(),
                        core=core or CytogeneticCore(),
                        readouts=("comet_alkaline", "growth", "damage_probe"),
                        duration_min=COMET_MIN, n_points=121)


def build_mn(core=None) -> VirtualAssay:
    return VirtualAssay(source=TabulatedSource(),
                        core=core or CytogeneticCore(),
                        readouts=("micronucleus_cbmn", "growth",
                                  "damage_probe"),
                        duration_min=MN_MIN, n_points=289)


# --------------------------------------------------------------------------
def processing_report(assay) -> None:
    """Where the damage sits over time, for a bulky-adduct former."""
    tr = assay.trajectory(demo("direct-acting bulky (4NQO-like)"), 1.0)
    print("Damage processing — bulky adduct former @ 1 uM, 4 h")
    print(f"  {'t/min':>6} {'bulky':>8} {'SSB':>8} {'DSB':>8} "
          f"{'acentric':>9}")
    for i in np.linspace(0, len(tr["t"]) - 1, 6).astype(int):
        print(f"  {tr['t'][i]:6.0f} {tr['adduct_bulky'][i]:8.3f} "
              f"{tr['ssb'][i]:8.3f} {tr['dsb'][i]:8.4f} "
              f"{tr['acentric'][i]:9.5f}")
    print("  the SSB pool is excision intermediate, not primary damage:\n"
          "  nothing cut the DNA, repair did.")


def series_report(assay, comp, protocol, extra=()) -> tuple:
    lo, hi = RANGES[comp.name]
    ser = dose_series(assay, comp, log_doses(lo, hi, 11), protocol=protocol)
    res = call_result(ser)
    print(f"\n{comp.name}  [{protocol.name}]")
    head = f"  {'dose/uM':>10} {'IR':>7} {'gate':>7} {'signal':>9}"
    for e in extra:
        head += f" {e[1]:>9}"
    print(head + "  valid")
    for r in ser["rows"]:
        w = assay.well(comp, r["dose_uM"])
        line = (f"  {r['dose_uM']:10.3f} {r['IR']:7.2f} "
                f"{r['growth_factor']:7.2f} {r['signal']:9.2f}")
        for e in extra:
            line += f" {w.get(e[0], float('nan')):9.2f}"
        print(line + f"  {str(r['valid']):5s}"
              + ("" if r["valid"] else "  <- gated"))
    print(f"  -> {res['verdict']} ({res['reason']}); "
          f"max IR(valid) = {res['max_IR_valid']:.2f}")
    return ser, res


def mechanism_table(comet, mn) -> dict:
    """The payoff: name the mechanism, not just the call.

    Step 2 established that a transcriptional reporter cannot distinguish a
    clastogen from an aneugen -- both just raise p53.  Two observables fix
    that: comet answers "was the DNA broken", centromere status answers "was
    a whole chromosome lost".
    """
    print("\nMechanism classification (same upstream, two instruments)")
    print(f"  {'compound':34s} {'comet':>8} {'MN':>8} {'%MN C+':>8}  "
          f"mechanism")
    out = {}
    for comp in PANEL:
        lo, hi = RANGES[comp.name]
        doses = log_doses(lo, hi, 11)
        c = call_result(dose_series(comet, comp, doses, protocol=COMET))
        m_ser = dose_series(mn, comp, doses, protocol=MN_CBMN)
        m = call_result(m_ser)
        # centromere composition at the most-induced valid MN well
        valid = [r for r in m_ser["rows"] if r["valid"] and r["dose_uM"] > 0]
        best = max(valid, key=lambda r: r["IR"]) if valid else None
        cpos = (mn.well(comp, best["dose_uM"])["mn_centromere_pos_pct"]
                if best else float("nan"))

        # The two observables evidence two components INDEPENDENTLY; they do
        # not select between them.  Written as an either/or, this classifier
        # called a compound that is plainly comet-positive "aneugenic" purely
        # because most of its micronuclei carried centromeres, and the
        # clastogenic half vanished from the verdict.  Real chemicals are not
        # obliged to use one channel.
        comet_pos = c["verdict"] == "POSITIVE"
        mn_pos = m["verdict"] == "POSITIVE"
        aneugenic = mn_pos and cpos >= CPOS_CALL          # whole chromosomes lost
        clastogenic = mn_pos and comet_pos                # DNA demonstrably broken
        if aneugenic and clastogenic:
            mech = f"MIXED ({cpos:.0f}% C+ micronuclei, comet-positive)"
        elif aneugenic:
            mech = "ANEUGENIC (whole chromosomes)"
        elif clastogenic:
            mech = "CLASTOGENIC (strand breakage)"
        elif mn_pos:
            mech = "clastogenic, breaks below comet resolution"
        elif comet_pos:
            mech = "breaks without chromosome damage"
        else:
            mech = "negative"
        print(f"  {comp.name:34s} {c['verdict'][:4]:>8} {m['verdict'][:4]:>8} "
              f"{cpos:8.1f}  {mech}")
        out[comp.name] = {"comet": c["verdict"], "mn": m["verdict"],
                          "pct_cpos": cpos, "mechanism": mech,
                          "aneugenic": aneugenic, "clastogenic": clastogenic,
                          "mn_series": m_ser}
    return out


def repair_experiment() -> dict:
    """Does faster excision repair make a bulky compound look MORE damaged?

    Not a curve fit -- a manipulation of the model.  Re-run one dose of one
    compound across a range of NER rates and read the comet signal.  If the
    comet's bulky-adduct signal really is a repair intermediate, the signal
    must rise with repair capacity, which is the opposite of the intuition
    that repair protects.
    """
    print("\nExperiment: comet signal vs excision-repair rate "
          "(bulky compound, 1 uM)")
    print(f"  {'k_ner':>8} {'bulky left':>11} {'SSB':>8} {'%tail':>8}")
    rows = []
    for k in (0.0, 0.005, 0.02, 0.08, 0.20):
        core = CytogeneticCore(k_ner=k)
        a = build_comet(core)
        w = a.well(demo("direct-acting bulky (4NQO-like)"), 1.0)
        tr = a.trajectory(demo("direct-acting bulky (4NQO-like)"), 1.0)
        rows.append({"k_ner": k, "pct_tail": w["pct_tail_dna"],
                     "ssb": float(tr["ssb"][-1]),
                     "bulky": float(tr["adduct_bulky"][-1])})
        print(f"  {k:8.3f} {rows[-1]['bulky']:11.3f} {rows[-1]['ssb']:8.3f} "
              f"{rows[-1]['pct_tail']:8.2f}")
    mono = all(b["pct_tail"] >= a["pct_tail"] - 1e-9
               for a, b in zip(rows, rows[1:]))
    print(f"  -> %tail rises with repair capacity: {mono}.  A repair-"
          "proficient cell\n     looks MORE damaged in the comet than a "
          "repair-deficient one, because\n     the assay is scoring the "
          "incisions, not the adducts.")
    return {"rows": rows, "monotone": mono}


def turnover_experiment() -> dict:
    """Why does the micronucleus dose-response turn over?

    A micronucleus needs a mitosis, so anything that stops division
    suppresses the endpoint.  Two candidates do that here: p53 -> p21 arrest,
    and outright lethality.  Disabling each in turn attributes the turnover
    instead of assuming it -- the first version of this asserted p53 arrest
    was the cause and the decomposition said otherwise.
    """
    comp = demo("direct-acting bulky (4NQO-like)")
    lo, hi = RANGES[comp.name]
    doses = log_doses(lo, hi, 11)
    base = CytogeneticCore()
    variants = {
        "normal": CytogeneticCore(),
        "no p53 arrest": CytogeneticCore(
            p53=dc_replace(base.p53, K_arrest=1e6)),
        "no lethality": CytogeneticCore(K_tox=1e9),
    }
    series = {k: [r for r in dose_series(build_mn(v), comp, doses,
                                         protocol=MN_CBMN)["rows"]
                  if r["dose_uM"] > 0]
              for k, v in variants.items()}

    print("\nExperiment: what suppresses the micronucleus endpoint at high dose?")
    print(f"  {'dose/uM':>10}" + "".join(f"{k:>16}" for k in variants))
    for i in range(len(series["normal"])):
        print(f"  {series['normal'][i]['dose_uM']:10.3f}"
              + "".join(f"{series[k][i]['signal']:16.1f}" for k in variants))

    def peak(rows):
        t = max(rows, key=lambda r: r["signal"])
        return t["signal"], t["dose_uM"], rows[-1]["signal"]

    out = {k: peak(v) for k, v in series.items()}
    for k, (pk, dose, last) in out.items():
        print(f"  {k:>16}: peak {pk:7.1f} at {dose:6.2f} uM, "
              f"top dose {last:7.1f}")

    pk_n, _, last_n = out["normal"]
    drop_n = 1.0 - last_n / pk_n
    drop_a = 1.0 - out["no p53 arrest"][2] / out["no p53 arrest"][0]
    drop_l = 1.0 - out["no lethality"][2] / out["no lethality"][0]
    print(f"  fall from peak to top dose: {100*drop_n:.0f}% normally, "
          f"{100*drop_a:.0f}% without p53 arrest, "
          f"{100*drop_l:.0f}% without lethality")
    print("  -> the turnover is a CYTOTOXICITY artefact, not a p53 effect: "
          "removing\n     p53 arrest barely changes it, removing lethality "
          "abolishes it.  Reading\n     the high-dose decline as 'less "
          "genotoxic' would be exactly backwards,\n     which is what the "
          "protocol's cytostasis limit exists to prevent.")
    return {"series": series, "peaks": out, "drop_normal": drop_n,
            "drop_no_arrest": drop_a, "drop_no_lethality": drop_l}


def checks(mech, repair, turnover) -> list:
    out = []

    def chk(label, ok, detail):
        out.append((label, bool(ok), detail))
        print(f"  [{'PASS' if ok else 'FAIL'}] {label}: {detail}")

    print("\nStructural checks")

    an = mech["aneugen (colchicine-like)"]
    chk("aneugen is comet-negative at every dose (nothing is broken)",
        an["comet"] == "NEGATIVE",
        f"comet verdict = {an['comet']}")

    chk("aneugen's micronuclei are centromere-POSITIVE",
        an["mn"] == "POSITIVE" and an["pct_cpos"] >= 50,
        f"MN {an['mn']}, {an['pct_cpos']:.1f}% centromere-positive")

    cl = mech["direct clastogen (bleomycin-like)"]
    chk("clastogen's micronuclei are centromere-NEGATIVE",
        cl["mn"] == "POSITIVE" and cl["pct_cpos"] < 10,
        f"MN {cl['mn']}, {cl['pct_cpos']:.1f}% centromere-positive")

    chk("the two mechanisms are told apart, which endpoint #2 could not",
        an["mechanism"].startswith("ANEUGENIC")
        and cl["mechanism"].startswith("CLASTOGENIC"),
        f"aneugen -> {an['mechanism']}; clastogen -> {cl['mechanism']}")

    mx = mech["mixed clastogen/aneugen"]
    chk("a two-mechanism compound is reported as MIXED, not forced to one",
        mx["aneugenic"] and mx["clastogenic"]
        and not (an["clastogenic"] or cl["aneugenic"]),
        f"mixed -> {mx['mechanism']}; the pure aneugen and pure clastogen "
        f"each still get exactly one component")

    nt = mech["non-genotoxic cytotoxicant"]
    chk("pure cytotoxicant is negative in both endpoints",
        nt["comet"] != "POSITIVE" and nt["mn"] != "POSITIVE",
        f"comet {nt['comet']}, MN {nt['mn']}")

    chk("comet scores excision intermediates, not adducts",
        repair["monotone"] and repair["rows"][0]["pct_tail"]
        < repair["rows"][-1]["pct_tail"],
        f"%tail {repair['rows'][0]['pct_tail']:.2f} at k_ner=0 -> "
        f"{repair['rows'][-1]['pct_tail']:.2f} at k_ner=0.20")

    # Assert the attribution, not just the shape: the turnover must survive
    # removing p53 arrest and must vanish when lethality is removed.
    chk("MN turnover is a cytotoxicity artefact, not a p53-arrest effect",
        turnover["drop_normal"] > 0.5
        and turnover["drop_no_arrest"] > 0.5
        and turnover["drop_no_lethality"] < 0.1,
        f"fall from peak: {100*turnover['drop_normal']:.0f}% normally, "
        f"{100*turnover['drop_no_arrest']:.0f}% without p53 arrest, "
        f"{100*turnover['drop_no_lethality']:.0f}% without lethality")

    return out


def figure(comet, mn, mech, repair, turnover):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"\n(figure skipped: {type(e).__name__}: {e})")
        return None

    fig, ax = plt.subplots(1, 4, figsize=(19, 4.3))

    a = ax[0]
    tr = comet.trajectory(demo("direct-acting bulky (4NQO-like)"), 1.0)
    a.plot(tr["t"], tr["adduct_bulky"], label="bulky adduct")
    a.plot(tr["t"], tr["ssb"], label="SSB (excision gap)")
    a.plot(tr["t"], tr["dsb"] * 20, label="DSB x20")
    a.set_xlabel("time / min"); a.set_ylabel("lesions / cell")
    a.set_title("damage processing, 4 h"); a.legend(fontsize=8)

    a = ax[1]
    for comp in PANEL:
        lo, hi = RANGES[comp.name]
        ser = dose_series(comet, comp, log_doses(lo, hi, 11), protocol=COMET)
        rows = ser["rows"][1:]
        a.plot([r["dose_uM"] for r in rows], [r["signal"] for r in rows],
               "-o", ms=3, label=comp.name.split(" (")[0])
    a.set_xscale("log"); a.set_xlabel("dose / uM"); a.set_ylabel("% tail DNA")
    a.set_title("comet (4 h)"); a.legend(fontsize=7)

    a = ax[2]
    for comp in PANEL:
        rows = mech[comp.name]["mn_series"]["rows"][1:]
        a.plot([r["dose_uM"] for r in rows], [r["signal"] for r in rows],
               "-o", ms=3, label=comp.name.split(" (")[0])
    a.set_xscale("log"); a.set_xlabel("dose / uM")
    a.set_ylabel("MN per 1000 BN cells")
    a.set_title("micronucleus (36 h)"); a.legend(fontsize=7)

    a = ax[3]
    for label, style in (("normal", "-o"), ("no p53 arrest", "--s"),
                         ("no lethality", ":^")):
        rows = turnover["series"][label]
        a.plot([r["dose_uM"] for r in rows], [r["signal"] for r in rows],
               style, ms=4, label=label)
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel("dose / uM")
    a.set_ylabel("MN per 1000 BN cells")
    a.set_title("turnover = cytotoxicity, not p53"); a.legend(fontsize=8)

    fig.suptitle("Virtual comet + micronucleus — one cytogenetic core, two "
                 "instruments (illustrative parameters, not fitted)",
                 fontsize=11)
    fig.tight_layout()
    os.makedirs(FIGDIR, exist_ok=True)
    path = os.path.join(FIGDIR, "genotox_comet_mn.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    print(f"\nfigure -> {path}")
    return path


def report() -> dict:
    print("=" * 70)
    print("GENOTOX STEP 3 — virtual comet + micronucleus (one core)")
    print("=" * 70)
    comet, mn = build_comet(), build_mn()
    print(f"comet : {comet.description}  [{COMET_MIN/60:.0f} h]")
    print(f"MN    : {mn.description}  [{MN_MIN/60:.0f} h]")
    print("Parameters are illustrative (H); absolute %tail and MN "
          "frequencies are not predictions.\n")

    processing_report(comet)
    series_report(comet, demo("direct-acting bulky (4NQO-like)"), COMET)
    series_report(mn, demo("direct-acting bulky (4NQO-like)"), MN_CBMN,
                  extra=(("mn_centromere_pos_pct", "%C+"), ("CBPI", "CBPI")))
    series_report(mn, demo("aneugen (colchicine-like)"), MN_CBMN,
                  extra=(("mn_centromere_pos_pct", "%C+"), ("CBPI", "CBPI")))

    mech = mechanism_table(comet, mn)
    repair = repair_experiment()
    turnover = turnover_experiment()
    ck = checks(mech, repair, turnover)
    path = figure(comet, mn, mech, repair, turnover)

    n_ok = sum(1 for _, ok, _ in ck if ok)
    print(f"\n{n_ok}/{len(ck)} structural checks passed.")
    return {"mechanism": mech, "repair": repair, "turnover": turnover,
            "checks": ck, "figure": path}


if __name__ == "__main__":
    report()
