"""
Benchmark runner — score a structure-reading upstream on reference chemicals.

Validates every benchmark structure, shows which alerts fire on what, then
scores the structural-alert source on three endpoints and reports where it
fails and why.

    python3 -m genotox.run_benchmark
"""
from __future__ import annotations

import numpy as np

from .alerts import StructuralAlertSource
from .assay import VirtualAssay
from .benchmark import REFERENCES, score_source, validate
from .core import SOSCore
from .cytogenetic import CytogeneticCore
from .damage import Compound, DamageFlux
from .doseresponse import (COMET, MN_CBMN, UMU, call_result, dose_series,
                           log_doses)

DOSES = log_doses(0.1, 300.0, 9)


def umu_assay(source):
    return VirtualAssay(source=source, core=SOSCore())


def comet_assay(source):
    return VirtualAssay(source=source, core=CytogeneticCore(),
                        readouts=("comet_alkaline", "growth", "damage_probe"),
                        duration_min=4 * 60, n_points=49)


def mn_assay(source):
    return VirtualAssay(source=source, core=CytogeneticCore(),
                        readouts=("micronucleus_cbmn", "growth",
                                  "damage_probe"),
                        duration_min=36 * 60, n_points=73)


def structures_section() -> dict:
    v = validate()
    print("Benchmark structures")
    print(f"  {v['n']} reference compounds; formulas checked against RDKit")
    if v["problems"]:
        for name, why in v["problems"]:
            print(f"  PROBLEM  {name}: {why}")
    else:
        print("  every structure parses and matches its stated formula")
    return v


def alerts_section(source) -> dict:
    print("\nWhich rules fire on what")
    print(f"  {'compound':28s} {'role':34s} rules fired")
    out = {}
    for r in REFERENCES:
        fired = [n for n, _ in source.explain(r.smiles)["fired"]]
        out[r.name] = fired
        print(f"  {r.name:28s} {r.role[:33]:34s} "
              f"{', '.join(fired) if fired else '(none)'}")
    return out


def channel_section(source) -> None:
    print("\nChannels emitted, versus the channels a correct source should hit")
    print(f"  {'compound':28s} {'emitted':46s} expected")
    for r in REFERENCES:
        f = source.flux_from_smiles(r.smiles)
        got = ", ".join(f"{k}:{v:.2f}" for k, v in f.nonzero().items()) or "-"
        want = ", ".join(r.channels) or "(none)"
        flag = "" if set(f.nonzero()) or not r.channels else "   <- missed"
        print(f"  {r.name:28s} {got[:45]:46s} {want}{flag}")


def score_section(source) -> dict:
    print("\nScores")
    print(f"  {'endpoint':10s} {'n':>3} {'TP':>3} {'FP':>3} {'TN':>3} "
          f"{'FN':>3} {'sens':>6} {'spec':>6} {'acc':>6}  skipped  "
          f"artifact-class declined")
    out = {}
    for label, build, proto in (("umu", umu_assay, UMU),
                                ("comet", comet_assay, COMET),
                                ("mn", mn_assay, MN_CBMN)):
        sc = score_source(source, build, proto, label, DOSES)
        out[label] = sc
        print(f"  {label:10s} {sc.n:3d} {sc.tp:3d} {sc.fp:3d} {sc.tn:3d} "
              f"{sc.fn:3d} {sc.sensitivity:6.2f} {sc.specificity:6.2f} "
              f"{sc.accuracy:6.2f}  {sc.skipped:7d}  "
              f"{sc.artifact_declined}/{sc.artifact_total}")
    print("  'artifact-class declined' counts in-vitro positives attributed "
          "to cytotoxicity\n  where the source correctly fired no damage "
          "channel.  Reproducing those\n  positives needs a cytotoxicity "
          "prediction, which no source here makes.")
    return out


def failure_section(source, scores) -> dict:
    """Name every miss and say which layer is responsible.

    A score with no error analysis is a number; the point of a benchmark is
    which kind of compound it fails on.
    """
    print("\nEvery miss, and where it belongs")
    misses = []
    for label, build, proto in (("umu", umu_assay, UMU),
                                ("comet", comet_assay, COMET),
                                ("mn", mn_assay, MN_CBMN)):
        for r in REFERENCES:
            want = r.expected.get(label, "?")
            if want == "?" or (r.artifact and want == "+"):
                continue
            comp = Compound(name=r.name, per_uM=DamageFlux(),
                            smiles=r.smiles,
                            direct_fraction=0.02 if r.needs_s9 else 1.0,
                            s9_fraction=0.95 if r.needs_s9 else 0.0)
            res = call_result(dose_series(build(source), comp, DOSES,
                                          s9=r.needs_s9, protocol=proto))
            got = res["verdict"]
            hit = (got == "POSITIVE") == (want == "+")
            if not hit:
                fired = [n for n, _ in source.explain(r.smiles)["fired"]]
                misses.append({"endpoint": label, "compound": r.name,
                               "want": want, "got": got, "fired": fired})
    for m in misses:
        why = ("no rule fired — alerts describe reactivity, and this is a "
               "target-binding property" if not m["fired"]
               else f"fired {', '.join(m['fired'])}")
        print(f"  {m['endpoint']:6s} {m['compound']:26s} want {m['want']} "
              f"got {m['got']:12s} {why}")
    return {"misses": misses}


def checks(v, fired, scores, failures) -> list:
    out = []

    def chk(label, ok, detail):
        out.append((label, bool(ok), detail))
        print(f"  [{'PASS' if ok else 'FAIL'}] {label}: {detail}")

    print("\nStructural checks")

    chk("every benchmark structure parses and matches its formula",
        not v["problems"], f"{v['n']} compounds, {len(v['problems'])} problems")

    chk("the species-split topo channels keep a gyrase poison out of the "
        "mammalian endpoints",
        "quinolone-3-carboxylate" in fired["ciprofloxacin"]
        and all(m["compound"] != "ciprofloxacin" or m["endpoint"] == "umu"
                for m in failures["misses"]),
        "ciprofloxacin fires topo_bacterial only; comet and MN stay negative")

    chk("the aneugen chemotype routes to the aneugenic channel",
        "benzimidazole carbamate" in fired["carbendazim"]
        and "benzimidazole carbamate" in fired["nocodazole"],
        "carbendazim and nocodazole both fire it")

    chk("alerts miss target-binding aneugens, and the run says so",
        not fired["colchicine"]
        and any(m["compound"] == "colchicine" for m in failures["misses"]),
        "colchicine fires nothing and is reported as a miss — a SMARTS list "
        "cannot see tubulin binding")

    neg = [r.name for r in REFERENCES if r.role.startswith("negative")
           or r.role.startswith("misleading")]
    clean = [n for n in neg if not fired[n]]
    chk("misleading positives and negatives mostly draw no alert",
        len(clean) >= len(neg) - 1,
        f"{len(clean)} of {len(neg)} draw nothing: {', '.join(clean)}")

    chk("specificity is not achieved by predicting nothing",
        scores["umu"].tp > 0 and scores["mn"].tp > 0,
        f"umu TP={scores['umu'].tp}, MN TP={scores['mn'].tp}")

    return out


def report() -> dict:
    print("=" * 74)
    print("GENOTOX BENCHMARK — a structure-reading upstream on reference "
          "chemicals")
    print("=" * 74)
    print("Labels are published reference classifications for control "
          "chemicals, not\nmeasurements made here. These compounds are "
          "famous, so a good score is a\nfloor the upstream must clear, "
          "never a validation. Metabolic activation is\nsupplied from the "
          "reference metadata, not predicted.\n")

    source = StructuralAlertSource()
    print(f"source: {source.provenance}\n")
    v = structures_section()
    fired = alerts_section(source)
    channel_section(source)
    scores = score_section(source)
    failures = failure_section(source, scores)
    ck = checks(v, fired, scores, failures)
    n_ok = sum(1 for _, ok, _ in ck if ok)
    print(f"\n{n_ok}/{len(ck)} structural checks passed.")
    return {"validate": v, "fired": fired, "scores": scores,
            "failures": failures, "checks": ck}


if __name__ == "__main__":
    report()
