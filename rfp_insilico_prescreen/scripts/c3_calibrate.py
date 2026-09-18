#!/usr/bin/env python3
"""
C3 - Does NetMHCpan work on THESE TWO alleles?

A pre-screen that tells a sponsor "do not pay for that tetramer" has to earn it.
Class I prediction is good on average, but the average is carried by well-studied
alleles (A*02:01 and friends). A*32:01 is not one of those. So the benchmark from
C1 is scored and the predictor is measured per allele.

Reports, per allele: ROC AUC, average precision, and - the number that actually
matters for this decision - what a call is worth at a chosen %Rank cut, in both
directions:

  PPV   if the pre-screen says "synthesise", how often is that peptide a real
        binder (i.e. how much tetramer budget is wasted)
  NPV   if the pre-screen says "skip", how often would that have been a real
        binder thrown away (i.e. how much science is lost)

For a go/no-go on an expensive reagent, NPV at the skip threshold is the
governing quantity, because a missed binder is unrecoverable while a wasted
synthesis is only money.

Output: results/c3_calibration.json, results/c3_benchmark_scored.tsv
"""
import csv
import json
import math
import os
import random
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import results_path, iedb_mhci  # noqa: E402

ENDPOINT = "https://tools-cluster-interface.iedb.org/tools_api/mhci/"
EL, BA = "netmhcpan_el", "netmhcpan_ba"
CHUNK = 900              # peptides per request
CUTS = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0]   # wide cuts included on purpose: skipping is
                                          # the irreversible move, so the rule needs room
                                          # to be permissive before it gives up
BOOT = 2000
SEED = 20260918
HELD_OUT_FRAC = 0.2      # reserved for C4, so the pre-screen is not scored on its own tuning set
MIN_SENSITIVITY = 0.95   # a skipped complex is unrecoverable; keep 95% of real binders in BUILD

# IEDB class I data is overwhelmingly eluted ligands, which are positive by
# construction: A*32:01 comes back ~96% positive. PPV and NPV computed at that
# prevalence describe IEDB's collection policy, not this decision. They are
# therefore recomputed at a stated prior for "a designed peptide proposed for a
# chosen allele actually binds it", and swept so the reader can move it.
ASSUMED_PREVALENCE = 0.30
PREVALENCE_SWEEP = [0.1, 0.2, 0.3, 0.5]


def tsv(name):
    with open(results_path(name)) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def auc(scores, labels):
    """Mann-Whitney with mid-ranks. Higher score = more likely positive."""
    pairs = sorted(zip(scores, labels))
    ranks, i, n = [0.0] * len(pairs), 0, len(pairs)
    while i < n:
        j = i
        while j + 1 < n and pairs[j + 1][0] == pairs[i][0]:
            j += 1
        r = (i + j) / 2.0 + 1
        for k in range(i, j + 1):
            ranks[k] = r
        i = j + 1
    P = sum(l for _, l in pairs)
    N = n - P
    if P == 0 or N == 0:
        return None
    rp = sum(r for r, (_, l) in zip(ranks, pairs) if l == 1)
    return (rp - P * (P + 1) / 2.0) / (P * N)


def average_precision(scores, labels):
    order = sorted(range(len(scores)), key=lambda i: -scores[i])
    tp = 0
    total_pos = sum(labels)
    if not total_pos:
        return None
    s = 0.0
    for rank, i in enumerate(order, 1):
        if labels[i] == 1:
            tp += 1
            s += tp / rank
    return s / total_pos


def predictive_values(sens, spec, prev):
    """PPV/NPV at a stated prior, instead of the benchmark's own prevalence."""
    tp = sens * prev
    fn = (1 - sens) * prev
    tn = spec * (1 - prev)
    fp = (1 - spec) * (1 - prev)
    ppv = tp / (tp + fp) if (tp + fp) else 0.0
    npv = tn / (tn + fn) if (tn + fn) else 0.0
    return round(ppv, 4), round(npv, 4)


def operating_point(ranks, labels, cut):
    tp = sum(1 for r, l in zip(ranks, labels) if r < cut and l == 1)
    fp = sum(1 for r, l in zip(ranks, labels) if r < cut and l == 0)
    fn = sum(1 for r, l in zip(ranks, labels) if r >= cut and l == 1)
    tn = sum(1 for r, l in zip(ranks, labels) if r >= cut and l == 0)
    d = lambda a, b: a / b if b else 0.0            # noqa: E731
    sens, spec = d(tp, tp + fn), d(tn, tn + fp)
    ppv, npv = predictive_values(sens, spec, ASSUMED_PREVALENCE)
    return {"cut": cut, "tp": tp, "fp": fp, "fn": fn, "tn": tn,
            "sensitivity": round(sens, 4), "specificity": round(spec, 4),
            "ppv_synthesise_at_assumed_prev": ppv,
            "npv_skip_at_assumed_prev": npv,
            "ppv_npv_by_prevalence": {
                str(p): list(predictive_values(sens, spec, p))
                for p in PREVALENCE_SWEEP},
            "n_flagged_in_benchmark": tp + fp}


def boot_auc_ci(scores, labels, n=BOOT):
    rng = random.Random(SEED)
    idx = list(range(len(labels)))
    vals = []
    for _ in range(n):
        s = [rng.choice(idx) for _ in idx]
        a = auc([scores[i] for i in s], [labels[i] for i in s])
        if a is not None:
            vals.append(a)
    vals.sort()
    if not vals:
        return None
    lo = vals[int(0.025 * len(vals))]
    hi = vals[int(0.975 * len(vals)) - 1]
    return [round(lo, 4), round(hi, 4)]


def score(peptides, allele):
    """%Rank per peptide for both heads. IEDB takes a peptide list directly for
    class I, so no concatenation trick is needed here - lengths are fixed."""
    out = {}
    peps = sorted(peptides)
    for i in range(0, len(peps), CHUNK):
        block = peps[i:i + CHUNK]
        for method, key in ((EL, "el_rank"), (BA, "ba_rank")):
            rows = iedb_mhci(ENDPOINT, method, block, [allele])
            for r in rows:
                p = r["peptide"]
                col = "percentile_rank" if "percentile_rank" in r else "rank"
                out.setdefault(p, {})[key] = float(r[col])
        print(f"    {allele}  {min(i+CHUNK, len(peps))}/{len(peps)}", flush=True)
    return out


def main():
    bench = tsv("c1_benchmark.tsv")
    by_allele = defaultdict(list)
    for r in bench:
        by_allele[r["allele"]].append((r["peptide"], int(r["label"])))

    cache = results_path("c3_benchmark_scored.tsv")
    done = {}
    if os.path.exists(cache):
        for r in tsv("c3_benchmark_scored.tsv"):
            done[(r["allele"], r["peptide"])] = r
        print(f"resuming: {len(done)} peptides already scored")

    rows = []
    for allele, items in by_allele.items():
        need = [p for p, _ in items if (allele, p) not in done]
        print(f"{allele}: {len(items)} peptides, {len(need)} to score")
        got = score(need, allele) if need else {}
        for pep, lab in items:
            if (allele, pep) in done:
                rows.append(done[(allele, pep)])
                continue
            g = got.get(pep)
            if not g:
                continue
            rows.append({"allele": allele, "peptide": pep, "label": lab,
                         "el_rank": g.get("el_rank"), "ba_rank": g.get("ba_rank")})
        with open(cache, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["allele", "peptide", "label",
                                              "el_rank", "ba_rank"],
                               delimiter="\t")
            w.writeheader()
            w.writerows(rows)

    # Hold out a slice per allele so C4's pre-screen is never scored on the
    # peptides that set its own threshold.
    rng = random.Random(SEED)
    held_out = []
    out = {"cuts": CUTS, "assumed_prevalence": ASSUMED_PREVALENCE,
           "prevalence_sweep": PREVALENCE_SWEEP,
           "min_sensitivity_for_skip_cut": MIN_SENSITIVITY,
           "held_out_frac": HELD_OUT_FRAC, "per_allele": {}}
    for allele in sorted(by_allele):
        sub = [r for r in rows if r["allele"] == allele and r["el_rank"] not in (None, "")]
        rng.shuffle(sub)
        n_hold = int(len(sub) * HELD_OUT_FRAC)
        hold, cal = sub[:n_hold], sub[n_hold:]
        held_out.extend([allele, r["peptide"]] for r in hold)

        labels = [int(r["label"]) for r in cal]
        el = [-float(r["el_rank"]) for r in cal]
        ranks = [float(r["el_rank"]) for r in cal]
        ba_pairs = [(-float(r["ba_rank"]), int(r["label"])) for r in cal
                    if r["ba_rank"] not in (None, "")]
        ops = [operating_point(ranks, labels, c) for c in CUTS]

        # The skip threshold. Sensitivity rises monotonically with the cut, so
        # "widest cut meeting a sensitivity floor" degenerates to "skip nothing"
        # and is worthless. The pre-screen only earns its keep by skipping, so
        # take the NARROWEST cut that still meets the floor: it skips the most
        # while losing no more than (1 - floor) of real binders.
        ok = [o for o in ops if o["sensitivity"] >= MIN_SENSITIVITY]
        met = bool(ok)
        chosen = min(ok, key=lambda o: o["cut"]) if met else max(
            ops, key=lambda o: o["sensitivity"])
        rec = {
            "n_calibration": len(cal), "n_held_out": len(hold),
            "n_positive": sum(labels), "n_negative": len(labels) - sum(labels),
            "positive_fraction": round(sum(labels) / len(labels), 4) if labels else None,
            "EL": {"auc": round(auc(el, labels) or 0, 4),
                   "auc_ci95": boot_auc_ci(el, labels),
                   "average_precision": round(average_precision(el, labels) or 0, 4)},
            "operating_points": ops,
            "recommended_skip_cut": chosen["cut"],
            "recommended_skip_npv": chosen["npv_skip_at_assumed_prev"],
            "recommended_skip_sensitivity": chosen["sensitivity"],
            "skip_cut_rule": (f"narrowest cut in {CUTS} retaining sensitivity "
                              f">= {MIN_SENSITIVITY}"),
            # What the rule is worth in practice: the share of a candidate pool
            # it would send to SKIP, at the stated prior.
            "expected_skip_fraction_at_assumed_prev": round(
                (1 - chosen["sensitivity"]) * ASSUMED_PREVALENCE
                + chosen["specificity"] * (1 - ASSUMED_PREVALENCE), 4),
            "sensitivity_floor_met": met,
            "skip_cut_warning": (None if met else
                                 f"no cut tested reaches sensitivity "
                                 f"{MIN_SENSITIVITY}; the widest available keeps "
                                 f"{chosen['sensitivity']:.2f}, so a skip on this "
                                 f"allele discards roughly "
                                 f"{(1-chosen['sensitivity'])*100:.0f}% of real binders"),
        }
        if ba_pairs:
            rec["BA"] = {"auc": round(auc([x for x, _ in ba_pairs],
                                          [y for _, y in ba_pairs]) or 0, 4)}
        out["per_allele"][allele] = rec
    out["held_out"] = held_out

    p = results_path("c3_calibration.json")
    with open(p, "w") as f:
        json.dump(out, f, indent=1)

    for allele, rec in out["per_allele"].items():
        ci = rec["EL"]["auc_ci95"]
        print(f"\n{allele}  calibration n={rec['n_calibration']} "
              f"({rec['n_positive']} pos / {rec['n_negative']} neg, "
              f"{rec['positive_fraction']*100:.0f}% positive)  held out {rec['n_held_out']}")
        print(f"  EL  ROC AUC {rec['EL']['auc']:.3f}"
              + (f"  95% CI [{ci[0]:.3f}, {ci[1]:.3f}]" if ci else ""))
        print(f"  {'%Rank cut':>10s} {'sens':>6s} {'spec':>6s} "
              f"{'PPV@30%':>8s} {'NPV@30%':>8s}")
        for op in rec["operating_points"]:
            print(f"  {op['cut']:>10.1f} {op['sensitivity']:6.3f} {op['specificity']:6.3f} "
                  f"{op['ppv_synthesise_at_assumed_prev']:8.3f} "
                  f"{op['npv_skip_at_assumed_prev']:8.3f}")
        print(f"  -> skip anything with %Rank >= {rec['recommended_skip_cut']}  "
              f"(NPV {rec['recommended_skip_npv']:.2f} at {ASSUMED_PREVALENCE:.0%} prior, "
              f"sensitivity {rec['recommended_skip_sensitivity']:.2f}, "
              f"skips {rec['expected_skip_fraction_at_assumed_prev']:.0%} of a candidate pool)")
        if rec.get("skip_cut_warning"):
            print(f"     WARNING: {rec['skip_cut_warning']}")

    print(f"\nwrote {p}")


if __name__ == "__main__":
    main()
