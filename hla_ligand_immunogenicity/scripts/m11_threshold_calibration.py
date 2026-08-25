#!/usr/bin/env python3
"""
M11 - Measure the decision rule against ground truth, and pick a better one.

The pipeline's scoring rule was chosen by argument: eluted-ligand scoring
over-calls, so require the binding-affinity head to agree. This module tests
that argument on the M10 benchmark of measured HLA-DR-restricted CD4 T-cell
responses, and reports what it actually buys.

Rules compared, all on the identical set of (peptide, allele) pairs:

  EL          eluted-ligand %Rank alone            (the industry default)
  BA          binding-affinity %Rank alone
  GEO         sqrt(EL x BA) - a continuous consensus, so the "does adding BA
              help" question has an answer that does not depend on a threshold
  MIN         min(EL, BA) - the permissive combination
  AND         EL < 1 and BA < 10 - the rule currently in the pipeline
              (a single operating point, so it gets sensitivity/specificity,
              not an AUC)

Statistics that keep the answer honest:

  * ROC AUC and average precision for every continuous rule.
  * A **paired bootstrap over 9-mer clusters**, not over pairs. Overlapping
    peptides from one study are one observation, and resampling pairs would
    shrink the confidence interval to whatever width the redundancy dictates.
    The reported interval is on the *difference* between two rules on the same
    resampled clusters, which is the only form that answers "is this rule
    better than that one".
  * Metrics stratified by self vs non-self source, because the label sets are
    confounded: if a rule's advantage disappears within stratum, it was
    detecting foreignness, not binding.
  * A prevalence-corrected PPV. The benchmark is roughly balanced; a real
    ligand scan is not, so raw PPV from a balanced set overstates what a flag
    means in practice.

One caveat governs how the absolute numbers may be read. NetMHCIIpan is
trained on IEDB binding-affinity and mass-spec eluted-ligand data. The labels
here are *T-cell assay outcomes*, a different endpoint the predictor was not
trained on - but many of these peptides also carry binding measurements in
IEDB, so partial training-set overlap is certain and cannot be excluded from
the outside. Treat every absolute AUC as an optimistic bound.

The *relative* comparison is far more robust: every rule uses the same two
possibly-leaky predictors on the same peptides, so leakage inflates them all
together and largely cancels in the difference. "Is the consensus rule better
than EL alone" is answerable here; "NetMHCIIpan achieves AUC X on epitope
prediction" is not.

Output: results/m11_calibration.json, results/m11_benchmark_scored.tsv
"""
import csv
import json
import math
import os
import random
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import load_config, results_path, iedb_mhcii  # noqa: E402

MAX_PER_ALLELE = 400      # keeps the IEDB load reasonable; sampling is seeded
SEED = 20260825
BOOTSTRAP = 800        # enough for a 95% interval to two decimals; the
                       # resample is O(n log n) per draw, so more is mostly heat
# Prevalence of true DR epitopes among 15-mers of an arbitrary protein. Used
# only to translate benchmark PPV into something meaningful for a real scan;
# 5% is a deliberately generous round number, and the effect of changing it is
# reported rather than hidden.
SCAN_PREVALENCE = 0.05


# --------------------------------------------------------------------- stats
def roc_auc(scores, labels):
    """Mann-Whitney U formulation, with ties handled by mid-ranks."""
    order = sorted(range(len(scores)), key=lambda i: scores[i])
    ranks = [0.0] * len(scores)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and scores[order[j + 1]] == scores[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    n1 = sum(labels)
    n0 = len(labels) - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    s1 = sum(r for r, l in zip(ranks, labels) if l == 1)
    return (s1 - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def average_precision(scores, labels):
    pairs = sorted(zip(scores, labels), key=lambda t: -t[0])
    tp = fp = 0
    total_pos = sum(labels)
    if not total_pos:
        return float("nan")
    ap, prev_recall = 0.0, 0.0
    for s, l in pairs:
        tp += l
        fp += 1 - l
        recall = tp / total_pos
        precision = tp / (tp + fp)
        ap += precision * (recall - prev_recall)
        prev_recall = recall
    return ap


def confusion(pred, labels):
    tp = sum(1 for p, l in zip(pred, labels) if p and l)
    fp = sum(1 for p, l in zip(pred, labels) if p and not l)
    fn = sum(1 for p, l in zip(pred, labels) if not p and l)
    tn = sum(1 for p, l in zip(pred, labels) if not p and not l)
    sens = tp / (tp + fn) if tp + fn else float("nan")
    spec = tn / (tn + fp) if tn + fp else float("nan")
    ppv = tp / (tp + fp) if tp + fp else float("nan")
    npv = tn / (tn + fn) if tn + fn else float("nan")
    den = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = ((tp * tn - fp * fn) / den) if den else float("nan")
    return {"tp": tp, "fp": fp, "fn": fn, "tn": tn,
            "sensitivity": round(sens, 4), "specificity": round(spec, 4),
            "ppv": round(ppv, 4), "npv": round(npv, 4), "mcc": round(mcc, 4)}


def ppv_at_prevalence(sens, spec, prev):
    if any(math.isnan(x) for x in (sens, spec)):
        return float("nan")
    num = sens * prev
    den = num + (1 - spec) * (1 - prev)
    return num / den if den else float("nan")


def cluster_bootstrap_auc_delta(rows, key_a, key_b, n=BOOTSTRAP, seed=SEED):
    """Paired bootstrap over clusters on AUC(a) - AUC(b)."""
    by_cluster = defaultdict(list)
    for r in rows:
        by_cluster[r["cluster"]].append(r)
    clusters = list(by_cluster)
    rng = random.Random(seed)
    deltas = []
    for _ in range(n):
        sample = []
        for _ in clusters:
            sample.extend(by_cluster[clusters[rng.randrange(len(clusters))]])
        labels = [r["label"] for r in sample]
        if not (0 < sum(labels) < len(labels)):
            continue
        a = roc_auc([r[key_a] for r in sample], labels)
        b = roc_auc([r[key_b] for r in sample], labels)
        if not (math.isnan(a) or math.isnan(b)):
            deltas.append(a - b)
    deltas.sort()
    if not deltas:
        return None
    lo = deltas[int(0.025 * len(deltas))]
    hi = deltas[min(int(0.975 * len(deltas)), len(deltas) - 1)]
    point = roc_auc([r[key_a] for r in rows], [r["label"] for r in rows]) - \
            roc_auc([r[key_b] for r in rows], [r["label"] for r in rows])
    return {"delta": round(point, 4), "ci95": [round(lo, 4), round(hi, 4)],
            "p_gt_0": round(sum(1 for d in deltas if d > 0) / len(deltas), 4)}


# ------------------------------------------------------------------- scoring
SPACER_AA = "G"
CHUNK_AA = 3000           # max residues per submitted record


def score_benchmark(cfg, bench, on_result=None):
    """
    Best EL and BA %Rank per (peptide, allele) from the IEDB predictors.

    Peptides are concatenated into a few long pseudo-proteins rather than
    submitted one FASTA record each, because the wall-clock cost of this
    endpoint tracks the number of *requests* far more than the amount of
    sequence in them: under load a single 60-residue request took 45 s, and
    20 peptides submitted as 20 records took 89 s. Whether that is per-record
    processing or queue latency was not isolated and does not matter here -
    either way, fewer and larger requests is the fix, and it turns a
    five-hour job into minutes.

    The concatenation is exact, not an approximation: NetMHCIIpan scores each
    k-mer independently of its surroundings when no context flag is set, so a
    frame lying wholly inside one peptide gets the score it would get alone.
    Peptides are separated by a spacer one residue shorter than the scan
    length, so every frame either lies wholly inside a peptide or touches the
    spacer, and spacer-touching frames are discarded by position.
    `verify_concatenation` checks this equivalence against standalone scoring
    before the run trusts it.

    `on_result` is called with each completed job's results so the caller can
    persist them as they arrive; at ~45 s per request a crash near the end of a
    two-hour run should not cost the whole run.
    """
    pcfg = cfg["prediction"]
    groups = defaultdict(set)          # (allele, scan_length) -> {peptide}
    for r in bench:
        groups[(r["allele"], min(15, len(r["peptide"])))].add(r["peptide"])

    jobs = []
    for (allele, length), peps in groups.items():
        spacer = SPACER_AA * (length - 1)
        cur, layout, pos = [], [], 1
        for pep in sorted(peps):
            if cur and pos + len(pep) - 1 > CHUNK_AA:
                jobs.append((allele, length, "".join(cur), layout))
                cur, layout, pos = [], [], 1
            if cur:
                cur.append(spacer)
                pos += len(spacer)
            cur.append(pep)
            layout.append((pep, pos, pos + len(pep) - 1))
            pos += len(pep)
        if cur:
            jobs.append((allele, length, "".join(cur), layout))
    print(f"{len(bench)} pairs -> {len(jobs)} concatenated IEDB calls", flush=True)

    out = defaultdict(dict)

    def run(job):
        allele, length, seq, layout = job
        res = {}
        for method, field in (("netmhciipan_el", "el"), ("netmhciipan_ba", "ba")):
            rows = iedb_mhcii(pcfg["endpoint"], method, {"cat": seq}, [allele], length)
            for r in rows:
                s0, e0 = int(r["start"]), int(r["end"])
                for pep, ps, pe in layout:
                    if ps <= s0 and e0 <= pe:          # frame wholly inside one peptide
                        key = (pep, allele)
                        v = float(r["rank"])
                        cur = res.setdefault(key, {}).get(field)
                        if cur is None or v < cur:
                            res[key][field] = v
                        break
        print(f"  ok {allele:18s} len{length} {len(layout):3d} peptides", flush=True)
        return res

    with ThreadPoolExecutor(max_workers=4) as ex:
        for res in ex.map(run, jobs):
            for k, v in res.items():
                out[k].update(v)
            if on_result:
                on_result(res)
    return out


def verify_concatenation(cfg, peptides, allele):
    """Score a handful of peptides standalone and concatenated; they must agree."""
    pcfg = cfg["prediction"]
    peptides = [p for p in peptides if len(p) >= 15][:8]
    if len(peptides) < 3:
        return True, "too few peptides to verify"
    solo = {}
    rows = iedb_mhcii(pcfg["endpoint"], "netmhciipan_el",
                      {p: p for p in peptides}, [allele], 15)
    for r in rows:
        v = float(r["rank"])
        if v < solo.get(r["id"], 9e9):
            solo[r["id"]] = v
    fake = [{"peptide": p, "allele": allele} for p in peptides]
    cat = score_benchmark(cfg, fake)
    bad = []
    for p in peptides:
        a, b = solo.get(p), cat.get((p, allele), {}).get("el")
        if a is None or b is None or abs(a - b) > 1e-6:
            bad.append((p, a, b))
    return not bad, (f"{len(peptides)} peptides agree exactly" if not bad
                     else f"MISMATCH: {bad[:3]}")


def main():
    cfg = load_config()
    with open(results_path("m10_benchmark.tsv")) as f:
        bench = [dict(r) for r in csv.DictReader(f, delimiter="\t")]
    for r in bench:
        r["label"] = int(r["label"])
        r["is_self"] = bool(int(r["is_self"]))
        r["cluster"] = int(r["cluster"])

    # seeded, per-allele, label-balanced subsample
    rng = random.Random(SEED)
    by_allele = defaultdict(list)
    for r in bench:
        by_allele[r["allele"]].append(r)
    sampled = []
    for allele, rows in sorted(by_allele.items()):
        pos = [r for r in rows if r["label"] == 1]
        neg = [r for r in rows if r["label"] == 0]
        cap = MAX_PER_ALLELE // 2
        rng.shuffle(pos)
        rng.shuffle(neg)
        sampled.extend(pos[:cap] + neg[:cap])
    print(f"benchmark {len(bench)} pairs -> {len(sampled)} sampled "
          f"({sum(r['label'] for r in sampled)} pos)")

    cache = results_path("m11_benchmark_scored.tsv")
    partial = results_path("m11_scores_partial.tsv")
    scores = {}
    for src in (cache, partial):
        if not os.path.exists(src):
            continue
        with open(src) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if r.get("el_rank") and r.get("ba_rank"):
                    scores[(r["peptide"], r["allele"])] = {
                        "el": float(r["el_rank"]), "ba": float(r["ba_rank"])}
        print(f"resuming: {len(scores)} pairs from {os.path.basename(src)}")
    todo = [r for r in sampled if (r["peptide"], r["allele"]) not in scores]
    if todo:
        ok, msg = verify_concatenation(
            cfg, sorted({r["peptide"] for r in todo}), todo[0]["allele"])
        print(f"concatenation check: {msg}", flush=True)
        if not ok:
            sys.exit("concatenated scoring does not reproduce standalone scoring")
        # append each job's results as they land, so an interrupted run resumes
        with open(partial, "a", newline="") as pf:
            pw = csv.writer(pf, delimiter="\t")
            if pf.tell() == 0:
                pw.writerow(["peptide", "allele", "el_rank", "ba_rank"])

            def persist(res):
                for (pep, allele), v in res.items():
                    if "el" in v and "ba" in v:
                        pw.writerow([pep, allele, v["el"], v["ba"]])
                pf.flush()

            scores.update(score_benchmark(cfg, todo, on_result=persist))

    rows = []
    for r in sampled:
        s = scores.get((r["peptide"], r["allele"]))
        if not s or "el" not in s or "ba" not in s:
            continue
        el, ba = max(s["el"], 1e-3), max(s["ba"], 1e-3)
        rows.append({**r, "el_rank": el, "ba_rank": ba,
                     # higher = more epitope-like, for every rule
                     "s_EL": -math.log10(el),
                     "s_BA": -math.log10(ba),
                     "s_GEO": -math.log10(math.sqrt(el * ba)),
                     "s_MIN": -math.log10(min(el, ba))})
    print(f"{len(rows)} pairs scored on both heads")

    with open(cache, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["peptide", "allele", "label", "is_self", "cluster",
                    "el_rank", "ba_rank"])
        for r in rows:
            w.writerow([r["peptide"], r["allele"], r["label"], int(r["is_self"]),
                        r["cluster"], r["el_rank"], r["ba_rank"]])

    labels = [r["label"] for r in rows]
    out = {"n_pairs": len(rows), "n_positive": sum(labels),
           "n_clusters": len({r["cluster"] for r in rows}),
           "rules": {}, "comparisons": {}, "operating_points": {},
           "stratified": {}, "scan_prevalence": SCAN_PREVALENCE}

    for name in ("EL", "BA", "GEO", "MIN"):
        k = "s_" + name
        out["rules"][name] = {
            "auc": round(roc_auc([r[k] for r in rows], labels), 4),
            "average_precision": round(average_precision([r[k] for r in rows], labels), 4),
        }

    for a, b in (("GEO", "EL"), ("BA", "EL"), ("MIN", "EL")):
        out["comparisons"][f"{a}_vs_{b}"] = cluster_bootstrap_auc_delta(
            rows, "s_" + a, "s_" + b)

    # ---- operating points -------------------------------------------------
    sb, wb = cfg["prediction"]["sb_rank"], cfg["prediction"]["wb_rank"]
    bac = cfg["prediction"]["ba_confirm_rank"]
    points = {
        f"AND (current): EL<{sb:g} and BA<{bac:g}":
            [r["el_rank"] < sb and r["ba_rank"] < bac for r in rows],
        f"EL<{sb:g} alone": [r["el_rank"] < sb for r in rows],
        f"EL<{wb:g} alone": [r["el_rank"] < wb for r in rows],
        f"BA<{bac:g} alone": [r["ba_rank"] < bac for r in rows],
    }
    for name, pred in points.items():
        c = confusion(pred, labels)
        c["ppv_at_scan_prevalence"] = round(
            ppv_at_prevalence(c["sensitivity"], c["specificity"], SCAN_PREVALENCE), 4)
        out["operating_points"][name] = c

    # ---- threshold sweep on the best continuous rule ----------------------
    best_rule = max(out["rules"], key=lambda k: out["rules"][k]["auc"])
    key = "s_" + best_rule
    # sweep on quantiles of the score rather than every distinct value: the
    # full sweep is O(n^2) and adds nothing but runtime at this sample size
    uniq = sorted({r[key] for r in rows})
    step = max(len(uniq) // 400, 1)
    cand = uniq[::step] + [uniq[-1]]
    sweep = []
    for t in cand:
        pred = [r[key] >= t for r in rows]
        c = confusion(pred, labels)
        c["threshold_score"] = round(t, 4)
        c["threshold_rank"] = round(10 ** (-t), 4)
        sweep.append(c)
    by_mcc = max(sweep, key=lambda c: (c["mcc"] if not math.isnan(c["mcc"]) else -9))
    at_spec = [c for c in sweep if not math.isnan(c["specificity"]) and c["specificity"] >= 0.95]
    at_sens = [c for c in sweep if not math.isnan(c["sensitivity"]) and c["sensitivity"] >= 0.80]
    out["best_continuous_rule"] = best_rule
    out["sweep"] = {
        "max_mcc": by_mcc,
        "highest_sensitivity_at_spec95": max(at_spec, key=lambda c: c["sensitivity"]) if at_spec else None,
        "highest_specificity_at_sens80": max(at_sens, key=lambda c: c["specificity"]) if at_sens else None,
    }
    for k in ("max_mcc", "highest_sensitivity_at_spec95", "highest_specificity_at_sens80"):
        c = out["sweep"][k]
        if c:
            c["ppv_at_scan_prevalence"] = round(
                ppv_at_prevalence(c["sensitivity"], c["specificity"], SCAN_PREVALENCE), 4)

    # ---- per-allele, to see whether one global threshold is defensible ----
    per_allele = {}
    for allele in sorted({r["allele"] for r in rows}):
        sub = [r for r in rows if r["allele"] == allele]
        ls = [r["label"] for r in sub]
        if len(sub) < 30 or not (0 < sum(ls) < len(ls)):
            continue
        key = "s_" + max(out["rules"], key=lambda k: out["rules"][k]["auc"])
        # threshold that maximises MCC within this allele, for comparison with
        # the global one; reported, not applied - per-allele tuning on this
        # much data would fit noise
        uq = sorted({r[key] for r in sub})
        best = max(({"t": t, **confusion([r[key] >= t for r in sub], ls)}
                    for t in uq[::max(len(uq) // 150, 1)]),
                   key=lambda c: (c["mcc"] if not math.isnan(c["mcc"]) else -9))
        per_allele[allele] = {
            "n": len(sub), "n_positive": sum(ls),
            "auc": round(roc_auc([r[key] for r in sub], ls), 3),
            "best_local_threshold_rank": round(10 ** (-best["t"]), 3),
            "best_local_mcc": best["mcc"],
        }
    out["per_allele"] = per_allele

    # ---- stratified by self / non-self ------------------------------------
    for stratum, subset in (("self", [r for r in rows if r["is_self"]]),
                            ("non_self", [r for r in rows if not r["is_self"]])):
        if len(subset) < 30 or len({r["label"] for r in subset}) < 2:
            out["stratified"][stratum] = {"n": len(subset), "note": "too few to evaluate"}
            continue
        ls = [r["label"] for r in subset]
        out["stratified"][stratum] = {
            "n": len(subset), "n_positive": sum(ls),
            **{name: round(roc_auc([r["s_" + name] for r in subset], ls), 4)
               for name in ("EL", "BA", "GEO", "MIN")},
        }

    with open(results_path("m11_calibration.json"), "w") as f:
        json.dump(out, f, indent=2)

    # ---- report -----------------------------------------------------------
    print(f"\n{out['n_pairs']} pairs ({out['n_positive']} positive) in "
          f"{out['n_clusters']} clusters\n")
    print(f"{'rule':6s} {'ROC AUC':>8s} {'avg prec':>9s}")
    for name, v in out["rules"].items():
        print(f"{name:6s} {v['auc']:8.3f} {v['average_precision']:9.3f}")
    print("\npaired bootstrap over clusters, AUC difference:")
    for k, v in out["comparisons"].items():
        if v:
            print(f"  {k:12s} {v['delta']:+.4f}  95% CI [{v['ci95'][0]:+.4f}, "
                  f"{v['ci95'][1]:+.4f}]  P(>0)={v['p_gt_0']:.3f}")
    print("\noperating points:")
    hdr = f"  {'rule':34s} {'sens':>6s} {'spec':>6s} {'PPV':>6s} {'MCC':>6s} {'PPV@scan':>9s}"
    print(hdr)
    for name, c in out["operating_points"].items():
        print(f"  {name:34s} {c['sensitivity']:6.3f} {c['specificity']:6.3f} "
              f"{c['ppv']:6.3f} {c['mcc']:6.3f} {c['ppv_at_scan_prevalence']:9.3f}")
    print(f"\nbest continuous rule: {best_rule}")
    for k in ("max_mcc", "highest_sensitivity_at_spec95", "highest_specificity_at_sens80"):
        c = out["sweep"][k]
        if c:
            print(f"  {k:32s} %Rank<{c['threshold_rank']:<7.3f} sens {c['sensitivity']:.3f} "
                  f"spec {c['specificity']:.3f} MCC {c['mcc']:.3f} "
                  f"PPV@scan {c['ppv_at_scan_prevalence']:.3f}")
    if out.get("per_allele"):
        aucs = [v["auc"] for v in out["per_allele"].values()]
        locs = [v["best_local_threshold_rank"] for v in out["per_allele"].values()]
        print(f"\nper-allele ({len(aucs)} alleles with >=30 pairs): "
              f"AUC {min(aucs):.3f}-{max(aucs):.3f}, "
              f"locally optimal %Rank cut {min(locs):.2f}-{max(locs):.2f}")

    print("\nstratified AUC:")
    for s, v in out["stratified"].items():
        if "note" in v:
            print(f"  {s:9s} n={v['n']}  {v['note']}")
        else:
            print(f"  {s:9s} n={v['n']:5d} ({v['n_positive']} pos)  "
                  + "  ".join(f"{k} {v[k]:.3f}" for k in ("EL", "BA", "GEO", "MIN")))


if __name__ == "__main__":
    main()
