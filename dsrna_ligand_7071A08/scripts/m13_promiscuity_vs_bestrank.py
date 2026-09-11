#!/usr/bin/env python3
"""
M13 - Does promiscuity add specificity over the best single-allele rank?

The pipeline does not flag peptides one DR molecule at a time. It scans a
ligand against the whole panel and reports a peptide as risk when several
molecules present it, weighting each by how common that molecule is. That
implies a claim nobody has tested: that *breadth across the panel* is a better
predictor of a real CD4 response than the single best %Rank the peptide
achieves anywhere.

It could easily be false. Breadth is a deterministic function of the same
per-allele scores, so it may add nothing but a smoothing of the best rank. If
it does add nothing, the honest simplification is to flag on best rank and use
population weighting only to describe who is affected, not to decide.

  label      a peptide is positive if IEDB holds a positive human DR-restricted
             T-cell assay for it on ANY molecule in the panel, negative if every
             record for it is negative
  predictors best_rank      min EL %Rank across all 25 DR molecules
             breadth_sb     how many molecules reach %Rank < 1
             breadth_wb     how many reach %Rank < 5
             pop_presenting weighted US/EU fraction carrying a presenting
                            molecule - the metric the pipeline actually scores on

THE CONFOUND, and why it is the whole difficulty here. A peptide's label
depends on how many DR molecules somebody happened to test it against: assayed
on five molecules it has five chances to come out positive, assayed on one it
has one. Peptides tested broadly are also, on average, the ones a study already
believed were promiscuous. So test count is correlated with the label *and*
plausibly with predicted breadth, and an uncorrected comparison would credit
breadth for detecting how much attention a peptide received.

Three things keep that from happening:
  * the AUC of test count *by itself* is reported first, as the size of the
    problem;
  * the comparison is repeated within each test-count stratum, where the count
    cannot vary; and
  * the primary analysis is restricted to peptides tested on exactly ONE
    molecule, where the label cannot depend on breadth of testing at all.
If breadth only wins in the pooled analysis and not in the single-allele
stratum, it was reading ascertainment, and the module says so.

Sampling takes at most one peptide per 9-mer cluster, so the rows are close to
independent before any bootstrap is applied.

Output: results/m13_promiscuity.json, results/m13_peptide_scored.tsv
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
from common import (load_config, results_path, data_path, iedb_mhcii,  # noqa: E402
                    CoverageModel, population_weights)
from m11_threshold_calibration import (roc_auc, average_precision, confusion,  # noqa: E402
                                       ppv_at_prevalence, SPACER_AA, CHUNK_AA)

MAX_PEPTIDES = 1200
SEED = 20260825
BOOTSTRAP = 800
SCAN_PREVALENCE = 0.05


def cluster_bootstrap_delta(rows, key_a, key_b, n=BOOTSTRAP, seed=SEED):
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
    if not deltas:
        return None
    deltas.sort()
    labels = [r["label"] for r in rows]
    point = (roc_auc([r[key_a] for r in rows], labels)
             - roc_auc([r[key_b] for r in rows], labels))
    return {"delta": round(point, 4),
            "ci95": [round(deltas[int(0.025 * len(deltas))], 4),
                     round(deltas[min(int(0.975 * len(deltas)), len(deltas) - 1)], 4)],
            "p_gt_0": round(sum(1 for d in deltas if d > 0) / len(deltas), 4)}


def score_panel(cfg, peptides, panel, cache_path):
    """
    EL %Rank for every (peptide, allele) over the whole panel.

    Only the eluted-ligand head: M11 established that the affinity head neither
    improves ranking nor removes false positives selectively, so paying for it
    again would double a long job for nothing.
    """
    pcfg = cfg["prediction"]
    have = defaultdict(dict)
    if os.path.exists(cache_path):
        with open(cache_path) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                have[r["peptide"]][r["allele"]] = float(r["el_rank"])
        print(f"resuming: {sum(len(v) for v in have.values())} peptide-allele scores cached")

    jobs = []
    for allele in panel:
        for length in (13, 14, 15):
            peps = sorted({p for p in peptides
                           if min(15, len(p)) == length
                           and allele not in have.get(p, {})})
            if not peps:
                continue
            spacer = SPACER_AA * (length - 1)
            cur, layout, pos = [], [], 1
            for pep in peps:
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
    if not jobs:
        return have
    print(f"{len(peptides)} peptides x {len(panel)} DR molecules "
          f"-> {len(jobs)} calls still to make", flush=True)

    def run(job):
        allele, length, seq, layout = job
        rows = iedb_mhcii(pcfg["endpoint"], "netmhciipan_el",
                          {"cat": seq}, [allele], length)
        res = {}
        for r in rows:
            s0, e0 = int(r["start"]), int(r["end"])
            for pep, ps, pe in layout:
                if ps <= s0 and e0 <= pe:
                    v = float(r["rank"])
                    if v < res.get((pep, allele), 9e9):
                        res[(pep, allele)] = v
                    break
        print(f"  ok {allele:18s} len{length} {len(layout):4d} peptides", flush=True)
        return res

    exists = os.path.exists(cache_path)
    with open(cache_path, "a", newline="") as cf:
        cw = csv.writer(cf, delimiter="\t")
        if not exists:
            cw.writerow(["peptide", "allele", "el_rank"])
        with ThreadPoolExecutor(max_workers=4) as ex:
            for res in ex.map(run, jobs):
                for (pep, allele), v in res.items():
                    have[pep][allele] = v
                    cw.writerow([pep, allele, v])
                cf.flush()
    return have


def main():
    cfg = load_config()
    sb = cfg["prediction"]["sb_rank"]
    wb = cfg["prediction"]["wb_rank"]
    with open(results_path("m2_panel_alleles.txt")) as f:
        panel = [l.strip() for l in f if l.strip()]
    with open(data_path("drb1_allele_frequencies.json")) as f:
        tables = json.load(f)
    coverage = CoverageModel(tables, population_weights(cfg)).weighted

    # ---- peptide-level labels --------------------------------------------
    pep = {}
    with open(results_path("m10_benchmark.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            d = pep.setdefault(r["peptide"], {
                "peptide": r["peptide"], "cluster": int(r["cluster"]),
                "is_self": bool(int(r["is_self"])), "pos": 0, "tested": set()})
            d["tested"].add(r["allele"])
            d["pos"] += int(r["label"])
    for d in pep.values():
        d["label"] = 1 if d["pos"] else 0
        d["n_tested"] = len(d["tested"])

    # one peptide per cluster, balanced, seeded
    rng = random.Random(SEED)
    by_cluster = defaultdict(list)
    for d in pep.values():
        by_cluster[d["cluster"]].append(d)
    reps = []
    for c in sorted(by_cluster):
        members = sorted(by_cluster[c], key=lambda d: d["peptide"])
        reps.append(rng.choice(members))
    pos = [d for d in reps if d["label"] == 1]
    neg = [d for d in reps if d["label"] == 0]
    rng.shuffle(pos)
    rng.shuffle(neg)
    cap = MAX_PEPTIDES // 2
    sample = pos[:cap] + neg[:cap]
    print(f"{len(pep)} distinct peptides in {len(by_cluster)} clusters "
          f"-> {len(sample)} sampled ({len(pos[:cap])} positive), one per cluster")

    scores = score_panel(cfg, [d["peptide"] for d in sample], panel,
                         results_path("m13_panel_scores.tsv"))

    rows = []
    for d in sample:
        s = scores.get(d["peptide"], {})
        if len(s) < len(panel):
            continue                      # incomplete panel; excluded, counted below
        ranks = {a: s[a] for a in panel}
        hits_sb = [a for a, v in ranks.items() if v < sb]
        hits_wb = [a for a, v in ranks.items() if v < wb]
        rows.append({
            **{k: d[k] for k in ("peptide", "cluster", "label", "is_self", "n_tested")},
            "best_rank": min(ranks.values()),
            "breadth_sb": len(hits_sb),
            "breadth_wb": len(hits_wb),
            "pop_presenting": coverage(hits_sb),
            "pop_presenting_wb": coverage(hits_wb),
            # higher = more epitope-like, for every predictor
            "s_best_rank": -math.log10(max(min(ranks.values()), 1e-3)),
            "s_breadth_sb": float(len(hits_sb)),
            "s_breadth_wb": float(len(hits_wb)),
            "s_pop_presenting": coverage(hits_sb),
            "s_pop_presenting_wb": coverage(hits_wb),
            "s_n_tested": float(d["n_tested"]),
        })
    dropped = len(sample) - len(rows)
    print(f"{len(rows)} peptides scored on the full panel"
          + (f" ({dropped} dropped for incomplete panel coverage)" if dropped else ""))

    with open(results_path("m13_peptide_scored.tsv"), "w", newline="") as f:
        cols = ["peptide", "cluster", "label", "is_self", "n_tested", "best_rank",
                "breadth_sb", "breadth_wb", "pop_presenting", "pop_presenting_wb"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    labels = [r["label"] for r in rows]
    preds = ("best_rank", "breadth_sb", "breadth_wb", "pop_presenting",
             "pop_presenting_wb", "n_tested")
    out = {
        "n_peptides": len(rows), "n_positive": sum(labels),
        "n_clusters": len({r["cluster"] for r in rows}),
        "scan_prevalence": SCAN_PREVALENCE,
        "predictors": {p: {
            "auc": round(roc_auc([r["s_" + p] for r in rows], labels), 4),
            "average_precision": round(
                average_precision([r["s_" + p] for r in rows], labels), 4)}
            for p in preds},
        "comparisons": {}, "operating_points": {}, "strata": {},
    }

    for p in ("breadth_sb", "breadth_wb", "pop_presenting", "pop_presenting_wb"):
        out["comparisons"][f"{p}_vs_best_rank"] = cluster_bootstrap_delta(
            rows, "s_" + p, "s_best_rank")

    # ---- operating points -------------------------------------------------
    points = {
        f"best %Rank < {sb:g} (any one molecule)": [r["best_rank"] < sb for r in rows],
        "breadth >= 2 molecules at %Rank < 1": [r["breadth_sb"] >= 2 for r in rows],
        "breadth >= 3 molecules at %Rank < 1": [r["breadth_sb"] >= 3 for r in rows],
        "breadth >= 4 molecules at %Rank < 1": [r["breadth_sb"] >= 4 for r in rows],
        "population presenting >= 20%": [r["pop_presenting"] >= 0.20 for r in rows],
        "population presenting >= 40%": [r["pop_presenting"] >= 0.40 for r in rows],
    }
    for name, pred in points.items():
        c = confusion(pred, labels)
        c["n_flagged"] = sum(pred)
        c["ppv_at_scan_prevalence"] = round(
            ppv_at_prevalence(c["sensitivity"], c["specificity"], SCAN_PREVALENCE), 4)
        out["operating_points"][name] = c

    # ---- the confound, three ways -----------------------------------------
    single = [r for r in rows if r["n_tested"] == 1]
    multi = [r for r in rows if r["n_tested"] > 1]
    out["confound"] = {
        "auc_of_test_count_alone": out["predictors"]["n_tested"]["auc"],
        "positive_rate_by_test_count": {},
        "single_allele_stratum": None,
        "by_test_count": {},
    }
    counts = defaultdict(lambda: [0, 0])
    for r in rows:
        b = min(r["n_tested"], 4)
        counts[b][0] += 1
        counts[b][1] += r["label"]
    for b in sorted(counts):
        n, k = counts[b]
        out["confound"]["positive_rate_by_test_count"][
            f"{b}{'+' if b == 4 else ''}"] = {"n": n, "positive_rate": round(k / n, 3)}

    for name, subset in (("single_allele", single), ("multi_allele", multi)):
        if len(subset) < 60 or len({r["label"] for r in subset}) < 2:
            out["strata"][name] = {"n": len(subset), "note": "too few to evaluate"}
            continue
        ls = [r["label"] for r in subset]
        entry = {"n": len(subset), "n_positive": sum(ls),
                 **{p: round(roc_auc([r["s_" + p] for r in subset], ls), 4)
                    for p in preds if p != "n_tested"}}
        entry["breadth_sb_vs_best_rank"] = cluster_bootstrap_delta(
            subset, "s_breadth_sb", "s_best_rank")
        entry["pop_presenting_vs_best_rank"] = cluster_bootstrap_delta(
            subset, "s_pop_presenting", "s_best_rank")
        out["strata"][name] = entry
    out["confound"]["single_allele_stratum"] = out["strata"].get("single_allele")

    with open(results_path("m13_promiscuity.json"), "w") as f:
        json.dump(out, f, indent=2)

    # ---- report -----------------------------------------------------------
    print(f"\n{out['n_peptides']} peptides ({out['n_positive']} positive) "
          f"in {out['n_clusters']} clusters\n")
    print(f"{'predictor':22s} {'ROC AUC':>8s} {'avg prec':>9s}")
    for p, v in out["predictors"].items():
        print(f"{p:22s} {v['auc']:8.3f} {v['average_precision']:9.3f}")
    print("\nvs best single-allele rank (paired bootstrap over clusters):")
    for k, v in out["comparisons"].items():
        if v:
            print(f"  {k:34s} {v['delta']:+.4f}  95% CI "
                  f"[{v['ci95'][0]:+.4f}, {v['ci95'][1]:+.4f}]  P(>0)={v['p_gt_0']:.3f}")
    print("\noperating points:")
    hdr = (f"  {'rule':38s} {'flagged':>8s} {'sens':>6s} {'spec':>6s} "
           f"{'PPV':>6s} {'MCC':>6s} {'PPV@scan':>9s}")
    print(hdr)
    for name, c in out["operating_points"].items():
        print(f"  {name:38s} {c['n_flagged']:8d} {c['sensitivity']:6.3f} "
              f"{c['specificity']:6.3f} {c['ppv']:6.3f} {c['mcc']:6.3f} "
              f"{c['ppv_at_scan_prevalence']:9.3f}")
    print(f"\nconfound: AUC of test count alone = "
          f"{out['confound']['auc_of_test_count_alone']:.3f}")
    for b, v in out["confound"]["positive_rate_by_test_count"].items():
        print(f"  tested on {b} molecule(s): n={v['n']:4d}  "
              f"positive rate {v['positive_rate']:.3f}")
    for name, v in out["strata"].items():
        if "note" in v:
            print(f"\n{name}: n={v['n']}  {v['note']}")
            continue
        print(f"\n{name} stratum (n={v['n']}, {v['n_positive']} positive):")
        print("  AUC  " + "  ".join(f"{p} {v[p]:.3f}" for p in preds if p != "n_tested"))
        for k in ("breadth_sb_vs_best_rank", "pop_presenting_vs_best_rank"):
            d = v.get(k)
            if d:
                print(f"  {k:30s} {d['delta']:+.4f}  95% CI "
                      f"[{d['ci95'][0]:+.4f}, {d['ci95'][1]:+.4f}]")


if __name__ == "__main__":
    main()
