#!/usr/bin/env python3
"""
M12 - Put an empirical bound on the tolerance discount.

M4 down-weights a predicted epitope whose 9-mer core is identical (weight 0) or
near-identical (weight 0.35) to a human proteome 9-mer. The 0.35 is a guess. It
is also load-bearing: it sets how much of a VHH framework's predicted content
survives into the score, and therefore where an antibody-derived ligand lands
against the Protein A benchmark.

The M10 benchmark can bound it. Among peptides the predictor calls binders,
compare how often self-like peptides actually elicit a T-cell response against
how often foreign ones do. That ratio is what the discount is trying to encode:

    discount = P(T-cell positive | binder, self-like)
             / P(T-cell positive | binder, foreign)

One bias makes this an **upper bound on the ratio, and therefore a floor on how
much discounting is warranted** - the estimate errs toward saying "self
peptides respond nearly as often", never the reverse. Self-derived peptides
enter IEDB overwhelmingly through autoimmunity, allergy and tumour-antigen
studies, which select for the self peptides that *do* respond. Self peptides
nobody responds to are not interesting enough to assay, so they are missing
from the denominator. A truly random self peptide is less likely to be an
epitope than this number says.

The module therefore reports the ratio, its bootstrap interval, and the
direction of the bias, and recommends a weight no larger than the estimate.

Output: results/m12_tolerance_weight.json
"""
import csv
import json
import math
import os
import random
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import load_config, results_path, data_path  # noqa: E402

BOOTSTRAP = 2000
SEED = 20260825


def iter_proteome(path):
    buf, name = [], None
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    yield "".join(buf)
                name, buf = line, []
            else:
                buf.append(line.strip())
    if name:
        yield "".join(buf)


def best_human_identity(peptides, proteome):
    """Max identity (9/9 or 8/9) of any 9-mer of each peptide to the human proteome."""
    cores = defaultdict(set)          # 9-mer -> peptides containing it
    for p in peptides:
        for i in range(len(p) - 8):
            cores[p[i:i + 9]].add(p)
    idx = defaultdict(set)
    for c in cores:
        for i in range(9):
            idx[c[:i] + "." + c[i + 1:]].add(c)
    print(f"  indexing {len(cores)} distinct 9-mers from {len(peptides)} peptides")

    best = defaultdict(int)
    n = 0
    for seq in iter_proteome(proteome):
        n += 1
        for i in range(len(seq) - 8):
            k = seq[i:i + 9]
            hits = set()
            for j in range(9):
                hits |= idx.get(k[:j] + "." + k[j + 1:], set())
            for c in hits:
                ident = sum(a == b for a, b in zip(c, k))
                for p in cores[c]:
                    if ident > best[p]:
                        best[p] = ident
    print(f"  screened {n} human proteins")
    return best


def ratio_ci(self_rows, foreign_rows, n=BOOTSTRAP, seed=SEED):
    """Bootstrap over clusters on the positive-rate ratio."""
    def by_cluster(rows):
        d = defaultdict(list)
        for r in rows:
            d[r["cluster"]].append(r)
        return d

    sc, fc = by_cluster(self_rows), by_cluster(foreign_rows)
    sk, fk = list(sc), list(fc)
    rng = random.Random(seed)
    out = []
    for _ in range(n):
        s = [x for _ in sk for x in sc[sk[rng.randrange(len(sk))]]]
        f = [x for _ in fk for x in fc[fk[rng.randrange(len(fk))]]]
        if not s or not f:
            continue
        ps = sum(r["label"] for r in s) / len(s)
        pf = sum(r["label"] for r in f) / len(f)
        if pf > 0:
            out.append(ps / pf)
    out.sort()
    if not out:
        return None
    return [round(out[int(0.025 * len(out))], 3),
            round(out[min(int(0.975 * len(out)), len(out) - 1)], 3)]


def main():
    cfg = load_config()
    scored_path = results_path("m11_benchmark_scored.tsv")
    if not os.path.exists(scored_path):
        sys.exit("run m11_threshold_calibration.py first")
    with open(scored_path) as f:
        rows = [dict(r) for r in csv.DictReader(f, delimiter="\t")]
    for r in rows:
        r["label"] = int(r["label"])
        r["cluster"] = int(r["cluster"])
        r["el_rank"] = float(r["el_rank"])
        r["ba_rank"] = float(r["ba_rank"])

    proteome = os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), cfg["tolerance_filter"]["proteome"])
    peptides = sorted({r["peptide"] for r in rows})
    print(f"screening {len(peptides)} benchmark peptides against the human proteome")
    best = best_human_identity(peptides, proteome)

    # Restrict to peptides the predictor would flag at all - the discount only
    # ever applies to predicted binders, so that is the population to measure on.
    sb = cfg["prediction"]["sb_rank"]
    binders = [r for r in rows if r["el_rank"] < sb]
    if len(binders) < 100:
        binders = [r for r in rows if r["el_rank"] < cfg["prediction"]["wb_rank"]]
        gate = f"EL %Rank < {cfg['prediction']['wb_rank']:g}"
    else:
        gate = f"EL %Rank < {sb:g}"

    out = {"gate": gate, "n_binders": len(binders), "strata": {}}
    for name, pred in (("exact_self_9mer", lambda p: best.get(p, 0) == 9),
                       ("near_self_9mer", lambda p: best.get(p, 0) == 8),
                       ("self_like_9_or_8", lambda p: best.get(p, 0) >= 8),
                       ("foreign", lambda p: best.get(p, 0) < 8)):
        sub = [r for r in binders if pred(r["peptide"])]
        out["strata"][name] = {
            "n": len(sub),
            "n_positive": sum(r["label"] for r in sub),
            "positive_rate": round(sum(r["label"] for r in sub) / len(sub), 4) if sub else None,
        }

    slf = [r for r in binders if best.get(r["peptide"], 0) >= 8]
    frn = [r for r in binders if best.get(r["peptide"], 0) < 8]
    est = None
    if slf and frn and out["strata"]["foreign"]["positive_rate"]:
        est = (out["strata"]["self_like_9_or_8"]["positive_rate"]
               / out["strata"]["foreign"]["positive_rate"])
        out["discount_estimate"] = round(est, 3)
        out["discount_ci95"] = ratio_ci(slf, frn)
    out["configured_discount_exact"] = cfg["tolerance_filter"]["discount_exact"]
    out["configured_discount_near"] = cfg["tolerance_filter"]["discount_tcrface"]
    out["bias_direction"] = (
        "Upper bound. Self-derived peptides reach IEDB through autoimmunity, allergy and "
        "tumour-antigen studies, which select for the self peptides that do respond; self "
        "peptides nobody responds to are rarely assayed and are missing from the denominator. "
        "The true positive rate for an arbitrary self peptide is lower than measured here, so "
        "the real discount is at least this strong.")
    out["recommendation"] = (
        f"Use a near-self weight no greater than {out.get('discount_estimate', 'n/a')}; "
        f"the configured {cfg['tolerance_filter']['discount_tcrface']} is "
        + ("consistent with the bound." if est is not None and
           cfg["tolerance_filter"]["discount_tcrface"] <= est
           else "ABOVE the bound and should be lowered.")
        if est is not None else "insufficient data to bound the weight.")

    with open(results_path("m12_tolerance_weight.json"), "w") as f:
        json.dump(out, f, indent=2)

    print(f"\ngate: {gate}   {len(binders)} predicted binders in the benchmark\n")
    print(f"{'stratum':22s} {'n':>6s} {'positive':>9s} {'rate':>7s}")
    for k, v in out["strata"].items():
        print(f"{k:22s} {v['n']:6d} {v['n_positive']:9d} "
              f"{(v['positive_rate'] if v['positive_rate'] is not None else float('nan')):7.3f}")
    if est is not None:
        ci = out["discount_ci95"]
        print(f"\nempirical discount (self-like / foreign positive rate): {est:.3f}"
              + (f"   95% CI [{ci[0]}, {ci[1]}]" if ci else ""))
        print(f"configured near-self weight: {cfg['tolerance_filter']['discount_tcrface']}")
    print(f"\n{out['bias_direction']}")
    print(f"\n{out['recommendation']}")


if __name__ == "__main__":
    main()
