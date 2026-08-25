#!/usr/bin/env python3
"""
M2 - HLA-DR panel design with a *computed* population-coverage guarantee.

The usual practice is to grab the DR subset of the IEDB class-II reference set
(15 molecules) and assert that it is "broadly representative". It is not:
measured against the IEDB allele-frequency tables that set reaches only ~90%
US-Caucasian / ~86% European DRB1 phenotypic coverage, short of a 95-98%
requirement.

This module instead:
  1. pulls the IEDB population-coverage frequency tables (the same data the
     IEDB web tool uses),
  2. builds a US composite population weighted by US Census 2020 shares plus a
     European population,
  3. greedily grows a DRB1 panel until the weighted coverage target is met,
  4. adds the four DRB3/4/5 molecules for second-DR-molecule breadth, and
  5. writes the coverage curve so the panel size is a defensible choice rather
     than a round number.

Coverage model: common.CoverageModel - single-locus Hardy-Weinberg over the
IEDB frequency tables, with the tool's renormalisation for populations whose
allele frequencies sum above 1. Verified against the IEDB CLI to two decimals.
"""
import json
import os
import sys
import tarfile
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from common import (load_config, results_path, data_path,  # noqa: E402
                    CoverageModel, population_weights)

POPCOV_URL = ("https://downloads.iedb.org/tools/population/3.0.2/"
              "IEDB_Population_Coverage-3.0.2.tar.gz")


def load_frequency_tables():
    """Return {population: {allele: allele_frequency}} for HLA-DRB1."""
    cache = data_path("drb1_allele_frequencies.json")
    if os.path.exists(cache):
        with open(cache) as f:
            return json.load(f)

    workdir = data_path("_iedb_popcov")
    os.makedirs(workdir, exist_ok=True)
    tgz = os.path.join(workdir, "popcov.tar.gz")
    if not os.path.exists(tgz):
        urllib.request.urlretrieve(POPCOV_URL, tgz)
    with tarfile.open(tgz) as t:
        t.extractall(workdir)

    pkg = os.path.join(workdir, "population_coverage")
    sys.path.insert(0, os.path.join(pkg, "deps", "population-coverage-pickle"))
    sys.path.insert(0, pkg)
    from population_coverage_pickle import population_coverage  # noqa: E402

    tables = {}
    for pop, loci in population_coverage["II"].items():
        if "HLA-DRB1" in loci:
            tables[pop] = {a: f for a, f in loci["HLA-DRB1"] if f > 0}
    with open(cache, "w") as f:
        json.dump(tables, f, indent=0, sort_keys=True)
    return tables


def greedy_panel(model, tables, target, min_alleles):
    """Grow a DRB1 panel by largest marginal weighted-coverage gain."""
    candidates = sorted({a for t in tables.values() for a in t})
    panel, curve = [], []
    while True:
        best, best_cov = None, model.weighted(panel)
        for a in candidates:
            if a in panel:
                continue
            c = model.weighted(panel + [a])
            if c > best_cov + 1e-12:
                best, best_cov = a, c
        if best is None:
            break
        panel.append(best)
        curve.append((len(panel), best, best_cov))
        if best_cov >= target and len(panel) >= min_alleles:
            break
    return panel, curve


def main():
    cfg = load_config()
    tables = load_frequency_tables()
    weights = population_weights(cfg)
    for pop in weights:
        if pop not in tables:
            sys.exit(f"population {pop!r} absent from the IEDB tables")
    model = CoverageModel(tables, weights)
    target = cfg["panel"]["coverage_target"]
    min_n = cfg["panel"]["min_alleles"]

    # The legacy panel this pipeline is replacing, for a like-for-like number.
    legacy = ["HLA-DRB1*01:01", "HLA-DRB1*03:01", "HLA-DRB1*04:01", "HLA-DRB1*04:05",
              "HLA-DRB1*07:01", "HLA-DRB1*08:02", "HLA-DRB1*09:01", "HLA-DRB1*11:01",
              "HLA-DRB1*12:01", "HLA-DRB1*13:02", "HLA-DRB1*15:01"]

    panel, curve = greedy_panel(model, tables, target, min_n)

    report_pops = (list(cfg["panel"]["us_weights"]) + [cfg["panel"]["eu_population"]]
                   + ["World"])
    out = {
        "weights": weights,
        "target": target,
        "drb1_panel": panel,
        "drb345_panel": cfg["panel"]["include_drb345"],
        "panel_size_drb1": len(panel),
        "panel_size_total": len(panel) + len(cfg["panel"]["include_drb345"]),
        "weighted_coverage": model.weighted(panel),
        "legacy_drb1_panel": legacy,
        "legacy_weighted_coverage": model.weighted(legacy),
        "per_population": {},
        "legacy_per_population": {},
        "curve": [{"n": n, "added": a, "weighted_coverage": c} for n, a, c in curve],
    }
    for p in report_pops:
        if p in tables:
            out["per_population"][p] = model.population(panel, p)
            out["legacy_per_population"][p] = model.population(legacy, p)

    with open(results_path("m2_panel.json"), "w") as f:
        json.dump(out, f, indent=2)

    with open(results_path("m2_panel_alleles.txt"), "w") as f:
        for a in panel + cfg["panel"]["include_drb345"]:
            f.write(a + "\n")

    print(f"DRB1 panel: {len(panel)} alleles "
          f"(+{len(cfg['panel']['include_drb345'])} DRB3/4/5 "
          f"= {out['panel_size_total']} DR molecules)")
    print(f"weighted US/EU DRB1 coverage: {out['weighted_coverage']*100:.2f}%  "
          f"(legacy 15-molecule panel: {out['legacy_weighted_coverage']*100:.2f}%)")
    print()
    print(f"{'population':30s} {'new':>8s} {'legacy':>8s}")
    for p in report_pops:
        if p in out["per_population"]:
            print(f"{p:30s} {out['per_population'][p]*100:7.2f}% "
                  f"{out['legacy_per_population'][p]*100:7.2f}%")
    print()
    print("greedy build order:")
    for row in out["curve"]:
        print(f"  {row['n']:2d}  {row['added']:18s} -> {row['weighted_coverage']*100:6.2f}%")


if __name__ == "__main__":
    main()
