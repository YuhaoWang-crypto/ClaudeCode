#!/usr/bin/env python3
"""
C5 - How hard is it to find the donors this RFP asks for?

The RFP asks for "human PBMC with the specified HLA antigens" (HLA-A*32:01 and
HLA-B*57:01) plus full HLA typing of every donor used. Nobody in the RFP says
how rare that is, and it drives the cost and the timeline more than any assay
choice does. This is pure arithmetic on published allele frequencies and it can
be done before a single peptide is ordered.

Frequencies come from the IEDB population-coverage tables (the same source the
class II panel work used), so they are traceable rather than remembered.

Phenotype frequency from allele frequency is single-locus Hardy-Weinberg:

    P(carries the allele) = 1 - (1 - f)^2

Two loci are combined as if independent. That is an approximation and it is
stated as one: A and B sit ~1.3 Mb apart on chromosome 6 and are in strong
linkage disequilibrium, so the true double-carrier rate must come from
haplotype tables. The independent estimate is reported as a planning figure
with that caveat attached, and the direction of the error is flagged.

Output: results/c5_donor_feasibility.json
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import data_path, results_path  # noqa: E402

ALLELES = {"HLA-A": "HLA-A*32:01", "HLA-B": "HLA-B*57:01"}
# US Census 2020 shares, the same weighting the class II work used.
US_WEIGHTS = {
    "United States Caucasoid": 0.577, "United States Hispanic": 0.189,
    "United States Black": 0.123, "United States Asian": 0.060,
    "United States Amerindian": 0.013, "United States Mestizo": 0.038,
}
US_EU_SPLIT = (0.5, 0.5)
EU_POP = "Europe"
# Donors needed carrying the allele(s), for the sample-size arithmetic.
WANT_DONORS = [5, 10, 20]


def load_class_i():
    pkg = os.path.abspath(data_path("_iedb_popcov", "population_coverage"))
    sys.path.insert(0, os.path.join(pkg, "deps", "population-coverage-pickle"))
    sys.path.insert(0, pkg)
    from population_coverage_pickle import population_coverage  # noqa: E402
    return population_coverage["I"]


def phenotype(f):
    """Hardy-Weinberg carrier frequency from allele frequency."""
    return 1.0 - (1.0 - f) ** 2


def weights():
    tot = sum(US_WEIGHTS.values())
    us_w, eu_w = US_EU_SPLIT
    w = {p: us_w * v / tot for p, v in US_WEIGHTS.items()}
    w[EU_POP] = w.get(EU_POP, 0.0) + eu_w
    return w


def screened_for(p, n):
    """Donors to type, in expectation, to end up with n carriers."""
    return None if p <= 0 else n / p


def main():
    tables = load_class_i()
    w = weights()
    missing = [p for p in w if p not in tables]
    if missing:
        raise SystemExit(f"population tables missing: {missing}")

    per_pop = {}
    for pop in w:
        loci = tables[pop]
        rec = {}
        for locus, allele in ALLELES.items():
            f = dict(loci.get(locus, [])).get(allele, 0.0)
            rec[allele] = {"allele_frequency": round(f, 5),
                           "carrier_frequency": round(phenotype(f), 5)}
        a = rec[ALLELES["HLA-A"]]["carrier_frequency"]
        b = rec[ALLELES["HLA-B"]]["carrier_frequency"]
        rec["either"] = round(1 - (1 - a) * (1 - b), 5)
        rec["both_independent"] = round(a * b, 5)
        per_pop[pop] = rec

    def wavg(get):
        return sum(w[p] * get(per_pop[p]) for p in w)

    a_name, b_name = ALLELES["HLA-A"], ALLELES["HLA-B"]
    comb = {
        a_name: round(wavg(lambda r: r[a_name]["carrier_frequency"]), 5),
        b_name: round(wavg(lambda r: r[b_name]["carrier_frequency"]), 5),
        "either": round(wavg(lambda r: r["either"]), 5),
        "both_independent": round(wavg(lambda r: r["both_independent"]), 5),
    }

    burden = {}
    for key in (a_name, b_name, "either", "both_independent"):
        p = comb[key]
        burden[key] = {str(n): (round(screened_for(p, n), 1)
                                if screened_for(p, n) else None)
                       for n in WANT_DONORS}

    out = {
        "alleles": ALLELES,
        "weighting": {"us_weights": US_WEIGHTS, "us_eu_split": list(US_EU_SPLIT),
                      "eu_population": EU_POP},
        "per_population": per_pop,
        "weighted_us_eu": comb,
        "donors_to_screen": burden,
        "caveats": [
            "Two loci combined as independent. HLA-A and HLA-B are ~1.3 Mb apart "
            "and in strong linkage disequilibrium, so the true double-carrier rate "
            "is not the product. B*57:01 in Europeans travels mostly on the 57.1 "
            "ancestral haplotype (A*01:01-B*57:01-C*06:02-DRB1*07:01), which "
            "carries A*01:01 rather than A*32:01 - so for THIS pair the "
            "independent product is more likely an over- than an under-estimate. "
            "Confirm against a haplotype table before costing the study.",
            "Registry-typed donors are not a random population sample; a "
            "commercial HLA-typed PBMC bank may already hold rare-allele donors, "
            "which changes the arithmetic from screening to availability.",
        ],
    }
    p = results_path("c5_donor_feasibility.json")
    with open(p, "w") as f:
        json.dump(out, f, indent=1)

    print(f"carrier frequency, weighted US/EU ({US_EU_SPLIT[0]:.0%} US : "
          f"{US_EU_SPLIT[1]:.0%} EU)\n")
    print(f"  {a_name:14s} {comb[a_name]*100:6.2f}%")
    print(f"  {b_name:14s} {comb[b_name]*100:6.2f}%")
    print(f"  {'either':14s} {comb['either']*100:6.2f}%")
    print(f"  {'both':14s} {comb['both_independent']*100:6.2f}%  "
          f"(independent approximation - see caveats)")
    print("\ndonors to HLA-type, in expectation, to obtain N carriers")
    print(f"  {'requirement':22s} " + "".join(f"N={n:<8d}" for n in WANT_DONORS))
    for key in (a_name, b_name, "either", "both_independent"):
        cells = "".join(f"{burden[key][str(n)]:<10.0f}" for n in WANT_DONORS)
        print(f"  {key:22s} {cells}")
    print(f"\nwrote {p}")


if __name__ == "__main__":
    main()
