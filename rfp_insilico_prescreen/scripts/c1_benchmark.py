#!/usr/bin/env python3
"""
C1 - Labelled benchmark for the two RFP alleles.

The RFP asks a CRO to build custom tetramers for HLA-A*32:01 and HLA-B*57:01.
Before an in-silico pre-screen can be allowed to say "don't synthesise that
complex", the predictor has to be shown to work ON THESE TWO ALLELES - not on
class I in general. A*32:01 in particular is not a well-studied allele, and a
predictor's average performance says nothing about a data-poor allele.

So: pull every IEDB MHC-ligand and binding record for these two alleles and
label each peptide.

  label 1   Positive in an MHC binding or ligand-elution assay
  label 0   Negative only
  excluded  both (counted and reported)

Restricted to 9-mers because that is what the RFP's uniblock peptides are and
what a tetramer is loaded with.

Output: results/c1_benchmark.tsv, results/c1_benchmark_summary.json
"""
import csv
import json
import os
import sys
import time
import urllib.error
import urllib.request
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import results_path  # noqa: E402

API = "https://query-api.iedb.org/mhc_search"
SELECT = ("linear_sequence,linear_sequence_length,mhc_allele_name,mhc_class,"
          "qualitative_measure,quantitative_measure,assay_names,structure_type,"
          "parent_source_antigen_source_org_name,pubmed_id")
ALLELES = ["HLA-A*32:01", "HLA-B*57:01"]
PAGE = 1000
LENGTH = 9
CANON = set("ACDEFGHIKLMNPQRSTVWY")
POSITIVE = {"Positive", "Positive-High", "Positive-Intermediate", "Positive-Low"}


def get(url, tries=4):
    """GET with retries; surfaces the response body, which is where this API
    explains a malformed query. A 4xx will not fix itself, so do not retry it."""
    last = None
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            body = ""
            try:
                body = e.read().decode()[:400]
            except Exception:
                pass
            last = f"{e} - {body}"
            if e.code < 500:
                break
            time.sleep(2 ** i)
        except Exception as e:
            last = e
            time.sleep(2 ** i)
    raise RuntimeError(f"{url} failed: {last}")


def fetch_allele(allele):
    """Every class I record IEDB holds for one allele.

    The '*' goes in raw: PostgREST rejects a percent-encoded '*' in an eq.
    filter. 'offset' requires 'order' or the API refuses the page.
    """
    rows, offset = [], 0
    while True:
        q = (f"{API}?select={SELECT}"
             f"&mhc_class=eq.I"
             f"&mhc_allele_name=eq.{allele}"
             f"&order=structure_id&limit={PAGE}&offset={offset}")
        page = get(q)
        rows.extend(page)
        if len(page) < PAGE:
            return rows
        offset += PAGE


def main():
    per_allele, labels = {}, {}
    raw_total = 0
    for allele in ALLELES:
        rows = fetch_allele(allele)
        raw_total += len(rows)
        seen = defaultdict(set)
        for r in rows:
            pep = (r.get("linear_sequence") or "").strip().upper()
            if len(pep) != LENGTH or set(pep) - CANON:
                continue
            q = r.get("qualitative_measure")
            if q is None:
                continue
            seen[pep].add("pos" if q in POSITIVE else "neg")
        pos = {p for p, s in seen.items() if s == {"pos"}}
        neg = {p for p, s in seen.items() if s == {"neg"}}
        amb = {p for p, s in seen.items() if len(s) > 1}
        per_allele[allele] = {
            "raw_records": len(rows), "peptides_9mer": len(seen),
            "positives": len(pos), "negatives": len(neg),
            "ambiguous_excluded": len(amb),
        }
        for p in pos:
            labels[(allele, p)] = 1
        for p in neg:
            labels[(allele, p)] = 0
        print(f"{allele:14s} {len(rows):6d} records -> {len(pos):5d} pos / "
              f"{len(neg):5d} neg 9-mers ({len(amb)} ambiguous dropped)")

    out = results_path("c1_benchmark.tsv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["allele", "peptide", "label"])
        for (allele, pep), lab in sorted(labels.items()):
            w.writerow([allele, pep, lab])

    summary = {"alleles": per_allele, "raw_records": raw_total,
               "labelled_pairs": len(labels), "length": LENGTH,
               "positive_terms": sorted(POSITIVE)}
    with open(results_path("c1_benchmark_summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    print(f"\nwrote {out}  ({len(labels)} labelled peptide-allele pairs)")


if __name__ == "__main__":
    main()
