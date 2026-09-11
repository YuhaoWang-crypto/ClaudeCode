#!/usr/bin/env python3
"""
M10 - Build a labelled benchmark of measured CD4 T-cell responses.

Everything upstream of this module is a *claim* about accuracy. The EL+BA
consensus rule was introduced on the argument that eluted-ligand scoring
over-calls - which is a reasonable argument and no evidence at all. Without
ground truth there is no way to know whether the rule raises specificity,
lowers sensitivity by more than it gains, or does nothing.

So: pull every HLA-DR-restricted T-cell assay result IEDB holds for the
molecules in the panel, in human hosts, and label each (peptide, allele) pair
by what was actually measured.

  label 1   at least one Positive T-cell assay for that peptide on that DR
            molecule
  label 0   only Negative assays for that pair
  excluded  pairs with both (genuinely contested; counted and reported)

Three biases are handled explicitly rather than ignored, because each one can
manufacture accuracy that does not exist:

  * **Redundancy.** Overlapping peptides from the same protein and study are
    not independent observations. Peptides sharing any 9-mer are clustered and
    the evaluation is reported per-cluster as well as per-pair.
  * **Source-organism confounding.** Positives skew toward pathogens and
    allergens, negatives toward self and toward peptides someone expected to
    fail. A predictor that merely recognises "foreign" would score well for the
    wrong reason, so the self/non-self split of each label is recorded and the
    metrics are reported stratified.
  * **Assay-context loss.** IEDB negatives include peptides that were never
    going to respond in that donor for reasons unrelated to binding. That
    inflates apparent specificity, so the negative set is reported with its
    provenance and never presented as a clean non-binder set.

Output: results/m10_benchmark.tsv, results/m10_benchmark_summary.json
"""
import csv
import json
import os
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import results_path, data_path  # noqa: E402

API = "https://query-api.iedb.org/tcell_search"
SELECT = ("linear_sequence,linear_sequence_length,mhc_restriction,"
          "mhc_allele_resolution,qualitative_measure,source_organism_name,"
          "parent_source_antigen_source_org_name,pubmed_id,structure_type")
PAGE = 1000
MIN_LEN, MAX_LEN = 13, 25
CANON = set("ACDEFGHIKLMNPQRSTVWY")


def get(url, tries=4):
    """GET with retries. Surfaces the response body on an HTTP error - this API
    explains what is wrong with a query there, and swallowing it costs a lot of
    time (the paging rule below was found exactly that way)."""
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
                break            # a malformed query will not fix itself
            time.sleep(2 ** i)
        except Exception as e:
            last = e
            time.sleep(2 ** i)
    raise RuntimeError(f"{url} failed: {last}")


def fetch_allele(allele):
    """Every human-host class II T-cell record for one DR molecule."""
    rows, offset = [], 0
    while True:
        q = (f"{API}?select={SELECT}"
             f"&mhc_class=eq.II"
             f"&host_organism_iri=eq.NCBITaxon:9606"
             f"&structure_type=eq.Linear%20peptide"
             # PostgREST rejects a percent-encoded '*' in an eq. filter, so the
             # allele name goes in raw; only the space in the type filter is encoded
             f"&mhc_restriction=eq.{allele}"
             # this API refuses an offset without an order, to keep paging stable
             f"&order=tcell_id&limit={PAGE}&offset={offset}")
        page = get(q)
        rows.extend(page)
        if len(page) < PAGE:
            return rows
        offset += PAGE


def is_self(row):
    org = (row.get("parent_source_antigen_source_org_name")
           or row.get("source_organism_name") or "")
    return "Homo sapiens" in org


def cluster_by_shared_9mer(peptides):
    """Union-find over peptides sharing any 9-mer; returns {peptide: cluster_id}."""
    parent = {p: p for p in peptides}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    index = defaultdict(list)
    for p in peptides:
        for i in range(len(p) - 8):
            index[p[i:i + 9]].append(p)
    for members in index.values():
        for m in members[1:]:
            union(members[0], m)
    roots = {}
    out = {}
    for p in peptides:
        r = find(p)
        out[p] = roots.setdefault(r, len(roots))
    return out


def main():
    with open(results_path("m2_panel_alleles.txt")) as f:
        panel = [l.strip() for l in f if l.strip()]

    pair = {}      # (peptide, allele) -> record
    n_raw = 0
    for allele in panel:
        rows = fetch_allele(allele)
        n_raw += len(rows)
        kept = 0
        for r in rows:
            seq = (r.get("linear_sequence") or "").strip().upper()
            if not (MIN_LEN <= len(seq) <= MAX_LEN) or not set(seq) <= CANON:
                continue
            # allele-level restriction only - "DR" or "DRB1*04" cannot be
            # assigned to a specific molecule in the panel
            if r.get("mhc_restriction") != allele:
                continue
            key = (seq, allele)
            rec = pair.setdefault(key, {
                "peptide": seq, "allele": allele, "n_pos": 0, "n_neg": 0,
                "self_pos": 0, "self_neg": 0, "pmids": set(), "orgs": set()})
            q = (r.get("qualitative_measure") or "")
            slf = is_self(r)
            if q.startswith("Positive"):
                rec["n_pos"] += 1
                rec["self_pos"] += slf
            elif q == "Negative":
                rec["n_neg"] += 1
                rec["self_neg"] += slf
            else:
                continue
            if r.get("pubmed_id"):
                rec["pmids"].add(str(r["pubmed_id"]))
            org = (r.get("parent_source_antigen_source_org_name")
                   or r.get("source_organism_name") or "")
            if org:
                rec["orgs"].add(org)
            kept += 1
        print(f"  {allele:20s} {len(rows):6d} records -> {kept:5d} usable", flush=True)

    labelled, ambiguous = [], 0
    for rec in pair.values():
        if rec["n_pos"] and rec["n_neg"]:
            ambiguous += 1
            continue
        if not rec["n_pos"] and not rec["n_neg"]:
            continue
        label = 1 if rec["n_pos"] else 0
        rec["label"] = label
        rec["n_assays"] = rec["n_pos"] + rec["n_neg"]
        rec["is_self"] = bool(rec["self_pos"] or rec["self_neg"])
        labelled.append(rec)

    clusters = cluster_by_shared_9mer(sorted({r["peptide"] for r in labelled}))
    for r in labelled:
        r["cluster"] = clusters[r["peptide"]]

    labelled.sort(key=lambda r: (r["allele"], -r["label"], r["peptide"]))
    cols = ["peptide", "allele", "label", "n_pos", "n_neg", "n_assays",
            "is_self", "cluster", "n_studies", "organisms"]
    with open(results_path("m10_benchmark.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for r in labelled:
            w.writerow([r["peptide"], r["allele"], r["label"], r["n_pos"],
                        r["n_neg"], r["n_assays"], int(r["is_self"]),
                        r["cluster"], len(r["pmids"]),
                        "; ".join(sorted(r["orgs"]))[:200]])

    pos = [r for r in labelled if r["label"] == 1]
    neg = [r for r in labelled if r["label"] == 0]
    summary = {
        "raw_records": n_raw,
        "labelled_pairs": len(labelled),
        "positives": len(pos),
        "negatives": len(neg),
        "ambiguous_excluded": ambiguous,
        "distinct_peptides": len({r["peptide"] for r in labelled}),
        "clusters": len({r["cluster"] for r in labelled}),
        "alleles_with_data": len({r["allele"] for r in labelled}),
        "self_fraction_positives": round(sum(r["is_self"] for r in pos) / max(len(pos), 1), 3),
        "self_fraction_negatives": round(sum(r["is_self"] for r in neg) / max(len(neg), 1), 3),
        "per_allele": {a: {
            "pos": sum(1 for r in pos if r["allele"] == a),
            "neg": sum(1 for r in neg if r["allele"] == a)} for a in panel},
    }
    with open(results_path("m10_benchmark_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\n{n_raw} raw records -> {len(labelled)} labelled (peptide, allele) pairs")
    print(f"  positives {len(pos)}   negatives {len(neg)}   "
          f"ambiguous excluded {ambiguous}")
    print(f"  {summary['distinct_peptides']} distinct peptides in "
          f"{summary['clusters']} 9-mer-sharing clusters")
    print(f"  self-derived: {summary['self_fraction_positives']*100:.0f}% of positives, "
          f"{summary['self_fraction_negatives']*100:.0f}% of negatives")
    print("\nper allele (pos/neg):")
    for a in panel:
        s = summary["per_allele"][a]
        if s["pos"] or s["neg"]:
            print(f"  {a:20s} {s['pos']:5d} / {s['neg']:5d}")


if __name__ == "__main__":
    main()
