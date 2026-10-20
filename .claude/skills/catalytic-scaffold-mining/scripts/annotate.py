#!/usr/bin/env python3
"""Stage 5/6 support: add lineage and prior-characterisation to a candidate table.

Two questions decide whether a structural hit is interesting:

  1. **Where on the tree is it?**  A hit in a lineage that already has a
     characterised member of this family is confirmation; a hit in a lineage
     that has none is the actual lead.
  2. **Is it already known?**  A pipeline that "discovers" poly(A) polymerase
     has discovered nothing. Separating annotated enzymes from genuinely
     uncharacterised proteins is what turns a hit list into a shortlist.

Lineage comes from the UniProt taxonomy API, cached on disk. Characterisation
is keyword-matched against the target description, which is a blunt instrument
and is reported as such -- `unknown` means "no informative annotation", not
"novel".

    python3 annotate.py out/urt1_proteome.candidates.tsv \
        --out out/urt1_proteome.annotated.tsv
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys

# Lineage markers, most specific first; first match wins.
GROUPS = [
    ("Viridiplantae", "Viridiplantae"),
    ("Metazoa", "Metazoa"),
    ("Fungi", "Fungi"),
    ("Discoba", "Discoba"),              # trypanosomes, Leishmania
    ("Amoebozoa", "Amoebozoa"),
    ("Alveolata", "Alveolata"),          # apicomplexa, ciliates, dinoflagellates
    ("Stramenopiles", "Stramenopiles"),
    ("Rhizaria", "Rhizaria"),
    ("Haptista", "Haptista"),
    ("Cryptophyceae", "Cryptophyta"),
    ("Rhodophyta", "Rhodophyta"),
    ("Metamonada", "Metamonada"),
    ("Eukaryota", "other Eukaryota"),
    ("Bacteria", "Bacteria"),
    ("Archaea", "Archaea"),
    ("Viruses", "Viruses"),
]

# Known members of the polymerase-beta-like nucleotidyltransferase world.
# Order matters: more specific patterns first.
KNOWN_CLASSES = [
    ("cGAS / cGLR", ("cyclic gmp-amp synthase", "cgas", "cgas-like receptor",
                     "cyclic gmp-amp")),
    ("bacterial CD-NTase", ("cyclic dipyrimidine nucleotide synthase",
                            "cyclic dipurine nucleotide synthase",
                            "cd-ntase", "dnc v", "dncv")),
    ("OAS", ("2'-5'-oligoadenylate synthase", "oligoadenylate synthase")),
    ("CCA-adding enzyme", ("cca-adding", "trna nucleotidyltransferase",
                           "trna adenylyltransferase",
                           "trna cytidylyltransferase")),
    ("TUTase / uridylyltransferase", ("uridylyltransferase", "uridylyl transferase",
                                      "tutase", "poly(u) polymerase",
                                      "poly u polymerase")),
    ("poly(A) polymerase", ("poly(a) polymerase", "poly a polymerase",
                            "poly(a) rna polymerase", "polyadenylate polymerase",
                            "polynucleotide adenylyltransferase",
                            "pap/25a-associated", "pap-associated", "pap ")),
    ("TENT family", ("terminal nucleotidyltransferase", "tent",
                     "terminal uridylyl")),
    ("Mab-21 family", ("mab-21", "mab21")),
    ("antibiotic / aminoglycoside NTase", ("aminoglycoside", "lincosamide",
                                           "streptomycin adenylyl")),
    ("glycosyl/ribitol transferase", ("ribitol-5-phosphate", "fukutin")),
    ("kanamycin/other transferase", ("nucleotidyltransferase domain",)),
]

# Substrings that mean "this annotation tells you nothing about the chemistry".
# Checked only after every KNOWN_CLASSES pattern has failed.
UNINFORMATIVE = (
    "uncharacterized", "putative", "hypothetical", "duf", "domain-containing",
    "-like protein", "unknown", "predicted", "expressed protein",
)


def classify_description(desc):
    d = (desc or "").lower().strip()
    for label, pats in KNOWN_CLASSES:
        for p in pats:
            if p in d:
                return label
    if not d:
        return "unknown"
    if any(u in d for u in UNINFORMATIVE):
        return "unknown"
    return "other annotated"


def fetch_lineages(tax_ids, cache_path):
    cache = {}
    if os.path.exists(cache_path):
        cache = json.load(open(cache_path))
    todo = [t for t in tax_ids if t and str(t) not in cache]
    for i in range(0, len(todo), 40):
        chunk = todo[i:i + 40]
        query = "+OR+".join(f"(tax_id:{t})" for t in chunk)
        url = (f"https://rest.uniprot.org/taxonomy/search?query={query}"
               f"&size=500&fields=id,scientific_name,lineage")
        try:
            out = subprocess.run(["curl", "-sS", "--max-time", "90", url],
                                 capture_output=True, text=True, check=True)
            data = json.loads(out.stdout)
        except Exception as exc:                        # noqa: BLE001
            print(f"  taxonomy fetch failed for chunk {i}: {exc}", file=sys.stderr)
            continue
        for r in data.get("results", []):
            names = [x.get("scientificName") for x in r.get("lineage", [])]
            cache[str(r.get("taxonId"))] = names
        print(f"  taxonomy {min(i + 40, len(todo))}/{len(todo)}", flush=True)
    json.dump(cache, open(cache_path, "w"))
    return cache


def group_of(lineage):
    if not lineage:
        return "unknown"
    s = set(lineage)
    for marker, label in GROUPS:
        if marker in s:
            return label
    return "unknown"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--out", required=True)
    ap.add_argument("--cache", default="taxonomy_cache.json")
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.table), delimiter="\t"))
    tax_ids = sorted({r.get("tax_id") for r in rows if r.get("tax_id")})
    print(f"{len(rows)} rows, {len(tax_ids)} distinct taxa")
    cache = fetch_lineages(tax_ids, args.cache)

    for r in rows:
        lin = cache.get(str(r.get("tax_id")), [])
        r["lineage_group"] = group_of(lin)
        r["lineage"] = "; ".join(lin[1:8])
        r["known_class"] = classify_description(r.get("description"))

    fields = list(rows[0].keys()) if rows else []
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
