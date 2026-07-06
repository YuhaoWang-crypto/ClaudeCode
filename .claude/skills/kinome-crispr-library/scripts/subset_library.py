#!/usr/bin/env python3
"""Subset a published pooled-CRISPR library to a target gene set and report coverage.

Generalized from a run that subset the Broad GPP Brunello library to the human
kinome (UniProt pkinfam): 512 kinases -> 506 covered (495 direct + 11 via HGNC
alias), 6 true gaps. Works for any gene list against any library file whose
columns include a gene symbol, an sgRNA spacer, and (optionally) an on-target score.

Inputs
------
  --genes FILE        one gene symbol per line (target set). OR:
  --pkinfam           fetch UniProt pkinfam and use all human protein kinases.
  --library FILE      the pooled library (TSV/CSV). Column names auto-detected;
                      override with --gene-col/--spacer-col/--score-col.
  --resolve-aliases   for unmatched symbols, query HGNC for prev/alias symbols
                      and retry against the library (needs network).

Outputs (to --outdir, default .)
  selected_guides.tsv   gene, library_symbol, spacer, score, GC, source
  coverage_stats.json   machine-readable coverage summary
  + a human-readable summary to stdout
"""
import argparse, csv, json, os, sys, urllib.request, statistics as st
from collections import defaultdict, Counter

COMP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
def gc(s): return 100 * sum(c in "GC" for c in s) / len(s) if s else 0.0


def fetch_pkinfam():
    url = "https://www.uniprot.org/docs/pkinfam.txt"
    txt = urllib.request.urlopen(url, timeout=60).read().decode("utf-8", "replace")
    import re
    genes = {}
    for line in txt.splitlines():
        m = re.match(r"^(\S+)\s+\S+_HUMAN\s*\(\s*([A-Z0-9]+)\s*\)", line)
        if m and m.group(1) != "-":
            genes[m.group(1)] = m.group(2)
    return set(genes)


def sniff(path):
    return "\t" if path.endswith((".tsv", ".txt")) else ","


def detect(cols, wants):
    low = {c.lower(): c for c in cols}
    for w in wants:
        for lc, orig in low.items():
            if w in lc:
                return orig
    return None


def hgnc_aliases(sym):
    url = f"https://rest.genenames.org/fetch/symbol/{sym}"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    try:
        d = json.load(urllib.request.urlopen(req, timeout=30))
        docs = d["response"]["docs"]
        if not docs:
            return []
        return (docs[0].get("prev_symbol", []) or []) + (docs[0].get("alias_symbol", []) or [])
    except Exception:
        return []


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--genes")
    src.add_argument("--pkinfam", action="store_true")
    p.add_argument("--library", required=True)
    p.add_argument("--gene-col"); p.add_argument("--spacer-col"); p.add_argument("--score-col")
    p.add_argument("--guides-per-gene", type=int, default=0, help="0 = keep all available")
    p.add_argument("--resolve-aliases", action="store_true")
    p.add_argument("--outdir", default=".")
    a = p.parse_args()

    target = fetch_pkinfam() if a.pkinfam else set(l.strip() for l in open(a.genes) if l.strip())
    print(f"target genes: {len(target)}")

    lib = defaultdict(list)
    with open(a.library) as f:
        r = csv.DictReader(f, delimiter=sniff(a.library))
        gcol = a.gene_col or detect(r.fieldnames, ["gene symbol", "gene", "target"])
        scol = a.spacer_col or detect(r.fieldnames, ["sgrna target sequence", "spacer", "sequence", "sgrna"])
        ocol = a.score_col or detect(r.fieldnames, ["rule set", "score", "efficacy"])
        if not gcol or not scol:
            sys.exit(f"could not detect gene/spacer columns in {r.fieldnames}; pass --gene-col/--spacer-col")
        for row in r:
            g = (row[gcol] or "").strip()
            if not g:
                continue
            try:
                sc = float(row[ocol]) if ocol and row.get(ocol) not in (None, "") else float("nan")
            except ValueError:
                sc = float("nan")
            lib[g].append((row[scol].strip().upper(), sc))
    libset = set(lib)

    alias_map = {}
    if a.resolve_aliases:
        for s in sorted(target - libset):
            for cand in hgnc_aliases(s):
                if cand in libset:
                    alias_map[s] = cand
                    break

    os.makedirs(a.outdir, exist_ok=True)
    rows = []; direct = via_alias = 0; gaps = []; scores = []
    for sym in sorted(target):
        libsym = sym if sym in libset else alias_map.get(sym)
        if not libsym:
            gaps.append(sym); continue
        direct += sym in libset and libsym == sym
        via_alias += libsym != sym
        guides = sorted(lib[libsym], key=lambda t: -(t[1] if t[1] == t[1] else -1))
        if a.guides_per_gene:
            guides = guides[:a.guides_per_gene]
        for i, (sp, sc) in enumerate(guides, 1):
            rows.append(dict(gene=sym, library_symbol=libsym, rank=i, spacer=sp,
                             score=("" if sc != sc else f"{sc:.4f}"), GC=f"{gc(sp):.0f}",
                             source=("direct" if libsym == sym else "alias")))
            if sc == sc:
                scores.append(sc)

    with open(os.path.join(a.outdir, "selected_guides.tsv"), "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)

    covered = len(target) - len(gaps)
    stats = dict(target=len(target), covered=covered, direct=direct, via_alias=via_alias,
                 gaps=gaps, guides=len(rows),
                 median_score=round(st.median(scores), 3) if scores else None,
                 frac_high=round(sum(x >= 0.6 for x in scores) / len(scores), 3) if scores else None)
    json.dump(stats, open(os.path.join(a.outdir, "coverage_stats.json"), "w"), indent=2)
    print(f"covered: {covered}/{len(target)} ({100*covered/len(target):.1f}%) "
          f"= {direct} direct + {via_alias} via alias")
    print(f"gaps ({len(gaps)}): {', '.join(gaps) if gaps else '(none)'}")
    if scores:
        print(f"on-target median {stats['median_score']} | >=0.6: {100*stats['frac_high']:.0f}%")
    print(f"wrote {a.outdir}/selected_guides.tsv ({len(rows)} guides) + coverage_stats.json")


if __name__ == "__main__":
    main()
