#!/usr/bin/env python3
"""Design de-novo SpCas9 (NGG) spacers for genes from their Ensembl canonical CDS.

Generalized from the run that filled the 5 kinome genes missing from Brunello
(CDK8, GRK1, PRAG1, RSKR, STK38). For each gene symbol: resolve the Ensembl
canonical transcript, pull its CDS, scan both strands for NGG protospacers in an
early-exon window, apply the standard hard filters, and pick a spread of guides.

Filters (reject a spacer if any fail): poly-T (TTTT), GC outside 30-70%,
homopolymer run >=5, internal BsmBI/BbsI site. Cut sites forced >= 8 nt apart.

Usage:
  python3 design_guides_from_cds.py --genes gaps.txt --out gap_guides.tsv \
      --guides-per-gene 6 --window 0.05 0.65
"""
import argparse, csv, json, time, urllib.request

COMP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
def rc(s): return "".join(COMP[c] for c in reversed(s))
def gc(s): return 100 * sum(c in "GC" for c in s) / len(s) if s else 0.0
BAD = ("CGTCTC", "GAGACG", "GAAGAC", "GTCTTC")  # BsmBI + BbsI (+rc)


def ok(sp):
    if "TTTT" in sp: return False
    if not (30 <= gc(sp) <= 70): return False
    if any(c * 5 in sp for c in "ACGT"): return False
    if any(b in sp for b in BAD): return False
    return True


def ensembl(url):
    req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
    return json.load(urllib.request.urlopen(req, timeout=40))


def canonical_cds(sym):
    xr = ensembl(f"https://rest.ensembl.org/xrefs/symbol/homo_sapiens/{sym}?content-type=application/json")
    ensg = next((x["id"] for x in xr if x["id"].startswith("ENSG")), None)
    if not ensg:
        return None, None
    look = ensembl(f"https://rest.ensembl.org/lookup/id/{ensg}?expand=1;content-type=application/json")
    trs = look.get("Transcript", [])
    canon = next((t for t in trs if t.get("is_canonical") == 1), trs[0] if trs else None)
    if not canon:
        return None, None
    enst = canon["id"]
    d = ensembl(f"https://rest.ensembl.org/sequence/id/{enst}?type=cds;content-type=application/json")
    return enst, d["seq"].upper()


def design(sym, cds, enst, n, lo_f, hi_f):
    L = len(cds); lo, hi = int(L * lo_f), int(L * hi_f)
    cands = []
    for i in range(0, L - 23):                        # sense PAM
        if cds[i + 21:i + 23] == "GG":
            sp = cds[i:i + 20]
            if lo <= i + 17 <= hi and ok(sp):
                cands.append((i + 17, "+", sp, cds[i + 20:i + 23]))
    for j in range(0, L - 23):                        # antisense PAM (CCN on sense)
        if cds[j:j + 2] == "CC":
            sp = rc(cds[j + 3:j + 23])
            if len(sp) == 20 and lo <= j + 6 <= hi and ok(sp):
                cands.append((j + 6, "-", sp, rc(cds[j:j + 3])))
    cands = sorted(set(cands), key=lambda c: (-(1 - abs(gc(c[2]) - 50) / 50), c[0]))
    chosen, used = [], []
    for cut, strand, sp, pam in cands:
        if all(abs(cut - u) >= 8 for u in used):
            chosen.append((cut, strand, sp, pam)); used.append(cut)
        if len(chosen) == n:
            break
    return sorted(chosen, key=lambda c: c[0]), L, len(cands)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--genes", required=True, help="one gene symbol per line")
    p.add_argument("--out", required=True)
    p.add_argument("--guides-per-gene", type=int, default=6)
    p.add_argument("--window", nargs=2, type=float, default=[0.05, 0.65], metavar=("LO", "HI"))
    a = p.parse_args()
    genes = [l.strip() for l in open(a.genes) if l.strip()]
    rows = []
    for sym in genes:
        try:
            enst, cds = canonical_cds(sym)
        except Exception as e:
            print(f"{sym}: ERROR {e}"); continue
        if not cds:
            print(f"{sym}: no canonical CDS found"); continue
        chosen, L, ncand = design(sym, cds, enst, a.guides_per_gene, *a.window)
        print(f"{sym}: CDS {L} nt, {ncand} candidates, chose {len(chosen)}")
        for rk, (cut, strand, sp, pam) in enumerate(chosen, 1):
            rows.append(dict(gene=sym, ensembl_transcript=enst, rank=rk, spacer=sp, PAM=pam,
                             strand=strand, cut_pos_in_cds=cut, cds_frac=f"{cut/L:.2f}",
                             GC=f"{gc(sp):.0f}", rs2_score="NA", offtarget="pending",
                             source="de-novo (Ensembl CDS scan)"))
        time.sleep(0.3)
    with open(a.out, "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    print(f"wrote {a.out}: {len(rows)} guides for {len(genes)} genes")


if __name__ == "__main__":
    main()
