#!/usr/bin/env python3
"""
M14 - Domain attribution for a chimeric ligand.

The stock pipeline scores a ligand as one object. That is the right unit when
the ligand is one domain, and the wrong unit when it is a fusion of a novel
binding domain onto a scaffold that is *already* a qualified affinity ligand -
because then the headline fold-change over the scaffold benchmark is partly a
comparison of the molecule with a piece of itself.

7071-A08 is exactly that case: residues 73-130 are a Protein A Z-domain variant
87.9% identical to ProteinA_Z, the batch's calibration anchor. So this module
splits the test article's risk by domain and asks the question the split makes
possible:

  of the epitope content in the Z-domain half, how much is *shared* with the
  clinically qualified Protein A leachate, and how much is *novel* - created by
  the engineering that distinguishes this scaffold from the wild-type domain?

Method
  1. Assign each M5 epitope to a domain by its core start position.
  2. Report, per domain, the absolute contribution to whole-protein pIRS and the
     domain-local epitope density (pIRS computed per 100 aa of that domain).
  3. Align the scaffold domain to the benchmark, map every scaffold epitope core
     onto benchmark coordinates, and classify it:
       shared  - the benchmark carries a predicted epitope overlapping >= 9
                 residues of the same window
       novel   - it does not
  4. Mark which novel cores contain a residue that differs between scaffold and
     benchmark, i.e. which are attributable to the engineering rather than to
     the prediction being noisy at the margin.

Output: results/m14_domain_attribution.json, results/m14_domain_epitopes.tsv
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import load_config, read_fasta, data_path, results_path  # noqa: E402

MIN_OVERLAP = 9          # same rule M6 uses to call two epitope windows the same
MIN_HUMAN_ALN = 30       # shortest alignment allowed to count as a human homologue

# Domain boundaries are a property of the construct, not of the pipeline, so
# they are declared here with the evidence for each call.
DOMAINS = {
    "dsRNA_7071_A08": [
        ("dsRBD", 1, 72, "double-stranded-RNA-binding domain fold"),
        ("ProteinA_Z_variant", 73, 130, "Protein A Z-domain variant"),
        ("C_term_Cys", 131, 132, "terminal Cys for resin coupling"),
    ],
}
SCAFFOLD = {"dsRNA_7071_A08": ("ProteinA_Z_variant", "ProteinA_Z")}


def tsv(name):
    with open(results_path(name)) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def nearest_human(sub, top=3):
    """Best local alignments of a domain against human Swiss-Prot.

    Seeded on shared 5-mers so the 20k-protein sweep stays cheap; the claim this
    supports ("this domain is not a human sequence") is a claim about the BEST
    hit, so a seed that can only miss remote homologues is acceptable here and
    is stated as such in the report.
    """
    hs = read_fasta(data_path("human_sprot.fasta"))
    seeds = {sub[i:i + 5] for i in range(len(sub) - 4)}
    cand = [(sum(1 for i in range(len(v) - 4) if v[i:i + 5] in seeds), k, v)
            for k, v in hs.items()]
    cand = sorted((c for c in cand if c[0] >= 2), reverse=True)[:400]
    hits = []
    for _, k, v in cand:
        a = align(sub, v)
        ident = sum(1 for x, y in zip(a[0], a[1]) if x == y and x != "-")
        # Rank on alignment score, never on percent identity: a local aligner
        # will happily return a 6-residue 100% match, and sorting on percent
        # puts that above a genuine 60-residue domain-level homology.
        hits.append({"entry": k.split("|")[-1] if "|" in k else k,
                     "accession": k.split("|")[1] if k.count("|") >= 2 else k,
                     "score": float(a.score),
                     "identity": ident, "aligned_length": len(a[0]),
                     "identity_pct": round(100.0 * ident / len(a[0]), 1)})
    hits = [h for h in hits if h["aligned_length"] >= MIN_HUMAN_ALN]
    hits.sort(key=lambda h: -h["score"])
    return hits[:top], len(cand)


def align(a, b):
    """Local BLOSUM62 alignment -> (aligned_a, aligned_b, a_start, b_start)."""
    from Bio import Align
    from Bio.Align import substitution_matrices
    al = Align.PairwiseAligner()
    al.mode = "local"
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score = -11, -1
    best = al.align(a, b)[0]
    return best


def main():
    cfg = load_config()
    seqs = read_fasta(data_path("sequences.fasta"))
    eps = tsv("m5_epitopes.tsv")
    anchor = cfg["benchmarks"]["anchor_low"]

    out = {"min_overlap": MIN_OVERLAP, "anchor": anchor, "articles": {}}
    rows_out = []

    for sid, domains in DOMAINS.items():
        if sid not in seqs:
            continue
        seq = seqs[sid]
        L = len(seq)
        mine = [e for e in eps if e["id"] == sid]
        whole_pirs = 100.0 / L * sum(float(e["pop_presenting"]) *
                                     float(e["tolerance_weight"]) for e in mine)

        per_domain = []
        for name, lo, hi, why in domains:
            inside = [e for e in mine if lo <= int(e["pos"]) <= hi - 8]
            dl = hi - lo + 1
            contrib = 100.0 / L * sum(float(e["pop_presenting"]) *
                                      float(e["tolerance_weight"]) for e in inside)
            density = (100.0 / dl * sum(float(e["pop_presenting"]) *
                                        float(e["tolerance_weight"]) for e in inside))
            per_domain.append({
                "domain": name, "start": lo, "end": hi, "length": dl,
                "evidence": why,
                "n_epitopes": len(inside),
                "n_foreign": sum(1 for e in inside
                                 if e["tolerance_class"] == "foreign"),
                "max_promiscuity": max([int(e["n_sb_alleles"]) for e in inside],
                                       default=0),
                "pirs_contribution": round(contrib, 4),
                "share_of_pirs": round(contrib / whole_pirs, 4) if whole_pirs else 0.0,
                "pirs_density_per_100aa": round(density, 4),
            })

        art = {"length": L, "pIRS": round(whole_pirs, 4), "domains": per_domain}

        # Is the novel domain actually foreign? Computed, not asserted.
        scaffold_name = SCAFFOLD.get(sid, (None, None))[0]
        for dom in per_domain:
            if dom["domain"] == scaffold_name or dom["length"] < 30:
                continue
            hits, n = nearest_human(seq[dom["start"] - 1:dom["end"]])
            dom["nearest_human"] = hits
            dom["human_candidates_aligned"] = n

        # ---- scaffold vs benchmark -------------------------------------
        if sid in SCAFFOLD and SCAFFOLD[sid][1] in seqs:
            dname, bench = SCAFFOLD[sid]
            lo, hi = next((d[1], d[2]) for d in domains if d[0] == dname)
            sub = seq[lo - 1:hi]
            bseq = seqs[bench]
            a = align(sub, bseq)
            # index maps: scaffold position (1-based in full ligand) -> benchmark position
            s2b, subs = {}, []
            sa, sb = a.aligned                    # blocks of (start, end) pairs
            # a.aligned yields numpy ints; cast so the block indices stay
            # JSON-serialisable all the way through
            sa = [(int(x), int(y)) for x, y in sa]
            sb = [(int(x), int(y)) for x, y in sb]
            for (s0, s1), (b0, b1) in zip(sa, sb):
                for k in range(s1 - s0):
                    sp, bp = s0 + k, b0 + k
                    s2b[lo + sp] = bp + 1
                    if sub[sp] != bseq[bp]:
                        subs.append({"ligand_pos": lo + sp, "ligand_aa": sub[sp],
                                     "benchmark_pos": bp + 1, "benchmark_aa": bseq[bp]})
            ident = sum(1 for (s0, s1), (b0, b1) in zip(sa, sb)
                        for k in range(s1 - s0) if sub[s0 + k] == bseq[b0 + k])
            aligned_len = sum(s1 - s0 for s0, s1 in sa)

            bench_eps = [e for e in eps if e["id"] == bench]
            sub_pos = {p["ligand_pos"] for p in subs}
            calls = []
            for e in mine:
                p = int(e["pos"])
                if not (lo <= p <= hi - 8):
                    continue
                win = [s2b.get(p + k) for k in range(9)]
                mapped = [w for w in win if w]
                shared, partner = False, None
                if mapped:
                    w0, w1 = min(mapped), max(mapped)
                    for b in bench_eps:
                        bp = int(b["pos"])
                        if min(bp + 8, w1) - max(bp, w0) + 1 >= MIN_OVERLAP:
                            shared, partner = True, b
                            break
                touched = sorted(sub_pos & set(range(p, p + 9)))
                call = {
                    "core": e["core"], "ligand_pos": p,
                    "benchmark_pos": min(mapped) if mapped else None,
                    "best_el_rank": float(e["best_el_rank"]),
                    "n_sb_alleles": int(e["n_sb_alleles"]),
                    "pop_presenting": float(e["pop_presenting"]),
                    "tolerance_class": e["tolerance_class"],
                    "shared_with_benchmark": shared,
                    "benchmark_core": partner["core"] if partner else None,
                    "substituted_positions_in_core": touched,
                    "engineering_attributable": bool(touched) and not shared,
                }
                calls.append(call)
                rows_out.append({"id": sid, "domain": dname, **call})

            novel = [c for c in calls if not c["shared_with_benchmark"]]
            art["scaffold_vs_benchmark"] = {
                "domain": dname, "benchmark": bench,
                "identity": f"{ident}/{aligned_len}",
                "identity_pct": round(100.0 * ident / aligned_len, 1) if aligned_len else None,
                "n_substitutions": len(subs),
                "substitutions": subs,
                "n_scaffold_epitopes": len(calls),
                "n_shared": sum(1 for c in calls if c["shared_with_benchmark"]),
                "n_novel": len(novel),
                "n_novel_engineering_attributable":
                    sum(1 for c in novel if c["engineering_attributable"]),
                "novel_pirs_contribution": round(
                    100.0 / L * sum(c["pop_presenting"] for c in novel), 4),
                "epitopes": calls,
            }

        # epitopes outside any scaffold, for the tsv
        for e in mine:
            p = int(e["pos"])
            dn = next((d[0] for d in domains if d[1] <= p <= d[2] - 8), "inter-domain")
            if sid in SCAFFOLD and dn == SCAFFOLD[sid][0]:
                continue
            rows_out.append({"id": sid, "domain": dn, "core": e["core"],
                             "ligand_pos": p, "benchmark_pos": None,
                             "best_el_rank": float(e["best_el_rank"]),
                             "n_sb_alleles": int(e["n_sb_alleles"]),
                             "pop_presenting": float(e["pop_presenting"]),
                             "tolerance_class": e["tolerance_class"],
                             "shared_with_benchmark": None, "benchmark_core": None,
                             "substituted_positions_in_core": [],
                             "engineering_attributable": None})
        out["articles"][sid] = art

    p = results_path("m14_domain_attribution.json")
    with open(p, "w") as f:
        json.dump(out, f, indent=1)

    cols = ["id", "domain", "core", "ligand_pos", "benchmark_pos", "best_el_rank",
            "n_sb_alleles", "pop_presenting", "tolerance_class",
            "shared_with_benchmark", "benchmark_core",
            "substituted_positions_in_core", "engineering_attributable"]
    t = results_path("m14_domain_epitopes.tsv")
    with open(t, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in sorted(rows_out, key=lambda r: (r["id"], r["ligand_pos"])):
            w.writerow({c: r.get(c) for c in cols})

    # ---- console ----------------------------------------------------
    for sid, art in out["articles"].items():
        print(f"\n{sid}  {art['length']} aa  pIRS {art['pIRS']:.2f}")
        print(f"{'domain':22s} {'range':>9s} {'ep':>3s} {'fgn':>4s} {'prom':>5s} "
              f"{'pIRS':>6s} {'share':>6s} {'per100aa':>9s}")
        print("-" * 76)
        for d in art["domains"]:
            print(f"{d['domain']:22s} {str(d['start'])+'-'+str(d['end']):>9s} "
                  f"{d['n_epitopes']:3d} {d['n_foreign']:4d} {d['max_promiscuity']:5d} "
                  f"{d['pirs_contribution']:6.2f} {d['share_of_pirs']*100:5.1f}% "
                  f"{d['pirs_density_per_100aa']:9.2f}")
        s = art.get("scaffold_vs_benchmark")
        if s:
            print(f"\n{s['domain']} vs {s['benchmark']}: {s['identity']} "
                  f"({s['identity_pct']}%), {s['n_substitutions']} substitutions")
            print(f"  scaffold epitopes {s['n_scaffold_epitopes']}: "
                  f"{s['n_shared']} shared with the qualified benchmark, "
                  f"{s['n_novel']} novel "
                  f"({s['n_novel_engineering_attributable']} containing a substituted residue)")
            print(f"  novel epitope contribution to whole-ligand pIRS: "
                  f"{s['novel_pirs_contribution']:.2f}")
            for c in s["epitopes"]:
                tag = ("shared" if c["shared_with_benchmark"]
                       else ("NOVEL*" if c["engineering_attributable"] else "novel"))
                print(f"    {c['core']}  lig {c['ligand_pos']:3d}  "
                      f"bench {c['benchmark_pos']}  EL {c['best_el_rank']:.2f}  "
                      f"{c['n_sb_alleles']} DR  pop {c['pop_presenting']*100:4.1f}%  "
                      f"{tag}"
                      + (f"  subs@{c['substituted_positions_in_core']}"
                         if c["substituted_positions_in_core"] else ""))
    print(f"\nwrote {p}\nwrote {t}")


if __name__ == "__main__":
    main()
