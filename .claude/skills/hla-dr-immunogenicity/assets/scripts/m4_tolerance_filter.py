#!/usr/bin/env python3
"""
M4 - Self / pre-existing-tolerance filter, with an empirical null.

The single largest source of false positives in a naive HLA-DR screen of an
antibody-derived ligand: most predicted binders sit in *framework* regions
whose 9-mer cores are near-identical to human immunoglobulin V germline. Those
cores are seen by a repertoire that has been negatively selected against them,
so counting them as immunogenic risk inflates the score of every VHH, scFv and
Fab-derived ligand equally and destroys the ability to rank them.

Each predicted 9-mer core is compared against every 9-mer of UniProt
Swiss-Prot *Homo sapiens* and classified by best identity:

  9/9   exact_human_9mer          weight 0     - the peptide itself is self
  8/9   near_human_9mer           weight 0.35  - one substitution from self
  <=7/9 foreign                   weight 1.0

The proteome pass uses a one-mismatch index, so it resolves 9/9 and 8/9 and
reports everything else as "foreign" - it deliberately does not claim a 7/9 or
6/9 number it cannot compute.

The 8/9 cut is not a guess. A 5-of-9 "TCR-face" pattern - the obvious
JanusMatrix-style shortcut - matches the human proteome by chance several
times per query (~20^-5 x 1.1e7 9-mers), so it flags essentially everything
and carries no information. At 8/9 the chance expectation is <0.1 hits per
query. The module measures this directly rather than asserting it: the same
proteome pass scores a **shuffled-sequence null** (each ligand's residues
permuted, cores re-extracted) and reports the null hit rate next to the real
one. If the null rate is not far below the real rate, the filter is noise and
the run says so.

Whether the human protein hit is an immunoglobulin V germline gene is recorded
separately - for an antibody-derived ligand that is the mechanistically
meaningful category.

Output: results/m4_core_tolerance.tsv, results/m4_filter_validation.json
"""
import csv
import json
import os
import random
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import load_config, read_fasta, data_path, results_path  # noqa: E402

NULL_SEED = 20260825


def iter_proteome(path):
    """Yield (accession, entry_name, description, sequence)."""
    hdr, buf = None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if hdr:
                    yield hdr + ("".join(buf),)
                h = line[1:].strip()
                parts = h.split("|")
                acc = parts[1] if len(parts) > 2 else h.split()[0]
                rest = parts[2] if len(parts) > 2 else h
                name = rest.split()[0]
                desc = rest[len(name):].strip()
                hdr, buf = (acc, name, desc), []
            else:
                buf.append(line.strip())
    if hdr:
        yield hdr + ("".join(buf),)


def build_index(cores):
    """{masked_9mer: {core}} for all 9 mask positions - enables <=1 mismatch lookup."""
    idx = defaultdict(set)
    for c in cores:
        for i in range(9):
            idx[c[:i] + "." + c[i + 1:]].add(c)
    return idx


def is_germline_v(name, desc):
    d = (name + " " + desc).lower()
    return ("immunoglobulin heavy variable" in d
            or "immunoglobulin kappa variable" in d
            or "immunoglobulin lambda variable" in d)


def main():
    cfg = load_config()
    proteome = os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), cfg["tolerance_filter"]["proteome"])
    seqs = read_fasta(data_path("sequences.fasta"))

    # ---- real cores -------------------------------------------------------
    cores = defaultdict(set)
    with open(results_path("m3_binding_long.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["call_el"] in ("SB", "WB"):
                cores[row["core"]].add(row["id"])
    cores = {c: v for c, v in cores.items() if len(c) == 9 and set(c) <= set("ACDEFGHIKLMNPQRSTVWY")}

    # ---- shuffled-sequence null ------------------------------------------
    rng = random.Random(NULL_SEED)
    null_cores = set()
    for sid, seq in seqs.items():
        chars = list(seq)
        rng.shuffle(chars)
        sh = "".join(chars)
        for i in range(0, len(sh) - 8, 3):     # sample every 3rd frame
            null_cores.add(sh[i:i + 9])
    null_cores -= set(cores)
    print(f"{len(cores)} predicted cores, {len(null_cores)} shuffled-null cores")

    all_cores = set(cores) | null_cores
    idx = build_index(all_cores)

    best = {}       # core -> (identity, acc, name, desc, human_9mer)
    n_prot = n_res = 0
    for acc, name, desc, seq in iter_proteome(proteome):
        n_prot += 1
        n_res += len(seq)
        for i in range(len(seq) - 8):
            k = seq[i:i + 9]
            hits = set()
            for j in range(9):
                hits |= idx.get(k[:j] + "." + k[j + 1:], set())
            for c in hits:
                ident = sum(a == b for a, b in zip(c, k))
                if ident > best.get(c, (0,))[0]:
                    best[c] = (ident, acc, name, desc, k)
    print(f"screened {n_prot} human proteins / {n_res:,} residues")

    def classify(c):
        ident = best.get(c, (0,))[0]
        if ident == 9:
            return "exact_human_9mer", 0.0
        if ident == 8:
            return "near_human_9mer", cfg["tolerance_filter"]["discount_tcrface"]
        return "foreign", 1.0

    # ---- validation: real vs shuffled null --------------------------------
    def rate(cs, min_ident):
        return sum(1 for c in cs if best.get(c, (0,))[0] >= min_ident) / max(len(cs), 1)

    validation = {
        "n_real_cores": len(cores),
        "n_null_cores": len(null_cores),
        "real_hit_rate_9of9": round(rate(cores, 9), 4),
        "null_hit_rate_9of9": round(rate(null_cores, 9), 4),
        "real_hit_rate_8of9": round(rate(cores, 8), 4),
        "null_hit_rate_8of9": round(rate(null_cores, 8), 4),
    }
    enr = validation["real_hit_rate_8of9"] / max(validation["null_hit_rate_8of9"], 1e-9)
    validation["enrichment_8of9_real_over_null"] = round(min(enr, 9999.0), 2)
    validation["filter_informative"] = (validation["null_hit_rate_8of9"] < 0.05
                                        and enr > 3.0)
    with open(results_path("m4_filter_validation.json"), "w") as f:
        json.dump(validation, f, indent=2)

    # ---- output -----------------------------------------------------------
    counts = defaultdict(int)
    with open(results_path("m4_core_tolerance.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["core", "from_sequences", "tolerance_class", "identity_to_human",
                    "human_hit_protein", "human_hit_entry", "human_hit_9mer",
                    "hit_is_germline_V", "weight"])
        for c in sorted(cores):
            cls, wt = classify(c)
            ident, acc, name, desc, k = best.get(c, (0, "-", "-", "", "-"))
            counts[cls] += 1
            w.writerow([c, ",".join(sorted(cores[c])), cls, ident, acc, name, k,
                        is_germline_v(name, desc) if ident >= 8 else False, wt])

    print("\ntolerance classification of predicted cores")
    for cls in ("exact_human_9mer", "near_human_9mer", "foreign"):
        n = counts[cls]
        print(f"  {cls:20s} {n:4d}  ({100*n/len(cores):5.1f}%)")
    print("\nfilter validation (real cores vs shuffled-sequence null)")
    for k in ("9of9", "8of9"):
        print(f"  >={k}: real {validation['real_hit_rate_'+k]*100:5.1f}%   "
              f"null {validation['null_hit_rate_'+k]*100:5.1f}%")
    print(f"  enrichment at 8/9: {validation['enrichment_8of9_real_over_null']}x")
    print(f"  filter informative: {validation['filter_informative']}")


if __name__ == "__main__":
    main()
