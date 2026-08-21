"""Step 7: selected neoantigens -> one concatemeric, synthesis-ready mRNA CDS.

This is the part that turns a ranked list into an actual product, and it is
where a naive concatenation quietly goes wrong: fusing 34 minigenes head-to-tail
creates 33 brand-new junction sequences that were never in the patient's tumor.
Any strong binder created there is a decoy epitope -- it competes for the same
HLA molecules and for the same T-cell response, and it is not tumor-specific.

So the construct step does four things:
  1. cut a fixed-length minigene around each mutation (mutation centered)
  2. order the minigenes to *minimize* predicted junction binders (greedy + 2-opt
     over a junction-cost matrix scored with the patient's own HLA alleles)
  3. re-scan the final junctions and report anything that still binds
  4. reverse-translate to a human-codon-optimized CDS and run mRNA-level QC
"""

from __future__ import annotations

import itertools
import random
import re
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .presentation import predict_iedb
from .features import f_presentation

# --------------------------------------------------------------------------
# 1. minigenes
# --------------------------------------------------------------------------


def minigene(mut_protein: str, mut_pos0: int, length: int = 25) -> Tuple[str, int]:
    """Cut a `length`-mer centered on the mutated residue.

    Near a protein terminus the window slides instead of padding, so the
    minigene stays a real subsequence of the mutant protein. Returns
    (peptide, offset_of_mutation_within_peptide).
    """
    half = length // 2
    start = max(0, min(mut_pos0 - half, len(mut_protein) - length))
    start = max(0, start)
    seq = mut_protein[start:start + length]
    return seq, mut_pos0 - start


def build_minigenes(selected: pd.DataFrame, proteome: Dict[str, str],
                    length: int = 25) -> pd.DataFrame:
    """Uses gene + protein change to rebuild the mutant protein, then cuts."""
    from .peptides import mutant_protein

    rows = []
    for _, r in selected.iterrows():
        gene, pc = r.get("gene"), str(r.get("protein_change", ""))
        wt = proteome.get(gene)
        ref_aa, pos1, alt_aa = None, None, None
        m = re.match(r"^([A-Z])(\d+)([A-Z])$", pc.lstrip("p."))
        if m:
            ref_aa, pos1, alt_aa = m.group(1), int(m.group(2)), m.group(3)
        if wt and ref_aa:
            mut, status = mutant_protein(wt, ref_aa, pos1, alt_aa)
        else:
            mut, status = None, "cannot rebuild mutant protein"
        if mut is None:
            # fall back to the epitope itself so the slot is never silently lost
            rows.append(dict(slot=r.get("slot"), var_id=r.get("var_id"), gene=gene,
                             protein_change=pc, minigene=r.get("mut_peptide"),
                             mut_offset_in_minigene=r.get("mut_offset"),
                             note=f"short minigene: {status}"))
            continue
        seq, off = minigene(mut, pos1 - 1, length)
        rows.append(dict(slot=r.get("slot"), var_id=r.get("var_id"), gene=gene,
                         protein_change=pc, minigene=seq,
                         mut_offset_in_minigene=off, note=""))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# 2-3. junction handling
# --------------------------------------------------------------------------

def junction_peptides(a: str, b: str, linker: str = "", k: int = 9) -> List[str]:
    """All k-mers spanning the a|linker|b boundary (present in neither minigene)."""
    joined = a + linker + b
    lo = max(0, len(a) - k + 1)
    hi = min(len(joined) - k, len(a) + len(linker) - 1)
    out = []
    for s in range(lo, hi + 1):
        pep = joined[s:s + k]
        if len(pep) == k and pep not in a and pep not in b:
            out.append(pep)
    return out


def junction_cost_matrix(minigenes: Sequence[str], alleles: Sequence[str],
                         linker: str = "", k: int = 9,
                         backend: str = "iedb", batch_size: int = 500,
                         verbose: bool = True) -> Tuple[np.ndarray, pd.DataFrame]:
    """cost[i][j] = worst (strongest) predicted presentation of any k-mer created
    by placing minigene j directly after minigene i."""
    n = len(minigenes)
    pairs = [(i, j) for i in range(n) for j in range(n) if i != j]
    pep_of_pair: Dict[Tuple[int, int], List[str]] = {}
    all_peps: List[str] = []
    for i, j in pairs:
        ps = junction_peptides(minigenes[i], minigenes[j], linker, k)
        pep_of_pair[(i, j)] = ps
        all_peps.extend(ps)
    all_peps = list(dict.fromkeys(all_peps))
    if verbose:
        print(f"  junction scan: {n} minigenes -> {len(pairs)} orderings, "
              f"{len(all_peps)} unique {k}-mers x {len(alleles)} alleles")
    preds = predict_iedb(all_peps, alleles, k, "I", batch_size=batch_size)
    best = (preds.groupby("peptide")["percentile_rank"].min()
            if not preds.empty else pd.Series(dtype=float))
    cost = np.zeros((n, n))
    for (i, j), ps in pep_of_pair.items():
        if not ps:
            continue
        vals = [f_presentation(best.get(p, 100.0)) for p in ps]
        cost[i, j] = max(vals) if vals else 0.0
    return cost, preds


def order_minimizing_junctions(cost: np.ndarray, seed: int = 0,
                               restarts: int = 20) -> List[int]:
    """Open-path TSP over the junction cost matrix: greedy nearest-neighbour from
    several starts, then 2-opt until no improvement."""
    n = cost.shape[0]
    if n <= 2:
        return list(range(n))
    rng = random.Random(seed)

    def path_cost(order):
        return sum(cost[order[i], order[i + 1]] for i in range(len(order) - 1))

    best_order, best_cost = None, float("inf")
    starts = list(range(n)) if n <= restarts else rng.sample(range(n), restarts)
    for s in starts:
        order, unused = [s], set(range(n)) - {s}
        while unused:
            last = order[-1]
            nxt = min(unused, key=lambda j: cost[last, j])
            order.append(nxt)
            unused.discard(nxt)
        improved = True
        while improved:
            improved = False
            for i in range(1, n - 1):
                for j in range(i + 1, n):
                    cand = order[:i] + order[i:j + 1][::-1] + order[j + 1:]
                    if path_cost(cand) < path_cost(order) - 1e-12:
                        order, improved = cand, True
        c = path_cost(order)
        if c < best_cost:
            best_order, best_cost = order, c
    return best_order


def scan_final_junctions(ordered_minigenes: Sequence[str], alleles: Sequence[str],
                         linker: str = "", ks=(8, 9, 10, 11),
                         flag_rank: float = 0.5) -> pd.DataFrame:
    """Predict every peptide created at the 33 junctions of the final order."""
    rows = []
    for k in ks:
        peps, meta = [], []
        for i in range(len(ordered_minigenes) - 1):
            for p in junction_peptides(ordered_minigenes[i], ordered_minigenes[i + 1], linker, k):
                peps.append(p)
                meta.append((i, p))
        if not peps:
            continue
        pr = predict_iedb(peps, alleles, k, "I")
        if pr.empty:
            continue
        for i, p in meta:
            sub = pr[pr["peptide"] == p]
            if sub.empty:
                continue
            best = sub.loc[sub["percentile_rank"].idxmin()]
            rows.append({"junction": i + 1, "length": k, "peptide": p,
                         "allele": best["allele"],
                         "percentile_rank": float(best["percentile_rank"]),
                         "flagged": bool(best["percentile_rank"] <= flag_rank)})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# 4. reverse translation + mRNA QC
# --------------------------------------------------------------------------

HUMAN_CODON_USAGE = {
    "A": [("GCC", .40), ("GCT", .26), ("GCA", .23), ("GCG", .11)],
    "R": [("AGA", .20), ("AGG", .20), ("CGG", .21), ("CGC", .19), ("CGA", .11), ("CGT", .08)],
    "N": [("AAC", .54), ("AAT", .46)],
    "D": [("GAC", .54), ("GAT", .46)],
    "C": [("TGC", .55), ("TGT", .45)],
    "Q": [("CAG", .75), ("CAA", .25)],
    "E": [("GAG", .58), ("GAA", .42)],
    "G": [("GGC", .34), ("GGA", .25), ("GGG", .25), ("GGT", .16)],
    "H": [("CAC", .59), ("CAT", .41)],
    "I": [("ATC", .48), ("ATT", .36), ("ATA", .16)],
    "L": [("CTG", .41), ("CTC", .20), ("TTG", .13), ("CTT", .13), ("CTA", .07), ("TTA", .07)],
    "K": [("AAG", .58), ("AAA", .42)],
    "M": [("ATG", 1.0)],
    "F": [("TTC", .55), ("TTT", .45)],
    "P": [("CCC", .33), ("CCT", .28), ("CCA", .27), ("CCG", .11)],
    "S": [("AGC", .24), ("TCC", .22), ("TCT", .18), ("AGT", .15), ("TCA", .15), ("TCG", .06)],
    "T": [("ACC", .36), ("ACA", .28), ("ACT", .24), ("ACG", .12)],
    "W": [("TGG", 1.0)],
    "Y": [("TAC", .57), ("TAT", .43)],
    "V": [("GTG", .47), ("GTC", .24), ("GTT", .18), ("GTA", .11)],
    "*": [("TGA", .47), ("TAA", .30), ("TAG", .23)],
}

RESTRICTION = {"EcoRI": "GAATTC", "BamHI": "GGATCC", "NotI": "GCGGCCGC",
               "BsaI": "GGTCTC", "BsmBI": "CGTCTC", "XbaI": "TCTAGA",
               "HindIII": "AAGCTT", "SapI": "GCTCTTC", "NheI": "GCTAGC"}

# mRNA-specific liabilities
CRYPTIC_POLYA = "AATAAA"
SPLICE_DONOR = re.compile(r"GGTAAG|GGTGAG")
SPLICE_ACCEPTOR = re.compile(r"[CT]{10}[ACGT]{1,3}AG")


def codon_optimize(protein: str, seed: int = 0, avoid: Sequence[str] = (),
                   max_homopolymer: int = 6, gc_window: int = 50,
                   gc_target: Tuple[float, float] = (40.0, 70.0),
                   rounds: int = 400) -> Dict[str, object]:
    """Frequency-weighted reverse translation, then local repair.

    Repair loop: while a liability exists (restriction site, homopolymer run,
    out-of-range GC window, cryptic polyA signal), re-draw the codons that
    overlap it. Synonymous only -- the protein is guaranteed unchanged.
    """
    rng = random.Random(seed)
    avoid_seqs = [RESTRICTION[a] for a in avoid if a in RESTRICTION] + \
                 [a for a in avoid if a not in RESTRICTION and set(a) <= set("ACGT")]

    def draw(aa):
        opts = HUMAN_CODON_USAGE.get(aa)
        if not opts:
            return "NNN"
        r, acc = rng.random(), 0.0
        for c, f in opts:
            acc += f
            if r <= acc:
                return c
        return opts[0][0]

    codons = [draw(aa) for aa in protein]

    def dna():
        return "".join(codons)

    def liabilities(s):
        out = []
        for site in avoid_seqs:
            for m in re.finditer(f"(?={site})", s):
                out.append((m.start(), m.start() + len(site), f"site:{site}"))
        for m in re.finditer(r"(A{%d,}|C{%d,}|G{%d,}|T{%d,})" % ((max_homopolymer + 1,) * 4), s):
            out.append((m.start(), m.end(), f"homopolymer:{m.group(0)[0]}x{len(m.group(0))}"))
        for m in re.finditer(f"(?={CRYPTIC_POLYA})", s):
            out.append((m.start(), m.start() + 6, "cryptic_polyA"))
        for m in SPLICE_DONOR.finditer(s):
            out.append((m.start(), m.end(), "splice_donor_motif"))
        for i in range(0, max(1, len(s) - gc_window + 1), 10):
            win = s[i:i + gc_window]
            if len(win) < gc_window:
                break
            gc = 100.0 * (win.count("G") + win.count("C")) / len(win)
            if not (gc_target[0] <= gc <= gc_target[1]):
                out.append((i, i + gc_window, f"gc_window:{gc:.0f}%"))
        return out

    repairs = []
    for _ in range(rounds):
        s = dna()
        libs = liabilities(s)
        if not libs:
            break
        start, end, what = libs[0]
        ci, cj = start // 3, min(len(codons) - 1, end // 3)
        for c in range(ci, cj + 1):
            codons[c] = draw(protein[c])
        repairs.append(what)

    s = dna()
    return {"dna": s, "repairs": repairs, "remaining": liabilities(s)}


def translate(dna: str) -> str:
    table = {}
    for aa, opts in HUMAN_CODON_USAGE.items():
        for c, _ in opts:
            table[c] = aa
    # complete the standard code for codons not in the usage table
    bases = "TCAG"
    std = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
    for i, b1 in enumerate(bases):
        for j, b2 in enumerate(bases):
            for k, b3 in enumerate(bases):
                table.setdefault(b1 + b2 + b3, std[i * 16 + j * 4 + k])
    return "".join(table.get(dna[i:i + 3], "X") for i in range(0, len(dna) - 2, 3))


def mrna_qc(dna: str, protein: str, rules) -> Dict[str, object]:
    gc = 100.0 * (dna.count("G") + dna.count("C")) / max(1, len(dna))
    u = 100.0 * dna.count("T") / max(1, len(dna))     # T in cDNA == U in mRNA
    flags = []
    if len(dna) > rules.max_cds_nt:
        flags.append(f"CDS {len(dna)} nt exceeds payload budget {rules.max_cds_nt} nt")
    if translate(dna).rstrip("*") != protein.rstrip("*"):
        flags.append("translation does not match the designed protein")
    for name, site in RESTRICTION.items():
        if site in dna and name in rules.avoid_sites:
            flags.append(f"residual restriction site {name}")
    if CRYPTIC_POLYA in dna:
        flags.append("cryptic polyadenylation signal AATAAA present")
    for m in re.finditer(r"(A{7,}|C{7,}|G{7,}|T{7,})", dna):
        flags.append(f"homopolymer run {m.group(0)[0]}x{len(m.group(0))} at {m.start()}")
    return {"length_nt": len(dna), "length_aa": len(protein), "gc_percent": round(gc, 1),
            "uridine_percent": round(u, 1), "flags": flags, "pass": not flags}


def assemble(selected: pd.DataFrame, proteome: Dict[str, str], alleles: Sequence[str],
             rules, optimize_order: bool = True, seed: int = 0,
             verbose: bool = True) -> Dict[str, object]:
    """Full construct build. Returns minigenes, order, junction report, CDS, QC."""
    mg = build_minigenes(selected, proteome, rules.epitope_length)
    seqs = list(mg["minigene"])
    order = list(range(len(seqs)))
    cost, junction_preds = None, None
    if optimize_order and len(seqs) > 2:
        cost, junction_preds = junction_cost_matrix(seqs, alleles, rules.linker,
                                                    verbose=verbose)
        naive_cost = sum(cost[i, i + 1] for i in range(len(seqs) - 1))
        order = order_minimizing_junctions(cost, seed=seed)
        opt_cost = sum(cost[order[i], order[i + 1]] for i in range(len(order) - 1))
    else:
        naive_cost = opt_cost = float("nan")

    mg_ordered = mg.iloc[order].reset_index(drop=True)
    mg_ordered.insert(0, "position", range(1, len(mg_ordered) + 1))
    ordered_seqs = list(mg_ordered["minigene"])

    protein = rules.signal_peptide + rules.linker.join(ordered_seqs)
    opt = codon_optimize(protein, seed=seed, avoid=rules.avoid_sites,
                         max_homopolymer=rules.max_homopolymer,
                         gc_target=rules.gc_target)
    qc = mrna_qc(opt["dna"], protein, rules)

    junctions = scan_final_junctions(ordered_seqs, alleles, rules.linker,
                                     flag_rank=rules.junction_scan_rank) \
        if len(ordered_seqs) > 1 else pd.DataFrame()

    return {
        "minigenes": mg_ordered,
        "order": order,
        "junction_cost_naive": naive_cost,
        "junction_cost_optimized": opt_cost,
        "junction_scan": junctions,
        "protein": protein,
        "cds": opt["dna"],
        "codon_repairs": opt["repairs"],
        "qc": qc,
    }
