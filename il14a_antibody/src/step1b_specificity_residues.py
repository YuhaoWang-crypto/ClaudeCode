"""Step 1b - Specificity-determining residue (SDR) mapping.

Step 1 showed the core problem with this target: within the secreted species
(aa 176-546) the structurally confident surface IS the coiled-coil, and the
coiled-coil is what TXLNA shares with beta-/gamma-taxilin. Picking a window by
composite score alone either buys structural confidence at the cost of paralog
cross-reactivity, or buys specificity from a disordered tail.

The way out is not to pick a "more unique window" - it is to pick a designable
window that CONTAINS TXLNA-unique positions, and then require the paratope to
contact those positions. This script finds them, using real global alignments
of TXLNA against both paralogs (Bio.Align.PairwiseAligner, BLOSUM62) rather
than k-mer heuristics.

Output: per-residue table of {aligned TXLNB residue, aligned TXLNG residue,
unique?, rel_SASA, pLDDT} plus a re-ranking of candidate windows by
"exposed unique residue content".
"""
import json
import sys

from Bio import Align
from Bio.Align import substitution_matrices

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from common import DATA, PATENT_ZONE, RESULTS, SECRETED_END, SECRETED_START, load_taxilins
from step1_epitope_scan import per_residue_structure, liabilities, overlaps


def align_pair(a, b):
    """Global alignment; returns list mapping each index of `a` to a residue of
    `b` or None for a gap."""
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = "global"
    aln = aligner.align(a, b)[0]
    mapping = [None] * len(a)
    for (a_start, a_end), (b_start, b_end) in zip(*aln.aligned):
        for off in range(a_end - a_start):
            mapping[a_start + off] = b[b_start + off]
    return mapping, aln.score


def main():
    alpha, beta, gamma = load_taxilins()
    struct = per_residue_structure(DATA / "AF-P40222.pdb")
    map_b, score_b = align_pair(alpha, beta)
    map_g, score_g = align_pair(alpha, gamma)

    rows = []
    for i, aa in enumerate(alpha, start=1):
        if i not in struct:
            continue
        _, rel_sasa, plddt = struct[i]
        rb, rg = map_b[i - 1], map_g[i - 1]
        # "unique" = differs from the aligned residue in BOTH paralogs
        unique = (rb != aa) and (rg != aa)
        rows.append({
            "pos": i, "aa": aa, "txlnb": rb, "txlng": rg,
            "unique": bool(unique), "rel_sasa": round(rel_sasa, 3),
            "plddt": round(plddt, 1),
            # an SDR is only useful to a paratope if it is solvent-exposed
            "exposed_unique": bool(unique and rel_sasa >= 0.25),
        })

    # --- re-rank windows by exposed-unique content -----------------------
    by_pos = {r["pos"]: r for r in rows}
    windows = []
    for length in (14, 16, 18, 20):
        for start in range(SECRETED_START, SECRETED_END - length + 2):
            end = start + length - 1
            if overlaps(start, end, *PATENT_ZONE):
                continue
            w = [by_pos[p] for p in range(start, end + 1) if p in by_pos]
            if len(w) != length:
                continue
            n_exp_uniq = sum(1 for r in w if r["exposed_unique"])
            mean_plddt = sum(r["plddt"] for r in w) / length
            mean_sasa = sum(r["rel_sasa"] for r in w) / length
            frag = "".join(r["aa"] for r in w)
            lia = liabilities(frag)
            # designability gate: a structure-based paratope needs a confident
            # backbone. Below pLDDT 70 the window is modelled as disordered.
            designable = mean_plddt >= 70.0
            windows.append({
                "start": start, "end": end, "sequence": frag,
                "n_exposed_unique": int(n_exp_uniq),
                "exposed_unique_positions": [
                    f"{r['aa']}{r['pos']}" for r in w if r["exposed_unique"]],
                "mean_plddt": round(mean_plddt, 1),
                "mean_rel_sasa": round(mean_sasa, 3),
                "designable": bool(designable),
                "n_liabilities": len(lia),
                "liabilities": lia,
                # SDR density per residue, gated on designability
                "sdr_score": round(
                    (n_exp_uniq / length) * mean_sasa * (mean_plddt / 100.0), 4),
            })

    designable = [w for w in windows if w["designable"]]
    designable.sort(key=lambda w: -w["sdr_score"])
    picked, seen = [], []
    for w in designable:
        if all(not overlaps(w["start"], w["end"], p["start"], p["end"]) for p in seen):
            picked.append(w)
            seen.append(w)
        if len(picked) == 6:
            break

    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "step1b_sdr_map.json").write_text(json.dumps({
        "alignment_scores": {"TXLNA_vs_TXLNB": score_b, "TXLNA_vs_TXLNG": score_g},
        "n_exposed_unique_in_secreted": sum(
            1 for r in rows if r["exposed_unique"] and SECRETED_START <= r["pos"]),
        "per_residue": rows,
        "top_designable_sdr_windows": picked,
    }, indent=2))

    print("=" * 100)
    print("GLOBAL ALIGNMENT (BLOSUM62) of TXLNA against its paralogs")
    print("=" * 100)
    ident_b = sum(1 for i, a in enumerate(alpha) if map_b[i] == a) / len(alpha)
    ident_g = sum(1 for i, a in enumerate(alpha) if map_g[i] == a) / len(alpha)
    print(f"  TXLNA vs TXLNB : {ident_b*100:5.1f}% identity over full length")
    print(f"  TXLNA vs TXLNG : {ident_g*100:5.1f}% identity over full length")
    sec = [r for r in rows if r["pos"] >= SECRETED_START]
    print(f"  secreted species aa {SECRETED_START}-{SECRETED_END}: "
          f"{sum(1 for r in sec if r['unique'])} unique residues, "
          f"{sum(1 for r in sec if r['exposed_unique'])} of them solvent-exposed")
    print()
    print("=" * 100)
    print("TOP DESIGNABLE WINDOWS RANKED BY EXPOSED SPECIFICITY-DETERMINING RESIDUES")
    print("(designable = mean pLDDT >= 70, i.e. a paratope can be built against it)")
    print("=" * 100)
    print(f"{'rank':5} {'range':11} {'pLDDT':7} {'relSASA':8} {'nSDR':5} {'sdr_score':10} sequence")
    print("-" * 100)
    for i, w in enumerate(picked, 1):
        print(f"{i:<5} {str((w['start'], w['end'])):11} {w['mean_plddt']:<7.1f} "
              f"{w['mean_rel_sasa']:<8.3f} {w['n_exposed_unique']:<5} "
              f"{w['sdr_score']:<10.4f} {w['sequence']}")
        print(f"      SDRs: {', '.join(w['exposed_unique_positions'])}")
        if w["liabilities"]:
            print(f"      liabilities: {', '.join(w['liabilities'])}")
    print()
    print("wrote results/step1b_sdr_map.json")


if __name__ == "__main__":
    main()
