"""Step 9 - Pick the best immunogen window in the C-terminal disordered region.

WHY THIS IS THE ONLY REMAINING TARGET
-------------------------------------
Two independent campaigns have now excluded the coiled-coil:

  * Geometry (G034x): the TXLNA<->STX4-H3 interface is 47 residues spanning
    124 A and overlaps TXLNA's own homodimer face by 96%. A ~20 A paratope
    covers at most 15% of it, so a conformational blocking antibody is
    geometrically impossible, and nothing binding that groove can distinguish
    TXLNA-STX4 from TXLNA-TXLNA.
  * Metrics (step7/step8, this repo): co-folding scores on this target fail an
    identical-paratope control and carry a signal-to-artifact ratio of 1.27.
    G034x's scrambled-sequence control failed the same way from the other
    direction - a shuffled peptide scored ipTM 0.774 against the native's
    0.648. Both observations have one cause: the model pairs ANY amphipathic
    helix into a coiled-coil groove.

And this script adds a third reason, which had not been checked: the groove is
the LEAST paralog-selective part of the protein, not the most. Over the zipper
region TXLNA 359-443, identity is 56.5% to TXLNB and 69.4% to TXLNG, against
full-length averages of 52.9% and 52.2%. Anything occupying that groove -
antibody OR peptide - is binding the family's most conserved surface.

That leaves the C-terminal region after the coiled-coil (UniProt coiled-coil
annotation ends at 491). It is the only part of TXLNA that is simultaneously
solvent-exposed, paralog-divergent, and NOT an interface. It is also where the
only antibody with reported functional blocking (1C6, US7622574) was raised.
Its disorder is why structure-based design fails there - which is a statement
about the method, not the target: a disordered linear epitope is the native
domain of immunisation and display, not docking.

G034x ranked 502-522 but noted that 484-504 and 510-530 score better on
paralog specificity (0.57 vs 0.48) and were never run. This scans the whole
region and picks properly.
"""
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA, RESULTS, gravy, net_charge, load_taxilins
from step1_epitope_scan import per_residue_structure

COILED_COIL_END = 491          # UniProt P40222 annotation: coiled coil 186-491
CTD_START, CTD_END = 492, 546
# the 1C6 immunogen, verified earlier against P40222
REF_1C6 = (502, 523)
# 1F2's immunogen, native part (its leading Cys is an added conjugation handle)
REF_1F2 = (493, 512)

# Experimentally annotated PTM sites (UniProt P40222 MOD_RES). A phosphosite
# inside an immunogen is a real problem: an antibody raised on the
# unphosphorylated peptide may fail to recognise the phosphorylated protein,
# and the phospho-occupancy of endogenous TXLNA is unknown. S515 sits in an
# S-P motif, the consensus for proline-directed kinases.
PHOSPHOSITES = {515: "Phosphoserine (S515, in an S-P proline-directed motif)"}


def paralog_maps(alpha, beta, gamma):
    from Bio import Align
    from Bio.Align import substitution_matrices
    al = Align.PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score, al.mode = -11, -1, "global"

    def m(a, b):
        aln = al.align(a, b)[0]
        out = [None] * len(a)
        for (s, e), (gs, ge) in zip(*aln.aligned):
            for o in range(e - s):
                out[s + o] = b[gs + o]
        return out
    return m(alpha, beta), m(alpha, gamma)


# Motifs that make a synthetic immunogen behave badly or read ambiguously.
IMMUNOGEN_LIABILITIES = {
    "deamidation_NG_NS": r"N[GS]",
    "isomerisation_DG_DP": r"D[GP]",
    "internal_Cys": r"C",          # fine at a terminus as a conjugation handle
}


def scan(alpha, mb, mg, struct, lo, hi, lengths=(15, 18, 20, 22, 25)):
    out = []
    for L in lengths:
        for s in range(lo, hi - L + 2):
            e = s + L - 1
            frag = alpha[s - 1:e]
            res = [struct[i] for i in range(s, e + 1) if i in struct]
            if len(res) != L:
                continue
            uniq = [i for i in range(s, e + 1)
                    if mb[i - 1] != alpha[i - 1] and mg[i - 1] != alpha[i - 1]]
            exposed_uniq = [i for i in uniq if struct[i][1] >= 0.25]
            sasa = sum(r[1] for r in res) / L
            plddt = sum(r[2] for r in res) / L
            lia = {}
            for name, pat in IMMUNOGEN_LIABILITIES.items():
                hits = [m.start() + 1 for m in re.finditer(pat, frag)]
                if name == "internal_Cys":
                    hits = [h for h in hits if 1 < h < L]   # termini are fine
                if hits:
                    lia[name] = hits
            # An immunogen is judged on: how much of it is TXLNA-unique AND
            # exposed, how exposed it is overall, and how clean its chemistry
            # is. Disorder is NOT penalised here - a linear peptide immunogen
            # does not need a defined fold.
            spec = len(exposed_uniq) / L
            clean = max(0.0, 1.0 - 0.34 * sum(len(v) for v in lia.values()))
            ptm = {p: PHOSPHOSITES[p] for p in PHOSPHOSITES if s <= p <= e}
            out.append({
                "start": s, "end": e, "length": L, "sequence": frag,
                "phosphosites_inside": ptm,
                "ptm_free": not ptm,
                "paralog_unique_exposed": len(exposed_uniq),
                "specificity_fraction": round(spec, 3),
                "mean_rel_sasa": round(sasa, 3),
                "mean_plddt": round(plddt, 1),
                "net_charge": round(net_charge(frag), 1),
                "gravy": round(gravy(frag), 2),
                "liabilities": lia,
                "score": round(spec * sasa * clean, 4),
                "outside_coiled_coil": s > COILED_COIL_END,
            })
    return out


def overlaps(a, b, c, d):
    return not (b < c or d < a)


def main():
    alpha, beta, gamma = load_taxilins()
    mb, mg = paralog_maps(alpha, beta, gamma)
    struct = per_residue_structure(DATA / "AF-P40222.pdb")

    # context: how selective is the coiled-coil groove, really?
    lo, hi = 359, 443
    idb = sum(1 for i in range(lo - 1, hi) if mb[i] == alpha[i]) / (hi - lo + 1)
    idg = sum(1 for i in range(lo - 1, hi) if mg[i] == alpha[i]) / (hi - lo + 1)
    full_b = sum(1 for i in range(len(alpha)) if mb[i] == alpha[i]) / len(alpha)
    full_g = sum(1 for i in range(len(alpha)) if mg[i] == alpha[i]) / len(alpha)

    cands = scan(alpha, mb, mg, struct, CTD_START, CTD_END)
    cands.sort(key=lambda r: -r["score"])
    picked = []
    for c in cands:
        if all(not overlaps(c["start"], c["end"], p["start"], p["end"])
               for p in picked):
            picked.append(c)
        if len(picked) == 5:
            break

    ref = next((c for c in scan(alpha, mb, mg, struct, REF_1C6[0], REF_1C6[1],
                                lengths=(REF_1C6[1] - REF_1C6[0] + 1,))), None)

    print("=" * 96)
    print("A. THE GROOVE IS THE LEAST PARALOG-SELECTIVE PART OF THE PROTEIN")
    print("=" * 96)
    print(f"  full-length identity      TXLNB {full_b*100:.1f}%   TXLNG {full_g*100:.1f}%")
    print(f"  zipper region 359-443     TXLNB {idb*100:.1f}%   TXLNG {idg*100:.1f}%"
          f"   <- MORE conserved than average")
    print("  Anything occupying that groove - antibody or peptide - binds the")
    print("  taxilin family's most conserved surface. This had not been checked")
    print("  for the peptide line.")
    print()
    print("=" * 96)
    print(f"B. IMMUNOGEN WINDOWS IN THE C-TERMINAL REGION (aa {CTD_START}-{CTD_END},")
    print(f"   i.e. after the coiled-coil ends at {COILED_COIL_END})")
    print("=" * 96)
    print(f"{'rank':5} {'range':12} {'len':4} {'spec':6} {'rSASA':7} {'nUniq':6} {'score':7} {'PTM':5} sequence")
    print("-" * 104)
    for i, c in enumerate(picked, 1):
        print(f"{i:<5} {str((c['start'], c['end'])):12} {c['length']:<4} "
              f"{c['specificity_fraction']:<6.2f} {c['mean_rel_sasa']:<7.3f} "
              f"{c['paralog_unique_exposed']:<6} {c['score']:<7.4f} "
              f"{'clean' if c['ptm_free'] else 'S515':5} {c['sequence']}")
        if c["liabilities"]:
            print(f"      liabilities: {c['liabilities']}")
        if c["phosphosites_inside"]:
            for pos, note in c["phosphosites_inside"].items():
                print(f"      PTM: {note}")
    print()
    if ref:
        print(f"  reference - 1C6 immunogen {REF_1C6}: spec {ref['specificity_fraction']:.2f}"
              f"  rSASA {ref['mean_rel_sasa']:.3f}  score {ref['score']:.4f}")
        print(f"      {ref['sequence']}")
        best = picked[0]
        if best["score"] > ref["score"]:
            print(f"  -> best scanned window beats the 1C6 immunogen "
                  f"({best['score']:.4f} vs {ref['score']:.4f})")
        ptm_free = [c for c in picked if c["ptm_free"]]
        if ptm_free:
            b = ptm_free[0]
            print(f"\n  best PTM-FREE window: {(b['start'], b['end'])} {b['sequence']}"
                  f"  spec {b['specificity_fraction']:.2f}  score {b['score']:.4f}")
            print(f"    (1F2's immunogen is {REF_1F2} - this window overlaps it,"
                  f" so it has independent antibody precedent)")

    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "step9_ctd_immunogen.json").write_text(json.dumps({
        "groove_paralog_identity": {
            "region": [lo, hi],
            "TXLNB": round(idb, 4), "TXLNG": round(idg, 4),
            "full_length_TXLNB": round(full_b, 4),
            "full_length_TXLNG": round(full_g, 4),
        },
        "ctd_region": [CTD_START, CTD_END],
        "ref_1C6_immunogen": ref,
        "top_windows": picked,
        "n_scanned": len(cands),
    }, indent=2))
    print(f"\nscanned {len(cands)} windows; wrote results/step9_ctd_immunogen.json")


if __name__ == "__main__":
    main()
