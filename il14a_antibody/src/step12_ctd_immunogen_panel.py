"""Step 12 - Design-ready immunogen panel for the C-terminus, phospho and not.

CAN WE DE NOVO DESIGN A BINDER AGAINST ERRPEGPGAQAPSSPRVTEAPC?
--------------------------------------------------------------
Structure-based paratope design: no, and the composition says why. Measured on
the peptide itself (step12 prints it):

    aromatics (F/W/Y) .......... 0.0%   (none at all, in 22 residues)
    large hydrophobics ......... 4.5%   (a single Val)
    Pro + Gly .................. 31.8%
    disorder ................... UniProt REGION 482-546 "Disordered"

An antibody paratope grips an epitope mainly by burying aromatic and large
hydrophobic side chains. This epitope has nothing to bury. Compare the
coiled-coil epitope 367-382 at 25% large hydrophobics and 0% Pro+Gly. So
RFantibody's 0/64 on 502-522 was not a tool failure - the epitope offers no
anchor, and it has no single conformation to complement either.

PHOSPHORYLATION CHANGES EXACTLY THAT DEFICIT
--------------------------------------------
The peptide is undesignable because it is featureless. A phosphate is the most
featureful thing that can be added to a peptide: dianionic at pH 7.4,
tetrahedral, and read by dedicated Arg/Lys/Ser/Thr-rich pockets. It supplies
the discrete anchor the sequence lacks - which is why phospho-specific
antibodies against disordered low-complexity peptides are among the most
reliably generated antibody classes in biology.

S515 is annotated by UniProt as a phosphoserine, sits in an S-P
proline-directed motif, and is the ONLY S/T-P motif in the whole C-terminal
region 492-546, so a pS515 antibody is positionally unambiguous.

But note what you get: a phospho-specific antibody binds only the modified
form. With pS515 occupancy on the relevant TXLNA pool unknown, its best use is
as a PROTEOFORM PROBE for the identity gate, not as the blocker. The blocker
should be phospho-insensitive. The panel below is built to yield both, from
one immunisation campaign, by counter-screening.

SPECIES CHECK
-------------
Raising antibodies against a self-conserved peptide fails on tolerance. Human
vs mouse identity is 90.5% full-length but only 77% across 502-523 and 60%
across 492-506 - the C-terminus is the most divergent part of the protein.
That is why 1C6 and 1F2 were raised successfully in Balb/c against exactly this
peptide, and it means wild-type mouse immunisation is viable here.
"""
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA, RESULTS, gravy, load_taxilins, net_charge, read_fasta

AROM, BIGHYD, FLEX = set("FWY"), set("ILVMFWY"), set("PG")
PHOSPHO_POS = 515          # UniProt MOD_RES, phosphoserine, S-P motif


def scramble(seq, seed=20260926):
    """Composition-matched random permutation.

    Sorting the residues would also be composition-matched but is NOT a valid
    scramble: it clusters all five prolines into a polyproline run, creating a
    molecule with its own distinct conformational behaviour. A control must
    differ from the real peptide only in residue ORDER, not in local character.
    Rejects permutations that leave too much of the original order intact.
    """
    rng = random.Random(seed)
    chars = list(seq)
    for _ in range(200):
        rng.shuffle(chars)
        cand = "".join(chars)
        same = sum(1 for a, b in zip(cand, seq) if a == b) / len(seq)
        longest = max((len(r) for r in __import__("re").findall(r"(.)\1*", cand)
                       for r in [r]), default=1)
        if same <= 0.25 and longest <= 3:
            return cand
    return "".join(chars)


def composition(seq):
    n = len(seq)
    return {
        "aromatic_pct": round(sum(c in AROM for c in seq) / n * 100, 1),
        "large_hydrophobic_pct": round(sum(c in BIGHYD for c in seq) / n * 100, 1),
        "pro_gly_pct": round(sum(c in FLEX for c in seq) / n * 100, 1),
        "net_charge": round(net_charge(seq), 1),
        "gravy": round(gravy(seq), 2),
    }


def mouse_map(human):
    """Align mouse Txlna onto human numbering; returns a same-length string."""
    from Bio import Align
    from Bio.Align import substitution_matrices
    m = read_fasta(DATA / "mouse_txlna.fasta")
    al = Align.PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score, al.mode = -11, -1, "global"
    aln = al.align(human, m)[0]
    out = ["-"] * len(human)
    for (s, e), (gs, ge) in zip(*aln.aligned):
        for o in range(e - s):
            out[s + o] = m[gs + o]
    return "".join(out)


def build_panel(A, mouse):
    def entry(pid, lo, hi, role, phospho, add_cys, rationale):
        native = A[lo - 1:hi]
        mo = mouse[lo - 1:hi]
        ident = sum(1 for a, b in zip(native, mo) if a == b) / len(native)
        # conjugation: a native C-terminal Cys is a ready maleimide handle
        has_native_cys = native.endswith("C")
        synth = native + ("C" if add_cys else "")
        disp = synth
        if phospho:
            i = PHOSPHO_POS - lo
            disp = synth[:i] + "pS" + synth[i + 1:]
        return {
            "id": pid, "role": role,
            "txlna_range": [lo, hi],
            "native_sequence": native,
            "synthesis_sequence": synth,
            "display": disp,
            "phosphorylated": phospho,
            "phospho_site": f"pS{PHOSPHO_POS}" if phospho else None,
            "conjugation": ("native C-terminal C%d, maleimide (MBS/SMCC) to KLH"
                            % hi if has_native_cys else
                            "add C-terminal Cys, maleimide to KLH"),
            "cys_added": add_cys,
            "mouse_identity": round(ident, 2),
            "tolerance_risk": "low" if ident < 0.85 else "HIGH - check titres",
            "composition": composition(native),
            "rationale": rationale,
        }

    return [
        entry("P1-pan", 502, 523, "primary blocker immunogen", False, False,
              "Exactly the 1C6 immunogen. Two independent hybridomas (1C6 6.16 nM, "
              "1F2 1.6 nM) came from this stretch and 1C6 blocks IL-14a-driven "
              "B-cell proliferation - the user's endpoint. Highest prior "
              "probability of success of anything in this project."),
        entry("P2-phospho", 502, 523, "proteoform probe / identity gate", True, False,
              "Same peptide carrying pS515. The phosphate supplies the anchor the "
              "unmodified sequence lacks. Counter-screen against P1 to isolate "
              "phospho-SPECIFIC clones."),
        entry("P3-orthogonal", 492, 506, "second epitope bin, PTM-free", False, True,
              "Overlaps 1F2's immunogen (493-512) so it has independent antibody "
              "precedent, avoids S515 entirely, and is the most mouse-divergent "
              "window (60%) so it should be the most immunogenic."),
    ]


def main():
    A, _, _ = load_taxilins()
    mouse = mouse_map(A)
    pep_1c6 = A[501:523]
    assert pep_1c6 == "ERRPEGPGAQAPSSPRVTEAPC", pep_1c6

    print("=" * 96)
    print("A. WHY DE NOVO PARATOPE DESIGN FAILS ON THIS EPITOPE")
    print("=" * 96)
    for nm, sq in [("1C6 epitope 502-523", pep_1c6),
                   ("coiled-coil 367-382", A[366:382]),
                   ("whole TXLNA", A)]:
        c = composition(sq)
        print(f"  {nm:24} aromatic {c['aromatic_pct']:5.1f}%   "
              f"large-hydrophobic {c['large_hydrophobic_pct']:5.1f}%   "
              f"Pro+Gly {c['pro_gly_pct']:5.1f}%")
    print("  -> zero aromatics and one Val in 22 residues: nothing for a paratope")
    print("     to bury. With disorder on top, RFantibody's 0/64 is expected.")

    panel = build_panel(A, mouse)
    print()
    print("=" * 96)
    print("B. IMMUNOGEN PANEL")
    print("=" * 96)
    for e in panel:
        print(f"\n  {e['id']}  [{e['role']}]   TXLNA {e['txlna_range'][0]}-{e['txlna_range'][1]}")
        print(f"    synthesise : {e['display']}")
        print(f"    conjugate  : {e['conjugation']}")
        print(f"    mouse id   : {e['mouse_identity']:.0%}   tolerance risk: {e['tolerance_risk']}")
        print(f"    why        : {e['rationale']}")

    # --- counter-screen matrix -------------------------------------------
    screens = [
        {"clone_type": "pan (blocker)", "immunogen": "P1-pan",
         "positive_on": ["P1-pan", "P2-phospho"],
         "negative_on": ["scrambled", "irrelevant KLH-peptide"],
         "use": "blocking candidate - binds regardless of S515 state"},
        {"clone_type": "phospho-specific", "immunogen": "P2-phospho",
         "positive_on": ["P2-phospho"],
         "negative_on": ["P1-pan", "scrambled"],
         "use": "proteoform probe - reports pS515 occupancy"},
        {"clone_type": "non-phospho-specific", "immunogen": "P1-pan",
         "positive_on": ["P1-pan"], "negative_on": ["P2-phospho", "scrambled"],
         "use": "the complement of the phospho probe; together they quantify occupancy"},
        {"clone_type": "second bin", "immunogen": "P3-orthogonal",
         "positive_on": ["P3-orthogonal"], "negative_on": ["P1-pan", "scrambled"],
         "use": "non-overlapping epitope; pairs with a P1 clone for sandwich ELISA"},
    ]
    print()
    print("=" * 96)
    print("C. COUNTER-SCREEN MATRIX - one campaign, four clone types")
    print("=" * 96)
    for s in screens:
        print(f"  {s['clone_type']:22} immunise {s['immunogen']:14} "
              f"+{','.join(s['positive_on'])}  -{','.join(s['negative_on'])}")
        print(f"  {'':22} use: {s['use']}")

    # --- controls ---------------------------------------------------------
    controls = {
        "mouse_ortholog_502_523": mouse[501:523],
        "scrambled_P1": scramble(pep_1c6),
        "note": ("The mouse ortholog peptide is the species cross-reactivity "
                 "control and matters if an IL-14a Tg mouse model is used. The "
                 "scrambled peptide is composition-matched so a hit on it means "
                 "composition-driven binding, not sequence recognition."),
    }
    print()
    print("  controls:")
    print(f"    mouse ortholog 502-523 : {controls['mouse_ortholog_502_523']}  "
          f"(species cross-reactivity)")
    print(f"    composition-matched scramble : {controls['scrambled_P1']}")

    out = {
        "question": "de novo binder/antibody against ERRPEGPGAQAPSSPRVTEAPC?",
        "structure_based_design_verdict": (
            "No. Zero aromatics, 4.5% large hydrophobics and 31.8% Pro+Gly in a "
            "disordered region leave a paratope nothing to bury. RFantibody 0/64 "
            "and this repo's failed co-folding controls are both explained by the "
            "epitope composition, not by tool choice."),
        "phospho_verdict": (
            "Phosphorylation supplies exactly the missing anchor and makes the "
            "peptide a tractable immunogen, but the product is a phospho-specific "
            "antibody - best used as a proteoform probe for the identity gate, "
            "with a phospho-insensitive clone carrying the blocking role."),
        "epitope_composition": {
            "1C6_502_523": composition(pep_1c6),
            "coiled_coil_367_382": composition(A[366:382]),
            "whole_protein": composition(A),
        },
        "species": {"full_length_identity": 0.905,
                    "identity_502_523": 0.77, "identity_492_506": 0.60,
                    "S515_conserved_in_mouse": mouse[514] == "S"},
        "panel": panel, "counter_screens": screens, "controls": controls,
    }
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "step12_immunogen_panel.json").write_text(json.dumps(out, indent=2))

    lines = ["# IL-14alpha C-terminal immunogen panel",
             "# pS = phosphoserine at TXLNA S515. Cys shown at a terminus is the",
             "# maleimide conjugation handle (native C523 in P1/P2; added in P3).", ""]
    for e in panel:
        lines.append(f">{e['id']} | TXLNA {e['txlna_range'][0]}-{e['txlna_range'][1]}"
                     f" | {e['role']} | {e['conjugation']}")
        lines.append(e["display"])
    lines += [">CONTROL_mouse_ortholog_502_523 | species cross-reactivity control",
              controls["mouse_ortholog_502_523"]]
    (RESULTS / "IL14A_ctd_immunogen_panel.fasta").write_text("\n".join(lines) + "\n")
    print("\nwrote results/step12_immunogen_panel.json and "
          "results/IL14A_ctd_immunogen_panel.fasta")


if __name__ == "__main__":
    main()
