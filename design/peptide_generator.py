#!/usr/bin/env python3
"""
De novo / target-informed peptide binder generator for Tirzepatide.

Tirzepatide (39 aa, GIP/GLP-1 dual agonist) forms a central amphipathic
alpha-helix (~res 6-28). Its hydrophobic/aromatic stripe presents:
  F6, Y10, I12, L14, I17, F22, V23, W25, L26, I27
and polar/charged residues E3, D9, D15, K16, Q19, K20, Q24 face outward.

Strategy families (structure-INDEPENDENT; scored later by Boltz-2 co-fold):
  A. Complementary amphipathic alpha-helices (heptad abcdefg design) that
     pack a hydrophobic face against the target helix and salt-bridge K16/K20/D9/D15.
  C. Receptor-mimetic peptides: idealized segments inspired by the incretin
     receptor ECD groove that natively grips this peptide class.
  D. Disulfide-stapled / macrocyclic hairpins presenting a hydrophobic loop
     (protease-resistant, easy SPPS).

Binder-derived families (B, E) are generated in a second module once the
top existing binders have Boltz structures and parsed interfaces.

Deterministic (no RNG) so the library is reproducible.
"""
import json, itertools

TZP = "Y-Aib-EGTFTSDYSI-Aib-LDKIAQ-K*-AFVQWLIAGGPSSGAPPPS-NH2"

cands = []  # list of {id, seq, family, note}

def add(cid, seq, family, note):
    seq = seq.upper().replace(" ", "")
    assert set(seq) <= set("ACDEFGHIKLMNPQRSTVWY"), f"{cid}: bad residue in {seq}"
    cands.append({"id": cid, "seq": seq, "family": family, "note": note})

# ---------------------------------------------------------------------------
# Family A: complementary amphipathic alpha-helices
# heptad positions a,d (+e sometimes) = hydrophobic knobs into target groove;
# b,c,f,g = solubilizing polar/charged, with e/g carrying charges for salt bridges.
# We vary (i) the hydrophobic knob identity, (ii) charge pattern, (iii) length.
# ---------------------------------------------------------------------------
# knob palettes complementary to target aromatics/aliphatics
knob_sets = {
    "aromWL": ["W", "L", "Y", "L", "F", "I", "W"],   # aromatic-rich to stack F6/Y10/W25/F22
    "aliphLIV": ["L", "I", "V", "L", "I", "V", "L"],  # tight aliphatic packing
    "mixFM":  ["F", "L", "M", "I", "W", "L", "V"],
}
# outward polar/charge templates for b,c,e,f,g (index within heptad)
# heptad: a b c d e f g   -> a,d hydrophobic; e,g charged; b,c,f polar
polar_face = {
    "acidic": {"b": "S", "c": "E", "e": "E", "f": "T", "g": "K"},  # net acidic-ish, salt bridge to K16/K20
    "basic":  {"b": "T", "c": "K", "e": "K", "f": "S", "g": "E"},  # to D9/D15/E3
    "neutral":{"b": "S", "c": "Q", "e": "N", "f": "T", "g": "S"},
}
hept_order = ["a", "b", "c", "d", "e", "f", "g"]
lengths = [28, 35]  # ~4 and 5 heptads
idx = 0
for kname, knobs in knob_sets.items():
    for pname, pf in polar_face.items():
        for L in lengths:
            seq = []
            hi = 0  # heptad counter for knob rotation
            for i in range(L):
                pos = hept_order[i % 7]
                if pos in ("a", "d"):
                    seq.append(knobs[hi % len(knobs)])
                    if pos == "d":
                        hi += 1
                else:
                    seq.append(pf[pos])
            # cap with helix N-cap (S/T) and C-cap, add a Trp anchor near middle
            s = "".join(seq)
            # place a single W anchor at ~ position 1/3 if none nearby (aromatic hotspot)
            idx += 1
            add(f"deA_{kname}_{pname}_{L}", "S" + s + "G", "A_helix",
                f"amphipathic helix; knobs={kname}; face={pname}; len~{L}")

# a few explicit hand-tuned amphipathic helices (strong hydrophobic moment)
add("deA_hm1", "SEELWKELAKLAEEWLKKLAELGKKLWELLKQLA", "A_helix", "high hydrophobic-moment aromatic helix")
add("deA_hm2", "TKEIWDEIKKLWEEIKKLWDELKKLGEIWEKLKK", "A_helix", "charged-complementary aromatic helix")
add("deA_hm3", "SQELLKKFAEELLKKFAEWLLKKFAEELGKKFA",  "A_helix", "Phe-ladder amphipathic helix")
add("deA_hm4", "GEEWLARLKELAEKWLARLKELAEKWLARLKKLG", "A_helix", "heptad LxxLxxx aromatic anchor helix")

# ---------------------------------------------------------------------------
# Family C: incretin-receptor-ECD-mimetic peptides
# The GIPR/GLP-1R N-terminal ECD grips the incretin C-terminal helix via a
# hydrophobic groove (aromatic + Trp cage). We build idealized groove peptides.
# ---------------------------------------------------------------------------
recep_mimetics = [
    ("deC_gipr_ecd1", "CQEGWLYWEATNKTGTWQELEDCGRIWLPWSK", "GIPR-ECD-like Trp-cage groove w/ disulfide"),
    ("deC_gipr_ecd2", "SWGEFTLCDWLHSQNMTREYCELLWNILKAG",   "GIPR-ECD groove variant"),
    ("deC_glp1r_ecd1","CLDEYACWHYAESGKSWSREYCDLLWRLIEG",   "GLP-1R-ECD-like aromatic groove"),
    ("deC_glp1r_ecd2","TWSEYLDWLLQKNGSWQRECEKLWNTLREAG",   "GLP-1R-ECD groove variant"),
    ("deC_groove1",   "WEYLDKAWQELGKSWNEIKKLWDELGKAWEK",   "idealized aromatic groove helix-loop"),
    ("deC_groove2",   "FEYWDKLAQEWGKSFNEIKKFWDELGKAFEW",   "Phe/Trp groove variant"),
]
for cid, seq, note in recep_mimetics:
    add(cid, seq, "C_receptor_mimetic", note)

# ---------------------------------------------------------------------------
# Family D: disulfide-stapled / macrocyclic hairpins
# Cys...Cys clamp a beta-hairpin presenting a hydrophobic binding loop with a
# Trp/Phe hotspot; type I' turn (NG / pG) between strands. Protease resistant.
# ---------------------------------------------------------------------------
turns = ["NG", "PG", "DG"]
loops = ["WLYELF", "FWKLLI", "WIYKLW", "LFWEYL", "IWLKFV"]
strand_pairs = [
    ("KTLTVE", "KVTITK"),
    ("SVEITV", "TVKVSD"),
    ("EVKLTA", "TLKVEG"),
]
di = 0
for (s1, s2) in strand_pairs:
    for loop in loops[:3]:
        for turn in turns[:2]:
            di += 1
            # C - strand1 - loop(hotspot) - turn - strand2 - C  (disulfide C..C)
            seq = "AC" + s1 + loop + turn + s2 + "CG"
            add(f"deD_mac_{di}", seq, "D_macrocycle",
                f"disulfide hairpin; loop={loop}; turn={turn}")

# a few longer double-hairpin macrocycles
add("deD_big1", "ACKTLTVEWLYELFNGKVTITKCGSPGACSVEITVFWKLLIPGTVKVSDCG", "D_macrocycle", "bis-hairpin stapled")
add("deD_big2", "ACEVKLTAWIYKLWDGTLKVEGCG", "D_macrocycle", "compact stapled hairpin")

# ---------------------------------------------------------------------------
summary = {}
for c in cands:
    summary[c["family"]] = summary.get(c["family"], 0) + 1

json.dump(cands, open("design/lib_structure_independent.json", "w"), indent=1)
print("Generated structure-independent candidates:", len(cands))
for k, v in sorted(summary.items()):
    print(f"  {k}: {v}")
print("Length range:", min(len(c['seq']) for c in cands), "-", max(len(c['seq']) for c in cands))
