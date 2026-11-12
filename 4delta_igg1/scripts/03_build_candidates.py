#!/usr/bin/env python3
"""Task 2: generate the 16 grafted candidates (4 sites x 4 graft strategies).

Design driver (from script 01): in its bioactive conformation 4-Delta is a
LINEAR pharmacophore - a mobile RGD at the N-terminal end followed by a rigid
amphipathic helix - whose two termini sit 35.8 +/- 2.1 A apart. Every target
site is a CLOSED loop whose anchors are only 6.9-13.0 A apart. So 23-30 A of
slack has to be absorbed by junction design, and the four strategies are four
different answers to that one problem.

  S1 minimal     GS | 4D | SG                      naive insertion, no slack budget
  S2 spanmatch   GS | 4D | (GGGGS)x2               C-return linker sized to the slack
  S3 replace     native loop deleted, N-cap | 4D | GPG+(GGGGS)  max room, helix capped
  S4 coretrunc   GS | 4D(1-18) | (GGGGS)x2         rigid helix core only, least strain

Controls: WT scaffolds, free 4-Delta, and a CH3 C-terminal fusion that has no
closed-loop span penalty at all (topology benchmark / ceiling).
"""
import json

FOURD = "FETLRGDLRILSILRHQNLLKEL"          # CgA(39-61) E46L
FOURD_CORE = FOURD[:18]                      # F39..N56 - keeps the fully helical core
assert FOURD_CORE == "FETLRGDLRILSILRHQN", FOURD_CORE

# canonical human IgG1 CH3, EU 341-447
CH3 = ("GQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGS"
       "FFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK")
CH3_LO = 341
assert len(CH3) == 447 - CH3_LO + 1 == 107
def eu(i):
    """residue at EU position i"""
    return CH3[i - CH3_LO]
# sanity-check the EU frame against the loop definitions
assert eu(355) == "R" and eu(362) == "Q", "AB loop frame"
assert eu(383) == "S" and eu(391) == "Y", "CD loop frame"
assert eu(413) == "D" and eu(422) == "V", "EF loop frame"
assert eu(447) == "K"

# trastuzumab Fv (reference IgG1 scaffold; client sequences were not supplied)
VH = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTI"
      "SADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS")
VL = ("DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTD"
      "FTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK")
CDRH3 = "WGGDGFYAMDY"
assert VH.count(CDRH3) == 1
cdr3_start = VH.index(CDRH3)                 # 0-based; loop = VH[cdr3_start:+11]

G4S = "GGGGS"

def junctions(strategy, payload, native_loop):
    """Return (n_junction, payload, c_junction, keep_native_loop)."""
    if strategy == "S1_minimal":
        return "GS", payload, "SG", True
    if strategy == "S2_spanmatch":
        return "GS", payload, G4S * 2, True
    if strategy == "S3_replace":
        # native loop removed -> maximum room; D/P N-cap nucleates the 4D helix,
        # GPG breaks it before the acceptor beta-strand resumes
        return "GD", payload, "GPG" + G4S, False
    if strategy == "S4_coretrunc":
        return "GS", FOURD_CORE, G4S * 2, True
    raise ValueError(strategy)

STRATEGIES = ["S1_minimal", "S2_spanmatch", "S3_replace", "S4_coretrunc"]
SITES = {
    # name:      (loop_lo, loop_hi) in EU numbering, for the CH3 sites
    "CH3-AB": (355, 362),
    "CH3-CD": (383, 391),
    "CH3-EF": (413, 422),
}
ANCHOR_CA = {"CDRH3": 6.86, "CH3-AB": 8.47, "CH3-CD": 12.99, "CH3-EF": 9.78}
FOURD_SPAN = 35.75

candidates = []
cid = 0

for site in ["CDRH3", "CH3-AB", "CH3-CD", "CH3-EF"]:
    for strat in STRATEGIES:
        cid += 1
        name = f"C{cid:02d}_{site}_{strat.split('_')[0]}"
        if site == "CDRH3":
            native = CDRH3
            nj, pay, cj, keep = junctions(strat, FOURD, native)
            insert = nj + pay + cj
            if keep:
                new_loop = native[:4] + insert + native[4:]   # inserted mid-loop
                placement = "inserted after VH CDRH3 position 4 (apex), native loop retained"
            else:
                new_loop = insert
                placement = "native CDRH3 replaced in full"
            vh_var = VH[:cdr3_start] + new_loop + VH[cdr3_start + len(CDRH3):]
            chains = {"H": vh_var, "L": VL}
            module = "Fv"
            delta_len = len(vh_var) - len(VH)
            graft_lo = cdr3_start + (4 if keep else 0)       # 0-based in the variant
            graft_hi = graft_lo + len(insert) - 1
            loop_desc = f"VH CDRH3 {CDRH3}"
        else:
            lo, hi = SITES[site]
            native = "".join(eu(i) for i in range(lo, hi + 1))
            nj, pay, cj, keep = junctions(strat, FOURD, native)
            insert = nj + pay + cj
            pre = CH3[: lo - CH3_LO]
            post = CH3[hi - CH3_LO + 1 :]
            if keep:
                # insert at the loop apex, splitting the native loop in half
                h = len(native) // 2
                new_loop = native[:h] + insert + native[h:]
                placement = (f"inserted at apex of {site} loop (EU {lo}-{hi}), "
                             f"after EU {lo + h - 1}; native loop retained")
                graft_lo = len(pre) + h
            else:
                new_loop = insert
                placement = f"native {site} loop (EU {lo}-{hi}) replaced in full"
                graft_lo = len(pre)
            ch3_var = pre + new_loop + post
            chains = {"A": ch3_var, "B": ch3_var}
            module = "CH3 homodimer"
            delta_len = len(ch3_var) - len(CH3)
            graft_hi = graft_lo + len(insert) - 1
            loop_desc = f"CH3 {site} EU {lo}-{hi} {native}"

        payload_len = len(pay)
        c_linker_len = len(cj)
        # straight-chain reach of the C-return linker, ~3.0 A per residue for Gly/Ser
        reach = 3.0 * c_linker_len
        slack = FOURD_SPAN - ANCHOR_CA[site]
        candidates.append({
            "id": name,
            "site": site,
            "strategy": strat,
            "module": module,
            "native_loop": loop_desc,
            "placement": placement,
            "payload": pay,
            "payload_len": payload_len,
            "n_junction": nj,
            "c_junction": cj,
            "insert": insert,
            "insert_len": len(insert),
            "length_change_vs_wt": delta_len,
            "graft_index_range_0based": [graft_lo, graft_hi],
            "chains": chains,
            "geometry_budget": {
                "anchor_CA_A": ANCHOR_CA[site],
                "payload_end_to_end_A": FOURD_SPAN if payload_len == 23 else None,
                "slack_to_absorb_A": round(slack, 2),
                "c_linker_residues": c_linker_len,
                "c_linker_max_reach_A": round(reach, 1),
                "slack_covered": bool(reach >= slack),
            },
        })

# ----------------------------------------------------------------- controls
controls = [
    {"id": "WT_Fv", "site": "none", "strategy": "control", "module": "Fv",
     "chains": {"H": VH, "L": VL}, "note": "unmodified trastuzumab Fv baseline"},
    {"id": "WT_CH3", "site": "none", "strategy": "control", "module": "CH3 homodimer",
     "chains": {"A": CH3, "B": CH3}, "note": "unmodified human IgG1 CH3 homodimer baseline"},
    {"id": "REF_4delta_free", "site": "none", "strategy": "control", "module": "peptide",
     "chains": {"P": FOURD},
     "note": "free 4-Delta peptide; validates the pipeline against PDB 6R2X"},
    {"id": "BM_CH3_Cterm_fusion", "site": "CH3-Cterm", "strategy": "benchmark",
     "module": "CH3 homodimer", "chains": {"A": CH3 + G4S + FOURD, "B": CH3 + G4S + FOURD},
     "note": ("topology benchmark: 4-Delta appended to the free C-terminus (K447) via "
              "G4S. No closed-loop span penalty - shows the attainable ceiling for RGD "
              "presentation and is NOT one of the 16 requested designs."),
     "graft_index_range_0based": [len(CH3) + len(G4S), len(CH3) + len(G4S) + len(FOURD) - 1]},
]

out = {
    "scaffold": {
        "antibody": "trastuzumab (reference human IgG1; client HC/LC sequences not supplied)",
        "VH": VH, "VL": VL, "CH3_341_447": CH3,
        "assumption_note": ("All designs are expressed as EU-numbered edits so they transfer "
                            "unchanged to the client's own IgG1 once its sequences arrive."),
    },
    "peptide": {"4delta": FOURD, "4delta_core_1_18": FOURD_CORE,
                "bioactive_end_to_end_A": FOURD_SPAN},
    "strategies": {
        "S1_minimal": "GS + payload + SG; native loop kept. No slack budget - naive insertion.",
        "S2_spanmatch": "GS + payload + (GGGGS)x2; native loop kept. C-return linker sized to the slack.",
        "S3_replace": "GD + payload + GPG+(GGGGS); native loop deleted. Max room, helix capped, beta-strands protected.",
        "S4_coretrunc": "GS + 4delta(1-18) + (GGGGS)x2; native loop kept. Rigid helix core only.",
    },
    "candidates": candidates,
    "controls": controls,
}

with open("results/03_candidates.json", "w") as fh:
    json.dump(out, fh, indent=2)

print(f"{len(candidates)} candidates + {len(controls)} controls\n")
hdr = f"{'id':<22}{'site':<9}{'strategy':<14}{'ins':>4}{'dLen':>6}{'slack':>7}{'reach':>7}  ok"
print(hdr); print("-" * len(hdr))
for c in candidates:
    g = c["geometry_budget"]
    print(f"{c['id']:<22}{c['site']:<9}{c['strategy']:<14}{c['insert_len']:>4}"
          f"{c['length_change_vs_wt']:>6}{g['slack_to_absorb_A']:>7.1f}"
          f"{g['c_linker_max_reach_A']:>7.1f}  {'Y' if g['slack_covered'] else 'N'}")
print("\ncontrols:", ", ".join(c["id"] for c in controls))
print("wrote results/03_candidates.json")
