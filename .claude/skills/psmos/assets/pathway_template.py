"""
Template: add a new pathway to psmos/pathways.py.

Copy the Pathway below into pathways.py, fill in the blanks, register it in the
PATHWAYS dict, then run:  python3 -m psmos.run_all <Key>

Worked example sketched here: Wnt/β-catenin (regeneration ↔ fibrosis).
The whole orthologs → cds → paralogs → alphagenome_r → evo2 → scoring pipeline
is reused automatically; you only declare the biology.
"""

# --- 1. If the pathway needs organisms not already in SPECIES, add them there:
#
# "planaria": Species("Schmidtea mediterranea (planaria)", 79327, 0.76, 0.80, "flatworm"),
#   (name, NCBI taxon id, tract prior 0..1, thru prior 0..1, clade)
# If the organism is on Ensembl, also add it to cds.ENSEMBL_SPECIES so it gets
# CDS + Compara + Evo2; otherwise it is an honest CDS/Evo2 gap.

WNT = Pathway(
    key="Wnt",
    label="Wnt/β-catenin signalling",
    reference_species="human",
    families=[
        Family("wnt_ligand", "Wnt ligand", ["WNT1", "WNT3A", "WNT5A"]),          # extend
        Family("frizzled", "Frizzled receptor", ["FZD1", "FZD5", "FZD7"]),        # extend
        Family("lrp", "LRP5/6 co-receptor", ["LRP5", "LRP6"]),
        Family("dvl", "Dishevelled", ["DVL1", "DVL2", "DVL3"]),
        Family("destruction", "APC/AXIN/GSK3 destruction complex",
               ["APC", "AXIN1", "AXIN2", "GSK3B", "CSNK1A1"]),
        # GATE: β-catenin effector + TCF/LEF DNA-binding partner. No canonical
        # transcriptional output without both.
        Family("bcatenin", "β-catenin (Armadillo) effector", ["CTNNB1"], is_gate=True),
        Family("tcf_lef", "TCF/LEF transcription factor",
               ["TCF7", "TCF7L1", "TCF7L2", "LEF1"], is_gate=True),
    ],
    # Per-species gene symbols for the GATE families only, as fallback lists.
    # Use each species' own gene names; add true negative controls by querying
    # the human symbol so the search genuinely returns nothing.
    ortholog_seed={
        "bcatenin": {
            "human":     ["CTNNB1"],
            "mouse":     ["Ctnnb1"],
            "zebrafish": ["ctnnb1", "ctnnb2"],
            "fly":       ["arm"],            # armadillo
            "planaria":  ["ctnnb1", "Smed-bcatenin1"],
            # e.g. "yeast": ["CTNNB1"],  # negative control (expect zero hits)
        },
        "tcf_lef": {
            "human":     ["TCF7L2"],
            "mouse":     ["Tcf7l2"],
            "zebrafish": ["tcf7l2", "tcf7l1a"],
            "fly":       ["pan"],            # pangolin (TCF)
            "planaria":  ["tcf", "Smed-tcf"],
        },
    },
)

# --- 2. Register it:
# PATHWAYS = {p.key: p for p in [NOTCH, HIPPO, WNT]}

# --- 3. Curated baseline for scoring:
#   * Lighter Notch-style path: add a CURATED_WNT row-dict in scoring.py.
#   * PSMOS six-layer path: add G/D/N/R/E/X + arch-role fits per species in
#     scoring_psmos.py (mirror CURATED_HIPPO), plus a build_wnt_dashboard.py.
#
# Everything else (gate, CDS, Compara redundancy, AlphaGenome R for human/mouse,
# Evo2 constraint) is computed with no extra code.
