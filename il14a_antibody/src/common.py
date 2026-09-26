"""Shared constants and helpers for the IL-14alpha blocking-antibody campaign.

All sequence/structure facts here are traceable to a primary source:
  - TXLNA / IL-14alpha canonical sequence : UniProt P40222 (546 aa)
  - beta-taxilin                          : UniProt Q8N3L3 (684 aa)
  - gamma-taxilin                         : UniProt Q9NUQ3 (528 aa)
  - structure                             : AlphaFold DB AF-P40222-F1 model v6
"""
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data"
RESULTS = Path(__file__).resolve().parent.parent / "results"
FIGURES = Path(__file__).resolve().parent.parent / "figures"

# --- IL-14alpha biology -----------------------------------------------------
# The secreted cytokine IL-14alpha arises from an alternate translation start
# at Met176 of TXLNA, so ONLY residues 176-546 are present on the
# extracellular/secreted species that a systemic IgG can reach.
SECRETED_START = 176
SECRETED_END = 546

# Zone claimed by US7622574 (mAb 1C6/1F2). Kept free for freedom-to-operate.
PATENT_ZONE = (493, 523)

# Epitopes proposed by the earlier G034 pipeline, retained for comparison.
LEGACY_EPITOPES = {
    "G034-EPI-A": (1, 17),
    "G034-EPI-B": (58, 98),
    "G034-EPI-C": (358, 384),
    "G034-EPI-D": (425, 445),
}

# --- amino-acid scales ------------------------------------------------------
# Tien et al. 2013 (PLoS ONE 8:e80635) theoretical max side-chain+backbone SASA
MAX_ASA = {
    "A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0,
    "Q": 225.0, "E": 223.0, "G": 104.0, "H": 224.0, "I": 197.0,
    "L": 201.0, "K": 236.0, "M": 224.0, "F": 240.0, "P": 159.0,
    "S": 155.0, "T": 172.0, "W": 285.0, "Y": 263.0, "V": 174.0,
}

KD_HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5,
    "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
    "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
    "Y": -1.3, "V": 4.2,
}

CHARGE_PH74 = {"R": 1.0, "K": 1.0, "H": 0.1, "D": -1.0, "E": -1.0}


def read_fasta(path):
    """Return the single sequence in a FASTA file as an uppercase string."""
    lines = Path(path).read_text().strip().split("\n")
    return "".join(l.strip() for l in lines if not l.startswith(">")).upper()


def net_charge(seq, ph74=True):
    table = CHARGE_PH74 if ph74 else {"R": 1, "K": 1, "H": 0.5, "D": -1, "E": -1}
    return sum(table.get(c, 0.0) for c in seq)


def gravy(seq):
    return sum(KD_HYDROPATHY[c] for c in seq) / len(seq) if seq else 0.0


def load_taxilins():
    """Return (alpha, beta, gamma) taxilin sequences."""
    return (
        read_fasta(DATA / "P40222.fasta"),
        read_fasta(DATA / "Q8N3L3.fasta"),
        read_fasta(DATA / "Q9NUQ3.fasta"),
    )
