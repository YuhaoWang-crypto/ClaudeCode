"""Amino-acid reference data used by the TCR triage pipeline.

Every table here is a literature constant, not a fitted parameter. Sources are
named inline so that any number can be traced back.
"""

from __future__ import annotations

AA20 = "ACDEFGHIKLMNPQRSTVWY"

# Residue volumes in cubic angstroms. Zamyatnin, Annu Rev Biophys 1984.
VOLUME = {
    "A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5,
    "Q": 143.8, "E": 138.4, "G": 60.1, "H": 153.2, "I": 166.7,
    "L": 166.7, "K": 168.6, "M": 162.9, "F": 189.9, "P": 112.7,
    "S": 89.0, "T": 116.1, "W": 227.8, "Y": 193.6, "V": 140.0,
}

# Kyte-Doolittle hydropathy. Kyte & Doolittle, J Mol Biol 1982.
HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# Formal side-chain charge at pH 7.4. Histidine is given a partial value
# because its pKa (~6.0) leaves it mostly neutral but titratable.
CHARGE = {
    "D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0, "H": 0.1,
    **{a: 0.0 for a in AA20 if a not in "DEKRH"},
}

# Theoretical maximum side-chain + backbone SASA for Gly-X-Gly.
# Tien et al., PLoS ONE 2013 (theoretical column).
MAX_SASA = {
    "A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0,
    "Q": 225.0, "E": 223.0, "G": 104.0, "H": 224.0, "I": 197.0,
    "L": 201.0, "K": 236.0, "M": 224.0, "F": 240.0, "P": 159.0,
    "S": 155.0, "T": 172.0, "W": 285.0, "Y": 263.0, "V": 174.0,
}

# Count of side-chain hydrogen-bond donor + acceptor groups.
HBOND_CAPACITY = {
    "R": 3, "K": 1, "W": 1, "N": 2, "Q": 2, "H": 2, "S": 1,
    "T": 1, "Y": 1, "C": 1, "D": 2, "E": 2,
    **{a: 0 for a in AA20 if a not in "RKWNQHSTYCDE"},
}

AROMATIC = set("FWYH")

# Side-chain atoms that carry the formal charge, used for electrostatic
# complementarity. Backbone atoms are excluded on purpose.
CHARGED_SIDECHAIN_ATOMS = {
    "ASP": ("OD1", "OD2"),
    "GLU": ("OE1", "OE2"),
    "LYS": ("NZ",),
    "ARG": ("NE", "NH1", "NH2"),
    "HIS": ("ND1", "NE2"),
}

# Side-chain nitrogen/oxygen atoms available for hydrogen bonding.
POLAR_SIDECHAIN_ATOMS = {
    "ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2"),
    "ASN": ("OD1", "ND2"), "GLN": ("OE1", "NE2"),
    "SER": ("OG",), "THR": ("OG1",), "TYR": ("OH",),
    "LYS": ("NZ",), "ARG": ("NE", "NH1", "NH2"),
    "HIS": ("ND1", "NE2"), "TRP": ("NE1",), "CYS": ("SG",),
}

BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}

# Peptide anchor preferences for HLA-A*02:01 (P2 and P-Omega). Compiled from
# the SYFPEITHI/IEDB consensus motif; used only as a coarse pre-filter in the
# off-target scan, never as a substitute for NetMHCpan.
HLA_A0201_ANCHORS = {2: set("LMIVQAT"), -1: set("VLIAMT")}

# Sequence liability motifs relevant to a therapeutic TCR construct.
# Deamidation/isomerisation motifs: Robinson & Robinson, PNAS 2001.
LIABILITY_MOTIFS = {
    "N-glycosylation sequon": r"N[^P][ST]",
    "Asn deamidation (NG)": r"NG",
    "Asp isomerisation (DG)": r"DG",
    "Asp isomerisation (DP)": r"DP",
    "Unpaired cysteine": r"C",
    "Met oxidation": r"M",
    "Trp oxidation": r"W",
}
