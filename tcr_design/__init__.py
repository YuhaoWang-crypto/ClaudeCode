"""Structure-based triage for TCR affinity-maturation candidates.

This package measures a TCR-pMHC interface, enumerates single mutants over the
contacting CDR positions, gates them on specificity and developability risk,
screens the target epitope against the human proteome, and writes the input
files for the tools that can actually score the survivors.

It deliberately does not predict affinity. See `scan.py` for why.
"""

from .liability import apply_gates, assess_liabilities
from .offtarget import peptide_contact_profile, scan_proteome
from .scan import Mutation, scan_complex, score_mutation
from .structure import Complex, Residue, load_complex

__version__ = "0.1.0"

__all__ = [
    "Complex",
    "Mutation",
    "Residue",
    "apply_gates",
    "assess_liabilities",
    "load_complex",
    "peptide_contact_profile",
    "scan_complex",
    "scan_proteome",
    "score_mutation",
]
