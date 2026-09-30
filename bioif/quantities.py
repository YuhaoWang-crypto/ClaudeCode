"""
bioif.quantities -- a small shared vocabulary, plus the lossless coercions
between its members.

This file is deliberately short. The point is not to standardise all of
biology; it is that a *small* agreed vocabulary at the seams, with automatic
algebra between equivalent encodings, removes most of the glue code that
makes model swapping expensive.
"""
from __future__ import annotations

import math

from .adapter import Coercion
from .core import Quantity

# --- chemistry -------------------------------------------------------------
#: pIC50, pKi, pKd and pEC50 are DIFFERENT quantities. ChEMBL pools all four
#: into one `pchembl_value` column, which is convenient and is exactly the
#: implicit pooling the contract exists to catch: a functional IC50 and an
#: equilibrium Kd are not interchangeable, and converting between them needs
#: information (substrate concentration, Km, mechanism) that an affinity
#: record does not carry.
PIC50 = Quantity("target_affinity", "pIC50", "log10", "", (2.0, 12.0))
PKI = Quantity("target_affinity", "pKi", "log10", "", (2.0, 12.0))
PKD = Quantity("target_affinity", "pKd", "log10", "", (2.0, 12.0))
PEC50 = Quantity("target_affinity", "pEC50", "log10", "", (2.0, 12.0))

#: ChEMBL `standard_type` -> the quantity it actually is.
BY_STANDARD_TYPE = {"IC50": PIC50, "Ki": PKI, "Kd": PKD, "EC50": PEC50}
IC50_NM = Quantity("target_affinity", "nM", "linear", "", (1e-3, 1e9))
OCCUPANCY = Quantity("target_occupancy", "fraction", "linear",
                     "unbound target", (0.0, 1.0))

# --- the convergence node --------------------------------------------------
#: Both the chemical arm and the genetic arm are forced through this quantity.
#: That is the whole trick for cross-modality reuse: a model trained on
#: genetic perturbation can only score a compound if the compound's effect is
#: first expressed in the same currency the genetic model was trained in.
RESIDUAL_ACTIVITY = Quantity("residual_target_activity", "fraction", "linear",
                             "untreated wild-type", (0.0, 1.0))

# --- genetics / transcriptomics -------------------------------------------
DELTA_PSI = Quantity("delta_psi", "fraction", "linear",
                     "reference allele", (-1.0, 1.0))
RNA_LFC = Quantity("transcript_abundance", "log2_ratio", "log2",
                   "isogenic wild-type", (-12.0, 12.0))
PROTEIN_LFC = Quantity("protein_abundance", "log2_ratio", "log2",
                       "isogenic wild-type", (-12.0, 12.0))

# --- phenotype -------------------------------------------------------------
FITNESS = Quantity("fitness_effect", "gene_effect", "linear",
                   "non-targeting control", (-3.0, 1.0))


# --------------------------------------------------------------------------
# Lossless coercions. These are pure algebra: the registry may insert them
# without asking, and they contribute exactly zero variance.
# --------------------------------------------------------------------------

def _pic50_to_nm(x, _claim):
    return 10.0 ** (9.0 - x)


def _nm_to_pic50(x, _claim):
    return 9.0 - math.log10(max(x, 1e-12))


def _ic50_to_occupancy(x_nm, claim):
    """
    Hill occupancy at the dose carried in the claim's context (n = 1).

    Note this is a *coercion* only because the dose is an explicit input: it
    is arithmetic once you have both numbers. The adapter that uses it must
    declare `requires_context=('dose_uM',)` so a missing dose refuses rather
    than defaulting.
    """
    dose_nm = (claim.context.dose_uM or 0.0) * 1000.0
    return dose_nm / (dose_nm + x_nm) if (dose_nm + x_nm) > 0 else 0.0


PIC50_TO_IC50 = Coercion("pIC50->IC50nM", PIC50, IC50_NM, _pic50_to_nm)
IC50_TO_PIC50 = Coercion("IC50nM->pIC50", IC50_NM, PIC50, _nm_to_pic50)
IC50_TO_OCCUPANCY = Coercion("IC50nM->occupancy", IC50_NM, OCCUPANCY,
                             _ic50_to_occupancy)
IC50_TO_OCCUPANCY.requires_context = ("dose_uM",)

LOSSLESS = [PIC50_TO_IC50, IC50_TO_PIC50, IC50_TO_OCCUPANCY]
