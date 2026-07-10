"""perturbomics — combine drug, CRISPR and in-silico perturbation signatures.

The unifying object across every modality (L1000/CMap drug profiles, CRISPR
knockout/interference/activation screens, single-cell perturbation DGE, and
Geneformer in-silico perturbation) is a *signature*: a signed, ranked vector
over genes. Once every perturbation is a Signature you can:

  * score connectivity (signature similarity) between ANY two perturbations —
    drug<->drug, drug<->CRISPR, perturbation<->disease   (connectivity.py)
  * mine drug + CRISPR combinations that jointly reverse a disease state
    (combine.py)
  * compute those signatures rigorously from raw single-cell counts
    (pseudobulk.py)

Rigor discipline (inherited from this repo): every function docstring marks
whether its output is a deterministic computation (rigorous) or a biological
interpretation that needs experimental validation (hypothesis).
"""

from .signature import Signature
from .connectivity import (
    enrichment_score,
    weighted_connectivity_score,
    connectivity_matrix,
    normalized_connectivity,
)
from .combine import rank_reversers, best_combinations, combination_score
from .pseudobulk import pseudobulk, signature_from_pseudobulk

__all__ = [
    "Signature",
    "enrichment_score",
    "weighted_connectivity_score",
    "connectivity_matrix",
    "normalized_connectivity",
    "rank_reversers",
    "best_combinations",
    "combination_score",
    "pseudobulk",
    "signature_from_pseudobulk",
]
