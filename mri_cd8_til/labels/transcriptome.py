"""Transcriptomic immune labels for the large tier-B cohort (I-SPY2, n=717).

Three published, transparent signatures are implemented.  All are simple
rank/z-score aggregates rather than a fitted deconvolution, deliberately: a
fitted deconvolution (CIBERSORTx and friends) adds another set of parameters
that would themselves need to be kept out of the validation fold, and buys
little for a signature this well characterised.

* ``CD8_CORE``      -- CD8 lineage plus cytotoxic effectors.
* ``CYTOLYTIC``     -- Rooney et al. 2015 cytolytic activity (GZMA, PRF1),
  defined there as the geometric mean of the two.
* ``TIS``           -- the 18-gene Tumour Inflammation Signature (Ayers et al.
  2017), the T-cell-inflamed phenotype used for pembrolizumab response.

The labels these produce support the ``inflamed_vs_rest`` contrast only; see
``mri_cd8_til.labels.schema``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .schema import TRANSCRIPTOMIC, Immunophenotype, LabelProvenance, LabelTask

CD8_CORE = (
    "CD8A",
    "CD8B",
    "GZMA",
    "GZMB",
    "GZMK",
    "PRF1",
    "CD3D",
    "CD3E",
    "CD2",
    "IFNG",
)

CYTOLYTIC = ("GZMA", "PRF1")

#: Ayers et al., J Clin Invest 2017;127(8):2930-2940.
TIS = (
    "CCL5",
    "CD27",
    "CD274",
    "CD276",
    "CD8A",
    "CMKLR1",
    "CXCL9",
    "CXCR6",
    "HLA-DQA1",
    "HLA-DRB1",
    "HLA-E",
    "IDO1",
    "LAG3",
    "NKG7",
    "PDCD1LG2",
    "PSMB10",
    "STAT1",
    "TIGIT",
)

SIGNATURES: dict[str, tuple[str, ...]] = {
    "cd8_core": CD8_CORE,
    "cytolytic": CYTOLYTIC,
    "tis": TIS,
}


@dataclass(frozen=True)
class SignatureScores:
    """Per-case immune-signature scores."""

    case_ids: tuple[str, ...]
    scores: dict[str, np.ndarray]

    def __getitem__(self, name: str) -> np.ndarray:
        return self.scores[name]

    def as_matrix(self) -> tuple[tuple[str, ...], np.ndarray]:
        names = tuple(sorted(self.scores))
        return names, np.column_stack([self.scores[n] for n in names])


def score_signatures(
    expression: np.ndarray,
    gene_names: list[str] | tuple[str, ...],
    case_ids: list[str] | tuple[str, ...],
    signatures: dict[str, tuple[str, ...]] | None = None,
    already_log2: bool = True,
) -> SignatureScores:
    """Score immune signatures from a genes-by-cases expression matrix.

    Parameters
    ----------
    expression:
        ``(n_genes, n_cases)`` matrix.  I-SPY2's GSE194040 gene-level table is
        already log2 and ComBat-adjusted across its two array platforms, which
        is why ``already_log2`` defaults to True.
    gene_names:
        Row labels of ``expression``.

    Notes
    -----
    Genes are z-scored **across cases** before aggregation.  That makes the
    score a within-cohort relative measure: it is comparable across patients in
    one cohort and *not* comparable across cohorts without re-standardisation.
    Treating such a score as an absolute value is a common way to manufacture
    apparent cross-cohort agreement.
    """
    signatures = signatures or SIGNATURES
    matrix = np.asarray(expression, dtype=float)
    if not already_log2:
        matrix = np.log2(matrix + 1.0)
    index = {g.upper(): i for i, g in enumerate(gene_names)}

    scores: dict[str, np.ndarray] = {}
    for name, genes in signatures.items():
        rows = [index[g] for g in genes if g in index]
        if not rows:
            raise KeyError(f"signature {name!r}: none of its genes are in the matrix")
        sub = matrix[rows, :]
        z = _zscore_rows(sub)
        if name == "cytolytic":
            # Rooney et al. define cytolytic activity as a geometric mean of
            # the raw (non-z-scored) expression values.
            scores[name] = np.exp(np.log(np.clip(sub, 1e-9, None)).mean(axis=0))
        else:
            scores[name] = z.mean(axis=0)
    return SignatureScores(case_ids=tuple(case_ids), scores=scores)


def _zscore_rows(matrix: np.ndarray) -> np.ndarray:
    mu = matrix.mean(axis=1, keepdims=True)
    sd = matrix.std(axis=1, ddof=1, keepdims=True)
    sd[sd == 0] = 1.0
    return (matrix - mu) / sd


@dataclass(frozen=True)
class InflamedThreshold:
    """Cut-point separating T-cell-inflamed from non-inflamed tumours."""

    signature: str
    cut: float
    quantile: float

    @classmethod
    def fit(
        cls, scores: SignatureScores, signature: str = "tis", quantile: float = 0.60
    ) -> "InflamedThreshold":
        """Fit on **training cases only** -- see ``labels.spatial`` for why."""
        values = scores[signature]
        return cls(signature=signature, cut=float(np.quantile(values, quantile)), quantile=quantile)

    def apply(self, scores: SignatureScores) -> np.ndarray:
        """Binary ``inflamed_vs_rest`` labels for the ``RM-whole`` task."""
        return (scores[self.signature] >= self.cut).astype(int)


def pseudo_phenotype(scores: SignatureScores, threshold: InflamedThreshold) -> list[Immunophenotype]:
    """Three-class labels with the desert/excluded split explicitly refused.

    Non-inflamed cases are returned as ``DESERT`` purely as a placeholder.  The
    function exists so that calling code fails loudly if it tries to feed
    transcriptomic labels into the ``desert_vs_excluded`` task: it must first
    consult :func:`supported_tasks`, which reports that task as unsupported.
    """
    inflamed = threshold.apply(scores).astype(bool)
    return [Immunophenotype.INFLAMED if flag else Immunophenotype.DESERT for flag in inflamed]


def supported_tasks() -> tuple[LabelTask, ...]:
    return TRANSCRIPTOMIC.supports


def provenance() -> LabelProvenance:
    return TRANSCRIPTOMIC
