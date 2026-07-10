"""The common currency: a signed, ranked perturbation signature over genes.

A ``Signature`` is nothing more than a mapping gene -> score, where the sign
means direction (up = induced by the perturbation, down = repressed) and the
magnitude means strength/confidence. Every data source in this skill is coerced
into this one representation:

  * Drug Repurposing Hub / CMap L1000  -> moderated z-score per landmark gene
  * CRISPR screen (MAGeCK / DESeq2)     -> log2 fold-change (or -log10 p * sign)
  * single-cell perturbation DGE        -> DESeq2 stat / shrunken log2FC
  * Geneformer in-silico perturbation   -> cosine shift per gene (see geneformer.md)

Keeping one representation is what makes cross-modality connectivity possible.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd


@dataclass
class Signature:
    """A signed ranked gene signature.

    Parameters
    ----------
    scores : mapping gene -> float
        Signed score. Positive = up-regulated by the perturbation.
    name : str
        Human label (e.g. ``"vorinostat_10uM_MCF7"`` or ``"CRISPRi_TP53"``).
    modality : str
        Free-text tag, e.g. ``"drug"``, ``"crispr_ko"``, ``"crispr_i"``,
        ``"crispr_a"``, ``"scrna_dge"``, ``"geneformer"``, ``"disease"``.
    meta : dict
        Anything else worth carrying (cell line, dose, MOA, target, phase...).
    """

    scores: pd.Series
    name: str = "signature"
    modality: str = "unknown"
    meta: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.scores, pd.Series):
            self.scores = pd.Series(dict(self.scores), dtype=float)
        self.scores = self.scores.astype(float)
        self.scores.index = self.scores.index.astype(str)
        # collapse duplicate gene ids by mean, drop non-finite
        self.scores = self.scores[np.isfinite(self.scores.values)]
        if self.scores.index.has_duplicates:
            self.scores = self.scores.groupby(level=0).mean()

    # ---- constructors ----------------------------------------------------

    @classmethod
    def from_deseq2(
        cls,
        result: pd.DataFrame,
        stat_col: str = "stat",
        gene_col: str | None = None,
        **kw,
    ) -> "Signature":
        """Build from a DESeq2 / PyDESeq2 results table (RIGOROUS input).

        Uses the Wald ``stat`` column by default — it already folds the
        log2 fold-change and its uncertainty into one signed, ranked number,
        which is exactly what a signature wants. Fall back to
        ``log2FoldChange`` only if ``stat`` is absent.
        """
        df = result.copy()
        if gene_col is not None:
            df = df.set_index(gene_col)
        if stat_col not in df.columns:
            stat_col = "log2FoldChange"
        return cls(scores=df[stat_col], **kw)

    @classmethod
    def from_ranked(
        cls, up: Sequence[str], down: Sequence[str] = (), **kw
    ) -> "Signature":
        """Build from two gene lists (up / down). Scores are +1 / -1.

        Handy for CMap-style discrete query sets or published gene lists.
        """
        s = pd.Series(
            {g: 1.0 for g in up} | {g: -1.0 for g in down}, dtype=float
        )
        return cls(scores=s, **kw)

    @classmethod
    def from_l1000(
        cls, zscores: Mapping[str, float] | pd.Series, **kw
    ) -> "Signature":
        """Build from L1000 / CMap moderated z-scores (one value per gene)."""
        return cls(scores=pd.Series(dict(zscores), dtype=float), **kw)

    # ---- views -----------------------------------------------------------

    def top(self, k: int = 50, direction: str = "up") -> list[str]:
        """The k strongest genes. ``direction`` in {up, down, abs}."""
        if direction == "up":
            s = self.scores.sort_values(ascending=False)
        elif direction == "down":
            s = self.scores.sort_values(ascending=True)
        elif direction == "abs":
            s = self.scores.reindex(self.scores.abs().sort_values(
                ascending=False).index)
        else:
            raise ValueError(direction)
        return list(s.head(k).index)

    def query_sets(self, k: int = 50) -> tuple[list[str], list[str]]:
        """The (up, down) tag sets used as a CMap-style query."""
        return self.top(k, "up"), self.top(k, "down")

    def align(self, genes: Iterable[str], fill: float = 0.0) -> pd.Series:
        """Reindex onto a shared gene universe (0 for missing)."""
        return self.scores.reindex(list(genes)).fillna(fill)

    def ranked_index(self) -> pd.Series:
        """Genes sorted high->low; the ranked reference list for GSEA/WTCS."""
        return self.scores.sort_values(ascending=False)

    def __len__(self) -> int:
        return len(self.scores)

    def __repr__(self) -> str:
        return (f"Signature({self.name!r}, modality={self.modality!r}, "
                f"n_genes={len(self)})")
