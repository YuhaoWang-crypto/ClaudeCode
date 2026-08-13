"""The public interface: predict knockdown responses in a new cell line.

Everything else in this package is benchmarking machinery.  This is the part you
call.

    from virtualcell.predict import VirtualCell

    vc = VirtualCell.from_atlas()                    # trains on the four lines
    out = vc.predict("my_controls.h5ad", ["TP53", "MYC", "RPL5"])
    out.expression      # (n_knockdowns, n_genes) predicted profile, log1p CP10K
    out.effect          # (n_knockdowns, n_genes) change from control
    out.top_genes("TP53", 20)                        # ranked differential genes

**Input.** Non-targeting control cells of the cell line you care about -- an
``.h5ad``, an ``AnnData``, or a mean expression vector -- plus the gene symbols
to knock down.  Nothing else.  No perturbation data from that cell line is
needed or used; that is the point.

**Output.** A :class:`Prediction`: the post-perturbation profile, the effect
relative to control, and a ranked differential-expression list per knockdown.

**What it is good for, and what it is not.**  Measured on four held-out CRISPRi
cell lines, this recovers the component of a knockdown response that transfers
between contexts -- about 43% of effect variance -- and clears the Virtual Cell
Challenge baseline on all three of its headline metrics.  It does not predict
the context-specific remainder, and for a gene perturbed in none of the training
lines it can say which genes will respond but not distinguish that knockdown
from another.  See ``RESULTS.md``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .data import CellLine, load_all
from .model import ContextTransferModel, Hyper

# The switches cross-validation chose on the four-line benchmark, used when a
# caller does not supply their own.  Source lines only; no target line was ever
# involved in selecting these.
DEFAULT_HYPER = Hyper(temp=0.1, smooth=0.15, n_neighbors=5, shrink=0.0,
                      gamma=0.5, mod_clip=5.0, beta=0.8, use_global=0.0,
                      on_target=True, rank=80, rank_mix=0.5, unseen_k=25,
                      renorm=0.0)


@dataclass
class Prediction:
    """Predicted knockdown responses for one cell line."""

    knockdowns: np.ndarray     # (n_kd,) gene symbols that were silenced
    genes: np.ndarray          # (n_genes,) ENSG ids
    symbols: np.ndarray        # (n_genes,) gene symbols
    control: np.ndarray        # (n_genes,) the control profile it started from
    effect: np.ndarray         # (n_kd, n_genes) change from control, log1p CP10K
    seen: np.ndarray           # (n_kd,) was this knockdown in the training lines?

    @property
    def expression(self) -> np.ndarray:
        """Post-perturbation profile, log1p counts-per-10k."""
        return self.control[None, :] + self.effect

    def top_genes(self, knockdown: str, n: int = 20):
        """The genes this knockdown moves most, strongest first."""
        i = int(np.flatnonzero(self.knockdowns == knockdown)[0])
        order = np.argsort(-np.abs(self.effect[i]))[:n]
        return [(str(self.symbols[j]), float(self.effect[i, j])) for j in order]

    def to_frame(self):
        """Long-form table: one row per knockdown x gene.  Needs pandas."""
        import pandas as pd

        n_kd, n_gene = self.effect.shape
        return pd.DataFrame({
            "knockdown": np.repeat(self.knockdowns, n_gene),
            "gene": np.tile(self.symbols, n_kd),
            "ensembl_id": np.tile(self.genes, n_kd),
            "control": np.tile(self.control, n_kd),
            "predicted": self.expression.ravel(),
            "effect": self.effect.ravel(),
            "log2fc": self.effect.ravel() / np.log(2.0),
        })

    def __repr__(self) -> str:
        return (f"Prediction({self.knockdowns.size} knockdowns x "
                f"{self.genes.size} genes, {int(self.seen.sum())} seen in training)")


class VirtualCell:
    """A fitted virtual cell: ask it what a knockdown does in your cell line."""

    def __init__(self, sources: list[CellLine], symbols: np.ndarray,
                 genes: np.ndarray, hyper: Hyper | None = None):
        self._sources = sources
        self._symbols = symbols
        self._genes = genes
        self._hyper = hyper or DEFAULT_HYPER

    @classmethod
    def from_atlas(cls, exclude: str | None = None,
                   hyper: Hyper | None = None) -> "VirtualCell":
        """Train on the packaged four-line CRISPRi atlas.

        ``exclude`` drops one line by name, which is how to get an honest
        zero-shot estimate for a line the atlas already contains.
        """
        lines, genes, symbols = load_all()
        if exclude is not None:
            keep = [cl for cl in lines if cl.name.lower() != exclude.lower()]
            if len(keep) == len(lines):
                raise ValueError(f"no cell line named {exclude!r}; have "
                                 f"{[cl.name for cl in lines]}")
            lines = keep
        return cls(lines, symbols, genes, hyper)

    def control_profile(self, controls) -> np.ndarray:
        """Turn control cells into the mean log1p CP10K profile the model wants.

        Accepts a path to an ``.h5ad``, an ``AnnData``, or an array already on
        this gene axis.  Counts are normalised the same way the training data
        was, so a caller can hand over raw counts without thinking about it.
        """
        if isinstance(controls, np.ndarray):
            if controls.shape[-1] != self._genes.size:
                raise ValueError(f"expected {self._genes.size} genes, "
                                 f"got {controls.shape[-1]}")
            return controls if controls.ndim == 1 else controls.mean(axis=0)

        import anndata as ad

        adata = ad.read_h5ad(controls) if isinstance(controls, (str, Path)) \
            else controls

        # match on symbol first, then on ENSG, so either annotation works
        var = adata.var
        keys = None
        for col in ("gene_name", "gene_symbol", "symbol"):
            if col in var:
                keys = np.asarray(var[col]).astype(str)
                target = self._symbols
                break
        if keys is None:
            names = np.asarray(adata.var_names).astype(str)
            keys, target = names, (self._genes if names[0].startswith("ENSG")
                                   else self._symbols)

        lookup: dict[str, int] = {}
        for i, k in enumerate(keys):
            lookup.setdefault(k, i)
        take = np.array([lookup.get(t, -1) for t in target])
        found = take >= 0
        if found.sum() < 0.5 * target.size:
            raise ValueError(f"only matched {found.sum()} of {target.size} genes; "
                             "check that var carries gene symbols or ENSG ids")

        X = adata.X
        counts = np.asarray(X.sum(axis=0)).ravel()
        profile = np.zeros(self._genes.size)
        totals = counts.sum() or 1.0
        profile[found] = np.log1p(counts[take[found]] / totals * 1e4)
        print(f"matched {found.sum():,}/{target.size:,} genes from "
              f"{adata.n_obs:,} control cells")
        return profile

    def predict(self, controls, knockdowns) -> Prediction:
        """Predict what each knockdown does to this cell line."""
        control = self.control_profile(controls)
        kd = np.asarray(knockdowns, dtype=str)

        model = ContextTransferModel(self._hyper)
        model.fit(self._sources, control, self._symbols)
        effect = model.predict(kd)

        known = {n for s in self._sources for n in s.names}
        seen = np.array([k in known for k in kd])
        return Prediction(knockdowns=kd, genes=self._genes, symbols=self._symbols,
                          control=control, effect=effect, seen=seen)

    def __repr__(self) -> str:
        return (f"VirtualCell(trained on {[s.name for s in self._sources]}, "
                f"{self._genes.size} genes)")
