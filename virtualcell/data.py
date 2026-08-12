"""Harmonise CRISPRi Perturb-seq pseudobulk from four cell lines into one tensor.

Sources
-------
K562 / RPE1   Replogle et al. Cell 2022, essential-gene CRISPRi screens.
              Pseudobulk h5ad from Figshare+ 20029387.
HepG2/Jurkat  Nadig et al. Nat Genet 2025, matched essential-gene CRISPRi.
              Pseudobulked here from GEO GSE264667 by ``prep_nadig``.

All four use the same CRISPRi library design, the same targeted essential gene
set, and the same Weissman-lab processing conventions, so a gene knocked down in
one line is directly comparable to the same gene knocked down in another.  That
is what makes them usable as a stand-in for the Virtual Cell Challenge 2 task,
where a model must transfer perturbation effects into a cell line it has never
seen any perturbation data for.

Every cell line is reduced to three arrays on a shared gene axis:

    ctrl  (n_control_replicates, G)   non-targeting pseudobulk, per gem group
    pert  (n_perturbations, G)        one pseudobulk per knocked-down gene
    names (n_perturbations,)          the knocked-down gene symbol

expressed as ``log1p`` of counts-per-10k, which is the space the Virtual Cell
Challenge metrics operate in.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

DATA = Path("/home/user/vcc_data")
CONTROL_LABEL = "non-targeting"


@dataclass
class CellLine:
    """One biological context: control replicates plus its perturbation panel."""

    name: str
    ctrl: np.ndarray          # (n_ctrl, G) log1p CP10K
    pert: np.ndarray          # (n_pert, G) log1p CP10K
    names: np.ndarray         # (n_pert,) knocked-down gene symbol
    ncells: np.ndarray        # (n_pert,) cells behind each pseudobulk
    ctrl_ncells: np.ndarray   # (n_ctrl,)
    genes: np.ndarray         # (G,) ENSG ids
    symbols: np.ndarray       # (G,) gene symbols

    # Both derived arrays are hot paths in the hyperparameter search and each is
    # tens of megabytes, so they are computed once and kept.
    _mu: np.ndarray | None = field(default=None, repr=False, compare=False)
    _delta: np.ndarray | None = field(default=None, repr=False, compare=False)

    @property
    def mu(self) -> np.ndarray:
        """Control mean profile - the only thing a zero-shot model may see."""
        if self._mu is None:
            self._mu = self.ctrl.mean(axis=0)
        return self._mu

    @property
    def delta(self) -> np.ndarray:
        """Perturbation effect, ``pert - control mean``.

        This matches ``perturbation_effect`` in Arc's ``cell-eval``.
        """
        if self._delta is None:
            self._delta = self.pert - self.mu[None, :]
        return self._delta

    def __repr__(self) -> str:
        return (f"CellLine({self.name}: {self.pert.shape[0]} perturbations, "
                f"{self.ctrl.shape[0]} control replicates, {self.pert.shape[1]} genes)")


def _norm_log(counts: np.ndarray, target_sum: float = 1e4) -> np.ndarray:
    """Counts-per-10k then log1p, the space the challenge metrics use."""
    totals = counts.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1.0
    return np.log1p(counts / totals * target_sum)


def _collapse_to_gene(mat: np.ndarray, labels: np.ndarray,
                      weights: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Average rows sharing a label, weighting by cell count.

    Replogle index perturbations by ``gene_transcript`` (``ZC3H18_P1P2``), so a
    gene targeted by two guide sets appears twice.  Cross-line comparison needs
    one row per gene, so those rows are pooled by cell count.
    """
    uniq = np.array(sorted(set(labels)))
    out = np.zeros((uniq.size, mat.shape[1]), dtype=np.float64)
    tot = np.zeros(uniq.size)
    index = {u: i for i, u in enumerate(uniq)}
    for row, lab, w in zip(mat, labels, weights):
        i = index[lab]
        out[i] += row * w
        tot[i] += w
    out /= np.maximum(tot, 1e-9)[:, None]
    return out, uniq


def load_replogle(path: Path, name: str) -> CellLine:
    """Load a Replogle pseudobulk h5ad (raw mean-UMI-per-cell matrix)."""
    import anndata as ad

    a = ad.read_h5ad(path)
    counts = np.asarray(a.X, dtype=np.float64)
    obs_names = np.asarray(a.obs_names)
    is_ctrl = a.obs["core_control"].values.astype(bool)
    # '10023_ZC3H18_P1P2_ENSG00000158545' -> 'ZC3H18'
    target = np.array([s.split("_")[1] for s in obs_names])

    # num_cells_filtered is absent for the non-targeting rows outside the
    # designated core control set; those rows are still usable controls, so fall
    # back to the unfiltered count rather than dropping them.
    ncells = a.obs["num_cells_filtered"].values.astype(float)
    fallback = a.obs["num_cells_unfiltered"].values.astype(float)
    ncells = np.where(np.isfinite(ncells), ncells, fallback)

    # A handful of rows are non-targeting but not flagged core_control; treat any
    # row whose target is the control label as a control replicate.
    is_ctrl = is_ctrl | (target == CONTROL_LABEL)

    pert_counts, pert_names = _collapse_to_gene(
        counts[~is_ctrl], target[~is_ctrl], ncells[~is_ctrl])
    pert_n, _ = _collapse_to_gene(
        ncells[~is_ctrl][:, None], target[~is_ctrl], np.ones((~is_ctrl).sum()))

    return CellLine(
        name=name,
        ctrl=_norm_log(counts[is_ctrl]),
        pert=_norm_log(pert_counts),
        names=pert_names,
        ncells=pert_n.ravel(),
        ctrl_ncells=ncells[is_ctrl],
        genes=np.asarray(a.var_names),
        symbols=np.asarray(a.var["gene_name"].astype(str)),
    )


def load_nadig(path: Path, name: str) -> CellLine:
    """Load a pseudobulk ``.npz`` produced by ``prep_nadig``."""
    z = np.load(path, allow_pickle=True)
    return CellLine(
        name=name,
        ctrl=_norm_log(z["ctrl_mean"].astype(np.float64)),
        pert=_norm_log(z["pert_mean"].astype(np.float64)),
        names=z["pert_names"],
        ncells=z["pert_ncells"],
        ctrl_ncells=z["ctrl_ncells"],
        genes=z["gene_ids"],
        symbols=z["gene_names"],
    )


def load_all(min_lines: int = 4) -> tuple[list[CellLine], np.ndarray, np.ndarray]:
    """Load every available cell line and project onto a shared gene axis.

    Returns the cell lines, the shared ENSG ids, and the shared gene symbols.
    """
    lines: list[CellLine] = []
    rep = DATA / "replogle"
    pb = DATA / "pseudobulk"

    if (rep / "K562_essential_raw_bulk.h5ad").exists():
        lines.append(load_replogle(rep / "K562_essential_raw_bulk.h5ad", "K562"))
    if (rep / "rpe1_raw_bulk.h5ad").exists():
        lines.append(load_replogle(rep / "rpe1_raw_bulk.h5ad", "RPE1"))
    if (pb / "hepg2.npz").exists():
        lines.append(load_nadig(pb / "hepg2.npz", "HepG2"))
    if (pb / "jurkat.npz").exists():
        lines.append(load_nadig(pb / "jurkat.npz", "Jurkat"))

    if len(lines) < min_lines:
        warnings.warn(f"only {len(lines)} cell lines available (wanted {min_lines})")

    shared = set(lines[0].genes)
    for cl in lines[1:]:
        shared &= set(cl.genes)
    shared_genes = np.array(sorted(shared))

    symbols = None
    projected = []
    for cl in lines:
        idx = {g: i for i, g in enumerate(cl.genes)}
        take = np.array([idx[g] for g in shared_genes])
        if symbols is None:
            symbols = cl.symbols[take]
        projected.append(CellLine(
            name=cl.name,
            ctrl=cl.ctrl[:, take],
            pert=cl.pert[:, take],
            names=cl.names,
            ncells=cl.ncells,
            ctrl_ncells=cl.ctrl_ncells,
            genes=shared_genes,
            symbols=symbols,
        ))
    return projected, shared_genes, symbols


def shared_perturbations(lines: list[CellLine]) -> np.ndarray:
    """Gene knockdowns measured in every cell line."""
    common = set(lines[0].names)
    for cl in lines[1:]:
        common &= set(cl.names)
    return np.array(sorted(common))


if __name__ == "__main__":
    lines, genes, symbols = load_all(min_lines=1)
    for cl in lines:
        print(cl)
    print(f"\nshared genes: {genes.size}")
    common = shared_perturbations(lines)
    print(f"shared perturbations across {len(lines)} lines: {common.size}")
