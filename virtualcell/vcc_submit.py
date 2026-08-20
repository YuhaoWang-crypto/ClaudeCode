"""Turn this project's predictions into a Virtual Cell Challenge 2026 submission.

The challenge wants something this project does not natively produce.  Its model
predicts a *mean* effect per knockdown; a submission is **cells** — 400 of them
per knockdown, per context, in raw integer counts:

    3 contexts (A, B, C) x 300 knockdowns x 400 cells = 360,000 cells
    18,533 genes, raw counts, no control cells, one .vcc covering all contexts

The gap matters more than it looks.  ``cell-eval`` computes the differential
expression score by running a rank-sum test **on the cells you submit**, so 400
identical copies of a predicted mean have zero within-group variance and
manufacture significance rather than failing honestly.  A predictor of literally
nothing scores well that way.  So the sampler here is not cosmetic.

**How cells are made.**  Resample real control cells of that context and apply
the predicted shift multiplicatively, with a per-cell knockdown efficiency
drawn around 1.  Starting from real cells keeps the gene-gene correlation
structure and the library-size distribution that a from-scratch sampler would
have to invent, and the per-cell efficiency term supplies the spread that CRISPRi
actually shows.

Nothing here talks to the network.  Getting the official ``gene_names.csv``,
``pert_counts.csv`` and control cells needs an API key, which belongs in the
user's own terminal:

    vcc login --token-stdin        # paste the token there, never into a chat
    vcc datasets download controls -o vcc_2026_controls.zip
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
from scipy import sparse

CONTEXTS = ("A", "B", "C")
CELLS_PER_PERT = 400
PERT_COL = "target_gene"
CONTEXT_COL = "context"


def effect_to_factor(control_mean: np.ndarray, effect: np.ndarray) -> np.ndarray:
    """Convert a log1p-space effect into the multiplicative factor on counts.

    The naive conversion ``exp(effect)`` is wrong, and quietly so.  ``log1p`` is
    only logarithmic where expression is high; near zero it is linear, so the
    same multiplicative change produces a much smaller log1p delta on a lowly
    expressed gene.  Deriving the factor from the effect alone therefore
    under-delivers: asking for -2.00 realised -0.73 on a test panel.

    Anchoring on the control mean fixes it — go to linear space, apply the
    effect there, and take the ratio.
    """
    before = np.expm1(control_mean)
    # Expression cannot go below zero, so an effect larger than the gene's own
    # baseline is unsatisfiable and gets clamped at the floor rather than
    # producing a negative target.
    after = np.expm1(np.maximum(control_mean + effect, 0.0))
    return np.where(before > 1e-9, after / np.maximum(before, 1e-9), 1.0)


def sample_cells(control_counts: sparse.csr_matrix, effect: np.ndarray,
                 n_cells: int, rng: np.random.Generator,
                 control_mean: np.ndarray | None = None,
                 efficiency_sd: float = 0.25) -> sparse.csr_matrix:
    """Realise one knockdown as cells, starting from real control cells.

    ``effect`` is the predicted change in log1p CP10K space.  Applying it as a
    multiplicative factor on counts keeps every cell's library size and
    gene-gene structure intact, and drawing the strength per cell reproduces the
    cell-to-cell variation in knockdown efficiency that makes the submitted
    population disperse the way a real one does.
    """
    n_pool = control_counts.shape[0]
    take = rng.integers(0, n_pool, size=n_cells)
    block = control_counts[take].tocsr()

    if control_mean is None:
        # derive the control mean from the pool itself, in log1p CP10K
        totals = control_counts.sum()
        per_gene = np.asarray(control_counts.sum(axis=0)).ravel()
        control_mean = np.log1p(per_gene / max(totals, 1.0) * 1e4)
    factor = effect_to_factor(control_mean, effect)
    log_factor = np.log(np.clip(factor, 1e-6, 1e6))

    strength = rng.normal(1.0, efficiency_sd, size=n_cells).clip(0.0, 2.5)

    # Work on the CSR buffers directly: fancy-indexing a sparse matrix to get its
    # nonzeros costs a sort per call, and this runs once per knockdown per
    # context.
    values = block.data.astype(np.float64)
    cols = block.indices
    rows = np.repeat(np.arange(n_cells), np.diff(block.indptr))
    scaled = values * np.exp(log_factor[cols] * strength[rows])

    # Counts are integers; round stochastically so the mean is preserved rather
    # than biased downward by truncation.
    floor = np.floor(scaled)
    scaled = floor + (rng.random(scaled.size) < (scaled - floor))

    out = sparse.csr_matrix((scaled.astype(np.int32), cols.copy(),
                             block.indptr.copy()), shape=block.shape)
    out.eliminate_zeros()
    return out


def dedupe_genes(genes: np.ndarray) -> np.ndarray:
    """Column mask keeping the first occurrence of each gene symbol.

    The submission axis is keyed by symbol and ``vcc prep`` rejects a gene list
    with any repeat: *"Each gene must appear exactly once, in order."*  An
    ENSG-indexed matrix can carry two ids that share a symbol -- HSPA14 does
    here -- so the axis has to be collapsed before writing.
    """
    seen: set[str] = set()
    keep = np.zeros(genes.size, dtype=bool)
    for i, g in enumerate(genes):
        if g not in seen:
            seen.add(g)
            keep[i] = True
    return keep


def build(predictions: dict[str, tuple[np.ndarray, np.ndarray]],
          controls: dict[str, sparse.csr_matrix], genes: np.ndarray,
          out_path: Path, n_cells: int = CELLS_PER_PERT,
          seed: int = 0) -> Path:
    """Assemble the submission AnnData.

    ``predictions`` maps context -> (knockdown names, effects) and ``controls``
    maps context -> that context's real control cells on the same gene axis.
    """
    keep = dedupe_genes(genes)
    if not keep.all():
        dropped = genes.size - int(keep.sum())
        print(f"  dropping {dropped} duplicate gene symbol(s) so the axis is "
              f"unique, as prep requires")
        genes = genes[keep]
        cols = np.flatnonzero(keep)
        controls = {c: m[:, cols].tocsr() for c, m in controls.items()}
        predictions = {c: (n, e[:, cols]) for c, (n, e) in predictions.items()}

    rng = np.random.default_rng(seed)
    targets: list[str] = []
    contexts: list[str] = []

    # A full submission is ~1.26 billion stored entries -- 10 GB as CSR, and
    # roughly double that at the moment vstack concatenates.  That does not fit
    # in a normal machine's memory, so blocks stream straight to the HDF5 file
    # and only the row pointer is kept, which is 360k int64s.
    out_path.parent.mkdir(parents=True, exist_ok=True)
    indptr = [0]
    nnz = 0
    per_cell_totals: list[np.ndarray] = []

    with h5py.File(out_path, "w") as h5:
        grp = h5.create_group("X")
        data = grp.create_dataset("data", shape=(0,), maxshape=(None,),
                                  dtype=np.int32, chunks=(1 << 20,),
                                  compression="gzip", compression_opts=1)
        indices = grp.create_dataset("indices", shape=(0,), maxshape=(None,),
                                     dtype=np.int32, chunks=(1 << 20,),
                                     compression="gzip", compression_opts=1)

        for context in CONTEXTS:
            if context not in predictions:
                raise ValueError(f"context {context!r} missing; a submission "
                                 f"must cover all of {', '.join(CONTEXTS)}")
            names, effects = predictions[context]
            pool = controls[context]
            # one control mean per context, reused for every knockdown in it
            mu = np.log1p(np.asarray(pool.sum(axis=0)).ravel()
                          / max(pool.sum(), 1.0) * 1e4)
            print(f"  context {context}: {len(names)} knockdowns x {n_cells} "
                  f"cells from a pool of {pool.shape[0]:,} controls", flush=True)

            for k, (name, effect) in enumerate(zip(names, effects)):
                block = sample_cells(pool, effect, n_cells, rng, control_mean=mu)
                b = block.nnz
                data.resize((nnz + b,))
                indices.resize((nnz + b,))
                data[nnz:nnz + b] = block.data
                indices[nnz:nnz + b] = block.indices
                indptr.extend((block.indptr[1:] + nnz).tolist())
                nnz += b
                targets.extend([name] * n_cells)
                contexts.extend([context] * n_cells)
                per_cell_totals.append(np.asarray(block.sum(1)).ravel())
                if (k + 1) % 100 == 0:
                    print(f"    {k + 1}/{len(names)} ({nnz / 1e6:.0f}M entries)",
                          flush=True)

        n_obs = len(targets)
        grp.create_dataset("indptr", data=np.array(indptr, dtype=np.int64),
                           compression="gzip", compression_opts=1)
        grp.attrs["encoding-type"] = "csr_matrix"
        grp.attrs["encoding-version"] = "0.1.0"
        grp.attrs["shape"] = np.array([n_obs, genes.size], dtype=np.int64)

        _write_frame(h5, "obs", np.array([f"cell_{i}" for i in range(n_obs)]),
                     {PERT_COL: np.array(targets), CONTEXT_COL: np.array(contexts)})
        _write_frame(h5, "var", genes, {"gene_name": genes})
        h5.attrs["encoding-type"] = "anndata"
        h5.attrs["encoding-version"] = "0.1.0"

    totals = np.concatenate(per_cell_totals)
    print(f"  wrote {out_path} — {n_obs:,} cells x {genes.size:,} genes, "
          f"{nnz:,} stored entries ({nnz / n_obs:.0f} genes/cell), "
          f"{out_path.stat().st_size / 1e6:.0f} MB")
    print(f"  per-cell counts: median {np.median(totals):.0f}, "
          f"max {totals.max():.0f} (cap is 1,000,000)")
    return out_path


def _write_frame(h5, name: str, index: np.ndarray,
                 columns: dict[str, np.ndarray]) -> None:
    """Write an AnnData-format dataframe group, categoricals and all."""
    grp = h5.create_group(name)
    grp.attrs["encoding-type"] = "dataframe"
    grp.attrs["encoding-version"] = "0.2.0"
    grp.attrs["_index"] = "_index"
    grp.attrs["column-order"] = np.array(list(columns), dtype=object)

    dt = h5py.special_dtype(vlen=str)
    grp.create_dataset("_index", data=index.astype(object), dtype=dt,
                       compression="gzip", compression_opts=1)
    for key, values in columns.items():
        cats, codes = np.unique(values, return_inverse=True)
        sub = grp.create_group(key)
        sub.attrs["encoding-type"] = "categorical"
        sub.attrs["encoding-version"] = "0.2.0"
        sub.attrs["ordered"] = False
        sub.create_dataset("categories", data=cats.astype(object), dtype=dt)
        sub.create_dataset("codes", data=codes.astype(np.int32),
                           compression="gzip", compression_opts=1)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path,
                    default=Path("/home/user/vcc_submission/prediction.h5ad"))
    ap.add_argument("--n-perts", type=int, default=300)
    ap.add_argument("--n-cells", type=int, default=CELLS_PER_PERT)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    # Stand-in for the official bundle: this project's own four cell lines, three
    # of them cast as contexts A/B/C.  The real run swaps in the downloaded
    # controls and the official 300-knockdown list; the code path is identical.
    from .data import load_all, shared_perturbations
    from .predict import VirtualCell

    lines, gene_ids, symbols = load_all()
    common = shared_perturbations(lines)
    rng = np.random.default_rng(args.seed)
    perts = np.sort(rng.choice(common, min(args.n_perts, common.size),
                               replace=False))
    print(f"knockdowns: {perts.size}")

    import anndata

    sc_dir = Path("/home/user/vcc_data/singlecell")
    stand_ins = {"A": "rpe1", "B": "hepg2", "C": "jurkat"}
    predictions, controls = {}, {}
    for context, tag in stand_ins.items():
        pretty = {"rpe1": "RPE1", "hepg2": "HepG2", "jurkat": "Jurkat"}[tag]
        vc = VirtualCell.from_atlas(exclude=pretty)
        line = next(cl for cl in lines if cl.name == pretty)
        predictions[context] = (perts, vc.predict(line.mu, perts).effect)

        cells = anndata.read_h5ad(sc_dir / f"{tag}.h5ad")
        ctrl = cells[cells.obs["gene"].values == "non-targeting"]
        # project the control cells onto the model's gene axis
        lookup: dict[str, int] = {}
        for i, s in enumerate(cells.var["gene_name"].astype(str).values):
            lookup.setdefault(s, i)
        take = np.array([lookup.get(s, -1) for s in symbols])
        have = take >= 0
        M = sparse.csr_matrix(ctrl.X)[:, take[have]]
        full = sparse.lil_matrix((M.shape[0], symbols.size), dtype=np.int32)
        full[:, np.flatnonzero(have)] = M
        controls[context] = full.tocsr()
        print(f"  context {context} <- {pretty}: {M.shape[0]:,} control cells, "
              f"{have.sum():,}/{symbols.size:,} genes")

    build(predictions, controls, symbols, args.out, args.n_cells, args.seed)
    print("\nvalidate it with:")
    print(f"  vcc prep {args.out} -g gene_names.csv --perts pert_counts.csv "
          f"-o prediction.vcc")


if __name__ == "__main__":
    main()
