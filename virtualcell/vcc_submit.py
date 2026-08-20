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
import numpy as np
from scipy import sparse

CONTEXTS = ("A", "B", "C")
CELLS_PER_PERT = 400
PERT_COL = "target_gene"
CONTEXT_COL = "context"


def sample_cells(control_counts: sparse.csr_matrix, effect: np.ndarray,
                 n_cells: int, rng: np.random.Generator,
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
    block = control_counts[take].astype(np.float64).tolil().tocsr()

    # log1p-space effect -> multiplicative factor on counts
    factor = np.expm1(np.abs(effect)) + 1.0
    factor = np.where(effect >= 0, factor, 1.0 / factor)

    strength = rng.normal(1.0, efficiency_sd, size=n_cells).clip(0.0, 2.5)
    out = block.tolil(copy=True)
    data = block.tocsr()
    rows, cols = data.nonzero()
    vals = np.asarray(data[rows, cols]).ravel()
    scaled = vals * np.power(factor[cols], strength[rows])

    # Counts are integers; round stochastically so the mean is preserved rather
    # than biased downward by truncation.
    floor = np.floor(scaled)
    scaled = floor + (rng.random(scaled.size) < (scaled - floor))
    keep = scaled > 0
    out = sparse.csr_matrix(
        (scaled[keep].astype(np.int32), (rows[keep], cols[keep])),
        shape=block.shape, dtype=np.int32)
    return out


def build(predictions: dict[str, tuple[np.ndarray, np.ndarray]],
          controls: dict[str, sparse.csr_matrix], genes: np.ndarray,
          out_path: Path, n_cells: int = CELLS_PER_PERT,
          seed: int = 0) -> Path:
    """Assemble the submission AnnData.

    ``predictions`` maps context -> (knockdown names, effects) and ``controls``
    maps context -> that context's real control cells on the same gene axis.
    """
    rng = np.random.default_rng(seed)
    blocks, targets, contexts = [], [], []

    for context in CONTEXTS:
        if context not in predictions:
            raise ValueError(f"context {context!r} missing; a submission must "
                             f"cover all of {', '.join(CONTEXTS)}")
        names, effects = predictions[context]
        pool = controls[context]
        print(f"  context {context}: {len(names)} knockdowns x {n_cells} cells "
              f"from a pool of {pool.shape[0]:,} controls")
        for name, effect in zip(names, effects):
            blocks.append(sample_cells(pool, effect, n_cells, rng))
            targets.extend([name] * n_cells)
            contexts.extend([context] * n_cells)

    X = sparse.vstack(blocks).tocsr()
    adata = ad.AnnData(
        X=X,
        obs={PERT_COL: np.array(targets), CONTEXT_COL: np.array(contexts)},
        var={"gene_name": genes},
    )
    adata.var_names = genes
    adata.obs_names = [f"cell_{i}" for i in range(adata.n_obs)]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_path, compression="gzip")
    nnz = X.nnz
    print(f"  wrote {out_path} — {adata.n_obs:,} cells x {adata.n_vars:,} genes, "
          f"{nnz:,} stored entries ({nnz / X.shape[0]:.0f} genes/cell), "
          f"{out_path.stat().st_size / 1e6:.0f} MB")
    print(f"  per-cell counts: median {np.median(np.asarray(X.sum(1)).ravel()):.0f} "
          f"(cap is 1,000,000)")
    return out_path


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
