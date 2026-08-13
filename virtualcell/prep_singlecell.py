"""Stream a Perturb-seq h5ad and write a cell-subsampled copy.

Stack consumes *cells*, not pseudobulk: its in-context prediction takes perturbed
cells from a source population as context and control cells of the target
population as the query.  The pseudobulk this project runs on therefore cannot
be fed to it, and the raw files are 5.6-10.7 GB each -- more than fits alongside
everything else.

This keeps a bounded number of cells per perturbation, streaming in row chunks so
the full matrix is never resident, and writes a compact sparse h5ad.  At 20 cells
per perturbation a 9.4 GB source becomes a few hundred megabytes, which is small
enough to hold all four cell lines at once.

Handles both schemas in play: the Nadig GEO files store ``X`` dense with obs
categoricals as codes plus a ``__categories`` group; the Replogle files may store
``X`` sparse.  Both label the targeted gene in ``obs/gene`` and use
``non-targeting`` for controls.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
from scipy import sparse

CONTROL_LABEL = "non-targeting"
CHUNK = 20_000


def _decode(arr) -> np.ndarray:
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in arr])


def _column(h5: h5py.File, group: str, name: str) -> np.ndarray:
    """Read an obs/var column, rebuilding a categorical if it is stored as codes."""
    node = h5[f"{group}/{name}"]
    cats = f"{group}/__categories/{name}"
    if cats in h5:
        return _decode(h5[cats][:])[node[:]]
    values = node[:]
    if values.dtype.kind in "OS":
        return _decode(values)
    return values


def _is_sparse(h5: h5py.File) -> bool:
    return isinstance(h5["X"], h5py.Group)


def subsample(path: Path, out: Path, per_pert: int = 20, per_control: int = 2000,
              seed: int = 0) -> None:
    t0 = time.time()
    rng = np.random.default_rng(seed)

    with h5py.File(path, "r") as h5:
        target = _column(h5, "obs", "gene")
        gene_ids = _column(h5, "var", "gene_id") if "var/gene_id" in h5 \
            else _decode(h5[h5["var"].attrs["_index"]][:]) if "_index" in h5["var"].attrs \
            else _column(h5, "var", h5["var"].attrs.get("_index", "gene_id"))
        gene_names = _column(h5, "var", "gene_name")

        # choose the rows to keep before touching X
        keep: list[int] = []
        for name in np.unique(target):
            rows = np.flatnonzero(target == name)
            budget = per_control if name == CONTROL_LABEL else per_pert
            if rows.size > budget:
                rows = rng.choice(rows, budget, replace=False)
            keep.extend(rows.tolist())
        keep = np.sort(np.array(keep))
        print(f"  keeping {keep.size:,} of {target.size:,} cells "
              f"({np.unique(target).size:,} conditions)")

        n_genes = h5["X"].shape[1] if not _is_sparse(h5) else \
            h5["X"].attrs["shape"][1]
        blocks = []
        sparse_input = _is_sparse(h5)
        if sparse_input:
            matrix = ad.read_h5ad(path, backed="r").X          # lazily indexed
            for start in range(0, keep.size, CHUNK):
                sel = keep[start:start + CHUNK]
                blocks.append(sparse.csr_matrix(matrix[sel, :]))
                print(f"    {min(start + CHUNK, keep.size):,}/{keep.size:,} "
                      f"({time.time() - t0:.0f}s)", flush=True)
        else:
            # dense on disk: walk the file in order and pick out the wanted rows
            for start in range(0, target.size, CHUNK):
                stop = min(start + CHUNK, target.size)
                wanted = keep[(keep >= start) & (keep < stop)]
                if wanted.size == 0:
                    continue
                block = h5["X"][start:stop, :]
                blocks.append(sparse.csr_matrix(block[wanted - start]))
                del block
                print(f"    {stop:,}/{target.size:,} ({time.time() - t0:.0f}s)",
                      flush=True)

    X = sparse.vstack(blocks).tocsr()
    adata = ad.AnnData(
        X=X,
        obs={"gene": target[keep]},
        var={"gene_id": gene_ids, "gene_name": gene_names},
    )
    adata.obs_names = [f"cell_{i}" for i in range(adata.n_obs)]
    adata.var_names = gene_ids
    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out, compression="gzip")
    print(f"  wrote {out} ({out.stat().st_size / 1e6:.1f} MB, "
          f"{adata.n_obs:,} cells x {adata.n_vars:,} genes, {time.time() - t0:.0f}s)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("h5ad", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--per-pert", type=int, default=20)
    ap.add_argument("--per-control", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    print(f"subsampling {args.h5ad}")
    subsample(args.h5ad, args.out, args.per_pert, args.per_control, args.seed)


if __name__ == "__main__":
    main()
