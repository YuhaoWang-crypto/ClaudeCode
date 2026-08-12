"""Pseudobulk the Nadig et al. 2025 (GSE264667) HepG2 / Jurkat CRISPRi Perturb-seq h5ad files.

The GEO files are dense float32 cells x genes matrices of raw UMI counts (5.6 GB
HepG2, 9.4 GB Jurkat), which is too large to hold in memory alongside anything
else.  We stream them in row chunks and accumulate two sets of sums:

  * one row per targeted gene   -> the perturbation pseudobulk
  * one row per gem group, over non-targeting cells only -> control replicates

The control replicates are what give us an empirical null for differential
expression later on, mirroring the ``core_control`` rows that Replogle et al.
ship in their own pseudobulk files.

Output is a compressed ``.npz`` per cell line holding mean UMI counts per cell,
which matches the convention used by the Replogle pseudobulk files.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import h5py
import numpy as np

CONTROL_LABEL = "non-targeting"
CHUNK = 20_000


def _decode(arr) -> np.ndarray:
    """h5py returns bytes for object dtype; normalise to str."""
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in arr])


def _categorical(h5: h5py.File, group: str, column: str) -> np.ndarray:
    """Rebuild a pandas-style categorical stored as codes + ``__categories``."""
    codes = h5[f"{group}/{column}"][:]
    cats = _decode(h5[f"{group}/__categories/{column}"][:])
    return cats[codes]


def pseudobulk(path: Path, out: Path) -> None:
    t0 = time.time()
    with h5py.File(path, "r") as h5:
        n_cells, n_genes = h5["X"].shape
        target = _categorical(h5, "obs", "gene")
        gem = h5["obs/gem_group"][:]
        gene_ids = _decode(h5["var/gene_id"][:])
        gene_names = _categorical(h5, "var", "gene_name")

        # Perturbation groups: every targeted gene, controls excluded.
        perts = np.array(sorted(set(target) - {CONTROL_LABEL}))
        pert_index = {p: i for i, p in enumerate(perts)}

        # Control groups: one per gem group, non-targeting cells only.
        ctrl_gems = np.array(sorted(set(gem[target == CONTROL_LABEL].tolist())))
        ctrl_index = {g: i for i, g in enumerate(ctrl_gems)}

        pert_sum = np.zeros((perts.size, n_genes), dtype=np.float64)
        pert_n = np.zeros(perts.size, dtype=np.int64)
        ctrl_sum = np.zeros((ctrl_gems.size, n_genes), dtype=np.float64)
        ctrl_n = np.zeros(ctrl_gems.size, dtype=np.int64)

        # Precompute the destination row for every cell so the inner loop is a
        # single np.add.at per chunk rather than a python-level group-by.
        row = np.full(n_cells, -1, dtype=np.int64)
        is_ctrl = target == CONTROL_LABEL
        row[~is_ctrl] = [pert_index[t] for t in target[~is_ctrl]]
        crow = np.full(n_cells, -1, dtype=np.int64)
        crow[is_ctrl] = [ctrl_index[g] for g in gem[is_ctrl]]

        print(f"  {n_cells:,} cells x {n_genes:,} genes | "
              f"{perts.size:,} perturbations | {ctrl_gems.size} control gem groups")

        for start in range(0, n_cells, CHUNK):
            stop = min(start + CHUNK, n_cells)
            block = h5["X"][start:stop, :].astype(np.float64)

            sl_row = row[start:stop]
            m = sl_row >= 0
            if m.any():
                np.add.at(pert_sum, sl_row[m], block[m])
                np.add.at(pert_n, sl_row[m], 1)

            sl_crow = crow[start:stop]
            m = sl_crow >= 0
            if m.any():
                np.add.at(ctrl_sum, sl_crow[m], block[m])
                np.add.at(ctrl_n, sl_crow[m], 1)

            del block
            print(f"    {stop:,}/{n_cells:,} ({time.time() - t0:.0f}s)", flush=True)

    # Mean UMI per cell, matching the Replogle pseudobulk convention.
    pert_mean = pert_sum / np.maximum(pert_n, 1)[:, None]
    ctrl_mean = ctrl_sum / np.maximum(ctrl_n, 1)[:, None]

    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out,
        pert_mean=pert_mean.astype(np.float32),
        pert_names=perts,
        pert_ncells=pert_n,
        ctrl_mean=ctrl_mean.astype(np.float32),
        ctrl_gems=ctrl_gems,
        ctrl_ncells=ctrl_n,
        gene_ids=gene_ids,
        gene_names=gene_names,
    )
    print(f"  wrote {out}  ({out.stat().st_size / 1e6:.1f} MB, {time.time() - t0:.0f}s)")
    print(f"  median cells/perturbation: {np.median(pert_n):.0f} | "
          f"control cells: {ctrl_n.sum():,}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("h5ad", type=Path)
    ap.add_argument("out", type=Path)
    args = ap.parse_args()
    print(f"pseudobulking {args.h5ad}")
    pseudobulk(args.h5ad, args.out)


if __name__ == "__main__":
    main()
