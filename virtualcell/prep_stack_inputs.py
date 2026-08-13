"""Build Stack's in-context inputs for the leave-one-cell-line-out benchmark.

Stack predicts by analogy: give it a *context* population where the conditions
of interest were measured, and a *query* population of control cells, and it
generates what the query cells would look like under each condition.  Mapped
onto this benchmark:

    base-adata (context)  perturbed cells from the source cell lines, labelled
                          by knocked-down gene, including their non-targeting
                          cells so a control condition is generated too
    test-adata (query)    non-targeting control cells of the held-out line

which is exactly the Challenge's framing -- only the held-out line's controls
cross into prediction, and no perturbation measured in it is ever shown.

Stack aligns genes by symbol against its own 15,012-gene universe, so the two
files only need a shared symbol space, not identical indices.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse

SC = Path("/home/user/vcc_data/singlecell")
CONTROL = "non-targeting"


def load(name: str) -> ad.AnnData:
    a = ad.read_h5ad(SC / f"{name}.h5ad")
    a.var["gene_name"] = a.var["gene_name"].astype(str)
    return a


def build(target: str, sources: list[str], perts: np.ndarray, out_dir: Path,
          n_control: int = 500, seed: int = 0) -> tuple[Path, Path]:
    rng = np.random.default_rng(seed)
    tgt = load(target)
    srcs = [load(s) for s in sources]

    # shared gene symbols across every line involved
    shared = set(tgt.var["gene_name"])
    for s in srcs:
        shared &= set(s.var["gene_name"])
    symbols = np.array(sorted(shared))
    print(f"  shared gene symbols: {symbols.size}")

    def project(a: ad.AnnData, rows: np.ndarray) -> sparse.csr_matrix:
        # first occurrence of each symbol, since symbols can repeat
        lookup: dict[str, int] = {}
        for i, sym in enumerate(a.var["gene_name"].values):
            lookup.setdefault(sym, i)
        cols = np.array([lookup[s] for s in symbols])
        return sparse.csr_matrix(a.X[rows][:, cols])

    # context: source-line cells for the requested knockdowns, plus controls
    wanted = set(perts) | {CONTROL}
    blocks, labels = [], []
    for name, s in zip(sources, srcs):
        keep = np.flatnonzero(np.isin(s.obs["gene"].values, list(wanted)))
        blocks.append(project(s, keep))
        labels.append(s.obs["gene"].values[keep])
        print(f"  context from {name}: {keep.size:,} cells")
    context = ad.AnnData(
        X=sparse.vstack(blocks).tocsr(),
        obs={"gene": np.concatenate(labels)},
        var={"gene_name": symbols},
    )
    context.var_names = symbols

    # query: control cells of the held-out line only
    ctrl_rows = np.flatnonzero(tgt.obs["gene"].values == CONTROL)
    if ctrl_rows.size > n_control:
        ctrl_rows = np.sort(rng.choice(ctrl_rows, n_control, replace=False))
    test = ad.AnnData(
        X=project(tgt, ctrl_rows),
        obs={"gene": np.array([CONTROL] * ctrl_rows.size)},
        var={"gene_name": symbols},
    )
    test.var_names = symbols
    print(f"  query controls from {target}: {ctrl_rows.size:,} cells")

    out_dir.mkdir(parents=True, exist_ok=True)
    ctx_path = out_dir / f"context_{target}.h5ad"
    test_path = out_dir / f"test_{target}.h5ad"
    context.obs_names = [f"c{i}" for i in range(context.n_obs)]
    test.obs_names = [f"q{i}" for i in range(test.n_obs)]
    context.write_h5ad(ctx_path, compression="gzip")
    test.write_h5ad(test_path, compression="gzip")
    print(f"  wrote {ctx_path.name} ({ctx_path.stat().st_size/1e6:.0f} MB) "
          f"and {test_path.name} ({test_path.stat().st_size/1e6:.0f} MB)")
    return ctx_path, test_path


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target", required=True)
    ap.add_argument("--n-perts", type=int, default=100)
    ap.add_argument("--n-control", type=int, default=500)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out-dir", type=Path, default=Path("/home/user/stack_bench"))
    args = ap.parse_args()

    lines = ["k562", "rpe1", "hepg2", "jurkat"]
    sources = [l for l in lines if l != args.target]

    # knockdowns measured in every line, so the comparison is like-for-like
    common = None
    for l in lines:
        names = set(np.unique(load(l).obs["gene"].values)) - {CONTROL}
        common = names if common is None else (common & names)
    common = np.array(sorted(common))
    rng = np.random.default_rng(args.seed)
    perts = np.sort(rng.choice(common, min(args.n_perts, common.size),
                               replace=False))
    print(f"target={args.target} sources={sources}")
    print(f"  shared knockdowns: {common.size}, using {perts.size}")

    build(args.target, sources, perts, args.out_dir, args.n_control, args.seed)
    np.savetxt(args.out_dir / f"perts_{args.target}.txt", perts, fmt="%s")
    print(f"  wrote perts_{args.target}.txt")


if __name__ == "__main__":
    main()
