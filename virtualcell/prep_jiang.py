"""Jiang et al. 2025 as a six-context estimate of per-gene transferability.

The challenge lists this screen among its training data, and the obvious hope --
that it adds source contexts covering the official panel -- does not survive the
coverage check: of its 218 curated signalling regulators, **9 overlap the
official 300, and the genome-wide K562 screen already covers all 9**.  It is a
literature-curated pathway screen (STAT/SMAD/IRF/MAPK cores), so it and the VCC
panel are two different slices of "non-essential regulator" with little in
common.  It also stimulates every sample with a cytokine for 24 h, so even the
overlapping knockdowns are measured in a different cell state.

None of that blocks the use it *is* good for.  The model's per-gene
transferability prior asks a question about the genes being **watched**, not the
genes being **knocked down**: does this gene's response survive a change of
context?  Answering it needs many contexts, not the target panel's
perturbations.  This screen supplies six cell lines from six tissues -- A549
lung, BXPC3 pancreas, HAP1, HT29 colon, K562, MCF7 breast -- four of them
epithelial, which is what the official contexts B and C are and what the
four-line atlas has none of.

What is parsed is the published DE archive, not the 20 GB of Seurat objects:
271 files of ``gene_ID`` plus log2FC, Mixscale beta and p per cell line, giving
a (perturbation x 6 x ~7,000) effect tensor for 324 MB.

    python -m virtualcell.prep_jiang

writes ``jiang_transferability.npz``.
"""

from __future__ import annotations

import argparse
import zipfile
from pathlib import Path

import numpy as np

from .data import DATA

ARCHIVE = DATA / "jiang" / "DE_results_all_pathway.zip"
OUT = DATA / "pseudobulk" / "jiang_transferability.npz"
LINES = ("A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7")


def load_tensor(archive: Path = ARCHIVE, verbose: bool = True
                ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(log2fc, perturbations, genes)`` with log2fc of shape (P, 6, G).

    The axis is the **union** over files, padded with NaN.  Each DE table tests
    only the genes expressed in that pathway's samples, so the sets differ
    substantially between perturbations -- intersecting all 271 leaves a single
    gene, which is why this takes the union and lets the downstream count filter
    decide what is usable.
    """
    z = zipfile.ZipFile(archive)
    files = sorted(n for n in z.namelist() if n.endswith("_DE_results.txt"))

    union: set[str] = set()
    per_file: list[tuple[str, list[str], np.ndarray]] = []
    for i, name in enumerate(files):
        with z.open(name) as fh:
            head = fh.readline().decode().split()
            cols = [head.index(f"log2FC_{c}") - 1 for c in LINES]
            genes, vals = [], []
            for line in fh:
                f = line.decode().split()
                genes.append(f[0])
                vals.append([np.nan if f[c + 1] == "NA" else float(f[c + 1])
                             for c in cols])
        per_file.append((name, genes, np.asarray(vals)))
        union |= set(genes)
        if verbose and (i + 1) % 50 == 0:
            print(f"  read {i + 1}/{len(files)} files "
                  f"({len(union):,} genes seen so far)", flush=True)

    axis = np.array(sorted(union))
    index = {g: i for i, g in enumerate(axis)}
    perts, rows = [], []
    for name, genes, arr in per_file:
        padded = np.full((axis.size, len(LINES)), np.nan)
        take = np.array([index[g] for g in genes])
        padded[take] = arr
        rows.append(padded)
        # 'Parse_IFNG/RUNX1_IFNG_pathway_DE_results.txt' -> 'RUNX1|IFNG'
        base = Path(name).name.replace("_pathway_DE_results.txt", "")
        pathway = Path(name).parent.name.replace("Parse_", "")
        perts.append(f"{base[:-(len(pathway) + 1)]}|{pathway}")

    tensor = np.stack(rows).transpose(0, 2, 1)      # (P, 6, G)
    return tensor, np.array(perts), axis


def transferability(tensor: np.ndarray, min_lines: int = 4
                    ) -> tuple[np.ndarray, np.ndarray]:
    """Per-gene fraction of response variance surviving the across-line mean.

    Same definition the four-line atlas uses, so the two are interchangeable in
    :attr:`Hyper.gene_w`.  Genes measured in fewer than ``min_lines`` cell lines
    for a given perturbation contribute nothing to that perturbation's row.

    Returns ``(fraction, n_usable_perturbations_per_gene)``.  With ``k`` lines a
    pure-noise gene still scores about ``1/k``, so six lines put the floor near
    0.17 against the four-line atlas's 0.33 -- a wider usable range, which is
    the point of using this screen at all.
    """
    ok = np.isfinite(tensor).sum(axis=1) >= min_lines      # (P, G)
    t = np.where(np.isfinite(tensor), tensor, np.nan)

    with np.errstate(invalid="ignore"):
        mean_over_lines = np.nanmean(t, axis=1)            # (P, G)
    mean_over_lines = np.where(ok, mean_over_lines, np.nan)

    shared = np.nanvar(mean_over_lines, axis=0)            # (G,)
    pooled = np.nanvar(np.where(ok[:, None, :], t, np.nan).reshape(-1, t.shape[2]),
                       axis=0)
    frac = shared / np.maximum(pooled, 1e-12)
    return np.clip(frac, 0.0, 1.0), ok.sum(axis=0)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--archive", type=Path, default=ARCHIVE)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()

    print(f"reading {args.archive}")
    tensor, perts, genes = load_tensor(args.archive)
    print(f"  {tensor.shape[0]} perturbations x {tensor.shape[1]} cell lines "
          f"x {tensor.shape[2]:,} genes")
    print(f"  measured fraction: {np.isfinite(tensor).mean():.1%}")

    frac, n = transferability(tensor)
    good = n >= 50
    print(f"\ntransferability over {len(LINES)} lines "
          f"({', '.join(LINES)}), {good.sum():,} genes with >=50 usable "
          f"perturbations:")
    print(f"  median {np.median(frac[good]):.3f}, "
          f"IQR {np.percentile(frac[good], 25):.3f}-"
          f"{np.percentile(frac[good], 75):.3f}, "
          f"range {frac[good].min():.3f}-{frac[good].max():.3f}")
    print(f"  pure-noise floor at {len(LINES)} lines is ~{1/len(LINES):.2f}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.out, genes=genes, fraction=frac, n_perturbations=n,
                        lines=np.array(LINES), perturbations=perts)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
