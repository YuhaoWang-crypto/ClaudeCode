"""Export the compact data package.

The raw inputs are 1.1 GB of Figshare pseudobulk plus 15 GB of GEO single-cell
h5ad, which do not belong in a repository.  What does belong is everything
derived from them that is small, hard to recompute, and useful on its own:

``perturbations.csv``   one row per knockdown x cell line -- effect magnitude,
                        cross-line reproducibility, DE burden, on-target
                        knockdown efficiency
``genes.csv``           one row per gene -- baseline expression per line, how
                        much it responds, how context-specific that response is
``model_artifacts.npz`` the fitted per-fold quantities a prediction needs that
                        are not simply re-derivable: context weights, per-gene
                        knockdown efficiency, the response-program basis

Anyone can rebuild the full inputs with ``fetch_data.sh``; these tables let them
check the headline claims without downloading 16 GB first.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from . import metrics as M
from .benchmark import de_truth
from .data import load_all, shared_perturbations
from .model import ContextTransferModel, Hyper, SourceBank

OUT = Path("virtualcell/datapackage")


def perturbation_table(lines, common) -> str:
    """Per-knockdown, per-cell-line summary statistics."""
    col = {s: i for i, s in enumerate(lines[0].symbols)}
    rows = ["cell_line,knockdown,n_cells,effect_l2,relative_effect,"
            "n_de_genes,on_target_residual,mean_cross_line_r"]

    # cross-line correlation of each knockdown's effect, averaged over pairs
    eff = {}
    for cl in lines:
        idx = {n: i for i, n in enumerate(cl.names)}
        eff[cl.name] = cl.delta[np.array([idx[p] for p in common])]
    pair_r = np.zeros(common.size)
    pair_n = 0
    for i in range(len(lines)):
        for j in range(i + 1, len(lines)):
            a, b = eff[lines[i].name], eff[lines[j].name]
            a = a - a.mean(axis=1, keepdims=True)
            b = b - b.mean(axis=1, keepdims=True)
            num = (a * b).sum(axis=1)
            den = np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1) + 1e-9
            pair_r += num / den
            pair_n += 1
    pair_r /= max(pair_n, 1)

    for cl in lines:
        idx = {n: i for i, n in enumerate(cl.names)}
        take = np.array([idx[p] for p in common])
        truth = de_truth(cl, common)
        n_de = truth.sig.sum(axis=1)
        mag = np.linalg.norm(cl.delta[take], axis=1)
        rel = mag / np.median(mag)
        for k, p in enumerate(common):
            g = col.get(p)
            residual = ""
            if g is not None and cl.mu[g] >= 0.1:
                before = np.expm1(cl.mu[g])
                if before > 0:
                    residual = f"{min(np.expm1(cl.pert[take[k], g]) / before, 1.0):.4f}"
            rows.append(f"{cl.name},{p},{cl.ncells[take[k]]:.0f},{mag[k]:.4f},"
                        f"{rel[k]:.4f},{n_de[k]},{residual},{pair_r[k]:.4f}")
    return "\n".join(rows) + "\n"


def gene_table(lines, genes, symbols, common) -> str:
    """Per-gene baseline expression and responsiveness, per cell line."""
    header = ["ensembl_id", "symbol"]
    header += [f"baseline_{cl.name}" for cl in lines]
    header += [f"response_sd_{cl.name}" for cl in lines]
    header += ["shared_response_fraction"]
    rows = [",".join(header)]

    eff = []
    for cl in lines:
        idx = {n: i for i, n in enumerate(cl.names)}
        eff.append(cl.delta[np.array([idx[p] for p in common])])
    stack = np.stack(eff)                       # (n_lines, n_pert, n_genes)
    # how much of a gene's response variance survives averaging across lines
    shared = stack.mean(axis=0).var(axis=0)
    total = stack.reshape(-1, stack.shape[-1]).var(axis=0)
    frac = shared / np.maximum(total, 1e-12)

    for gi in range(genes.size):
        vals = [genes[gi], symbols[gi]]
        vals += [f"{cl.mu[gi]:.4f}" for cl in lines]
        vals += [f"{e[:, gi].std():.4f}" for e in eff]
        vals.append(f"{frac[gi]:.4f}")
        rows.append(",".join(vals))
    return "\n".join(rows) + "\n"


def model_artifacts(lines, symbols, hyper: dict[str, Hyper]) -> dict:
    """Per-fold fitted quantities that a prediction actually needs."""
    art: dict[str, np.ndarray] = {}
    for i, target in enumerate(lines):
        sources = [s for j, s in enumerate(lines) if j != i]
        hp = hyper.get(target.name, Hyper())
        model = ContextTransferModel(hp)
        weights = model.context_weights(sources, target.mu)
        bank = SourceBank.build(sources, weights)
        tag = target.name
        art[f"{tag}/source_names"] = np.array([s.name for s in sources])
        art[f"{tag}/context_weights"] = weights.astype(np.float32)
        art[f"{tag}/perturbations"] = bank.names
        art[f"{tag}/reliability"] = bank.reliability.astype(np.float32)
        art[f"{tag}/knockdown_residual"] = bank.kd.astype(np.float32)
        art[f"{tag}/global_effect"] = bank.global_effect.astype(np.float32)
        # the leading response programs, which are the compact, reusable part
        art[f"{tag}/program_basis"] = bank.basis[:50].astype(np.float32)
    return art


def main() -> None:
    lines, genes, symbols = load_all()
    common = shared_perturbations(lines)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "tables").mkdir(exist_ok=True)

    ctx = Path("results/context.json")
    hyper: dict[str, Hyper] = {}
    if ctx.exists():
        raw = json.loads(ctx.read_text())["hyper"]
        for name, h in raw.items():
            clean = {}
            for k, v in h.items():
                default = getattr(Hyper(), k)
                if isinstance(default, bool):
                    clean[k] = bool(v)
                elif isinstance(default, int):
                    clean[k] = int(float(v))
                else:
                    clean[k] = float(v)
            hyper[name] = Hyper(**clean)
        (OUT / "tables" / "hyperparameters.json").write_text(
            json.dumps(raw, indent=2, default=str))

    print("writing perturbations.csv ...")
    (OUT / "tables" / "perturbations.csv").write_text(
        perturbation_table(lines, common))
    print("writing genes.csv ...")
    (OUT / "tables" / "genes.csv").write_text(
        gene_table(lines, genes, symbols, common))
    print("writing model_artifacts.npz ...")
    np.savez_compressed(OUT / "model_artifacts.npz",
                        **model_artifacts(lines, symbols, hyper),
                        genes=genes, symbols=symbols,
                        shared_perturbations=common)

    for f in sorted(OUT.rglob("*")):
        if f.is_file():
            print(f"  {f.stat().st_size / 1e6:8.2f} MB  {f}")


if __name__ == "__main__":
    main()
