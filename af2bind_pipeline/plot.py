"""Optional figures for an AF2BIND run. Needs matplotlib; everything else does not.

    python -m af2bind_pipeline.plot results/6w70 --out figures/af2bind_6w70.png

Draws the p(bind) profile along the sequence, marking true ligand-contact
residues when the input structure contained a ligand, plus the per-bait
activation heat map for the top-scoring residues.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from . import core, structure, weights


def _load_run(run_dir: Path):
    data = np.load(run_dir / "features.npz", allow_pickle=False)
    report = json.loads((run_dir / "report.json").read_text())
    meta = json.loads(str(data["meta"]))
    seeds = tuple(report.get("seeds", [0]))
    heads = weights.load_heads(meta.get("mask_sidechains", True), seeds)
    pred = core.predict(data["features"], heads)
    return data, report, pred


def make_figure(run_dir: str | Path, out_path: str | Path, top_n: int = 15):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    run_dir = Path(run_dir)
    data, report, pred = _load_run(run_dir)
    p_bind = pred["p_bind"]
    chains = [str(c) for c in data["chain"]]
    resis = [int(r) for r in data["resi"]]
    keys = list(zip(chains, resis))

    positives: set = set()
    pdb_path = Path(report["target"])
    if pdb_path.exists():
        residues = structure.parse_pdb(pdb_path)
        lig = structure.ligands(residues)
        target_chains = {c for c, _ in keys}
        if lig and len(target_chains) == 1:
            lig = structure.assign_ligands_to_chain(
                residues, lig, next(iter(target_chains))
            )
        if lig:
            positives = structure.contact_residues(
                structure.protein_residues(residues), lig, 5.0
            ) & set(keys)

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(12, 7), gridspec_kw={"height_ratios": [1, 1.4]}
    )

    x = np.arange(len(p_bind))
    ax1.plot(x, p_bind, lw=1.0, color="#333333", zorder=2)
    ax1.fill_between(x, 0, p_bind, color="#7799cc", alpha=0.35, zorder=1)
    if positives:
        mask = np.array([k in positives for k in keys])
        ax1.scatter(
            x[mask], p_bind[mask], s=26, color="#cc3311", zorder=3,
            label=f"ligand contact within 5 A (n={int(mask.sum())})",
        )
        ax1.legend(loc="upper right", frameon=False, fontsize=9)
    ax1.axhline(0.5, ls="--", lw=0.8, color="grey")
    ax1.set_xlim(0, len(p_bind) - 1)
    ax1.set_ylim(0, 1)
    ax1.set_xlabel("residue (sequential index)")
    ax1.set_ylabel("p(bind)")
    title = f"AF2BIND: {pdb_path.name}, chain {report['meta'].get('chain_requested','?')}"
    if "ligand_validation" in report:
        v = report["ligand_validation"]
        title += (
            f"   ROC-AUC {v['roc_auc']}   AP {v['average_precision']}"
            f"   precision@15 {v['top_n']['15']['precision_at_15']}"
        )
    ax1.set_title(title, fontsize=11)

    order = np.argsort(-p_bind)[:top_n]
    aa_matrix, aa_labels = core.blosum_reorder(pred["p_bind_aa"])
    m = aa_matrix[order].T
    vmax = float(np.abs(m).max()) or 1.0
    im = ax2.imshow(
        m, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax, interpolation="nearest"
    )
    ax2.set_yticks(range(len(aa_labels)))
    ax2.set_yticklabels(aa_labels, fontsize=8)
    ax2.set_xticks(range(len(order)))
    ax2.set_xticklabels(
        [f"{data['resn'][j]}{resis[j]}" for j in order], rotation=90, fontsize=8
    )
    ax2.set_ylabel("bait amino acid")
    ax2.set_title(
        "per-bait logit contribution for the top residues "
        "(red favours binding, each column sums to that residue's total logit)",
        fontsize=10,
    )
    fig.colorbar(im, ax=ax2, pad=0.01, fraction=0.02)

    fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="af2bind-plot")
    ap.add_argument("run_dir", help="an output directory from af2bind_pipeline.run")
    ap.add_argument("--out", default=None)
    ap.add_argument("--top-n", type=int, default=15)
    args = ap.parse_args(argv)
    out = args.out or str(Path(args.run_dir) / "af2bind.png")
    print(f"[af2bind] wrote {make_figure(args.run_dir, out, args.top_n)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
