"""Demo 1 — predict knockdown responses in a cell line, from controls only.

    python demos/demo1_predict.py

This is the everyday use: you have control cells of some cell line and a list of
genes you are thinking of silencing, and you want to know what happens before
spending the reagents.

Here the "new" cell line is K562, held out of training entirely so the numbers
are honest — the model sees K562's control profile and nothing else from it.
Substitute your own ``.h5ad`` of control cells to use it for real.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # repo root

import numpy as np

from virtualcell.predict import VirtualCell

# Knockdowns spanning the range of effect strengths, so the output is not all
# one flavour: two ribosome/spliceosome genes that gut a cell, a proteasome
# subunit, a transcription factor, and a metabolic gene that barely registers.
KNOCKDOWNS = ["RPL5", "SF3B1", "PSMA1", "MYC", "ACOT12"]


def main() -> None:
    # exclude="K562" makes K562 genuinely unseen; drop it to train on all four
    vc = VirtualCell.from_atlas(exclude="K562")
    print(vc)

    # Stand-in for your own control cells: K562's measured control profile.
    from virtualcell.data import load_all
    lines, _, _ = load_all()
    k562_controls = next(cl for cl in lines if cl.name == "K562").mu

    pred = vc.predict(k562_controls, KNOCKDOWNS)
    print(pred, "\n")

    print(f"{'knockdown':<10}{'seen?':<8}{'|effect|':>10}{'genes moved':>14}"
          f"{'on-target':>12}")
    print("-" * 54)
    sym_to_col = {s: i for i, s in enumerate(pred.symbols)}
    for i, kd in enumerate(pred.knockdowns):
        e = pred.effect[i]
        on = e[sym_to_col[kd]] if kd in sym_to_col else np.nan
        print(f"{kd:<10}{'yes' if pred.seen[i] else 'NO':<8}"
              f"{np.linalg.norm(e):>10.2f}{int((np.abs(e) > 0.1).sum()):>14}"
              f"{on:>12.2f}")

    print("\nTop genes moved by SF3B1 knockdown (splicing factor):")
    for gene, effect in pred.top_genes("SF3B1", 12):
        print(f"   {gene:<12}{effect:+.3f}")

    out = Path("demos/out"); out.mkdir(parents=True, exist_ok=True)
    frame = pred.to_frame()
    frame.to_csv(out / "demo1_predictions.csv.gz", index=False,
                 compression="gzip")
    print(f"\nfull table -> {out/'demo1_predictions.csv.gz'} "
          f"({len(frame):,} rows: {pred.knockdowns.size} knockdowns x "
          f"{pred.genes.size} genes)")
    print("columns: knockdown, gene, ensembl_id, control, predicted, effect, log2fc")


if __name__ == "__main__":
    main()
