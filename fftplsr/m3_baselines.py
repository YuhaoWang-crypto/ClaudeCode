"""M3 -- is FFT-PLSR worth it? Held-out comparison against simpler models.

The paper reports FFT-PLSR's prospective accuracy but never what a simpler model
would have scored on the same split, so "the FFT encoding works" is untested
against "almost anything works here". This module runs that test on the paper's
own splits, train on what was known then, score on what was measured after.

Why the comparison is sharp rather than pedantic: in a recombination task only a
dozen positions ever change, so a 227-bin spectrum of a 454-residue protein is a
deterministic function of at most 12 bits. The spectrum cannot encode anything the
bits do not, and the bits are directly regressable.

Two metrics, because they answer different questions:

R2 / Pearson r
    How well does the model describe the whole held-out set? What a paper reports.
top-k mean fitness
    If you had ordered this model's k best untested variants, what would you have
    measured? What a design budget actually buys. Reported next to the oracle
    (the k genuinely best) and chance (the held-out mean).

Run: ``python3 -m fftplsr.m3_baselines``
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import datasets
from .baselines import ALL_BASELINES
from .encode import MutationEncoder
from .model import fit_predict, r2, screen_indices
from .variants import mutated_positions

TOP_K = 8


def _split(train_table, held_table, parent, label):
    return {
        "label": label,
        "parent": parent,
        "train_labels": train_table["Variants"].tolist(),
        "train_y": train_table["Fitness"].to_numpy(dtype=float),
        "held_labels": held_table["Variants"].tolist(),
        "held_y": held_table["Fitness"].to_numpy(dtype=float),
    }


def splits() -> list[dict]:
    """The three prospective splits available in the released data."""
    out = []

    t1, parent = datasets.trainset(1)
    t2, _ = datasets.trainset(2)
    out.append(
        _split(t1, t2[~t2["Variants"].isin(t1["Variants"])], parent, "R1: 13 singles -> 25 combos")
    )

    panel, _ = datasets.measured_panel("IFRS")
    out.append(
        _split(
            t2,
            panel[~panel["Variants"].isin(t2["Variants"])],
            parent,
            "R2: 38 -> 64 later IFRS variants",
        )
    )

    t4, com1 = datasets.trainset(4)
    panel_c, _ = datasets.measured_panel("Com1-IFRS")
    held_c = panel_c[~panel_c["Variants"].isin(t4["Variants"])]
    if len(held_c) > 2:
        out.append(_split(t4, held_c, com1, f"R3: 120 -> {len(held_c)} later Com1 variants"))
    return out


def _top_k_mean(pred, truth, k=TOP_K) -> float:
    order = np.argsort(-np.asarray(pred))
    return float(np.asarray(truth)[order[:k]].mean())


def evaluate_split(split, n_rounds=1, cv=None, verbose=True) -> pd.DataFrame:
    train_labels, train_y = split["train_labels"], split["train_y"]
    held_labels, held_y = split["held_labels"], split["held_y"]
    parent = split["parent"]

    positions = sorted(set(mutated_positions(train_labels)) | set(mutated_positions(held_labels)))
    encoder = MutationEncoder(parent, positions=positions)

    selection = screen_indices(
        lambda codes: encoder.encode(train_labels, codes),
        train_y,
        n_rounds=n_rounds,
        cv=cv,
        verbose=False,
    )
    pred = fit_predict(
        encoder.encode(train_labels, selection.indices),
        train_y,
        encoder.encode(held_labels, selection.indices),
        selection.n_components,
    )

    rows = [
        {
            "model": f"FFT-PLSR [{selection.label}]",
            "heldR2": r2(held_y, pred),
            "pearson": float(np.corrcoef(held_y, pred)[0, 1]),
            f"top{TOP_K}": _top_k_mean(pred, held_y),
        }
    ]
    for cls in ALL_BASELINES:
        try:
            mdl = cls().fit(train_labels, train_y)
            p = np.asarray(mdl.predict(held_labels), dtype=float)
            if not np.isfinite(p).all():
                raise ValueError("non-finite predictions")
            rows.append(
                {
                    "model": cls.name,
                    "heldR2": r2(held_y, p),
                    "pearson": float(np.corrcoef(held_y, p)[0, 1]) if p.std() > 0 else float("nan"),
                    f"top{TOP_K}": _top_k_mean(p, held_y),
                }
            )
        except Exception as exc:  # a baseline that cannot run is reported, not hidden
            rows.append({"model": cls.name, "heldR2": float("nan"), "pearson": float("nan"),
                         f"top{TOP_K}": float("nan"), "note": str(exc)[:40]})

    table = pd.DataFrame(rows).sort_values("heldR2", ascending=False, ignore_index=True)
    table.insert(0, "split", split["label"])
    return table


def main() -> None:
    pd.set_option("display.width", 220)
    all_tables = []
    for split in splits():
        print(f"\n{'=' * 90}\n{split['label']}  "
              f"(train n={len(split['train_labels'])}, held-out n={len(split['held_labels'])})\n{'=' * 90}")
        table = evaluate_split(split)
        held_y = split["held_y"]
        oracle = float(np.sort(held_y)[-TOP_K:].mean())
        chance = float(held_y.mean())
        print(table.drop(columns=["split"]).to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        print(f"  reference: oracle top{TOP_K} = {oracle:.3f}   "
              f"chance (held-out mean) = {chance:.3f}   best measured = {held_y.max():.3f}")
        all_tables.append(table)

    combined = pd.concat(all_tables, ignore_index=True)
    print(f"\n{'=' * 90}\nSUMMARY: held-out R2 by model and split\n{'=' * 90}")
    pivot = combined.assign(
        family=lambda d: np.where(d["model"].str.startswith("FFT-PLSR"), "FFT-PLSR", d["model"])
    ).pivot_table(index="family", columns="split", values="heldR2")
    print(pivot.to_string(float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
