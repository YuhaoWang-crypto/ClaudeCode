"""M1 -- reproduce the paper's two prospective FFT-PLSR rounds on the PylRS data.

This is the module that decides whether the method has actually been internalised
or merely paraphrased. It re-runs the paper's own rounds and checks three things,
in increasing strictness:

1. **Descriptor selection.** Does screening 566 AAindex entries by out-of-fold MSE
   land on the same winner the authors reported in `output/` (OOBM850103 for round
   1, RADA880104 for round 2), with the same `n_components` and cvMSE?
2. **Prospective accuracy.** Trained only on what was known at the time, what is
   R2 on variants measured *afterwards*? Paper: 0.843 (round 1), 0.835 (round 2).
3. **Did it find the winner?** Does Com1-IFRS (D2N/T56P/R61K/H62Y/H63Y) rank in
   the top handful of a 4096-variant space, i.e. would this pipeline have bought
   the same experiment?

Run: ``python3 -m fftplsr.m1_reproduce``
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.model_selection import KFold

from . import datasets
from .encode import MutationEncoder
from .model import fit_predict, r2, screen_indices
from .variants import enumerate_combinations, mutated_positions

#: Com1-IFRS, the round-2 winner, read off the two released sequences rather than
#: from the text: 7 substitutions, not the 5 a quick read of the paper suggests.
COM1 = "D2N/V31I/T56P/R61K/H62Y/T122S/S193R"

#: What the authors' own `output/` screen tables report, for a direct comparison.
PAPER_OUTPUTS = {
    "round1": {"index": "OOBM850103", "n_components": 10, "cvMSE": 0.35942376545456894},
    "round2": {"index": "RADA880104", "n_components": 4, "cvMSE": 0.32209445034217343},
}
PAPER_PROSPECTIVE_R2 = {"round1": 0.843, "round2": 0.835}


def _fit_round(parent, train_table, n_rounds=1, cv=None, space=None, verbose=True):
    labels = train_table["Variants"].tolist()
    y = train_table["Fitness"].to_numpy(dtype=float)
    query = list(space) if space is not None else []
    positions = sorted(set(mutated_positions(labels)) | set(mutated_positions(query)))
    encoder = MutationEncoder(parent, positions=positions)
    selection = screen_indices(
        lambda codes: encoder.encode(labels, codes), y, n_rounds=n_rounds, cv=cv, verbose=verbose
    )
    return encoder, selection, labels, y


def round1(verbose=True) -> dict:
    """Train on 13 variants (parent + 12 singles); test on the 25 measured later."""
    parent = datasets.ifrs()
    train, _ = datasets.trainset(1)
    later, _ = datasets.trainset(2)

    held = later[~later["Variants"].isin(train["Variants"])].reset_index(drop=True)
    space = enumerate_combinations(parent, datasets.ROUND12_SITES, min_order=2)

    if verbose:
        print(f"\n=== Round 1: train {len(train)} -> hold out {len(held)} ===")
    encoder, selection, labels, y = _fit_round(parent, train, space=space, verbose=verbose)

    X_train = encoder.encode(labels, selection.indices)
    pred_held = fit_predict(
        X_train, y, encoder.encode(held["Variants"], selection.indices), selection.n_components
    )
    held_r2 = r2(held["Fitness"], pred_held)

    scored = pd.DataFrame(
        {
            "variant": space,
            "predicted": fit_predict(
                X_train, y, encoder.encode(space, selection.indices), selection.n_components
            ),
        }
    ).sort_values("predicted", ascending=False, ignore_index=True)

    return {
        "name": "round1",
        "selection": selection,
        "held_out": held.assign(predicted=pred_held),
        "held_r2": held_r2,
        "space": scored,
        "com1_rank": int(scored.index[scored["variant"] == COM1][0]) + 1,
        "n_space": len(space),
    }


def round2(verbose=True) -> dict:
    """Train on 38 variants (singles + measured combos); test on everything measured after."""
    parent = datasets.ifrs()
    train, _ = datasets.trainset(2)
    panel, _ = datasets.measured_panel("IFRS")

    held = panel[~panel["Variants"].isin(train["Variants"])].reset_index(drop=True)
    space = enumerate_combinations(parent, datasets.ROUND12_SITES, min_order=2)

    if verbose:
        print(f"\n=== Round 2: train {len(train)} -> hold out {len(held)} ===")
    encoder, selection, labels, y = _fit_round(parent, train, space=space, verbose=verbose)

    X_train = encoder.encode(labels, selection.indices)
    pred_held = fit_predict(
        X_train, y, encoder.encode(held["Variants"], selection.indices), selection.n_components
    )
    scored = pd.DataFrame(
        {
            "variant": space,
            "predicted": fit_predict(
                X_train, y, encoder.encode(space, selection.indices), selection.n_components
            ),
        }
    ).sort_values("predicted", ascending=False, ignore_index=True)

    return {
        "name": "round2",
        "selection": selection,
        "held_out": held.assign(predicted=pred_held),
        "held_r2": r2(held["Fitness"], pred_held),
        "space": scored,
        "com1_rank": int(scored.index[scored["variant"] == COM1][0]) + 1,
        "n_space": len(space),
    }


def round3_com1_singles(verbose=True) -> dict:
    """Train on 96 Com1-background singles with 3 greedy descriptors (paper's flag=3, cv=10).

    `shuffle=False` reproduces the unshuffled KFold that scikit-learn's
    GridSearchCV used in the original notebook.
    """
    parent = datasets.com1()
    train, _ = datasets.trainset(3)
    if verbose:
        print(f"\n=== Round 3 (site triage): train {len(train)}, 3 descriptors, 10-fold ===")
    encoder, selection, labels, y = _fit_round(
        parent, train, n_rounds=3, cv=KFold(n_splits=10, shuffle=False), verbose=verbose
    )
    return {"name": "round3", "selection": selection}


def _report_row(res) -> dict:
    sel = res["selection"]
    paper = PAPER_OUTPUTS.get(res["name"], {})
    return {
        "round": res["name"],
        "index picked": sel.label,
        "paper index": paper.get("index", "-"),
        "match": "yes" if sel.label == paper.get("index") else "no",
        "n_comp": sel.n_components,
        "cvMSE": round(sel.cv_mse, 4),
        "paper cvMSE": round(paper["cvMSE"], 4) if "cvMSE" in paper else float("nan"),
        "held-out R2": round(res["held_r2"], 3),
        "paper R2": PAPER_PROSPECTIVE_R2.get(res["name"], float("nan")),
        "n held out": len(res["held_out"]),
        "Com1 rank": f"{res['com1_rank']}/{res['n_space']}",
    }


def main() -> None:
    pd.set_option("display.width", 200)
    results = [round1(), round2()]
    table = pd.DataFrame([_report_row(r) for r in results])

    print("\n" + "=" * 100)
    print("REPRODUCTION OF THE PAPER'S PROSPECTIVE ROUNDS")
    print("=" * 100)
    print(table.to_string(index=False))

    for res in results:
        held = res["held_out"].sort_values("Fitness", ascending=False)
        rho = np.corrcoef(held["Fitness"], held["predicted"])[0, 1]
        top = held.head(5)
        print(f"\n--- {res['name']}: held-out set (n={len(held)}), Pearson r = {rho:.3f} ---")
        print(top.to_string(index=False))
        print("  best 5 predicted:")
        print(res["space"].head(5).to_string(index=False))

    sel3 = round3_com1_singles()["selection"]
    print(f"\nRound 3 descriptors: {sel3.label} (paper: QIAN880114_OOBM770105_QIAN880125)")
    print(f"  cvMSE={sel3.cv_mse:.4f} cvR2={sel3.cv_r2:.3f} n_components={sel3.n_components}")


if __name__ == "__main__":
    main()
