"""Reproduce the published PylRS-libY SSM predictions, and assert they are exact.

The upstream Colab notebook (`4.prediction/reference_data/
reproduce_pylrsY_prediction.ipynb`) ships its own stored output. This script
re-runs the same inference and checks the numbers against that output to six
decimal places, so "it ran" and "it ran correctly" are different exit codes.

Run under the Python 3.8 / AutoGluon 0.4.0 environment from ENVIRONMENT.md:

    venv38/bin/python -W ignore pylrs/run_prediction.py --repo /path/to/PylRS-libY

✅ Verified 2026-10-06: all 20 published values reproduce exactly.
"""

from __future__ import annotations

import argparse
import pathlib
import pickle
import sys

#: The axial-conformer rows the notebook printed, as (row index, P(recognise), id, rank).
PUBLISHED_AXIAL = [
    (153, 0.263559, "306A", 1.0),
    (164, 0.240706, "306S", 7.0),
    (169, 0.223808, "306V", 20.0),
]

#: The equatorial-conformer rows the notebook printed.
PUBLISHED_EQUATORIAL = [
    (200, 0.248343, "348K", 4.0),
    (205, 0.229698, "348C", 13.0),
    (210, 0.225479, "346W", 16.0),
    (226, 0.226075, "309V", 15.0),
    (227, 0.247506, "309W", 5.0),
    (233, 0.238118, "346K", 9.0),
    (238, 0.227939, "302Y", 14.0),
    (240, 0.224450, "302V", 19.0),
    (282, 0.225438, "346D", 17.0),
    (321, 0.253156, "305P", 2.0),
    (324, 0.252993, "306A", 3.0),
    (325, 0.235091, "306G", 10.0),
    (329, 0.224905, "306K", 18.0),
    (330, 0.240698, "306I", 8.0),
    (336, 0.229785, "306R", 12.0),
    (340, 0.230716, "306V", 11.0),
    (341, 0.242327, "306T", 6.0),
]

MODEL = "AutogluonModels/ag-20220629_024330/"
TASKS = {
    "TCOY": ("reference_data/B1_TCOY_SSM_evo256.pickle", "reference_data/B1_ssm_tags.pickle"),
    "TCOC": ("reference_data/106_TCOC_SSM_evo256.pickle", "reference_data/106_ssm_tags.pickle"),
}
#: The library is 171 single mutants scored twice, once per TCO ring conformer.
N_PER_CONFORMER = 171


def score(repo: pathlib.Path, task: str):
    """Return a DataFrame of P(recognise) per SSM variant, as the notebook computes it."""
    import pandas as pd
    from autogluon.tabular import TabularPredictor

    feature_file, tag_file = TASKS[task]
    predictor = TabularPredictor.load(str(repo / "4.prediction" / MODEL))
    with open(repo / "4.prediction" / feature_file, "rb") as fh:
        X, _ = pickle.load(fh)
    with open(repo / "4.prediction" / tag_file, "rb") as fh:
        tags = pickle.load(fh)

    proba = predictor.predict_proba(pd.DataFrame(X))
    df = pd.concat([proba, pd.Series(tags, name="id").reset_index(drop=True)], axis=1)
    df["rank"] = df[True].rank(axis=0, method="min", ascending=False)
    df["conformer"] = ["ax"] * N_PER_CONFORMER + ["eq"] * (len(df) - N_PER_CONFORMER)
    return df


def check(df, published, lo, hi, label) -> list[str]:
    """Compare the rows the notebook printed against this run; return failures."""
    failures = []
    window = df[lo:hi]
    shown = window[window["rank"] <= 20]
    if len(shown) != len(published):
        failures.append(
            f"{label}: expected {len(published)} rows with rank<=20, got {len(shown)}"
        )
    for row, proba, ident, rank in published:
        if row not in df.index:
            failures.append(f"{label}: row {row} missing")
            continue
        got_p, got_id, got_rank = df.at[row, True], df.at[row, "id"], df.at[row, "rank"]
        if abs(got_p - proba) > 5e-7:
            failures.append(f"{label} row {row}: P {got_p:.6f} != published {proba:.6f}")
        if str(got_id) != ident:
            failures.append(f"{label} row {row}: id {got_id} != published {ident}")
        if got_rank != rank:
            failures.append(f"{label} row {row}: rank {got_rank} != published {rank}")
    return failures


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo", type=pathlib.Path, required=True,
                    help="checkout of github.com/wendao/PylRS-libY")
    ap.add_argument("--out", type=pathlib.Path, default=pathlib.Path("."),
                    help="where to write the ranked CSVs")
    args = ap.parse_args(argv)

    if not (args.repo / "4.prediction" / MODEL).exists():
        raise SystemExit(f"{args.repo} does not look like a PylRS-libY checkout")

    failures = []
    for task in TASKS:
        df = score(args.repo, task)
        proba = df[True]
        print(f"\n=== {task}: {len(df)} SSM variants "
              f"({N_PER_CONFORMER} axial + {len(df) - N_PER_CONFORMER} equatorial) ===")
        print(f"  P(recognise): min={proba.min():.4f} max={proba.max():.4f} "
              f"mean={proba.mean():.4f}; variants with P>0.5: {(proba > 0.5).sum()}")
        # Rename before itertuples: pandas turns the bool column labels into
        # positional names (_1, _2), which silently reads back the wrong one.
        ranked = (
            df[["id", "conformer", True, "rank"]]
            .rename(columns={True: "P_recognise"})
            .sort_values("P_recognise", ascending=False)
        )
        print("  top 12: " + ", ".join(
            f"{r.id}/{r.conformer}({r.P_recognise:.3f})" for r in ranked.head(12).itertuples()))

        out = args.out / f"{task}_ranked.csv"
        ranked.to_csv(out, index=False)
        print(f"  wrote {out}")

        if task == "TCOY":
            failures += check(df, PUBLISHED_AXIAL, 0, N_PER_CONFORMER, "axial")
            failures += check(df, PUBLISHED_EQUATORIAL, N_PER_CONFORMER, len(df), "equatorial")

    print("\n" + "=" * 70)
    if failures:
        print("REPRODUCTION FAILED:")
        for f in failures:
            print("  -", f)
        return 1
    n = len(PUBLISHED_AXIAL) + len(PUBLISHED_EQUATORIAL)
    print(f"REPRODUCTION EXACT: all {n} published TCOY values match to 6 decimals.")
    print("Note no variant exceeds P=0.5 -- this model is a ranker, not a calibrated")
    print("classifier. Its usable output is the order, not the probability.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
