"""Stage 7 - can the panel actually call selectivity, or only potency?

A profiling panel's headline output is a per-compound activity *profile*, and its
most-used derived quantity is the gap between two targets. That gap is a
difference of two predictions, and a difference can be far noisier than either
term - two models each at Spearman 0.8 can still produce a selectivity ranking
worth nothing if their errors are correlated in the wrong way.

Nothing in a single-target benchmark tests this. The panel contains two paralog
pairs that make it directly testable against measured data:

  HDAC1 / HDAC6    same family, same fluorogenic deacetylation readout
  PTP1B / PTPN11   same family, both phosphatases

For compounds with measured activity on both members of a pair, this compares
the **measured** difference against the **predicted** difference, using
scaffold-split out-of-fold predictions on both sides - so neither model has seen
the compound it is being scored on.

The comparison is made three ways, because they can disagree and each answers a
different question:

* Spearman on the signed difference - does predicted selectivity rank measured
  selectivity?
* Sign agreement - does the panel at least get the direction right?
* The same, restricted to compounds whose measured gap exceeds the two targets'
  combined noise floor - because a measured gap smaller than assay noise is not
  a selectivity fact to be reproduced.

If a pair fails here, the selectivity column in the panel output is not
trustworthy for that pair, and the report says so rather than shipping the number.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from p02_benchmark import cv_predict, featurise  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"

PAIRS = [("HDAC1", "HDAC6"), ("PTP1B", "PTPN11")]
N_SEEDS = 3


def oof_for(name: str, model_name: str) -> pd.DataFrame:
    """Scaffold-split out-of-fold predictions, averaged over fold partitions."""
    df = pd.read_csv(DATA / f"train_{name}.csv")
    y = df["pAffinity"].to_numpy(dtype=float)
    groups = df["scaffold"].to_numpy()
    fp, desc = featurise(df["std_smiles"].tolist())
    X = np.hstack([fp, desc]).astype(np.float32)
    preds = np.vstack(
        [cv_predict(X, y, groups, model_name, s, True) for s in range(N_SEEDS)]
    )
    return pd.DataFrame(
        {
            "molecule_chembl_id": df["molecule_chembl_id"],
            "measured": y,
            "oof_pred": preds.mean(axis=0),
            "scaffold": groups,
        }
    )


def main() -> None:
    bench = {b["name"]: b for b in json.loads((RESULTS / "benchmark.json").read_text())}
    prov = {
        t["name"]: t
        for t in json.loads((RESULTS / "panel_provenance.json").read_text())["targets"]
    }

    out = []
    cache: dict[str, pd.DataFrame] = {}
    for a, b in PAIRS:
        if a not in bench or b not in bench:
            print(f"{a}/{b}: one of the pair is not modellable, skipped")
            continue
        for t in (a, b):
            if t not in cache:
                print(f"computing out-of-fold predictions for {t} "
                      f"({bench[t]['best_model']}) ...", flush=True)
                cache[t] = oof_for(t, bench[t]["best_model"])

        m = cache[a].merge(cache[b], on="molecule_chembl_id", suffixes=(f"_{a}", f"_{b}"))
        if len(m) < 30:
            print(f"{a}/{b}: only {len(m)} shared compounds, too few")
            continue

        d_meas = m[f"measured_{a}"] - m[f"measured_{b}"]
        d_pred = m[f"oof_pred_{a}"] - m[f"oof_pred_{b}"]

        # A measured gap smaller than the two assays' combined noise is not a
        # selectivity fact; requiring it to be reproduced would be unfair in both
        # directions, so the strict subset is reported separately.
        na = prov[a]["noise_floor"]["noise_floor_log"] or 0.5
        nb = prov[b]["noise_floor"]["noise_floor_log"] or 0.5
        combined = float(np.hypot(na, nb))
        strong = d_meas.abs() >= combined

        def report(mask, label) -> dict:
            dm, dp = d_meas[mask], d_pred[mask]
            if len(dm) < 20:
                return {"subset": label, "n": int(len(dm)), "note": "too few"}
            nonzero = dm != 0
            return {
                "subset": label,
                "n": int(len(dm)),
                "spearman_delta": round(float(spearmanr(dm, dp)[0]), 3),
                "sign_agreement": round(
                    float((np.sign(dm[nonzero]) == np.sign(dp[nonzero])).mean()), 3
                ),
                "median_abs_error_of_delta": round(float((dp - dm).abs().median()), 3),
                "measured_delta_sd": round(float(dm.std()), 3),
            }

        res = {
            "pair": f"{a}/{b}",
            "n_shared_compounds": int(len(m)),
            "combined_noise_floor_log": round(combined, 3),
            "per_target_scaffold_spearman": {
                a: bench[a]["models"][bench[a]["best_model"]]["scaffold_spearman_mean"],
                b: bench[b]["models"][bench[b]["best_model"]]["scaffold_spearman_mean"],
            },
            "all": report(pd.Series(True, index=m.index), "all shared"),
            "gap_above_noise": report(strong, f"|measured gap| >= {combined:.2f} log"),
        }
        out.append(res)

        print(f"\n{a} vs {b}: {len(m)} shared compounds, combined noise {combined:.2f} log")
        for key in ("all", "gap_above_noise"):
            r = res[key]
            if "spearman_delta" not in r:
                print(f"  {r['subset']:<28s} n={r['n']:<5d} {r.get('note')}")
                continue
            print(
                f"  {r['subset']:<28s} n={r['n']:<5d} "
                f"rho(delta) {r['spearman_delta']:+.3f}  "
                f"sign agreement {r['sign_agreement']:.0%}  "
                f"median |err| {r['median_abs_error_of_delta']:.2f} log "
                f"(measured sd {r['measured_delta_sd']:.2f})"
            )
        strict = res["gap_above_noise"]
        if "spearman_delta" in strict:
            verdict = (
                "selectivity is predictable for this pair"
                if strict["spearman_delta"] >= 0.4 and strict["sign_agreement"] >= 0.7
                else "selectivity is NOT reliably predictable for this pair"
            )
            res["verdict"] = verdict
            print(f"  -> {verdict}")

    (RESULTS / "selectivity_validation.json").write_text(json.dumps(out, indent=2))
    print(f"\nwrote {RESULTS / 'selectivity_validation.json'}")


if __name__ == "__main__":
    main()
