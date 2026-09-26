"""Stage 2 - per-target QSAR benchmark, scaffold split as the headline.

One model is selected per target, by scaffold-split Spearman, from a fingerprint
ridge / random forest / gradient boosting bake-off. Whatever wins, wins - if the
linear model beats the trees on a target, that is what gets used and reported.

Every number that goes in the report is guarded:

* **Scaffold-grouped CV is the headline.** Random-split CV is computed alongside
  purely to quantify how much analogue leakage inflates the same pipeline. In
  the CDK2 ensemble work the leakage factor reached 4.5x, and the arm that
  benefited most from it was the one that failed its permutation null honestly.
  A pipeline that only reports random-split numbers cannot tell those apart.

* **A label-permutation null per target**, through the identical scaffold-split
  pipeline. A target is only called modellable if its real score clears its own
  null, not zero.

* **The assay noise floor sits next to RMSE.** p01 measured the median
  replicate range per target. An RMSE near or below that is fitting assay noise,
  not chemistry, and the report says so per target rather than globally.

* **Repeated over N_SEEDS fold partitions**, so every score carries a spread and
  targets are comparable with their uncertainty visible.

The selected model is refitted on all data and pickled for the screening stage,
together with the training fingerprints - the applicability domain in p04 is
computed against those, so it has to be the exact same feature definition.
"""

from __future__ import annotations

import argparse
import json
import pickle
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, rdFingerprintGenerator
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import GroupKFold, KFold

import lightgbm as lgb

RDLogger.DisableLog("rdApp.*")
warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
MODELS = ROOT / "models"

FP_BITS = 2048
FP_RADIUS = 2
N_FOLDS = 5
N_SEEDS = 3
N_PERMUTATIONS = 10
N_JOBS = 4

DESCRIPTORS = {
    "MolWt": Descriptors.MolWt,
    "MolLogP": Crippen.MolLogP,
    "TPSA": Descriptors.TPSA,
    "NumHDonors": Descriptors.NumHDonors,
    "NumHAcceptors": Descriptors.NumHAcceptors,
    "NumRotatableBonds": Descriptors.NumRotatableBonds,
    "NumAromaticRings": Descriptors.NumAromaticRings,
    "FractionCSP3": Descriptors.FractionCSP3,
    "HeavyAtomCount": Descriptors.HeavyAtomCount,
    "RingCount": Descriptors.RingCount,
    "NHOHCount": Descriptors.NHOHCount,
    "NOCount": Descriptors.NOCount,
}

_gen = rdFingerprintGenerator.GetMorganGenerator(radius=FP_RADIUS, fpSize=FP_BITS)


def featurise(smiles: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """Return (ECFP4 bit matrix as uint8, descriptor matrix as float32)."""
    fps, desc = [], []
    for smi in smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            fps.append(np.zeros(FP_BITS, dtype=np.uint8))
            desc.append(np.zeros(len(DESCRIPTORS), dtype=np.float32))
            continue
        fps.append(_gen.GetFingerprintAsNumPy(mol).astype(np.uint8))
        desc.append(np.array([f(mol) for f in DESCRIPTORS.values()], dtype=np.float32))
    return np.vstack(fps), np.vstack(desc)


def make_models(seed: int) -> dict:
    return {
        "ridge": RidgeCV(alphas=np.logspace(-2, 3, 12)),
        "rf": RandomForestRegressor(
            n_estimators=300,
            min_samples_leaf=2,
            max_features="sqrt",
            n_jobs=N_JOBS,
            random_state=seed,
        ),
        "lgbm": lgb.LGBMRegressor(
            n_estimators=600,
            learning_rate=0.05,
            num_leaves=63,
            min_child_samples=10,
            subsample=0.8,
            subsample_freq=1,
            colsample_bytree=0.3,
            reg_lambda=1.0,
            n_jobs=N_JOBS,
            random_state=seed,
            verbose=-1,
        ),
    }


def cv_predict(X, y, groups, model_name, seed, grouped=True):
    """Out-of-fold predictions under one fold partition."""
    oof = np.full(len(y), np.nan)
    rng = np.random.default_rng(seed)
    if grouped:
        # GroupKFold is deterministic; permuting the group labels per seed gives
        # genuinely different scaffold partitions across repeats.
        uniq = np.unique(groups)
        perm = rng.permutation(len(uniq))
        remap = dict(zip(uniq, perm))
        splits = GroupKFold(n_splits=N_FOLDS).split(
            X, y, groups=np.array([remap[g] for g in groups])
        )
    else:
        splits = KFold(n_splits=N_FOLDS, shuffle=True, random_state=seed).split(X)
    for tr, te in splits:
        m = make_models(seed)[model_name]
        m.fit(X[tr], y[tr])
        oof[te] = m.predict(X[te])
    return oof


def metrics(y, pred) -> dict:
    ok = np.isfinite(pred)
    y, pred = y[ok], pred[ok]
    ss_res = float(((y - pred) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    return {
        "spearman": float(spearmanr(y, pred)[0]),
        "rmse": float(np.sqrt(((y - pred) ** 2).mean())),
        "mae": float(np.abs(y - pred).mean()),
        "r2": float(1 - ss_res / ss_tot) if ss_tot > 0 else np.nan,
    }


def benchmark_target(name: str, prov: dict, args) -> dict | None:
    path = DATA / f"train_{name}.csv"
    if not path.exists():
        return None
    df = pd.read_csv(path)
    y = df["pAffinity"].to_numpy(dtype=float)
    groups = df["scaffold"].to_numpy()
    fp, desc = featurise(df["std_smiles"].tolist())
    X = np.hstack([fp, desc]).astype(np.float32)

    noise = prov["noise_floor"]["noise_floor_log"]
    print(
        f"\n{name}: {len(df)} compounds, {df['scaffold'].nunique()} scaffolds, "
        f"X={X.shape}, noise floor {noise} log"
    )

    per_model = {}
    for mname in ("ridge", "rf", "lgbm"):
        t0 = time.time()
        sc = [metrics(y, cv_predict(X, y, groups, mname, s, True)) for s in range(args.seeds)]
        rnd = metrics(y, cv_predict(X, y, groups, mname, 0, False))
        agg = {
            f"scaffold_{k}_mean": round(float(np.mean([s[k] for s in sc])), 4)
            for k in sc[0]
        } | {
            f"scaffold_{k}_sd": round(float(np.std([s[k] for s in sc])), 4) for k in sc[0]
        } | {f"random_{k}": round(v, 4) for k, v in rnd.items()}
        agg["fit_seconds"] = round(time.time() - t0, 1)
        per_model[mname] = agg
        print(
            f"  {mname:<6s} scaffold rho {agg['scaffold_spearman_mean']:+.3f}"
            f"+-{agg['scaffold_spearman_sd']:.3f}  RMSE {agg['scaffold_rmse_mean']:.3f}  "
            f"| random rho {agg['random_spearman']:+.3f}  ({agg['fit_seconds']:.0f}s)"
        )

    best = max(per_model, key=lambda m: per_model[m]["scaffold_spearman_mean"])

    # --- permutation null for the selected model ---------------------------
    null = []
    rng = np.random.default_rng(0)
    for p in range(args.perms):
        yp = rng.permutation(y)
        null.append(metrics(yp, cv_predict(X, yp, groups, best, p, True))["spearman"])
    null_stats = {
        "mean": round(float(np.mean(null)), 4),
        "sd": round(float(np.std(null)), 4),
        "p95": round(float(np.percentile(null, 95)), 4),
        "max": round(float(np.max(null)), 4),
        "n_permutations": len(null),
    }
    real = per_model[best]["scaffold_spearman_mean"]
    clears = bool(real > null_stats["p95"])
    print(
        f"  best={best}  rho {real:+.3f} vs null p95 {null_stats['p95']:+.3f} -> "
        f"{'CLEARS null' if clears else 'DOES NOT clear null'}"
    )

    # --- per-assay-group signal -------------------------------------------
    oof = cv_predict(X, y, groups, best, 0, True)
    by_group = {}
    for g, sub in df.groupby("assay_group"):
        idx = sub.index.to_numpy()
        if len(idx) >= 30:
            by_group[g] = {
                "n": int(len(idx)),
                "spearman": round(float(spearmanr(y[idx], oof[idx])[0]), 4),
            }

    # --- refit on everything and persist for screening --------------------
    MODELS.mkdir(parents=True, exist_ok=True)
    final = make_models(0)[best]
    final.fit(X, y)

    # A random forest is always persisted alongside, even when it is not the
    # selected model, purely as the uncertainty estimator: per-tree spread is a
    # usable prediction interval and neither ridge nor LightGBM gives one here.
    # The point prediction always comes from the selected model; the two are kept
    # separate in the output so a wide interval is never mistaken for a different
    # prediction.
    unc = make_models(0)["rf"]
    if best == "rf":
        unc = final
    else:
        unc.fit(X, y)

    with open(MODELS / f"{name}.pkl", "wb") as fh:
        pickle.dump(
            {
                "target": name,
                "model_name": best,
                "model": final,
                "uncertainty_rf": unc,
                "train_fp": fp,  # for the applicability domain in p04
                "train_smiles": df["std_smiles"].tolist(),
                "train_pAffinity": y,
                "train_scaffolds": groups,
                "fp_bits": FP_BITS,
                "fp_radius": FP_RADIUS,
                "descriptor_names": list(DESCRIPTORS),
            },
            fh,
        )

    rmse = per_model[best]["scaffold_rmse_mean"]
    return {
        "name": name,
        "n_compounds": int(len(df)),
        "n_scaffolds": int(df["scaffold"].nunique()),
        "best_model": best,
        "models": per_model,
        "permutation_null": null_stats,
        "clears_null": clears,
        "leakage_factor_random_over_scaffold": (
            round(per_model[best]["random_spearman"] / real, 2) if real > 0 else None
        ),
        "noise_floor_log": noise,
        "rmse_vs_noise_floor": (
            round(rmse / noise, 2) if noise else None
        ),
        "per_assay_group": by_group,
        "modellable": bool(clears and real >= 0.3),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=N_SEEDS)
    ap.add_argument("--perms", type=int, default=N_PERMUTATIONS)
    ap.add_argument("--only", default="", help="comma-separated target subset")
    args = ap.parse_args()

    prov = json.loads((RESULTS / "panel_provenance.json").read_text())
    by_name = {t["name"]: t for t in prov["targets"]}
    wanted = (
        [t.strip() for t in args.only.split(",") if t.strip()]
        if args.only
        else [t["name"] for t in prov["targets"] if t["status"] == "ok"]
    )

    out = []
    for name in wanted:
        res = benchmark_target(name, by_name[name], args)
        if res:
            out.append(res)

    print(f"\n{'=' * 84}")
    print(f"{'target':<10s} {'model':<6s} {'scaffold rho':>16s} {'RMSE':>6s} "
          f"{'noise':>6s} {'leak':>5s} {'null p95':>9s}  verdict")
    for r in sorted(out, key=lambda r: -r["models"][r["best_model"]]["scaffold_spearman_mean"]):
        m = r["models"][r["best_model"]]
        print(
            f"{r['name']:<10s} {r['best_model']:<6s} "
            f"{m['scaffold_spearman_mean']:+.3f}+-{m['scaffold_spearman_sd']:.3f}   "
            f"{m['scaffold_rmse_mean']:>6.3f} {str(r['noise_floor_log']):>6s} "
            f"{str(r['leakage_factor_random_over_scaffold']):>5s} "
            f"{r['permutation_null']['p95']:>+9.3f}  "
            f"{'modellable' if r['modellable'] else ('clears null, weak' if r['clears_null'] else 'NOT modellable')}"
        )
    (RESULTS / "benchmark.json").write_text(json.dumps(out, indent=2))
    print(f"\nwrote {RESULTS / 'benchmark.json'}  ({len(out)} targets)")


if __name__ == "__main__":
    main()
