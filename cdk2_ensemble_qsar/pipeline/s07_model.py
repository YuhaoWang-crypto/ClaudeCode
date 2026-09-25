"""Stage 7 - does the structural ensemble add anything over 2D chemistry?

The benchmark is built so that the answer can come back "no". Arms:

  score_only_<slice>  docking score as a ranker, zero parameters. The honest
                      "what does docking alone buy you" number.
  score_only_best     best (most negative) score across slices, still zero
                      parameters.
  2d                  ECFP4 + physicochemical descriptors. The baseline that
                      3D features have to beat to be worth their compute.
  score_1             single-slice docking score, fitted.
  score_ens           all slice scores + min/mean/sd.
  ifp_1               per-residue interaction fingerprint, reference slice.
  ifp_ens             interaction fingerprint pooled over slices (min/mean/sd
                      per channel). The sd channel is the ensemble-specific
                      signal: it encodes how conformation-sensitive a ligand's
                      contact pattern is, which no single structure can say.
  2d_plus_ifp_ens     fusion - the only arm that can show 3D adds *incremental*
                      information rather than merely correlating with 2D.

Protocol:

* Scaffold-grouped 5-fold CV is the headline. Random 5-fold is run alongside
  purely to show how much it inflates the numbers - analogue series split
  across train and test leak, and a random split reports that leakage as skill.
* Repeated over N_SEEDS different fold partitions; every number is mean +- sd
  across seeds, so arms are compared with their spread visible.
* A label-permutation null is run through the identical scaffold-split
  pipeline. An arm is only interesting if it clears its own null, not zero.
* Ranking metrics (EF5%, BEDROC) are reported next to regression metrics
  because a screening cascade cares about the top of the list, not global R2.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Crippen, Descriptors, QED, rdFingerprintGenerator
from scipy.stats import pearsonr, spearmanr
from sklearn.model_selection import GroupKFold, KFold

import lightgbm as lgb

RDLogger.DisableLog("rdApp.*")
warnings.filterwarnings("ignore", category=UserWarning)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
STRUCT = ROOT / "structures"

N_FOLDS = 5
N_SEEDS = 5
N_PERMUTATIONS = 20
ACTIVE_THRESHOLD = 7.0  # pIC50 >= 7 (100 nM) counts as an active for EF/BEDROC
BEDROC_ALPHA = 20.0
FP_BITS = 2048
FP_RADIUS = 2


# --------------------------------------------------------------------------
# metrics
# --------------------------------------------------------------------------
def enrichment_factor(y_true: np.ndarray, y_pred: np.ndarray, frac: float) -> float:
    """EF at the top `frac` of the ranked list."""
    actives = y_true >= ACTIVE_THRESHOLD
    base = actives.mean()
    if base == 0:
        return np.nan
    n_top = max(1, int(round(frac * len(y_true))))
    top = np.argsort(-y_pred)[:n_top]
    return float(actives[top].mean() / base)


def bedroc(y_true: np.ndarray, y_pred: np.ndarray, alpha: float = BEDROC_ALPHA) -> float:
    """Truskett/Bajorath BEDROC - early-recognition-weighted ranking score."""
    actives = (y_true >= ACTIVE_THRESHOLD).astype(int)
    n, n_act = len(actives), int(actives.sum())
    if n_act in (0, n):
        return np.nan
    order = np.argsort(-y_pred)
    ranks = np.where(actives[order] == 1)[0] + 1
    ra = n_act / n
    rie_num = np.exp(-alpha * ranks / n).sum() / n_act
    rie_den = (
        (1 - np.exp(-alpha)) / (n * (np.exp(alpha / n) - 1))
    )
    rie = rie_num / rie_den
    factor = ra * np.sinh(alpha / 2) / (np.cosh(alpha / 2) - np.cosh(alpha / 2 - alpha * ra))
    return float(rie * ra * alpha / (1 - np.exp(-alpha)) / factor) if factor else np.nan


def score_all(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    ok = np.isfinite(y_pred)
    if ok.sum() < 10:
        return {k: np.nan for k in ("pearson_r", "spearman_rho", "rmse", "r2", "ef5", "ef1", "bedroc")}
    yt, yp = y_true[ok], y_pred[ok]
    ss_res = ((yt - yp) ** 2).sum()
    ss_tot = ((yt - yt.mean()) ** 2).sum()
    return {
        "pearson_r": float(pearsonr(yt, yp)[0]),
        "spearman_rho": float(spearmanr(yt, yp)[0]),
        "rmse": float(np.sqrt(((yt - yp) ** 2).mean())),
        "r2": float(1 - ss_res / ss_tot) if ss_tot > 0 else np.nan,
        "ef5": enrichment_factor(yt, yp, 0.05),
        "ef1": enrichment_factor(yt, yp, 0.01),
        "bedroc": bedroc(yt, yp),
    }


# --------------------------------------------------------------------------
# feature blocks
# --------------------------------------------------------------------------
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
    "QED": QED.qed,
}


def build_2d(smiles: list[str]) -> tuple[np.ndarray, list[str]]:
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=FP_RADIUS, fpSize=FP_BITS)
    rows, names = [], None
    for smi in smiles:
        mol = Chem.MolFromSmiles(smi)
        fp = np.array(gen.GetFingerprintAsNumPy(mol), dtype=float)
        desc = [f(mol) for f in DESCRIPTORS.values()]
        rows.append(np.concatenate([fp, desc]))
        if names is None:
            names = [f"ecfp4_{i}" for i in range(FP_BITS)] + list(DESCRIPTORS)
    return np.vstack(rows), names


def pool_ifp(long: pd.DataFrame, ligand_ids: list[str], slices: list[str]):
    """min / mean / sd of every IFP channel across the conformational slices."""
    chan = [
        c for c in long.columns
        if c not in ("ligand_id", "receptor", "score") and not c.endswith("_missing")
    ]
    wide = long.pivot_table(index="ligand_id", columns="receptor", values=chan)
    wide = wide.reindex(ligand_ids)

    blocks, names = [], []
    for c in chan:
        sub = wide[c][slices].to_numpy(dtype=float)
        blocks.extend([np.nanmin(sub, axis=1), np.nanmean(sub, axis=1), np.nanstd(sub, axis=1)])
        names.extend([f"{c}__min", f"{c}__mean", f"{c}__sd"])
    return np.vstack(blocks).T, names


def single_ifp(long: pd.DataFrame, ligand_ids: list[str], slice_id: str):
    sub = long[long["receptor"] == slice_id].set_index("ligand_id")
    chan = [
        c for c in sub.columns
        if c not in ("receptor", "score") and not c.endswith("_missing")
    ]
    return sub.reindex(ligand_ids)[chan].to_numpy(dtype=float), list(chan)


# --------------------------------------------------------------------------
# evaluation
# --------------------------------------------------------------------------
def cv_predict(X, y, groups, seed, grouped=True):
    """Out-of-fold predictions under one fold partition."""
    n = len(y)
    oof = np.full(n, np.nan)
    rng = np.random.default_rng(seed)

    if grouped:
        # GroupKFold is deterministic, so shuffle the group labels per seed to
        # get genuinely different scaffold partitions across repeats.
        uniq = np.unique(groups)
        perm = rng.permutation(len(uniq))
        remap = {g: perm[i] for i, g in enumerate(uniq)}
        shuffled = np.array([remap[g] for g in groups])
        splitter = GroupKFold(n_splits=N_FOLDS).split(X, y, groups=shuffled)
    else:
        splitter = KFold(n_splits=N_FOLDS, shuffle=True, random_state=seed).split(X)

    for tr, te in splitter:
        model = lgb.LGBMRegressor(
            n_estimators=400,
            learning_rate=0.05,
            num_leaves=31,
            min_child_samples=10,
            subsample=0.8,
            subsample_freq=1,
            colsample_bytree=0.6,
            reg_lambda=1.0,
            random_state=seed,
            n_jobs=2,
            verbose=-1,
        )
        Xtr = np.nan_to_num(X[tr], nan=0.0, posinf=0.0, neginf=0.0)
        Xte = np.nan_to_num(X[te], nan=0.0, posinf=0.0, neginf=0.0)
        model.fit(Xtr, y[tr])
        oof[te] = model.predict(Xte)
    return oof


def evaluate(name, X, y, groups, grouped=True):
    per_seed = []
    for seed in range(N_SEEDS):
        oof = cv_predict(X, y, groups, seed, grouped=grouped)
        per_seed.append(score_all(y, oof))
    df = pd.DataFrame(per_seed)
    row = {"arm": name, "n_features": X.shape[1]}
    for m in df.columns:
        row[f"{m}_mean"] = float(df[m].mean())
        row[f"{m}_sd"] = float(df[m].std())
    return row


def main() -> None:
    manifest = json.loads((STRUCT / "ensemble_manifest.json").read_text())
    slices = [s["pdb_id"] for s in manifest["slices"]]
    reference = manifest["reference"]

    lig = pd.read_csv(DATA / "screening_set.csv")
    long = pd.read_csv(RESULTS / "ifp_long.csv.gz")
    scores = pd.read_csv(RESULTS / "docking_scores.csv", index_col=0)

    # Keep only ligands that docked successfully into every slice, so all arms
    # are compared on exactly the same compounds.
    complete = scores.dropna(axis=0, how="any")
    complete = complete[[c for c in slices if c in complete.columns]].dropna()
    lig = lig[lig["ligand_id"].isin(complete.index)].reset_index(drop=True)
    ligand_ids = lig["ligand_id"].tolist()
    print(
        f"{len(ligand_ids)} ligands with a pose in all {len(slices)} slices "
        f"(of {len(scores)} docked)"
    )

    y = lig["pchembl"].to_numpy(dtype=float)
    groups = lig["scaffold"].to_numpy()
    actives = (y >= ACTIVE_THRESHOLD)
    print(
        f"actives (pIC50 >= {ACTIVE_THRESHOLD}): {actives.sum()} "
        f"({actives.mean():.1%});  scaffolds: {len(set(groups))}"
    )

    S = complete.reindex(ligand_ids)[slices].to_numpy(dtype=float)

    # ---- zero-parameter baselines -------------------------------------
    print("\n--- zero-parameter ranking baselines (no model fitted) ---")
    zero_rows = []
    for j, sl in enumerate(slices):
        m = score_all(y, -S[:, j])  # more negative score = better binder
        zero_rows.append({"arm": f"score_only_{sl}", "n_features": 1, **{f"{k}_mean": v for k, v in m.items()}})
        print(
            f"  {sl:>6s}  spearman {m['spearman_rho']:+.3f}  "
            f"EF5% {m['ef5']:.2f}  BEDROC {m['bedroc']:.3f}"
        )
    m = score_all(y, -S.min(axis=1))
    zero_rows.append({"arm": "score_only_best", "n_features": 1, **{f"{k}_mean": v for k, v in m.items()}})
    print(f"  {'best':>6s}  spearman {m['spearman_rho']:+.3f}  EF5% {m['ef5']:.2f}  BEDROC {m['bedroc']:.3f}")

    # ---- feature blocks -------------------------------------------------
    print("\nbuilding feature blocks ...")
    X2d, n2d = build_2d(lig["smiles"].tolist())
    Xifp1, nifp1 = single_ifp(long, ligand_ids, reference)
    XifpN, nifpN = pool_ifp(long, ligand_ids, slices)
    Xs1 = S[:, [slices.index(reference)]]
    XsN = np.hstack([S, S.min(1, keepdims=True), S.mean(1, keepdims=True), S.std(1, keepdims=True)])
    print(
        f"  2d={X2d.shape}  ifp_1={Xifp1.shape}  ifp_ens={XifpN.shape}  "
        f"score_ens={XsN.shape}"
    )

    arms = {
        "2d": X2d,
        "score_1": Xs1,
        "score_ens": XsN,
        "ifp_1": Xifp1,
        "ifp_ens": XifpN,
        "ifp_ens_plus_score": np.hstack([XifpN, XsN]),
        "2d_plus_ifp_ens": np.hstack([X2d, XifpN, XsN]),
    }

    print(f"\n--- scaffold-split {N_FOLDS}-fold CV, {N_SEEDS} partitions ---")
    rows = []
    for name, X in arms.items():
        row = evaluate(name, X, y, groups, grouped=True)
        rows.append(row)
        print(
            f"  {name:<22s} rho {row['spearman_rho_mean']:+.3f}+-{row['spearman_rho_sd']:.3f}  "
            f"R2 {row['r2_mean']:+.3f}+-{row['r2_sd']:.3f}  "
            f"RMSE {row['rmse_mean']:.3f}  EF5% {row['ef5_mean']:.2f}  "
            f"BEDROC {row['bedroc_mean']:.3f}"
        )

    print(f"\n--- random-split {N_FOLDS}-fold CV (leakage reference, NOT the headline) ---")
    rand_rows = []
    for name, X in arms.items():
        row = evaluate(name + "__randomsplit", X, y, groups, grouped=False)
        rand_rows.append(row)
        print(
            f"  {name:<22s} rho {row['spearman_rho_mean']:+.3f}+-{row['spearman_rho_sd']:.3f}  "
            f"R2 {row['r2_mean']:+.3f}"
        )

    # ---- label-permutation null ----------------------------------------
    print(f"\n--- label-permutation null ({N_PERMUTATIONS} permutations, scaffold split) ---")
    null = {}
    rng = np.random.default_rng(0)
    for name in ("2d", "ifp_ens"):
        vals = []
        for p in range(N_PERMUTATIONS):
            yp = rng.permutation(y)
            oof = cv_predict(arms[name], yp, groups, seed=p, grouped=True)
            vals.append(score_all(yp, oof)["spearman_rho"])
        null[name] = {
            "mean": float(np.mean(vals)),
            "sd": float(np.std(vals)),
            "p95": float(np.percentile(vals, 95)),
            "max": float(np.max(vals)),
        }
        print(
            f"  {name:<10s} null rho mean {null[name]['mean']:+.3f} "
            f"sd {null[name]['sd']:.3f}  95th pct {null[name]['p95']:+.3f}"
        )

    out = pd.DataFrame(zero_rows + rows + rand_rows)
    out.to_csv(RESULTS / "model_comparison.csv", index=False)
    (RESULTS / "permutation_null.json").write_text(json.dumps(null, indent=2))
    print(f"\nwrote {RESULTS / 'model_comparison.csv'}")


if __name__ == "__main__":
    main()
