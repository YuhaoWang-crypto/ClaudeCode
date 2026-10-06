"""Virtual screen for pAzF-accepting MjTyrRS variants, with its own calibration.

The model is the multiclass attribution classifier from `baseline.py`: trained on
all 76 measured clones to predict which of 13 ncAAs a variant was selected for,
then read as P(pAzF) for candidates it has never seen.

Three parts, in the order that makes the output interpretable:

1. **Calibration first.** Leave one known pAzF clone out, retrain, and score it
   against a large pool of random library members. Its percentile tells you what
   the screen is worth *before* you read the pick-list. A screen whose own known
   positives land mid-pool is not a screen.

2. **Saturation single mutants** at the hotspot positions, as the reference
   workflow specifies. ⚠️ Read these as building blocks, not candidates: every
   published pAzF synthetase carries 4-5 substitutions, and the wild type is not
   a pAzF enzyme, so a single mutant is not expected to work. What this ranking
   says is which single substitutions carry pAzF character.

3. **Combinatorial pick-list**, which is the part you would actually order:
   recombinations of the residues observed across the seven published pAzF clones
   at their five positions, scored and diversity-filtered, with the known clones
   removed.

Run (needs fair-esm + torch):  python3 pylrs/screen_pazf.py
"""

from __future__ import annotations

import argparse
import itertools
import pathlib
import sys
import warnings

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")
sys.path.insert(0, str(pathlib.Path(__file__).parent))

from baseline import esm_embeddings, onehot_features  # noqa: E402
from tyrrs_dataset import apply_mutations, wild_type  # noqa: E402

AA20 = "ACDEFGHIKLMNPQRSTVWY"
TARGET = "pAzF"
#: Positions the pAzF campaign itself randomised.
PAZF_POSITIONS = [32, 107, 158, 159, 162]
#: Hotspots across all MjTyrRS campaigns, for the saturation scan.
HOTSPOTS = [32, 65, 67, 70, 107, 108, 109, 155, 158, 159, 162]
ESM_LARGE = "esm2_t33_650M_UR50D"


def featurise(mutation_sets, wt, model_name, wt_emb=None):
    seqs = [apply_mutations(wt, m) for m in mutation_sets]
    emb = esm_embeddings(seqs, model_name)
    if wt_emb is None:
        wt_emb = esm_embeddings([wt], model_name)[0]
    return emb - wt_emb, wt_emb


def fit_score(X_train, y_train, X_query, target):
    model = make_pipeline(StandardScaler(),
                          LogisticRegression(max_iter=5000, C=1.0,
                                             class_weight="balanced"))
    model.fit(X_train, y_train)
    proba = model.predict_proba(X_query)
    col = list(model.classes_).index(target)
    return proba[:, col]


def random_library(wt, positions, n, rng, exclude):
    """Random members of the campaign's own randomised library."""
    out = []
    while len(out) < n:
        muts = []
        for pos in positions:
            aa = AA20[rng.integers(len(AA20))]
            if aa != wt[pos - 1]:
                muts.append(f"{wt[pos - 1]}{pos}{aa}")
        key = "/".join(muts)
        if muts and key not in exclude:
            exclude.add(key)
            out.append(key)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--variants", type=pathlib.Path,
                    default=pathlib.Path("pylrs/data/tyrrs_variants.csv"))
    ap.add_argument("--campaigns", type=pathlib.Path,
                    default=pathlib.Path("pylrs/data/literature_campaigns.csv"))
    ap.add_argument("--out", type=pathlib.Path, default=pathlib.Path("pylrs/data"))
    ap.add_argument("--pool", type=int, default=2000, help="calibration pool size")
    ap.add_argument("--pick", type=int, default=20)
    ap.add_argument("--quota", type=int, default=2,
                    help="max picks sharing one motif in the spread list")
    ap.add_argument("--motif-positions", default="32,158,162",
                    help="positions whose residues define a motif")
    ap.add_argument("--no-rerank", action="store_true",
                    help="skip the slow ESM-2 650M re-rank of the shortlist")
    args = ap.parse_args(argv)

    rng = np.random.default_rng(0)
    wt = wild_type(args.campaigns)
    frame = pd.read_csv(args.variants).fillna({"mutations": ""})
    pos = frame[(frame["label"] == 1) & ~frame.duplicated(subset=["mutations"])]
    train_muts = pos["mutations"].tolist()
    y = pos["ncAA"].to_numpy()
    print(f"training clones: {len(pos)} over {pos['ncAA'].nunique()} ncAAs "
          f"({(y == TARGET).sum()} are {TARGET})")

    # ---- candidate sets -----------------------------------------------------
    ssm = [f"{wt[p - 1]}{p}{a}" for p in HOTSPOTS for a in AA20 if a != wt[p - 1]]

    observed = {p: set() for p in PAZF_POSITIONS}
    for muts in pos.loc[pos["ncAA"] == TARGET, "mutations"]:
        for mut in muts.split("/"):
            p = int(mut[1:-1])
            if p in observed:
                observed[p].add(mut[-1])
    for p in PAZF_POSITIONS:
        observed[p].add(wt[p - 1])          # allow keeping the wild-type residue
    combos = []
    for choice in itertools.product(*(sorted(observed[p]) for p in PAZF_POSITIONS)):
        muts = [f"{wt[p - 1]}{p}{a}" for p, a in zip(PAZF_POSITIONS, choice)
                if a != wt[p - 1]]
        if len(muts) >= 2:
            combos.append("/".join(muts))
    combos = sorted(set(combos) - set(train_muts))
    print(f"candidates: {len(ssm)} saturation singles, {len(combos)} recombinations "
          f"of residues seen in the {TARGET} clones")

    calib_pool = random_library(wt, PAZF_POSITIONS, args.pool, rng,
                                set(train_muts) | set(combos))

    # ---- features -----------------------------------------------------------
    # One-hot over the hotspot positions, not ESM, for the bulk screen. Measured
    # on this dataset (baseline.py): one-hot reaches 0.621 top-1 attribution
    # against 0.672 for ESM-2 650M -- a real but modest gap -- while 650M costs
    # ~4.6 s/sequence on CPU here, i.e. hours for the ~5k candidates. The
    # shortlist is re-ranked with 650M afterwards, where the cost is affordable.
    everything = train_muts + ssm + combos + calib_pool
    X_all = onehot_features(pd.DataFrame({"mutations": everything}), HOTSPOTS, wt)
    cut1, cut2, cut3 = (len(train_muts), len(train_muts) + len(ssm),
                        len(train_muts) + len(ssm) + len(combos))
    X_train, X_ssm, X_combo, X_pool = (X_all[:cut1], X_all[cut1:cut2],
                                       X_all[cut2:cut3], X_all[cut3:])

    # ---- 1. calibration -----------------------------------------------------
    print("\n" + "=" * 72)
    print(f"1. CALIBRATION -- where do known {TARGET} clones rank among "
          f"{args.pool} random library members?")
    print("=" * 72)
    target_idx = np.flatnonzero(y == TARGET)
    rows = []
    for i in target_idx:
        keep = np.ones(len(y), bool)
        keep[i] = False
        scores = fit_score(X_train[keep], y[keep],
                           np.vstack([X_train[i:i + 1], X_pool]), TARGET)
        rank = int((scores[1:] >= scores[0]).sum()) + 1
        rows.append({"clone": pos["clone"].iloc[i],
                     "mutations": train_muts[i],
                     "rank": rank, "of": args.pool + 1,
                     "percentile": round(100 * rank / (args.pool + 1), 2)})
    calib = pd.DataFrame(rows).sort_values("rank")
    print(calib.to_string(index=False))
    median = calib["percentile"].median()
    print(f"\n  median percentile of a held-out {TARGET} clone: {median:.2f}% "
          f"(chance = 50%)")
    print(f"  enrichment vs chance: {50 / median:.0f}x" if median > 0 else "")
    top10 = (calib["percentile"] <= 10).mean()
    print(f"  {top10:.0%} of held-out clones land in the top 10% of the pool")

    # ---- 2 and 3. scoring ---------------------------------------------------
    ssm_scores = fit_score(X_train, y, X_ssm, TARGET)
    combo_scores = fit_score(X_train, y, X_combo, TARGET)

    ssm_table = (pd.DataFrame({"mutation": ssm, "P_pAzF": ssm_scores})
                 .sort_values("P_pAzF", ascending=False, ignore_index=True))
    combo_table = (pd.DataFrame({"mutations": combos, "P_pAzF": combo_scores})
                   .assign(n_mut=lambda d: d["mutations"].str.count("/") + 1)
                   .sort_values("P_pAzF", ascending=False, ignore_index=True))

    print("\n" + "=" * 72)
    print("2. SATURATION SINGLE MUTANTS at the hotspots (⚠️ building blocks, "
          "not candidates)")
    print("=" * 72)
    print(ssm_table.head(15).to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    print("\n" + "=" * 72)
    print(f"3. PICK-LIST -- top {args.pick} recombinations, known clones excluded")
    print("=" * 72)
    picks = combo_table.head(args.pick)
    print(picks.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    # The calibration above shows the model reproduces the dominant pAzF motif and
    # misses the structurally distinct clones. A top-N list therefore elaborates
    # one family. This second list forces coverage: at most `quota` picks may
    # share the same residue choice at the decisive positions.
    print("\n" + "-" * 72)
    print(f"3b. SPREAD PICK-LIST -- same scores, but <= {args.quota} picks per motif")
    print(f"    (motif = residues at {args.motif_positions})")
    print("-" * 72)
    motif_pos = [int(p) for p in args.motif_positions.split(",")]
    seen: dict[tuple, int] = {}
    chosen = []
    for row in combo_table.itertuples():
        muts = {int(m[1:-1]): m[-1] for m in row.mutations.split("/")}
        motif = tuple(muts.get(p, wt[p - 1]) for p in motif_pos)
        if seen.get(motif, 0) >= args.quota:
            continue
        seen[motif] = seen.get(motif, 0) + 1
        chosen.append(row.Index)
        if len(chosen) == args.pick:
            break
    spread = combo_table.loc[chosen].reset_index(drop=True)
    spread.insert(0, "motif", ["".join(
        {int(m[1:-1]): m[-1] for m in r.split("/")}.get(p, wt[p - 1])
        for p in motif_pos) for r in spread["mutations"]])
    print(spread.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    def motifs_of(frame):
        return {"".join({int(m[1:-1]): m[-1] for m in r.split("/")}.get(p, wt[p - 1])
                        for p in motif_pos) for r in frame["mutations"]}

    print(f"\n  distinct motifs: {len(motifs_of(picks))} in the top-{args.pick} list, "
          f"{len(motifs_of(spread))} in the spread list")

    args.out.mkdir(parents=True, exist_ok=True)
    spread.to_csv(args.out / "pazf_picklist_spread.csv", index=False)

    # ---- 4. re-rank the shortlist with the stronger model -------------------
    if args.no_rerank:
        picks.to_csv(args.out / "pazf_picklist.csv", index=False)
        ssm_table.to_csv(args.out / "pazf_ssm_ranked.csv", index=False)
        combo_table.to_csv(args.out / "pazf_recombination_ranked.csv", index=False)
        calib.to_csv(args.out / "pazf_calibration.csv", index=False)
        print(f"\n-> {args.out} (re-rank skipped)")
        return 0

    short = picks["mutations"].tolist()
    print(f"\nre-ranking the {len(short)} picks with {ESM_LARGE} "
          f"({len(train_muts) + len(short)} sequences)...")
    Xb, _ = featurise(train_muts + short, wt, ESM_LARGE)
    big = fit_score(Xb[:len(train_muts)], y, Xb[len(train_muts):], TARGET)
    picks = picks.assign(
        P_pAzF_esm650=big,
        rank_onehot=np.arange(1, len(picks) + 1),
        rank_esm650=(-pd.Series(big)).rank(method="min").astype(int).to_numpy())
    print(picks.sort_values("rank_esm650").to_string(
        index=False, float_format=lambda v: f"{v:.4f}"))
    agree = float(np.corrcoef(picks["rank_onehot"], picks["rank_esm650"])[0, 1])
    print(f"\n  rank correlation between the two feature sets on the shortlist: "
          f"{agree:.2f}")

    args.out.mkdir(parents=True, exist_ok=True)
    picks.to_csv(args.out / "pazf_picklist.csv", index=False)
    ssm_table.to_csv(args.out / "pazf_ssm_ranked.csv", index=False)
    combo_table.to_csv(args.out / "pazf_recombination_ranked.csv", index=False)
    calib.to_csv(args.out / "pazf_calibration.csv", index=False)
    print(f"\n-> {args.out}/pazf_{{ssm,recombination}}_ranked.csv, pazf_calibration.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
