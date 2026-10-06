"""Sequence-only baseline on the MjTyrRS dataset, with the confounds controlled.

Two tasks, because they answer different questions and only one of them avoids
assumed negatives.

Task A -- ncAA attribution (decoy-free)
    Given a selected clone, which ncAA was it selected for? Leave-one-clone-out
    over the 62 measured positives. Uses **only measured data** -- no assumed
    negatives anywhere. If sequence carries ncAA-specific information, this beats
    chance; if it does not, sequence alone can only supply a promiscuity prior.

Task B -- selected vs library (needs decoys)
    Leave-one-ncAA-out: rank the held-out ncAA's positives above random members of
    its own randomised library. ⚠️ The decoys are unselected, not measured
    negatives, so this measures "looks like something selection would keep", not
    "recognises this ncAA". Reported second, and not as the headline.

⚠️ The confound that decides whether Task A means anything
----------------------------------------------------------
Each campaign randomised a *different* set of positions (pBpa 5, HQ-Ala 9, 2-NPA
10). So a model can often name the ncAA from **which positions carry a mutation**
-- a fact about library design, not about substrate recognition. Both tasks are
therefore run twice:

    free        all mutated positions          <- inflated by library design
    controlled  only positions randomised in EVERY campaign (32, 158, 162)

The gap between them is the result, not a nuisance.

Features, all sequence-only, as the reference workflow specifies:

    esm2        mean-pooled ESM-2 embedding of the full 306-residue variant
    esm2-delta  that embedding minus the wild type's (isolates the mutations)
    pssm        log-odds of the clone's residues under the OTHER ncAAs' positives
    onehot      indicator of (position, residue)

Run:  python3 pylrs/baseline.py --variants pylrs/data/tyrrs_variants.csv
"""

from __future__ import annotations

import argparse
import pathlib
import warnings

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")

AA20 = "ACDEFGHIKLMNPQRSTVWY"
ESM_MODEL = "esm2_t12_35M_UR50D"
RNG = 0


# ----------------------------------------------------------------- feature sets

def parse_mutations(text: str) -> dict[int, str]:
    out = {}
    for mut in filter(None, str(text).split("/")):
        out[int(mut[1:-1])] = mut[-1]
    return out


def onehot_features(frame: pd.DataFrame, positions: list[int], wt: str) -> np.ndarray:
    """Indicator of the residue at each allowed position (wild type included)."""
    index = {(p, a): i for i, (p, a) in enumerate((p, a) for p in positions for a in AA20)}
    X = np.zeros((len(frame), len(index)))
    for row, text in enumerate(frame["mutations"]):
        muts = parse_mutations(text)
        for pos in positions:
            X[row, index[(pos, muts.get(pos, wt[pos - 1]))]] = 1.0
    return X


def pssm_scores(train: pd.DataFrame, query: pd.DataFrame,
                positions: list[int], wt: str) -> np.ndarray:
    """Log-odds of each query clone's residues under the training positives.

    The dumb control: "selected clones look alike". If ESM cannot beat this, the
    language model is contributing nothing over residue composition.
    """
    counts = {p: {a: 0.5 for a in AA20} for p in positions}   # Laplace prior
    pos_train = train[train["label"] == 1]
    for text in pos_train["mutations"]:
        muts = parse_mutations(text)
        for pos in positions:
            counts[pos][muts.get(pos, wt[pos - 1])] += 1.0
    logodds = {p: {a: np.log(c / sum(counts[p].values())) for a, c in counts[p].items()}
               for p in positions}

    scores = np.zeros((len(query), len(positions)))
    for row, text in enumerate(query["mutations"]):
        muts = parse_mutations(text)
        for col, pos in enumerate(positions):
            scores[row, col] = logodds[pos][muts.get(pos, wt[pos - 1])]
    return scores


def esm_embeddings(sequences: list[str], model_name: str = ESM_MODEL) -> np.ndarray:
    """Mean-pooled final-layer ESM-2 representation per sequence."""
    import esm
    import torch

    model, alphabet = getattr(esm.pretrained, model_name)()
    model.eval()
    converter = alphabet.get_batch_converter()
    layer = model.num_layers

    out = []
    with torch.no_grad():
        for start in range(0, len(sequences), 8):
            batch = [(f"s{i}", s) for i, s in enumerate(sequences[start:start + 8])]
            _, _, tokens = converter(batch)
            reps = model(tokens, repr_layers=[layer])["representations"][layer]
            for i, (_, seq) in enumerate(batch):
                out.append(reps[i, 1:len(seq) + 1].mean(0).numpy())
    return np.stack(out)


# ----------------------------------------------------------------------- tasks

def task_attribution(frame: pd.DataFrame, feats: dict[str, np.ndarray],
                     label: str, keep: np.ndarray = None) -> list[dict]:
    """Leave-one-clone-out: name the ncAA a selected clone was selected for.

    `keep` restricts to one row per distinct mutation set. Several campaigns
    report the same sequence under different clone names (1BF6 and 1BF9 are
    identical, as are four of the boronate clones); holding one out while its
    twin stays in training makes the task trivial.
    """
    pos = frame["label"].to_numpy() == 1
    if keep is not None:
        pos = pos & keep
    y = frame.loc[pos, "ncAA"].to_numpy()
    counts = pd.Series(y).value_counts(normalize=True)
    chance = float((counts ** 2).sum())     # stratified-random top-1 accuracy

    rows = [{"task": "attribution", "scope": label, "features": "chance (stratified)",
             "top1": chance, "mrr": float("nan"), "n": len(y)}]

    for name, X_all in feats.items():
        X = X_all[pos]
        hits, rr = [], []
        for i in range(len(y)):
            train = np.ones(len(y), bool)
            train[i] = False
            if len(set(y[train])) < 2:
                continue
            model = make_pipeline(StandardScaler(),
                                  LogisticRegression(max_iter=5000, C=1.0,
                                                     class_weight="balanced"))
            model.fit(X[train], y[train])
            proba = model.predict_proba(X[i:i + 1])[0]
            order = np.argsort(-proba)
            ranked = list(model.classes_[order])
            hits.append(ranked[0] == y[i])
            rr.append(1.0 / (ranked.index(y[i]) + 1) if y[i] in ranked else 0.0)
        rows.append({"task": "attribution", "scope": label, "features": name,
                     "top1": float(np.mean(hits)), "mrr": float(np.mean(rr)),
                     "n": len(hits)})
    return rows


def task_selected_vs_library(frame: pd.DataFrame, feats: dict[str, np.ndarray],
                             label: str) -> list[dict]:
    """Leave-one-ncAA-out: rank real positives above their library decoys."""
    y = frame["label"].to_numpy()
    ncaa = frame["ncAA"].to_numpy()
    rows = []
    for name, X in feats.items():
        aucs, aps = [], []
        for held in sorted(set(ncaa[y == 1])):
            test = ncaa == held
            train = ~test
            if y[train].sum() < 3 or not (0 < y[test].sum() < test.sum()):
                continue
            model = make_pipeline(StandardScaler(),
                                  LogisticRegression(max_iter=5000, C=0.1,
                                                     class_weight="balanced"))
            model.fit(X[train], y[train])
            proba = model.predict_proba(X[test])[:, 1]
            aucs.append(roc_auc_score(y[test], proba))
            aps.append(average_precision_score(y[test], proba))
        rows.append({"task": "selected-vs-library", "scope": label, "features": name,
                     "mean_auc": float(np.mean(aucs)), "mean_ap": float(np.mean(aps)),
                     "n_folds": len(aucs)})
    return rows


# ------------------------------------------------------------------------ main

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--variants", type=pathlib.Path,
                    default=pathlib.Path("pylrs/data/tyrrs_variants.csv"))
    ap.add_argument("--campaigns", type=pathlib.Path,
                    default=pathlib.Path("pylrs/data/literature_campaigns.csv"))
    ap.add_argument("--esm-model", default=ESM_MODEL,
                    help="PLM used for all rows (decoys included)")
    ap.add_argument("--esm-model-large", default="esm2_t33_650M_UR50D",
                    help="larger PLM, used for the 62 positives in Task A only")
    args = ap.parse_args(argv)

    import sys
    sys.path.insert(0, str(pathlib.Path(__file__).parent))
    from tyrrs_dataset import apply_mutations, wild_type

    wt = wild_type(args.campaigns)
    frame = pd.read_csv(args.variants).fillna({"mutations": ""})

    # Positions randomised in EVERY campaign -- the confound-free subset.
    per_campaign = [set(int(p) for p in s.split(","))
                    for s in frame["positions"].unique()]
    common = sorted(set.intersection(*per_campaign))
    allpos = sorted(set.union(*per_campaign))
    print(f"MjTyrRS {len(wt)} aa | {int((frame['label'] == 1).sum())} positives, "
          f"{int((frame['label'] == 0).sum())} decoys, {frame['ncAA'].nunique()} ncAAs")
    print(f"positions randomised in every campaign (controlled scope): {common}")
    print(f"union of randomised positions (free scope):               {allpos}\n")

    # One row per distinct mutation set, so identical clones cannot leak.
    first = ~frame.duplicated(subset=["mutations"])
    keep = (frame["label"].to_numpy() == 1) & first.to_numpy()
    print(f"positives after collapsing identical mutation sets: {int(keep.sum())} "
          f"of {int((frame['label'] == 1).sum())}\n")

    seqs = [apply_mutations(wt, m) for m in frame["mutations"]]
    print(f"computing ESM-2 embeddings ({args.esm_model}) for {len(frame)} sequences...")
    emb = esm_embeddings(seqs, args.esm_model)
    wt_emb = esm_embeddings([wt], args.esm_model)[0]
    print(f"  embedding dim {emb.shape[1]}")

    # A larger PLM on the 62 positives only -- Task A needs no decoys, so this
    # stays cheap while testing whether model size changes the conclusion.
    big = np.zeros((len(frame), 1))
    if args.esm_model_large:
        idx = np.flatnonzero(keep)
        print(f"computing {args.esm_model_large} embeddings for {len(idx)} positives...")
        sub = esm_embeddings([seqs[i] for i in idx], args.esm_model_large)
        wt_big = esm_embeddings([wt], args.esm_model_large)[0]
        big = np.zeros((len(frame), sub.shape[1]))
        big[idx] = sub - wt_big
        print(f"  embedding dim {sub.shape[1]}")
    print()

    results = []
    for label, positions in [("free", allpos), ("controlled", common)]:
        feats = {
            "esm2": emb,
            "esm2-delta": emb - wt_emb,
            "pssm": pssm_scores(frame, frame, positions, wt),
            "onehot": onehot_features(frame, positions, wt),
        }
        if label == "controlled":
            # ESM sees the whole sequence, so it cannot be restricted to a position
            # subset; dropping it here keeps the controlled comparison honest.
            feats = {k: v for k, v in feats.items() if not k.startswith("esm2")}

        # The large PLM is embedded for positives only, so decoy rows would be
        # all-zero and trivially separable. It is valid for Task A, which uses no
        # decoys, and must never reach Task B.
        feats_attr = dict(feats)
        if label == "free" and args.esm_model_large and big.shape[1] > 1:
            feats_attr[f"{args.esm_model_large.split('_')[1]}-delta"] = big

        results += task_attribution(frame, feats_attr, label, keep=first.to_numpy())
        results += task_selected_vs_library(frame, feats, label)

    table = pd.DataFrame(results)
    attribution = table[table["task"] == "attribution"].dropna(axis=1, how="all")
    selection = table[table["task"] == "selected-vs-library"].dropna(axis=1, how="all")

    print("=" * 74)
    print("TASK A -- ncAA attribution (decoy-free, measured data only)")
    print("=" * 74)
    print(attribution.drop(columns=["task"]).to_string(
        index=False, float_format=lambda v: f"{v:.3f}"))

    print("\n" + "=" * 74)
    print("TASK B -- selected vs library decoys  (⚠️ negatives assumed, not measured)")
    print("=" * 74)
    print(selection.drop(columns=["task"]).to_string(
        index=False, float_format=lambda v: f"{v:.3f}"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
