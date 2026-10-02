"""
bioif.real.qsar -- the C->F and C->G nodes, with label-conditional conformal.

Two models, both on measured public labels:

  C->F   SMILES -> P(SR-p53 reporter active)      Tox21, 6,774 labelled
  C->G   SMILES -> P(Ames positive)               Hansen N=6512

Three things here are not the default thing to do, and each is deliberate:

1. **Scaffold split, not random split.** Bemis-Murcko scaffolds are assigned
   whole to train or test. A random split on a dataset this congeneric
   reports a number you cannot reproduce on new chemistry (§1A axis 10).

2. **Label-conditional (Mondrian) conformal, not marginal.** SR-p53 has 6.2%
   positives. A marginal 90% conformal classifier can hit its guarantee by
   predicting "inactive" for everything, which is useless. Calibrating
   separately per class forces the guarantee to hold for actives too, and
   the price shows up as prediction sets that say {0,1} -- "I don't know" --
   rather than as a silent failure on the minority class.

3. **Prediction SETS, not scores.** The output of a conformalised classifier
   is one of {0}, {1}, {0,1} (abstain) or {} (both labels rejected: an
   out-of-distribution signal). That maps directly onto the contract's right
   to refuse.
"""
from __future__ import annotations

import collections
import math
import statistics
from dataclasses import dataclass, field

import numpy as np

from . import tox


# --------------------------------------------------------------------------
# Features and splits
# --------------------------------------------------------------------------

def featurize(smiles: list[str], n_bits: int = 2048, radius: int = 2):
    """Morgan count fingerprints. Rows for unparseable SMILES are all-zero."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdFingerprintGenerator
    RDLogger.DisableLog("rdApp.*")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
    X = np.zeros((len(smiles), n_bits), dtype=np.float32)
    ok = np.zeros(len(smiles), dtype=bool)
    for i, s in enumerate(smiles):
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        X[i] = np.asarray(gen.GetCountFingerprintAsNumPy(m), dtype=np.float32)
        ok[i] = True
    return X, ok


def scaffold(smiles: str) -> str:
    """Bemis-Murcko scaffold SMILES; '' when it cannot be computed."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem.Scaffolds import MurckoScaffold
    RDLogger.DisableLog("rdApp.*")
    try:
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return ""
        return MurckoScaffold.MurckoScaffoldSmiles(mol=m, includeChirality=False)
    except Exception:
        return ""


def scaffold_split(smiles: list[str], fracs=(0.4, 0.3, 0.3), seed: int = 0):
    """
    Assign whole scaffolds to train / calibrate / test.

    Largest scaffold groups are dealt out first to keep the sizes close to
    the requested fractions; the shuffle only breaks ties among equal-sized
    groups, so the split is stable.
    """
    groups: dict[str, list[int]] = collections.defaultdict(list)
    for i, s in enumerate(smiles):
        groups[scaffold(s) or f"__singleton_{i}"].append(i)
    order = sorted(groups.values(), key=lambda g: (-len(g), smiles[g[0]]))
    rng = np.random.default_rng(seed)
    targets = [f * len(smiles) for f in fracs]
    buckets: list[list[int]] = [[], [], []]
    for g in order:
        deficits = [targets[k] - len(buckets[k]) for k in range(3)]
        k = int(np.argmax(deficits)) if max(deficits) > 0 else \
            int(rng.integers(3))
        buckets[k] += g
    return buckets


# --------------------------------------------------------------------------
# Model + label-conditional conformal
# --------------------------------------------------------------------------

@dataclass
class ConformalClassifier:
    model: object
    q: dict[int, float]                 # class -> nonconformity threshold
    alpha: float
    n_cal: dict[int, int]
    label: str = ""

    def proba(self, X) -> np.ndarray:
        return self.model.predict_proba(X)[:, 1]

    def predict_set(self, X) -> list[frozenset]:
        """{0}, {1}, {0,1} (abstain) or frozenset() (both rejected)."""
        p1 = self.proba(X)
        out = []
        for p in p1:
            s = set()
            if (1.0 - (1.0 - p)) <= self.q[0] + 1e-12:   # s(x,0) = 1 - P(0)
                s.add(0)
            if (1.0 - p) <= self.q[1] + 1e-12:           # s(x,1) = 1 - P(1)
                s.add(1)
            out.append(frozenset(s))
        return out

    def attainable(self, y: int) -> bool:
        return math.isfinite(self.q.get(y, math.inf))


def _conformal_q(scores: list[float], alpha: float) -> float:
    n = len(scores)
    if n == 0:
        return math.inf
    k = math.ceil((n + 1) * (1.0 - alpha))
    return math.inf if k > n else sorted(scores)[k - 1]


def fit_conformal_classifier(Xtr, ytr, Xca, yca, alpha: float = 0.1,
                             label: str = "", C: float = 1.0):
    from sklearn.linear_model import LogisticRegression
    m = LogisticRegression(max_iter=4000, C=C, class_weight="balanced")
    m.fit(Xtr, ytr)
    p = m.predict_proba(Xca)
    q, n_cal = {}, {}
    for cls in (0, 1):
        idx = [i for i, y in enumerate(yca) if y == cls]
        scores = [1.0 - p[i, cls] for i in idx]      # nonconformity
        q[cls] = _conformal_q(scores, alpha)
        n_cal[cls] = len(idx)
    return ConformalClassifier(m, q, alpha, n_cal, label)


@dataclass
class Eval:
    label: str
    n_test: int
    prevalence: float
    auroc: float
    ap: float
    cov: dict[int, float] = field(default_factory=dict)
    set_sizes: dict[int, float] = field(default_factory=dict)


def evaluate(cc: ConformalClassifier, X, y) -> Eval:
    from sklearn.metrics import average_precision_score, roc_auc_score
    y = np.asarray(y)
    p = cc.proba(X)
    sets = cc.predict_set(X)
    cov = {}
    for cls in (0, 1):
        idx = np.where(y == cls)[0]
        cov[cls] = float(np.mean([cls in sets[i] for i in idx])) if len(idx) \
            else float("nan")
    sizes = collections.Counter(len(s) for s in sets)
    return Eval(cc.label, len(y), float(y.mean()),
                float(roc_auc_score(y, p)) if len(set(y)) > 1 else float("nan"),
                float(average_precision_score(y, p)) if len(set(y)) > 1
                else float("nan"),
                cov, {k: v / len(sets) for k, v in sorted(sizes.items())})


# --------------------------------------------------------------------------
# Datasets
# --------------------------------------------------------------------------

def p53_dataset():
    """Tox21 SR-p53: labelled compounds only."""
    rows = [r for r in tox.load_tox21() if r[tox.P53] not in ("", "NA")]
    smi = [r["smiles"] for r in rows]
    y = np.array([int(float(r[tox.P53])) for r in rows])
    return rows, smi, y


def ames_dataset():
    rows = tox.load_ames()
    smi = [r["smiles"] for r in rows]
    y = np.array([int(float(r["ames"])) for r in rows])
    return rows, smi, y


FOLD_FRACS = (0.4, 0.3, 0.3)


def fold_id(dataset: str, fracs=FOLD_FRACS, seed: int = 0) -> str:
    """
    A name for the exact evaluation fold, so two adapters' measured metrics
    can be compared only when they really are comparable.
    """
    return f"tox21/{dataset}/scaffold/{'-'.join(str(f) for f in fracs)}/seed{seed}"


def build(dataset: str, alpha: float = 0.1, seed: int = 0, fracs=FOLD_FRACS):
    """Featurize, scaffold-split, fit, conformalise, evaluate."""
    rows, smi, y = p53_dataset() if dataset == "p53" else ames_dataset()
    X, ok = featurize(smi)
    keep = np.where(ok)[0]
    smi = [smi[i] for i in keep]
    X, y = X[keep], y[keep]
    rows = [rows[i] for i in keep]
    tr, ca, te = scaffold_split(smi, fracs=fracs, seed=seed)
    cc = fit_conformal_classifier(X[tr], y[tr], X[ca], y[ca], alpha,
                                  label=dataset)
    return {"cc": cc, "eval": evaluate(cc, X[te], y[te]),
            "X": X, "y": y, "rows": rows, "smiles": smi,
            "splits": (tr, ca, te), "fold_id": fold_id(dataset, fracs, seed)}
