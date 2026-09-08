"""The AF2BIND head: AlphaFold2 pair features -> per-residue p(bind).

Pure numpy. Everything here is a transcription of the reference implementation's
`af2bind()` function (github.com/sokrypton/af2bind, MIT), restructured so the
feature-extraction step and the logistic head are separable and cacheable.

Shapes, for a target of length L with the standard 20 bait residues:

    pair representation   (L + 20, L + 20, 128)
    pair_A  target->bait  (L, 20, 128)
    pair_B  bait->target  (L, 20, 128)     [transposed back to target-major]
    features x            (L, 2 * 20 * 128) == (L, 5120)
    logit_aa              (L, 20)          per-bait logit contribution
    p_bind                (L,)             sigmoid(logit_aa.sum(-1))
"""

from __future__ import annotations

import numpy as np

#: One residue of each amino acid type, in the order the head expects. These are
#: the 20 "bait" residues appended to the target as separate single-residue
#: chains; their pair-representation coupling to the target is the signal.
BAIT_SEQ = "ACDEFGHIKLMNPQRSTVWY"
N_BAIT = len(BAIT_SEQ)
PAIR_DIM = 128
N_FEATURES = 2 * N_BAIT * PAIR_DIM  # 5120

#: Order used by the paper's activation heat map: BLOSUM-style grouping, which
#: puts chemically similar baits next to each other.
BLOSUM_ORDER = "CSTAGPDEQNHRKMILVWYF"


def features_from_pair(pair: np.ndarray, n_bait: int = N_BAIT) -> np.ndarray:
    """Extract the (L, 5120) AF2BIND feature matrix from an AF2 pair tensor.

    `pair` is ``outputs["representations"]["pair"]`` of a prediction on
    target + `n_bait` appended bait residues, i.e. shape (L + n_bait,
    L + n_bait, 128). The last `n_bait` rows/columns are the baits.
    """
    pair = np.asarray(pair)
    if pair.ndim != 3 or pair.shape[0] != pair.shape[1]:
        raise ValueError(f"expected a square pair tensor, got {pair.shape}")
    if pair.shape[-1] != PAIR_DIM:
        raise ValueError(f"expected {PAIR_DIM} pair channels, got {pair.shape[-1]}")
    n_target = pair.shape[0] - n_bait
    if n_target < 1:
        raise ValueError(
            f"pair tensor of size {pair.shape[0]} has no target left after "
            f"removing {n_bait} baits"
        )
    pair_a = pair[:-n_bait, -n_bait:]                 # (L, n_bait, 128)
    pair_b = pair[-n_bait:, :-n_bait].swapaxes(0, 1)  # (L, n_bait, 128)
    x = np.concatenate(
        [pair_a.reshape(n_target, -1), pair_b.reshape(n_target, -1)], axis=-1
    )
    return np.asarray(x, dtype=np.float32)


def logit_aa(x: np.ndarray, head: dict) -> np.ndarray:
    """Per-bait-amino-acid logit contributions, shape (L, 20).

    Summing over the 20 baits gives the binding logit; the individual terms are
    the "activation" profile the paper correlates with ligand chemistry.
    """
    x = np.asarray(x, dtype=np.float32)
    if x.shape[-1] != N_FEATURES:
        raise ValueError(f"expected {N_FEATURES} features, got {x.shape[-1]}")
    z = (x - head["mean"]) / head["std"]
    z = z * head["w"][:, 0] + head["b"] / x.shape[-1]
    # (L, 2, 20, 128) -> sum over the A/B direction and the 128 channels
    return z.reshape(z.shape[0], 2, N_BAIT, -1).sum((1, 3))


def sigmoid(z: np.ndarray) -> np.ndarray:
    """Numerically stable logistic function."""
    z = np.asarray(z, dtype=np.float64)
    out = np.empty_like(z)
    pos, neg = z >= 0, z < 0
    out[pos] = 1.0 / (1.0 + np.exp(-z[pos]))
    e = np.exp(z[neg])
    out[neg] = e / (1.0 + e)
    return out


def predict(x: np.ndarray, heads: list[dict]) -> dict[str, np.ndarray]:
    """Score features with one or more heads.

    With a single head this reproduces the reference implementation exactly.
    With several (the 10 released folds), logits are averaged *before* the
    sigmoid, which is the standard way to ensemble logistic models and keeps the
    per-bait decomposition additive.
    """
    if not heads:
        raise ValueError("no heads supplied")
    per_seed = np.stack([logit_aa(x, h) for h in heads])  # (S, L, 20)
    mean_aa = per_seed.mean(0)
    logit = mean_aa.sum(-1)
    return {
        "p_bind": sigmoid(logit),
        "logit": logit,
        "p_bind_aa": mean_aa,
        "p_bind_per_seed": sigmoid(per_seed.sum(-1)),  # (S, L)
    }


def blosum_reorder(p_bind_aa: np.ndarray) -> tuple[np.ndarray, list[str]]:
    """Reorder the (L, 20) bait axis from A..Y into the paper's BLOSUM grouping."""
    idx = np.array([BAIT_SEQ.index(c) for c in BLOSUM_ORDER])
    return np.asarray(p_bind_aa)[:, idx], list(BLOSUM_ORDER)
