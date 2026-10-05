"""AAindex + FFT sequence encoding ("protein spectrum").

The encoding, following Cadet et al.'s Innov'SAR and the FP4COM implementation:

1. map every residue to a scalar through an AAindex entry      -> x[0..L-1]
2. subtract the mean of x                                      (kills the DC term)
3. discrete Fourier transform                                  F[j] = sum_k x_k e^(-2i.pi.j.k/L)
4. take the magnitude spectrum and divide by its maximum       -> amplitudes in [0, 1]
5. keep the first floor(L/2) bins                              (the spectrum is symmetric)

Several AAindex entries are concatenated along the feature axis.

Two code paths produce *bit-identical* features:

`encode_sequences`
    The reference path. One FFT per sequence; works on arbitrary sequences.

`MutationEncoder`
    The fast path, for the case this method is actually used in: many variants of
    one parent differing at a handful of positions. Because steps 1-3 are linear
    in the residue values, the spectrum of a variant is an exact rank-k update of
    the parent's spectrum,

        F_var[j] = F_parent[j] + sum_p d_p e^(-2i.pi.j.p/L) - (sum_p d_p) [j == 0]

    with d_p the change in descriptor value at position p. That turns an
    O(n . L log L) encode into one O(n . k . L) matrix product, which is what makes
    exhaustive spaces (2^20 variants) and 566-index screens affordable.
    `test_encode.py` asserts the two paths agree to floating-point tolerance.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np

from .variants import parse_variant

__all__ = [
    "aaindex_values",
    "all_indices",
    "available_indices",
    "incomplete_indices",
    "IncompleteIndexError",
    "spectrum",
    "encode_sequences",
    "MutationEncoder",
]


@lru_cache(maxsize=None)
def _aaindex1():
    try:
        from aaindex import aaindex1
    except ImportError as exc:  # pragma: no cover - dependency is declared
        raise ImportError(
            "the `aaindex` package is required for AAindex encodings: pip install aaindex"
        ) from exc
    return aaindex1


AA20 = "ACDEFGHIKLMNPQRSTVWY"


@lru_cache(maxsize=None)
def all_indices() -> tuple[str, ...]:
    """Accession codes of every AAindex1 entry (566 in aaindex 1.0.5 and 1.3.2)."""
    return tuple(_aaindex1().record_codes())


@lru_cache(maxsize=None)
def incomplete_indices() -> tuple[str, ...]:
    """Entries missing a value for at least one of the 20 standard residues.

    These cannot be used as encodings. The `aaindex` package changed here: all 566
    entries carry full values in 1.0.5, while 1.3.2 returns ``None`` for 13 of them
    (the AVBF000101-109 series, GUYH850103, ROSM880104, ROSM880105, YANJ020101).

    ⚠️ Reproducibility consequence: the paper's round-3 model used the triple
    ``AVBF000109_JUNJ780101_JUKT750101``, and AVBF000109 is one of the affected
    entries. That model therefore cannot be rebuilt from current `aaindex` data at
    all -- not approximately, but not at all -- without pinning `aaindex==1.0.5`
    or supplying the missing value by hand.
    """
    return tuple(
        code
        for code in all_indices()
        if any(aaindex_values(code).get(aa) is None for aa in AA20)
    )


@lru_cache(maxsize=None)
def available_indices(complete_only: bool = True) -> tuple[str, ...]:
    """Accession codes usable as encodings.

    `complete_only` (the default) excludes entries with missing residue values,
    which would otherwise fail mid-screen and be silently dropped from the
    candidate pool. Pass ``False`` to get all 566.
    """
    if not complete_only:
        return all_indices()
    skip = set(incomplete_indices())
    return tuple(code for code in all_indices() if code not in skip)


@lru_cache(maxsize=None)
def aaindex_values(code: str) -> dict[str, float]:
    """The residue values of one AAindex entry, keyed by one-letter code.

    Values may be ``None`` where the source record has none; `_lookup_table`
    rejects those rather than letting them become NaN features.
    """
    try:
        record = _aaindex1()[code]
    except Exception as exc:
        raise KeyError(f"unknown AAindex accession {code!r}") from exc
    return dict(record.values)


class IncompleteIndexError(ValueError):
    """Raised when an AAindex entry lacks a value for some standard residue."""


def _lookup_table(code: str) -> dict[str, float]:
    table = aaindex_values(code)
    missing = [aa for aa in AA20 if table.get(aa) is None]
    if missing:
        raise IncompleteIndexError(
            f"AAindex entry {code} has no value for {', '.join(missing)}; it cannot be used "
            "as an encoding. Your `aaindex` package version omits values that aaindex 1.0.5 "
            f"supplied for {len(incomplete_indices())} entries -- pin aaindex==1.0.5 to use them."
        )
    return table


def _numeric(seq: str, table: dict[str, float]) -> np.ndarray:
    try:
        return np.array([table[aa] for aa in seq], dtype=np.float64)
    except KeyError as exc:
        raise KeyError(
            f"residue {exc.args[0]!r} has no value in this AAindex entry "
            "(non-standard residues are not supported)"
        ) from exc


def spectrum(x: np.ndarray, pad_to: int | None = None, drop_dc: bool = True) -> np.ndarray:
    """Normalized half magnitude spectrum of a residue-value vector.

    `pad_to` zero-pads *after* mean-centering, so that sequences of different
    length yield feature vectors of equal width.

    `drop_dc` discards bin 0. **Keep this on.** Step 2 subtracts the mean, which
    sets the DC term to exactly zero, so bin 0 contains nothing but floating-point
    round-off (~1e-16). That is harmless until the regressor standardizes its
    columns -- `PLSRegression(scale=True)` divides each feature by its standard
    deviation, and for bin 0 that deviation *is* the round-off, so the column gets
    amplified to O(1) and PLS fits noise through it. The effect is not academic:
    on the paper's round-1 training set, including bin 0 makes the cross-validated
    error depend on whether the features were stored as float32 or float64
    (cvMSE 0.359 vs 0.900 at 10 components) and changes which of the 566 AAindex
    entries "wins". With bin 0 dropped the two agree to machine precision.
    See `m2_dc_artifact.py` for the demonstration.
    """
    x = np.asarray(x, dtype=np.float64)
    x = x - x.mean()
    if pad_to is not None:
        if pad_to < len(x):
            raise ValueError(f"pad_to={pad_to} is shorter than the sequence ({len(x)})")
        x = np.concatenate([x, np.zeros(pad_to - len(x))])
    n = len(x)
    mags = np.abs(np.fft.rfft(x))  # bins 0..n//2; the rest of the spectrum mirrors these
    peak = mags.max()
    if peak == 0.0:
        # A constant descriptor (every residue identical) carries no signal.
        return np.zeros(n // 2 - (1 if drop_dc else 0))
    return (mags / peak)[1 if drop_dc else 0 : n // 2]


def encode_sequences(
    sequences, indices, pad_to: int | None = None, drop_dc: bool = True
) -> np.ndarray:
    """Reference encoder: ``(n_sequences, n_indices * n_bins)`` feature matrix."""
    sequences = list(sequences)
    indices = [indices] if isinstance(indices, str) else list(indices)
    if not indices:
        raise ValueError("at least one AAindex accession is required")
    if not sequences:
        raise ValueError("no sequences to encode")

    blocks = []
    for code in indices:
        table = _lookup_table(code)
        blocks.append(
            np.stack(
                [
                    spectrum(_numeric(seq, table), pad_to=pad_to, drop_dc=drop_dc)
                    for seq in sequences
                ]
            )
        )
    return np.concatenate(blocks, axis=1)


class MutationEncoder:
    """Fast exact encoder for variants of a fixed parent sequence.

    Parameters
    ----------
    parent:
        The parent (background) sequence every variant is expressed against.
    positions:
        1-based positions that may be mutated. Defaults to every position, which
        is still correct but gives up the speedup. Passing the handful of
        positions actually in play is what makes the fast path fast.
    """

    def __init__(self, parent: str, positions=None):
        self.parent = parent
        self.length = len(parent)
        self.n_bins = self.length // 2
        if positions is None:
            positions = range(1, self.length + 1)
        self.positions = sorted(set(int(p) for p in positions))
        for pos in self.positions:
            if not 1 <= pos <= self.length:
                raise ValueError(f"position {pos} is outside the parent (length {self.length})")
        self._slot = {pos: i for i, pos in enumerate(self.positions)}

        # Phase matrix: rows = mutable positions, columns = retained FFT bins.
        bins = np.arange(self.n_bins + 1)  # includes Nyquist, needed for the max
        idx = np.array([p - 1 for p in self.positions], dtype=np.float64)
        self._phase = np.exp(-2j * np.pi * np.outer(idx, bins) / self.length)

        self._cache: dict[str, tuple[np.ndarray, np.ndarray, dict[str, int]]] = {}

    def _parent_spectrum(self, code: str) -> tuple[np.ndarray, np.ndarray, dict[str, int]]:
        """``(parent rFFT bins 0..n//2, per-position descriptor deltas, residue->column)``."""
        if code not in self._cache:
            table = _lookup_table(code)
            x = _numeric(self.parent, table)
            f_parent = np.fft.rfft(x - x.mean())
            # delta[i, a] = value(residue a) - value(parent residue at positions[i])
            order = sorted(table)
            lut = np.array([table[a] for a in order])
            col = {a: j for j, a in enumerate(order)}
            base = np.array([table[self.parent[p - 1]] for p in self.positions])
            self._cache[code] = (f_parent, lut[None, :] - base[:, None], col)
        return self._cache[code]

    def _picks(self, variants) -> np.ndarray:
        """``(n_variants, n_positions)`` matrix of substituted residues, as ord() codes; -1 = parent."""
        picks = np.full((len(variants), len(self.positions)), -1, dtype=np.int64)
        for row, name in enumerate(variants):
            for sub in parse_variant(name):
                slot = self._slot.get(sub.pos)
                if slot is None:
                    raise ValueError(
                        f"variant {name!r} mutates position {sub.pos}, which was not declared "
                        "mutable for this encoder"
                    )
                if self.parent[sub.index] != sub.wt:
                    raise ValueError(
                        f"{sub} in {name!r} disagrees with the parent residue "
                        f"{self.parent[sub.index]}{sub.pos}"
                    )
                picks[row, slot] = ord(sub.mut)
        return picks

    def encode(self, variants, indices, drop_dc: bool = True) -> np.ndarray:
        """Encode variant *labels* (not sequences) into the same features as `encode_sequences`.

        `drop_dc` defaults to True; see `spectrum` for why bin 0 must go.
        """
        variants = list(variants)
        indices = [indices] if isinstance(indices, str) else list(indices)
        if not variants:
            raise ValueError("no variants to encode")
        if not indices:
            raise ValueError("at least one AAindex accession is required")

        picks = self._picks(variants)
        blocks = []
        for code in indices:
            f_parent, delta, col = self._parent_spectrum(code)
            # Per-variant descriptor change at each mutable position.
            d = np.zeros(picks.shape, dtype=np.float64)
            rows, cols = np.nonzero(picks >= 0)
            if len(rows):
                letters = [chr(c) for c in picks[rows, cols]]
                try:
                    aa_cols = np.array([col[a] for a in letters])
                except KeyError as exc:
                    raise KeyError(f"residue {exc.args[0]!r} is not in AAindex entry {code}") from exc
                d[rows, cols] = delta[cols, aa_cols]

            f = f_parent[None, :] + d @ self._phase
            f[:, 0] -= d.sum(axis=1)  # mean-centering shifts only the DC bin
            mags = np.abs(f)
            peak = mags.max(axis=1, keepdims=True)
            peak[peak == 0.0] = 1.0
            blocks.append((mags / peak)[:, 1 if drop_dc else 0 : self.n_bins])
        return np.concatenate(blocks, axis=1)
