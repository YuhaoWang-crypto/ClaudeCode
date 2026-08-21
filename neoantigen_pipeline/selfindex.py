"""A memory-sane index of every k-mer in the human proteome.

Two questions need answering millions of times:

  "does this mutant peptide already exist somewhere in the self proteome?"
      -> the tumor-specificity gate
  "which self peptide is exactly one substitution away from this one?"
      -> the wild-type counterpart, used to mine validated neoepitopes

A Python `set` of 11.4 million 9-mer strings costs the better part of a
gigabyte and makes the second question slow. Encoding each k-mer as a base-20
integer in a sorted int64 array costs ~90 MB, answers membership by binary
search, and answers the one-mismatch question by generating the 19*k
substitution variants of the query and looking each one up.
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Optional

import numpy as np

AA = "ACDEFGHIKLMNPQRSTVWY"
_CODE = np.full(256, -1, dtype=np.int64)
for _i, _a in enumerate(AA):
    _CODE[ord(_a)] = _i


def _encode_peptide(pep: str) -> Optional[int]:
    v = 0
    for c in pep:
        code = _CODE[ord(c)] if ord(c) < 256 else -1
        if code < 0:
            return None
        v = v * 20 + int(code)
    return v


def _decode(value: int, k: int) -> str:
    out = []
    for _ in range(k):
        out.append(AA[value % 20])
        value //= 20
    return "".join(reversed(out))


class SelfKmerIndex:
    """Sorted int64 array of every valid k-mer in a proteome."""

    def __init__(self, proteome: Dict[str, str], k: int):
        self.k = k
        joined = "*".join(proteome.values())
        raw = np.frombuffer(joined.encode("ascii", "replace"), dtype=np.uint8)
        codes = _CODE[raw]
        n = len(codes) - k + 1
        if n <= 0:
            self.values = np.empty(0, dtype=np.int64)
            return
        vals = np.zeros(n, dtype=np.int64)
        valid = np.ones(n, dtype=bool)
        for j in range(k):
            win = codes[j:j + n]
            vals = vals * 20 + np.where(win >= 0, win, 0)
            valid &= win >= 0
        self.values = np.unique(vals[valid])

    def __len__(self) -> int:
        return int(len(self.values))

    def __contains__(self, pep: str) -> bool:
        if len(pep) != self.k:
            return False
        v = _encode_peptide(pep)
        if v is None:
            return False
        i = np.searchsorted(self.values, v)
        return bool(i < len(self.values) and self.values[i] == v)

    def contains_many(self, peptides: Iterable[str]) -> np.ndarray:
        peps = list(peptides)
        enc = np.array([_encode_peptide(p) if len(p) == self.k else None
                        for p in peps], dtype=object)
        out = np.zeros(len(peps), dtype=bool)
        ok = np.array([e is not None for e in enc])
        if not ok.any():
            return out
        vals = np.array([int(e) for e in enc[ok]], dtype=np.int64)
        idx = np.searchsorted(self.values, vals)
        idx_clipped = np.clip(idx, 0, max(0, len(self.values) - 1))
        hit = (idx < len(self.values)) & (self.values[idx_clipped] == vals)
        out[ok] = hit
        return out

    def one_mismatch(self, pep: str) -> List[str]:
        """Every self k-mer differing from `pep` at exactly one position.

        Returns [] when the peptide itself occurs in the proteome (then it is
        self, not a neoepitope) -- the caller wants that distinction.
        """
        if len(pep) != self.k or pep in self:
            return []
        base = _encode_peptide(pep)
        if base is None:
            return []
        cands: List[int] = []
        pw = [20 ** (self.k - 1 - i) for i in range(self.k)]
        for i, c in enumerate(pep):
            ci = int(_CODE[ord(c)]) if ord(c) < 256 else -1
            if ci < 0:
                return []
            for a in range(20):
                if a == ci:
                    continue
                cands.append(base + (a - ci) * pw[i])
        arr = np.array(cands, dtype=np.int64)
        idx = np.searchsorted(self.values, arr)
        idx_c = np.clip(idx, 0, max(0, len(self.values) - 1))
        hit = (idx < len(self.values)) & (self.values[idx_c] == arr)
        return [_decode(int(v), self.k) for v in arr[hit]]


def build_indexes(proteome: Dict[str, str], lengths: Iterable[int]) -> Dict[int, SelfKmerIndex]:
    return {int(k): SelfKmerIndex(proteome, int(k)) for k in sorted(set(int(x) for x in lengths))}
