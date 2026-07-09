"""Minimal FASTA parsing/validation (stdlib only)."""

from __future__ import annotations

import os
from typing import Dict, List, Tuple

_AA = set("ACDEFGHIKLMNPQRSTVWYXBZJUO*")


def parse_fasta(path: str) -> List[Tuple[str, str]]:
    """Return list of (id, sequence). `id` is the first token after '>'."""
    records: List[Tuple[str, str]] = []
    seq_id, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if seq_id is not None:
                    records.append((seq_id, "".join(chunks)))
                seq_id = line[1:].split()[0] or f"seq{len(records)+1}"
                chunks = []
            else:
                chunks.append(line.strip())
        if seq_id is not None:
            records.append((seq_id, "".join(chunks)))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records


def validate(records: List[Tuple[str, str]]) -> List[str]:
    """Return a list of human-readable warnings (empty == clean)."""
    warnings: List[str] = []
    seen = set()
    for sid, seq in records:
        if sid in seen:
            warnings.append(f"duplicate sequence id: {sid}")
        seen.add(sid)
        if not seq:
            warnings.append(f"{sid}: empty sequence")
            continue
        bad = sorted({c for c in seq.upper() if c not in _AA})
        if bad:
            warnings.append(f"{sid}: non-standard residues {bad}")
    return warnings


def ids(records: List[Tuple[str, str]]) -> List[str]:
    return [r[0] for r in records]


def write_single(record: Tuple[str, str], path: str) -> str:
    sid, seq = record
    with open(path, "w") as fh:
        fh.write(f">{sid}\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")
    return path
