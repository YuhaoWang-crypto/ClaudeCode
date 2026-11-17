r"""Sequence-identity clustering with mmseqs2 when available (optimisation #7).

mmseqs2 easy-cluster gives real sequence-identity clustering at scale; without
the binary we fall back to the dependency-free k-mer Jaccard clusterer in
``dedup.py``. Same interface either way: sequences -> integer cluster labels.

    # install: conda install -c bioconda mmseqs2   (or a static binary on PATH)
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from .dedup import cluster_sequences as _jaccard_cluster


def has_mmseqs() -> bool:
    return shutil.which("mmseqs") is not None


def cluster(seqs: list[str], min_seq_id: float = 0.4, coverage: float = 0.8,
            method: str = "auto") -> np.ndarray:
    """Return cluster labels (len n). method: 'auto' | 'mmseqs' | 'jaccard'."""
    if not seqs:
        return np.array([], dtype=int)
    use_mmseqs = method == "mmseqs" or (method == "auto" and has_mmseqs())
    if use_mmseqs:
        try:
            return _mmseqs_cluster(seqs, min_seq_id, coverage)
        except Exception:
            pass   # fall through to the proxy
    # Map an identity threshold to the Jaccard-proxy threshold (rough monotone).
    return _jaccard_cluster(seqs, threshold=max(0.2, min_seq_id * 0.6))


def _mmseqs_cluster(seqs: list[str], min_seq_id: float, coverage: float) -> np.ndarray:
    with tempfile.TemporaryDirectory(prefix="mmseqs_") as d:
        dd = Path(d)
        fasta = dd / "in.fasta"
        fasta.write_text("".join(f">{i}\n{s}\n" for i, s in enumerate(seqs)))
        subprocess.run(
            ["mmseqs", "easy-cluster", str(fasta), str(dd / "res"), str(dd / "tmp"),
             "--min-seq-id", str(min_seq_id), "-c", str(coverage), "--cov-mode", "1",
             "-v", "0"],
            check=True, capture_output=True)
        # res_cluster.tsv: <rep_id>\t<member_id> per line.
        tsv = (dd / "res_cluster.tsv").read_text().strip().splitlines()
        rep_of: dict[int, int] = {}
        for line in tsv:
            rep, member = line.split("\t")
            rep_of[int(member)] = int(rep)
        reps = {}
        labels = np.empty(len(seqs), dtype=int)
        for i in range(len(seqs)):
            r = rep_of.get(i, i)
            labels[i] = reps.setdefault(r, len(reps))
        return labels
