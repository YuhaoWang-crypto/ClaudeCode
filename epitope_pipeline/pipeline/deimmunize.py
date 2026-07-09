"""Structure-preserving de-immunization loop (ProteinMPNN + epitope re-scoring).

Idea ("epitope masking" by inverse folding)
--------------------------------------------
1. Predict epitopes on the antigen (B-cell: BepiPred; T-cell: NetMHCpan).
2. Give ProteinMPNN the backbone (PDB) and *fix every residue except the ones
   inside predicted epitopes* (and any functional/active-site residues you name).
   ProteinMPNN then redesigns only the epitope-region surface residues with
   sequences that still fold into the same backbone.
3. Re-score each designed variant's epitope load with an amino-acid-level
   predictor and keep the variants whose epitope load dropped while ProteinMPNN's
   own score (sequence<->structure compatibility, a foldability proxy) stayed good.

Why the re-scorer is NOT EDEN
-----------------------------
EDEN requires a *native nucleotide CDS* and explicitly rejects amino-acid /
back-translated / synthetic sequences as out-of-distribution. ProteinMPNN emits
amino-acid sequences, so EDEN cannot validly score them. The right oracle here is
an AA-level epitope model: BepiPred (B-cell) or NetMHCpan %Rank (T-cell). This
module therefore takes a pluggable `epitope_score` callable and defaults to the
DTU tools already wired into this package.

Nothing here replaces wet-lab validation: ProteinMPNN preserves fold
compatibility, not necessarily function, stability, or true immunogenicity.
"""

from __future__ import annotations

import json
import os
import subprocess
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence


# --------------------------------------------------------------------------- #
# Designed-variant record
# --------------------------------------------------------------------------- #
@dataclass
class Variant:
    name: str
    sequence: str
    mpnn_score: Optional[float] = None      # ProteinMPNN global score (lower=better fit)
    seq_recovery: Optional[float] = None    # fraction identical to wild type
    epitope_score: Optional[float] = None   # from the AA-level oracle (lower=less immunogenic)
    n_mutations: Optional[int] = None
    extra: dict = field(default_factory=dict)

    def as_row(self) -> dict:
        return {
            "name": self.name, "sequence": self.sequence,
            "mpnn_score": self.mpnn_score, "seq_recovery": self.seq_recovery,
            "epitope_score": self.epitope_score, "n_mutations": self.n_mutations,
        }


# --------------------------------------------------------------------------- #
# Fixing scheme: everything fixed EXCEPT epitope residues (+ never-touch stays fixed)
# --------------------------------------------------------------------------- #
def positions_to_redesign(seq_len: int, epitope_positions: Sequence[int],
                          keep_fixed: Optional[Sequence[int]] = None) -> List[int]:
    """Return the 1-indexed positions ProteinMPNN is allowed to vary.

    Only epitope residues are opened up; explicit `keep_fixed` (e.g. catalytic
    residues) are removed from the redesign set even if flagged as epitope.
    """
    keep = set(keep_fixed or [])
    design = sorted({p for p in epitope_positions if 1 <= p <= seq_len} - keep)
    return design


def fixed_positions_dict(pdb_id: str, chain: str, seq_len: int,
                         redesign_positions: Sequence[int]) -> dict:
    """Build ProteinMPNN --fixed_positions JSON: list the positions to KEEP fixed
    (all positions except the redesign set)."""
    design = set(redesign_positions)
    fixed = [p for p in range(1, seq_len + 1) if p not in design]
    return {pdb_id: {chain: fixed}}


# --------------------------------------------------------------------------- #
# ProteinMPNN backend (official CLI: dauparas/ProteinMPNN protein_mpnn_run.py)
# --------------------------------------------------------------------------- #
def run_proteinmpnn_cli(pdb_path: str, fixed_positions: dict, out_dir: str,
                        num_seq: int = 16, sampling_temp: float = 0.1,
                        binary: str = "protein_mpnn_run.py",
                        timeout: int = 1800) -> str:
    """Run the official ProteinMPNN CLI and return the output FASTA path.

    Requires the ProteinMPNN repo on PATH (MIT-licensed, CPU-friendly). Raises
    if the binary is missing so the caller can fall back / report cleanly.
    """
    import shutil
    exe = shutil.which(binary)
    if not exe:
        raise FileNotFoundError(
            f"'{binary}' not on PATH. Clone github.com/dauparas/ProteinMPNN and "
            f"expose protein_mpnn_run.py, or pass a custom `binary`.")
    os.makedirs(out_dir, exist_ok=True)
    fp_path = os.path.join(out_dir, "fixed_positions.jsonl")
    with open(fp_path, "w") as fh:
        json.dump(fixed_positions, fh)
    cmd = [exe, "--pdb_path", pdb_path, "--fixed_positions_jsonl", fp_path,
           "--out_folder", out_dir, "--num_seq_per_target", str(num_seq),
           "--sampling_temp", str(sampling_temp)]
    subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=timeout)
    seqs_dir = os.path.join(out_dir, "seqs")
    fastas = [os.path.join(seqs_dir, f) for f in os.listdir(seqs_dir)
              if f.endswith(".fa")] if os.path.isdir(seqs_dir) else []
    if not fastas:
        raise RuntimeError("ProteinMPNN produced no seqs/*.fa output")
    return fastas[0]


def parse_mpnn_fasta(path: str) -> List[Variant]:
    """Parse ProteinMPNN's output FASTA.

    Header format: `>...,  score=1.23, ... seq_recovery=0.45`. The first record
    is the wild-type reconstruction; subsequent ones are the samples.
    """
    variants: List[Variant] = []
    header, seq = None, []
    idx = 0

    def flush():
        nonlocal header, seq, idx
        if header is None:
            return
        s = "".join(seq)
        meta = {}
        for tok in header.replace(">", "").split(","):
            if "=" in tok:
                k, v = tok.split("=", 1)
                try:
                    meta[k.strip()] = float(v.strip())
                except ValueError:
                    meta[k.strip()] = v.strip()
        variants.append(Variant(
            name=f"design_{idx}" if idx else "wildtype", sequence=s,
            mpnn_score=meta.get("score"),
            seq_recovery=meta.get("seq_recovery"),
            extra=meta))
        idx += 1

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                header, seq = line, []
            elif line:
                seq.append(line.strip())
    flush()
    return variants


def hamming(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b) if x != y) + abs(len(a) - len(b))


# --------------------------------------------------------------------------- #
# The de-immunization loop
# --------------------------------------------------------------------------- #
def deimmunize(
    wildtype_seq: str,
    epitope_positions: Sequence[int],
    redesign_fn: Callable[[List[int]], List[Variant]],
    epitope_score: Callable[[str], float],
    keep_fixed: Optional[Sequence[int]] = None,
    max_seq_recovery_drop: float = 0.35,
    max_mpnn_score: Optional[float] = None,
) -> Dict:
    """Run one round of structure-preserving de-immunization.

    Parameters
    ----------
    wildtype_seq : original antigen amino-acid sequence.
    epitope_positions : 1-indexed positions predicted to be epitope residues.
    redesign_fn : callable given the redesign positions, returns candidate
        Variants (this is where ProteinMPNN plugs in; injected for testability).
    epitope_score : AA-level oracle; lower == less immunogenic (BepiPred / NetMHCpan).
    keep_fixed : positions never to mutate (active site, disulfides, ...).
    max_seq_recovery_drop : reject variants that diverge too far from wild type.
    max_mpnn_score : optional ceiling on ProteinMPNN score (foldability guard).

    Returns a dict with the wild-type baseline, all scored variants, and the
    accepted (lower-immunogenicity, structure-preserving) subset, ranked.
    """
    design_pos = positions_to_redesign(len(wildtype_seq), epitope_positions, keep_fixed)
    wt_epitope = epitope_score(wildtype_seq)

    candidates = [v for v in redesign_fn(design_pos)
                  if v.sequence and v.name != "wildtype"]
    for v in candidates:
        v.epitope_score = epitope_score(v.sequence)
        v.n_mutations = hamming(v.sequence, wildtype_seq)
        if v.seq_recovery is None and len(v.sequence) == len(wildtype_seq):
            v.seq_recovery = 1.0 - v.n_mutations / len(wildtype_seq)

    def acceptable(v: Variant) -> bool:
        if v.epitope_score is None or v.epitope_score >= wt_epitope:
            return False  # must actually reduce immunogenicity
        if v.seq_recovery is not None and v.seq_recovery < (1.0 - max_seq_recovery_drop):
            return False  # drifted too far from the natural protein
        if max_mpnn_score is not None and v.mpnn_score is not None \
                and v.mpnn_score > max_mpnn_score:
            return False  # poor sequence<->structure fit
        return True

    accepted = sorted((v for v in candidates if acceptable(v)),
                      key=lambda v: (v.epitope_score, v.mpnn_score or 0.0))

    return {
        "wildtype_epitope_score": wt_epitope,
        "redesign_positions": design_pos,
        "n_candidates": len(candidates),
        "n_accepted": len(accepted),
        "best": accepted[0].as_row() if accepted else None,
        "accepted": [v.as_row() for v in accepted],
        "all": [v.as_row() for v in candidates],
    }
