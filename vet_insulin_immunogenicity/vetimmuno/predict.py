"""MHC-II binding backends, and the %rank problem for custom molecules.

Two backends, with very different standing:

``NetMHCIIpanBackend``   the real predictor, run in custom-molecule mode
                         (``-mhcfsa``) against DLA/FLA sequences. Requires a
                         licensed local install. This is what a real assessment
                         would use -- subject to the applicability-domain
                         warning from ``groove.py``.

``IllustrativeScorer``   a transparent physicochemical pocket-complementarity
                         score with NO training data behind it. It exists so the
                         pipeline, the controls and the reports can be exercised
                         end-to-end without a licensed binary. Its numbers are
                         plumbing, not evidence, and every artefact it produces
                         is labelled ILLUSTRATIVE.

The %rank problem
-----------------
In custom-molecule mode NetMHCIIpan cannot emit %Rank: rank reference
distributions only exist for its built-in alleles. ``BackgroundRank`` rebuilds
one the same way the tool does -- score a large set of background peptides with
the *same* molecule, then express any peptide's score as its percentile against
that molecule's own distribution. This makes scores comparable across molecules
whatever backend produced them, and it works identically for DLA and FLA.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from .epitope import CORE_LENGTH, Core, Peptide
from .groove import ALPHA_POCKETS, BETA_POCKETS, Molecule

# Swiss-Prot average amino-acid composition (UniProtKB/Swiss-Prot release
# statistics). Used to generate background peptides for rank calibration.
SWISSPROT_FREQ: Dict[str, float] = {
    "A": 0.0825, "Q": 0.0393, "L": 0.0965, "S": 0.0663, "R": 0.0553, "E": 0.0672,
    "K": 0.0580, "T": 0.0535, "N": 0.0406, "G": 0.0707, "M": 0.0241, "W": 0.0110,
    "D": 0.0546, "H": 0.0227, "F": 0.0386, "Y": 0.0292, "C": 0.0138, "I": 0.0591,
    "P": 0.0474, "V": 0.0686,
}

# Kyte-Doolittle hydropathy, formal charge at pH 7, side-chain volume (A^3).
HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
CHARGE = {aa: 0.0 for aa in HYDROPATHY}
CHARGE.update({"D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0, "H": 0.1})
VOLUME = {
    "A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5, "Q": 143.8,
    "E": 138.4, "G": 60.1, "H": 153.2, "I": 166.7, "L": 166.7, "K": 168.6,
    "M": 162.9, "F": 189.9, "P": 112.7, "S": 89.0, "T": 116.1, "W": 227.8,
    "Y": 193.6, "V": 140.0,
}

# Anchor residue index (0-based) within the 9-mer core for each pocket, and the
# relative weight of that pocket. P1 dominates MHC-II binding energetics.
ANCHORS: Dict[str, Tuple[int, float]] = {
    "P1": (0, 0.40), "P4": (3, 0.20), "P6": (5, 0.20), "P9": (8, 0.20)
}


@dataclass
class Hit:
    """One (peptide, register, molecule) scored triple."""

    chain: str
    peptide_start: int
    peptide: str
    core: str
    core_start: int
    molecule: str
    locus: str
    score: float
    rank: float = float("nan")
    model: str = ""


def _z(values: Iterable[float], value: float) -> float:
    arr = np.fromiter(values, dtype=float)
    sd = arr.std()
    return 0.0 if sd == 0 else (value - arr.mean()) / sd


_H_MEAN, _H_SD = float(np.mean(list(HYDROPATHY.values()))), float(np.std(list(HYDROPATHY.values())))
_V_MEAN, _V_SD = float(np.mean(list(VOLUME.values()))), float(np.std(list(VOLUME.values())))


def _zh(aa: str) -> float:
    return (HYDROPATHY.get(aa, 0.0) - _H_MEAN) / _H_SD


def _zv(aa: str) -> float:
    return (VOLUME.get(aa, 0.0) - _V_MEAN) / _V_SD


def _lut(fn) -> np.ndarray:
    table = np.zeros(256)
    for aa in HYDROPATHY:
        table[ord(aa)] = fn(aa)
    return table


_LUT_H = _lut(_zh)
_LUT_V = _lut(_zv)
_LUT_Q = _lut(lambda aa: CHARGE.get(aa, 0.0))


class IllustrativeScorer:
    """⚠️ ILLUSTRATIVE pocket-complementarity score -- not a trained predictor.

    For each pocket the anchor residue of the core is compared with the residues
    lining that pocket in the molecule's groove pseudosequence:

      hydrophobic anchor in a hydrophobic pocket   -> favourable
      charged anchor facing opposite charge        -> favourable
      bulky anchor in a pocket lined by bulky side chains -> penalised

    Every term is a textbook physicochemical rule. None of it is fitted to
    binding data, so the output ranks peptides plausibly and predicts nothing.
    Use it to exercise the pipeline; use NetMHCIIpan to make claims.

    One consequence is visible in the output and worth stating plainly: this
    scorer reads only the four anchor positions, so two cores that differ solely
    at a non-anchor position receive an identical score. A trained model does
    not behave that way. If two rows in a report share a score, that is the
    surrogate's blindness, not a finding.
    """

    name = "illustrative-pocket-surrogate"
    trained = False

    def __init__(self, steric_weight: float = 0.5):
        self.steric_weight = steric_weight
        self._pockets: Dict[str, Dict[str, Tuple[float, float, float]]] = {}

    def _lining(self, mol: Molecule, pocket: str) -> str:
        from .groove import ALPHA_CONTACT, BETA_CONTACT
        table = BETA_POCKETS if mol.chain == "beta" else ALPHA_POCKETS
        positions = table.get(pocket, ())
        contacts = BETA_CONTACT if mol.chain == "beta" else ALPHA_CONTACT
        idx = {p: i for i, p in enumerate(contacts)}
        return "".join(
            mol.pseudoseq[idx[p]] for p in positions
            if p in idx and mol.pseudoseq[idx[p]] != "-"
        )

    def pocket_params(self, mol: Molecule) -> Dict[str, Tuple[float, float, float]]:
        """Mean hydrophobicity / charge / volume of each pocket's lining."""
        cached = self._pockets.get(mol.name)
        if cached is not None:
            return cached
        params: Dict[str, Tuple[float, float, float]] = {}
        for pocket in ANCHORS:
            lining = self._lining(mol, pocket)
            if not lining:
                continue
            params[pocket] = (
                float(np.mean([_zh(x) for x in lining])),
                float(np.mean([CHARGE.get(x, 0.0) for x in lining])),
                float(np.mean([_zv(x) for x in lining])),
            )
        self._pockets[mol.name] = params
        return params

    def score_core(self, core: str, mol: Molecule) -> float:
        return float(self.score([core], mol)[0])

    def score(self, cores: Sequence[str], mol: Molecule) -> np.ndarray:
        """Vectorised over peptides -- the rank calibration scores 20k at a time."""
        if not cores:
            return np.zeros(0)
        n, width = len(cores), len(cores[0])
        chars = np.frombuffer("".join(cores).encode(), dtype=np.uint8).reshape(n, width)
        params = self.pocket_params(mol)
        total = np.zeros(n)
        for pocket, (anchor_idx, weight) in ANCHORS.items():
            if pocket not in params or anchor_idx >= width:
                continue
            h_p, q_p, v_p = params[pocket]
            col = chars[:, anchor_idx]
            fit = (
                _LUT_H[col] * h_p
                - _LUT_Q[col] * q_p
                - self.steric_weight * np.maximum(0.0, _LUT_V[col] + v_p)
            )
            total += weight * fit
        return total


class NetMHCIIpanBackend:
    """NetMHCIIpan-4.3 in custom-molecule mode. Requires a licensed local install."""

    name = "netmhciipan-4.3"
    trained = True

    def __init__(self, binary: Optional[str] = None, length: int = 15):
        self.binary = binary or os.environ.get("NETMHCIIPAN") or shutil.which("netMHCIIpan")
        self.length = length

    @property
    def available(self) -> bool:
        return bool(self.binary) and Path(self.binary).exists()

    def score(self, cores: Sequence[str], mol: Molecule) -> np.ndarray:
        if not self.available:
            raise RuntimeError(
                "netMHCIIpan binary not found. Set $NETMHCIIPAN, or run the "
                "generated netmhciipan/run_netmhciipan.sh on a licensed host and "
                "re-run with --backend netmhciipan --parse-only."
            )
        with_tmp = Path(os.environ.get("TMPDIR", "/tmp")) / f"vetimmuno_{os.getpid()}"
        with_tmp.mkdir(parents=True, exist_ok=True)
        pep = with_tmp / "peptides.txt"
        pep.write_text("\n".join(cores) + "\n")
        mhc = with_tmp / "molecule.fasta"
        mhc.write_text(f">{re.sub(r'[^A-Za-z0-9]', '_', mol.name)}\n{mol.sequence}\n")
        proc = subprocess.run(
            [self.binary, "-f", str(pep), "-inptype", "1", "-mhcfsa", str(mhc)],
            capture_output=True, text=True, check=True,
        )
        return parse_netmhciipan(proc.stdout, len(cores))


def parse_netmhciipan(text: str, expected: int) -> np.ndarray:
    """Pull the EL/binding score column out of NetMHCIIpan stdout.

    Tolerant by design: column layout differs between releases and between
    built-in-allele and custom-molecule mode, so the header row is read rather
    than assumed.
    """
    scores: List[float] = []
    columns: Optional[List[str]] = None
    score_idx: Optional[int] = None
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split()
        if "Peptide" in fields and score_idx is None:
            columns = fields
            for candidate in ("Score_EL", "Score", "1-log50k", "Score_BA"):
                if candidate in columns:
                    score_idx = columns.index(candidate)
                    break
            continue
        if score_idx is not None and len(fields) > score_idx:
            try:
                scores.append(float(fields[score_idx]))
            except ValueError:
                continue
    if len(scores) != expected:
        raise RuntimeError(
            f"parsed {len(scores)} scores from NetMHCIIpan output, expected {expected}"
        )
    return np.array(scores)


class BackgroundRank:
    """Per-molecule %rank from a background peptide distribution.

    This is the fix for "custom molecules have no %Rank". A molecule's raw score
    scale is arbitrary; its percentile against a fixed background of random
    peptides is not, and is comparable across molecules and across species.
    """

    def __init__(self, backend, n_background: int = 20000, seed: int = 20260826,
                 core_length: int = CORE_LENGTH):
        self.backend = backend
        self.seed = seed
        self.core_length = core_length
        self.n_background = n_background
        self._peptides = self._draw()
        self._dist: Dict[str, np.ndarray] = {}

    def _draw(self) -> List[str]:
        rng = np.random.default_rng(self.seed)
        letters = np.array(list(SWISSPROT_FREQ))
        probs = np.array([SWISSPROT_FREQ[a] for a in letters])
        probs = probs / probs.sum()
        draws = rng.choice(letters, size=(self.n_background, self.core_length), p=probs)
        return ["".join(row) for row in draws]

    def distribution(self, mol: Molecule) -> np.ndarray:
        if mol.name not in self._dist:
            self._dist[mol.name] = np.sort(self.backend.score(self._peptides, mol))
        return self._dist[mol.name]

    def rank(self, scores: np.ndarray, mol: Molecule) -> np.ndarray:
        """Percent of background peptides scoring at least as high (lower = stronger)."""
        dist = self.distribution(mol)
        above = len(dist) - np.searchsorted(dist, scores, side="left")
        return 100.0 * above / len(dist)


# Binder classification follows the NetMHCIIpan-4.x convention for MHC-II.
STRONG_RANK = 2.0
WEAK_RANK = 10.0


def classify(rank: float) -> str:
    if rank <= STRONG_RANK:
        return "strong"
    if rank <= WEAK_RANK:
        return "weak"
    return "none"


def score_peptides(peptides: Sequence[Peptide], molecules: Sequence[Molecule],
                   backend, ranker: BackgroundRank) -> List[Hit]:
    """Score every (register, molecule) pair and attach a background %rank."""
    cores: List[Core] = []
    owners: List[Peptide] = []
    for pep in peptides:
        for core in pep.cores:
            cores.append(core)
            owners.append(pep)
    if not cores:
        return []
    seqs = [c.sequence for c in cores]
    hits: List[Hit] = []
    for mol in molecules:
        scores = backend.score(seqs, mol)
        ranks = ranker.rank(scores, mol)
        for core, pep, s, r in zip(cores, owners, scores, ranks):
            hits.append(Hit(
                chain=core.chain, peptide_start=pep.start, peptide=pep.sequence,
                core=core.sequence, core_start=core.start, molecule=mol.name,
                locus=mol.locus, score=float(s), rank=float(r), model=backend.name,
            ))
    return hits
