"""Human-proteome cross-reactivity screen.

Why this step is not optional
-----------------------------
Affinity-enhanced TCRs have caused lethal off-target toxicity in the clinic.
The most instructive case is an affinity-matured MAGE-A3 TCR that cross-reacted
with a titin-derived peptide presented on cardiomyocytes, causing fatal cardiac
toxicity in two patients (Linette et al., Blood 2013; Cameron et al., Sci
Transl Med 2013). The cross-reactive peptide was not similar to the intended
target by eye; it was similar at the positions the TCR actually reads.

So the screen here does not compare whole peptides. It weights each peptide
position by how much contact the TCR makes with it in the crystal structure,
then looks for human peptides that match at those positions, after filtering
for peptides that could be presented by the relevant HLA allele at all.

Limits, stated plainly
----------------------
* This is a *similarity screen*, not a validated cross-reactivity predictor.
  It produces a prioritised shortlist for experimental work, nothing more.
* The HLA anchor filter is a coarse motif. Real presentation should be scored
  with NetMHCpan or an equivalent, and confirmed by mass spectrometry of the
  relevant tissue.
* Position-weighted similarity cannot see conformational effects: a peptide
  that matches on paper may present a different backbone in the groove.
* It screens the canonical reviewed proteome. It does not cover isoforms,
  post-translational modifications, non-canonical translation products,
  microbial peptides, or tissue-restricted splice variants.
* The field standard remains a combinatorial-peptide-library X-scan of the
  actual receptor, used to build the motif that is then searched against the
  proteome. Deriving the motif from one structure, as here, is cheaper and
  weaker.
"""

from __future__ import annotations

import gzip
import re
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from Bio.Align import substitution_matrices

from .constants import HLA_A0201_ANCHORS
from .structure import Complex, _count_close_atoms, _heavy_atoms

BLOSUM62 = substitution_matrices.load("BLOSUM62")

UNIPROT_STREAM = (
    "https://rest.uniprot.org/uniprotkb/stream"
    "?query=%28proteome%3A{proteome}%29%20AND%20%28reviewed%3A{reviewed}%29"
    "&format=fasta&compressed=true"
)
HUMAN_PROTEOME = "UP000005640"

HEADER = re.compile(r"^>(?:sp|tr)\|(?P<acc>[^|]+)\|(?P<entry>\S+)\s+(?P<name>.*?)(?:\s+OS=|$)")
GENE = re.compile(r"\bGN=(?P<gene>\S+)")


@dataclass
class ProteomeHit:
    """A human peptide that resembles the target at TCR-contacted positions."""

    accession: str
    entry: str
    gene: str
    protein: str
    position: int
    peptide: str
    similarity: float
    identical_contacts: int
    contact_positions_matched: str

    def as_row(self) -> dict:
        return {
            "peptide": self.peptide,
            "similarity": round(self.similarity, 4),
            "identical_at_contacts": self.identical_contacts,
            "contact_match": self.contact_positions_matched,
            "gene": self.gene,
            "accession": self.accession,
            "protein": self.protein,
            "position": self.position,
        }


def download_proteome(cache_dir: Path, reviewed: bool = True) -> Path:
    """Fetch the human proteome FASTA from UniProt, caching it."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    suffix = "reviewed" if reviewed else "all"
    path = cache_dir / f"human_{suffix}.fasta.gz"
    if path.exists() and path.stat().st_size > 1_000_000:
        return path
    url = UNIPROT_STREAM.format(
        proteome=HUMAN_PROTEOME, reviewed="true" if reviewed else "false"
    )
    with urllib.request.urlopen(url, timeout=600) as response:
        path.write_bytes(response.read())
    return path


def iter_proteome(path: Path):
    """Yield ``(accession, entry, gene, protein_name, sequence)`` per entry."""
    accession = entry = gene = name = None
    chunks: list[str] = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if accession:
                    yield accession, entry, gene, name, "".join(chunks)
                match = HEADER.match(line.strip())
                if match:
                    accession = match.group("acc")
                    entry = match.group("entry")
                    name = match.group("name")
                else:
                    accession, entry, name = line[1:].strip(), "?", "?"
                gene_match = GENE.search(line)
                gene = gene_match.group("gene") if gene_match else "-"
                chunks = []
            else:
                chunks.append(line.strip())
    if accession:
        yield accession, entry, gene, name, "".join(chunks)


def peptide_contact_profile(complex_obj: Complex, cutoff: float = 5.0) -> list[int]:
    """Count TCR heavy-atom contacts to each peptide position.

    This is the weighting that makes the screen meaningful: it is a direct
    measurement of which peptide residues the receptor actually reads.
    """
    structure_path = complex_obj.path
    from Bio.PDB import PDBParser

    model = next(iter(PDBParser(QUIET=True).get_structure("x", str(structure_path))))
    peptide_chain = complex_obj.chains["peptide"]
    tcr_chains = complex_obj.tcr_chain_ids()

    tcr_atoms = _heavy_atoms(
        [r for cid in tcr_chains for r in model[cid] if r.id[0] == " "]
    )
    profile: list[int] = []
    for residue in model[peptide_chain]:
        if residue.id[0] != " ":
            continue
        atoms = _heavy_atoms([residue])
        profile.append(_count_close_atoms(atoms, tcr_atoms, cutoff))
    return profile


def _similarity(
    target: str, candidate: str, weights: np.ndarray
) -> tuple[float, int, str]:
    """Weighted BLOSUM62 similarity restricted to the TCR-contacted positions."""
    total = 0.0
    best = 0.0
    identical = 0
    pattern = []
    for i, (a, b) in enumerate(zip(target, candidate)):
        weight = weights[i]
        if weight <= 0:
            pattern.append(".")
            continue
        try:
            score = float(BLOSUM62[a, b])
            self_score = float(BLOSUM62[a, a])
        except (KeyError, IndexError):
            return 0.0, 0, ""
        floor = -4.0
        total += weight * (score - floor)
        best += weight * (self_score - floor)
        if a == b:
            identical += 1
            pattern.append(a)
        elif score > 0:
            pattern.append(b.lower())
        else:
            pattern.append("-")
    if best <= 0:
        return 0.0, 0, ""
    return total / best, identical, "".join(pattern)


def scan_proteome(
    complex_obj: Complex,
    proteome_path: Path,
    anchors: dict[int, set[str]] | None = None,
    top_n: int = 40,
    min_similarity: float = 0.72,
    contact_weight_floor: int = 4,
    exclude_self: bool = True,
    require_identity_top_n: int = 1,
) -> tuple[list[ProteomeHit], dict]:
    """Search the proteome for peptides resembling the target at contact positions.

    Returns the ranked hits plus a dictionary describing the motif used, so the
    search is reproducible and auditable.
    """
    target = complex_obj.peptide
    length = len(target)
    profile = peptide_contact_profile(complex_obj)
    if len(profile) != length:
        profile = (profile + [0] * length)[:length]

    weights = np.array(
        [count if count >= contact_weight_floor else 0 for count in profile],
        dtype=float,
    )
    if weights.sum() == 0:
        weights = np.array(profile, dtype=float)

    anchors = HLA_A0201_ANCHORS if anchors is None else anchors
    # Anchor keys follow the immunology convention: positive values are 1-based
    # positions from the N terminus, negative values count back from the C
    # terminus, so -1 is the C-terminal anchor (P-omega) whatever the length.
    anchor_checks = [
        (position - 1 if position > 0 else length + position, allowed)
        for position, allowed in anchors.items()
    ]
    anchor_checks = [
        (index, allowed) for index, allowed in anchor_checks if 0 <= index < length
    ]

    # Combinatorial X-scan of a real receptor typically shows a few peptide
    # positions where substitution abolishes recognition. Those positions are
    # effectively invariant in any cross-reactive peptide, so requiring exact
    # identity at the most heavily contacted ones is what turns a broad
    # similarity list into a risk list.
    hot_positions = [
        int(i) for i in np.argsort(-weights)[:require_identity_top_n] if weights[i] > 0
    ]

    hits: list[ProteomeHit] = []
    windows_seen = 0
    windows_passing_anchors = 0

    for accession, entry, gene, name, sequence in iter_proteome(proteome_path):
        if len(sequence) < length:
            continue
        for start in range(len(sequence) - length + 1):
            windows_seen += 1
            candidate = sequence[start : start + length]
            ok = True
            for index, allowed in anchor_checks:
                if candidate[index] not in allowed:
                    ok = False
                    break
            if not ok:
                continue
            for index in hot_positions:
                if candidate[index] != target[index]:
                    ok = False
                    break
            if not ok:
                continue
            windows_passing_anchors += 1
            similarity, identical, pattern = _similarity(target, candidate, weights)
            if similarity < min_similarity:
                continue
            if exclude_self and candidate == target:
                continue
            hits.append(
                ProteomeHit(
                    accession=accession,
                    entry=entry,
                    gene=gene,
                    protein=name,
                    position=start + 1,
                    peptide=candidate,
                    similarity=similarity,
                    identical_contacts=identical,
                    contact_positions_matched=pattern,
                )
            )

    hits.sort(key=lambda h: (-h.similarity, -h.identical_contacts))
    # Keep one best hit per gene so a single repetitive protein cannot fill the
    # shortlist.
    deduplicated: list[ProteomeHit] = []
    seen_genes: set[str] = set()
    for hit in hits:
        key = (hit.gene, hit.peptide)
        if key in seen_genes:
            continue
        seen_genes.add(key)
        deduplicated.append(hit)
        if len(deduplicated) >= top_n:
            break

    motif = {
        "target_peptide": target,
        "tcr_contacts_per_position": profile,
        "weights_used": weights.tolist(),
        "weighted_positions": [
            f"P{i + 1}{target[i]}" for i in range(length) if weights[i] > 0
        ],
        "anchor_filter": {
            f"P{i + 1}": "".join(sorted(allowed)) for i, allowed in anchor_checks
        },
        "identity_required_at": [f"P{i + 1}{target[i]}" for i in hot_positions],
        "windows_scanned": windows_seen,
        "windows_passing_all_filters": windows_passing_anchors,
        "hits_above_threshold": len(hits),
        "min_similarity": min_similarity,
    }
    return deduplicated, motif
