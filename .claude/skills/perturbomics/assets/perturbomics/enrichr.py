"""Load REAL perturbation signatures from Enrichr gene-set libraries.

Enrichr (maayanlab.cloud) hosts curated perturbation libraries as flat,
tab-separated files — each line is ``term <TAB> desc <TAB> gene1 <TAB> gene2 ...``
These are genuine, downloadable perturbation signatures (no login), which makes
them the pragmatic real-data on-ramp for this skill when the full L1000 `.gctx`
compendium on clue.io is gated. Relevant libraries:

  * Disease_Perturbations_from_GEO_up / _down   (CREEDS disease DE signatures)
  * LINCS_L1000_Chem_Pert_up / _down            (drug perturbation signatures)
  * LINCS_L1000_CRISPR_KO_Consensus_Sigs        (CRISPR KO; terms "GENE Up/Down")
  * Single/ Gene_Perturbations_from_GEO_up / _down (genetic perturbations)

Every loader returns plain ``Signature`` objects (up/down -> ±1), so the rest
of the package (connectivity, combine) works unchanged.

RIGOR: the loaded gene lists and all downstream scores are ✅ deterministic.
The biology (a disease term's genes, a drug's mechanism) is only as good as the
underlying GEO/L1000 study — treat individual signatures as ⚠️ evidence, and
prefer consensus across replicates (see ``consensus_signature``).
"""

from __future__ import annotations

import os
import urllib.request
from collections import Counter
from typing import Callable

from .signature import Signature

_ENRICHR = ("https://maayanlab.cloud/Enrichr/geneSetLibrary"
            "?mode=text&libraryName={name}")


def load_library(name_or_path: str, cache_dir: str = ".") -> dict[str, list[str]]:
    """Return ``{term: [genes]}`` for an Enrichr library name OR a local path.

    If given a bare library name, downloads to ``cache_dir`` and caches. If
    given an existing path, parses it directly (offline).
    """
    if os.path.exists(name_or_path):
        path = name_or_path
    else:
        path = os.path.join(cache_dir, f"{name_or_path}.txt")
        if not os.path.exists(path):
            os.makedirs(cache_dir, exist_ok=True)
            urllib.request.urlretrieve(_ENRICHR.format(name=name_or_path), path)

    out: dict[str, list[str]] = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            term = parts[0]
            genes = [g.split(",")[0].strip().upper()
                     for g in parts[2:] if g.strip()]
            if genes:
                out[term] = genes
    return out


def paired_signatures(
    up: dict[str, list[str]], down: dict[str, list[str]],
    modality: str, name_fn: Callable[[str], str] | None = None,
    meta_fn: Callable[[str], dict] | None = None,
) -> list[Signature]:
    """Pair an _up and _down library by shared term -> two-tailed Signatures."""
    sigs = []
    for term in up:
        if term not in down:
            continue
        name = name_fn(term) if name_fn else term
        meta = meta_fn(term) if meta_fn else {"term": term}
        sigs.append(Signature.from_ranked(
            up=up[term], down=down[term], name=name, modality=modality,
            meta=meta))
    return sigs


def crispr_ko_signatures(
    consensus: dict[str, list[str]], modality: str = "crispr_ko"
) -> list[Signature]:
    """Parse LINCS_L1000_CRISPR_KO_Consensus_Sigs ("GENE Up"/"GENE Down") ->
    one two-tailed Signature per knocked-out gene."""
    up: dict[str, list[str]] = {}
    down: dict[str, list[str]] = {}
    for term, genes in consensus.items():
        toks = term.rsplit(" ", 1)
        if len(toks) != 2:
            continue
        gene, direction = toks[0], toks[1].lower()
        (up if direction == "up" else down)[gene] = genes
    sigs = []
    for gene in up:
        if gene not in down:
            continue
        sigs.append(Signature.from_ranked(
            up=up[gene], down=down[gene], name=f"KO_{gene}",
            modality=modality, meta={"target": gene}))
    return sigs


def consensus_signature(
    up: dict[str, list[str]], down: dict[str, list[str]],
    match: Callable[[str], bool], name: str, modality: str = "disease",
    min_recurrence: int = 2,
) -> Signature:
    """Aggregate all terms satisfying ``match`` into ONE consensus signature.

    Genes are scored by signed recurrence across the matching replicate terms
    (up count minus down count), keeping only genes seen at least
    ``min_recurrence`` times — a robust real disease/perturbation signature that
    down-weights single-study noise. ✅ deterministic; the consensus is the
    honest way to use CREEDS replicates instead of one arbitrary sample.
    """
    up_c: Counter = Counter()
    down_c: Counter = Counter()
    n = 0
    for term in up:
        if not match(term):
            continue
        n += 1
        up_c.update(up.get(term, []))
        down_c.update(down.get(term, []))
    scores: dict[str, float] = {}
    for g in set(up_c) | set(down_c):
        net = up_c.get(g, 0) - down_c.get(g, 0)
        if abs(net) >= min_recurrence:
            scores[g] = float(net)
    return Signature(scores, name=name, modality=modality,
                     meta={"n_replicates": n, "n_genes": len(scores)})
