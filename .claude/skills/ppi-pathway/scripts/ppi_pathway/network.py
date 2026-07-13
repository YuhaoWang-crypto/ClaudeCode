"""
network.py — load the bulk STRING / BioGRID files into a networkx graph.

Use this for genome-wide, offline work (the STRING API caps at ~a few thousand
identifiers per call). Call ``ppi_pathway.download`` first to cache the files.

STRING nodes are Ensembl protein ids like ``9606.ENSP00000269305``. We attach
the human-readable gene symbol from ``protein.info`` so you can query by symbol.
BioGRID rows already carry official gene symbols (columns 8/9 in tab3).
"""
from __future__ import annotations

import gzip
from pathlib import Path
from typing import Iterable

from . import HUMAN_TAXON, STRING_VERSION
from .download import data_dir


def _open(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


# --------------------------------------------------------------------------- #
# STRING
# --------------------------------------------------------------------------- #
def load_string_info(taxon: int = HUMAN_TAXON,
                     version: str = STRING_VERSION) -> dict[str, str]:
    """STRING protein id -> preferred gene symbol."""
    path = data_dir() / f"{taxon}.protein.info.v{version}.txt.gz"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} missing — run `python -m ppi_pathway.download --source string`")
    mapping: dict[str, str] = {}
    with _open(path) as fh:
        next(fh, None)  # header: #string_protein_id preferred_name protein_size annotation
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                mapping[parts[0]] = parts[1]
    return mapping


def load_string(taxon: int = HUMAN_TAXON, version: str = STRING_VERSION,
                min_score: int = 400, relabel_to_symbol: bool = True):
    """Return a networkx.Graph of STRING interactions with combined_score >=
    ``min_score`` (0-1000). Edge attribute ``score`` holds the combined score.

    ``min_score`` guidance: 150 low, 400 medium, 700 high, 900 highest.
    Filtering at load time keeps the human graph to a manageable size.
    """
    import networkx as nx

    links = data_dir() / f"{taxon}.protein.links.v{version}.txt.gz"
    if not links.exists():
        raise FileNotFoundError(
            f"{links} missing — run `python -m ppi_pathway.download --source string`")

    symbol = load_string_info(taxon, version) if relabel_to_symbol else {}
    g = nx.Graph()
    with _open(links) as fh:
        header = next(fh).split()
        # header: protein1 protein2 combined_score   (space-separated!)
        for line in fh:
            a, b, s = line.split()
            sc = int(s)
            if sc < min_score:
                continue
            na = symbol.get(a, a) if relabel_to_symbol else a
            nb = symbol.get(b, b) if relabel_to_symbol else b
            # keep the strongest edge if a symbol maps from two proteins
            if g.has_edge(na, nb):
                if sc > g[na][nb]["score"]:
                    g[na][nb]["score"] = sc
            else:
                g.add_edge(na, nb, score=sc)
    g.graph["source"] = f"STRING v{version}"
    g.graph["taxon"] = taxon
    g.graph["min_score"] = min_score
    return g


# --------------------------------------------------------------------------- #
# BioGRID (tab3)
# --------------------------------------------------------------------------- #
# tab3 column indices we care about (0-based). See the BioGRID tab3 spec.
_T3 = {
    "sym_a": 7,          # Official Symbol Interactor A
    "sym_b": 8,          # Official Symbol Interactor B
    "system": 11,        # Experimental System
    "system_type": 12,   # physical / genetic
    "throughput": 17,    # Low/High Throughput
    "organism_a": 35,    # Organism ID Interactor A
    "organism_b": 36,
}


def load_biogrid(path: str | Path | None = None, taxon: int = HUMAN_TAXON,
                 physical_only: bool = True, same_organism: bool = True):
    """Load a BioGRID tab3 file into a networkx.Graph keyed by gene symbol.

    Edge attributes: ``systems`` (set of experimental systems) and ``evidence``
    (count of supporting records). By default keeps only *physical* interactions
    within the same organism — the right choice for building a PPI network.
    """
    import networkx as nx

    if path is None:
        # find the extracted organism file in the cache
        hits = sorted(data_dir().glob("BIOGRID-ORGANISM-*.tab3.txt"))
        if not hits:
            raise FileNotFoundError(
                "No BioGRID tab3 file cached — run "
                "`python -m ppi_pathway.download --source biogrid`")
        path = hits[0]
    path = Path(path)

    g = nx.Graph()
    with _open(path) as fh:
        next(fh, None)  # header row
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) <= _T3["organism_b"]:
                continue
            if physical_only and c[_T3["system_type"]].lower() != "physical":
                continue
            if same_organism and (c[_T3["organism_a"]] != c[_T3["organism_b"]]):
                continue
            a, b = c[_T3["sym_a"]], c[_T3["sym_b"]]
            if not a or not b or a == "-" or b == "-" or a == b:
                continue
            sysname = c[_T3["system"]]
            if g.has_edge(a, b):
                g[a][b]["systems"].add(sysname)
                g[a][b]["evidence"] += 1
            else:
                g.add_edge(a, b, systems={sysname}, evidence=1)
    g.graph["source"] = f"BioGRID ({path.name})"
    g.graph["taxon"] = taxon
    return g


# --------------------------------------------------------------------------- #
# combine
# --------------------------------------------------------------------------- #
def consensus_graph(string_g, biogrid_g):
    """Overlay STRING and BioGRID; mark each edge with which source(s) support
    it (``in_string`` / ``in_biogrid``). Edges present in both are the most
    trustworthy for downstream analysis."""
    import networkx as nx

    g = nx.Graph()
    for u, v, d in string_g.edges(data=True):
        g.add_edge(u, v, score=d.get("score"), in_string=True, in_biogrid=False)
    for u, v, d in biogrid_g.edges(data=True):
        if g.has_edge(u, v):
            g[u][v]["in_biogrid"] = True
            g[u][v]["evidence"] = d.get("evidence")
        else:
            g.add_edge(u, v, in_string=False, in_biogrid=True,
                       evidence=d.get("evidence"))
    both = sum(1 for _, _, d in g.edges(data=True)
               if d.get("in_string") and d.get("in_biogrid"))
    g.graph["source"] = "STRING ∩∪ BioGRID"
    g.graph["edges_in_both"] = both
    return g
