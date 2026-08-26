"""Data acquisition layer.

Every external resource is fetched once and cached under ``data/cache``. All
downstream modules read from the cache, so a run is reproducible offline once
the cache is warm (``--offline`` refuses to hit the network at all).

Sources
-------
UniProtKB REST      reviewed insulin entries; mature A/B chains are taken from
                    the entry's own ``Peptide`` features, never hardcoded.
IPD-MHC (EBI)       ``MHC_prot.fasta`` -- the curated non-human MHC protein
                    release. Contains the canine DLA class II alleles.
NCBI E-utilities    feline MHC class II beta chains. IPD-MHC has *no* feline
                    entries, so the cat panel has to be assembled from GenBank.
IPD-IMGT/HLA        human DRB/DQA1/DQB1 protein alignments, used as the
                    reference numbering frame and as part of the NetMHCIIpan
                    training-space proxy.
"""

from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

PKG_ROOT = Path(__file__).resolve().parent.parent
CACHE = PKG_ROOT / "data" / "cache"

UNIPROT_INSULIN = {
    "human": "P01308",
    "dog": "P01321",
    "cat": "P06306",
    "bovine": "P01317",
    "pig": "P01315",
}

IPD_MHC_URL = "https://ftp.ebi.ac.uk/pub/databases/ipd/mhc/MHC_prot.fasta"
IMGT_HLA_BASE = "https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Mouse H-2 class II chains (I-A/I-E, several haplotypes). NetMHCIIpan-4.3's
# training data covers mouse, so these belong in the training-space proxy.
MOUSE_MHC2_BETA = ["P01921", "P06342", "P14483", "P04230", "P06345"]
MOUSE_MHC2_ALPHA = ["P01910"]

_OFFLINE = False


class OfflineError(RuntimeError):
    """Raised when a resource is missing from the cache in offline mode."""


def set_offline(flag: bool) -> None:
    global _OFFLINE
    _OFFLINE = flag


def _cached(name: str, fetch, binary: bool = False):
    """Return the cached bytes/str for ``name``, fetching once if absent."""
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / name
    if path.exists() and path.stat().st_size > 0:
        return path.read_bytes() if binary else path.read_text()
    if _OFFLINE:
        raise OfflineError(f"{name} not in cache and --offline was requested")
    payload = fetch()
    if binary:
        path.write_bytes(payload)
    else:
        path.write_text(payload)
    return payload


def _get(url: str, retries: int = 4) -> bytes:
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "vetimmuno/1.0"})
            with urllib.request.urlopen(req, timeout=300) as fh:
                return fh.read()
        except Exception as exc:  # network hiccup -> exponential backoff
            last = exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"failed to fetch {url}: {last}")


def _post(url: str, data: Dict[str, str], retries: int = 4) -> bytes:
    last = None
    body = urllib.parse.urlencode(data).encode()
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, data=body, headers={"User-Agent": "vetimmuno/1.0"})
            with urllib.request.urlopen(req, timeout=300) as fh:
                return fh.read()
        except Exception as exc:
            last = exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"failed to POST {url}: {last}")


# --------------------------------------------------------------------------
# FASTA
# --------------------------------------------------------------------------

def parse_fasta(text: str) -> List[Tuple[str, str]]:
    """Parse FASTA into ``(header, sequence)`` pairs, preserving order."""
    records: List[Tuple[str, List[str]]] = []
    for line in text.splitlines():
        if line.startswith(">"):
            records.append((line[1:].strip(), []))
        elif records and line.strip():
            records[-1][1].append(line.strip())
    return [(h, "".join(s)) for h, s in records]


def write_fasta(path: Path, records: Iterable[Tuple[str, str]], width: int = 60) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


# --------------------------------------------------------------------------
# Insulin
# --------------------------------------------------------------------------

def uniprot_entry(accession: str) -> dict:
    def fetch() -> str:
        url = (
            f"https://rest.uniprot.org/uniprotkb/{accession}.json"
            "?fields=accession,id,protein_name,organism_name,ft_peptide,ft_chain,sequence"
        )
        return _get(url).decode()

    return json.loads(_cached(f"uniprot_{accession}.json", fetch))


def insulin_chains(species: str) -> Dict[str, object]:
    """Mature insulin A and B chains for ``species``, straight from UniProt features.

    The chain boundaries come from the entry's ``Peptide`` features so that the
    result is traceable to the database record rather than to this file.
    """
    accession = UNIPROT_INSULIN[species]
    entry = uniprot_entry(accession)
    seq = entry["sequence"]["value"]
    chains: Dict[str, str] = {}
    for feat in entry.get("features", []):
        desc = feat.get("description") or ""
        start = feat["location"]["start"]["value"]
        end = feat["location"]["end"]["value"]
        if "A chain" in desc:
            chains["A"] = seq[start - 1:end]
        elif "B chain" in desc:
            chains["B"] = seq[start - 1:end]
    missing = {"A", "B"} - set(chains)
    if missing:
        raise RuntimeError(f"{accession}: no mature-chain feature for {sorted(missing)}")
    return {
        "species": species,
        "accession": accession,
        "entry_name": entry.get("uniProtkbId"),
        "organism": entry["organism"]["scientificName"],
        "A": chains["A"],
        "B": chains["B"],
    }


def all_insulins() -> Dict[str, Dict[str, object]]:
    return {sp: insulin_chains(sp) for sp in UNIPROT_INSULIN}


# --------------------------------------------------------------------------
# IPD-MHC (canine DLA + bovine BoLA)
# --------------------------------------------------------------------------

def ipd_mhc_records() -> List[Tuple[str, str]]:
    text = _cached("MHC_prot.fasta", lambda: _get(IPD_MHC_URL).decode())
    return parse_fasta(text)


def ipd_alleles(prefixes: Iterable[str]) -> Dict[str, str]:
    """``{allele_name: protein_sequence}`` for every allele matching a prefix.

    IPD-MHC headers look like ``IPD-MHC:DLA04814 DLA-DQA1*012:01:2 81 bp``;
    field 2 is the official allele name.
    """
    prefixes = tuple(prefixes)
    out: Dict[str, str] = {}
    for header, seq in ipd_mhc_records():
        parts = header.split()
        if len(parts) < 2:
            continue
        name = parts[1]
        if name.startswith(prefixes):
            out[name] = seq.replace("*", "").replace("X", "")
    return out


# --------------------------------------------------------------------------
# NCBI (feline FLA)
# --------------------------------------------------------------------------

FELINE_QUERY = (
    '(Felis catus[Organism]) AND '
    '(FLA-DRB[All Fields] OR DRB[Title] OR "class II antigen beta"[Title] OR '
    '"DR beta"[Title] OR "DRB beta"[Title] OR "histocompatibility antigen"[Title])'
)


def ncbi_feline_class2() -> List[Tuple[str, str]]:
    """Feline MHC class II protein records from GenBank/RefSeq.

    IPD-MHC ships no feline sequences, so this is the only route to an FLA
    panel -- and the reason the cat panel is weaker than the dog one.
    """

    def fetch() -> str:
        url = (
            f"{EUTILS}/esearch.fcgi?db=protein&retmax=500&retmode=json&term="
            + urllib.parse.quote(FELINE_QUERY)
        )
        ids = json.loads(_get(url).decode())["esearchresult"]["idlist"]
        if not ids:
            raise RuntimeError("NCBI returned no feline class II records")
        chunks = []
        for i in range(0, len(ids), 200):
            payload = {
                "db": "protein",
                "rettype": "fasta",
                "retmode": "text",
                "id": ",".join(ids[i:i + 200]),
            }
            chunks.append(_post(f"{EUTILS}/efetch.fcgi", payload).decode())
            time.sleep(0.4)  # stay under the E-utilities rate limit
        return "".join(chunks)

    return parse_fasta(_cached("feline_class2.fasta", fetch))


# --------------------------------------------------------------------------
# Human / mouse reference space
# --------------------------------------------------------------------------

def imgt_hla(locus: str) -> Dict[str, str]:
    """Human class II protein sequences from IPD-IMGT/HLA (``DRB``/``DQA1``/``DQB1``)."""
    fname = f"{locus}_prot.fasta"
    text = _cached(f"imgt_{fname}", lambda: _get(f"{IMGT_HLA_BASE}/{fname}").decode())
    out: Dict[str, str] = {}
    for header, seq in parse_fasta(text):
        parts = header.split()
        name = parts[1] if len(parts) > 1 else parts[0]
        out[name] = seq.replace("*", "").replace("X", "")
    return out


def mouse_mhc2(chain: str = "beta") -> Dict[str, str]:
    accs = MOUSE_MHC2_BETA if chain == "beta" else MOUSE_MHC2_ALPHA
    out: Dict[str, str] = {}
    for acc in accs:
        entry = uniprot_entry(acc)
        out[f"H-2:{acc}"] = entry["sequence"]["value"]
    return out
