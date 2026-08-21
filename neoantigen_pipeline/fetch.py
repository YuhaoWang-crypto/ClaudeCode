"""Public-data fetchers, all cached on disk.

Sources (all open access, no token):
  cBioPortal REST      somatic mutations + RNA expression for a TCGA sample
  UniProt REST         human reference proteome (canonical, reviewed)
  IEDB query-API       experimentally positive T-cell epitopes (TCR prior + benchmark)
  IEDB tools-cluster   NetMHCpan / NetMHCIIpan predictions (see presentation.py)
"""

from __future__ import annotations

import gzip
import json
import os
import time
from typing import Dict, List, Optional

import requests

from .config import CACHE_DIR

CBIO = "https://www.cbioportal.org/api"
UNIPROT = "https://rest.uniprot.org"
IEDB_QUERY = "https://query-api.iedb.org"

_UA = {"User-Agent": "neoantigen-selection-skill/1.0 (research demo)"}


def _cache_path(name: str) -> str:
    os.makedirs(CACHE_DIR, exist_ok=True)
    return os.path.join(CACHE_DIR, name)


def _cached_json(name: str, producer, force: bool = False):
    path = _cache_path(name)
    if os.path.exists(path) and not force:
        with open(path) as fh:
            return json.load(fh)
    obj = producer()
    with open(path, "w") as fh:
        json.dump(obj, fh)
    return obj


def _post(url: str, payload, tries: int = 4, timeout: int = 300):
    last = None
    for i in range(tries):
        try:
            r = requests.post(url, json=payload, headers=_UA, timeout=timeout)
            r.raise_for_status()
            return r.json()
        except Exception as exc:            # noqa: BLE001 - network retry
            last = exc
            time.sleep(2 ** i)
    raise RuntimeError(f"POST {url} failed: {last}")


def _get(url: str, tries: int = 4, timeout: int = 300, **kw):
    last = None
    for i in range(tries):
        try:
            r = requests.get(url, headers=_UA, timeout=timeout, **kw)
            r.raise_for_status()
            return r
        except Exception as exc:            # noqa: BLE001 - network retry
            last = exc
            time.sleep(2 ** i)
    raise RuntimeError(f"GET {url} failed: {last}")


# --------------------------------------------------------------------------
# cBioPortal
# --------------------------------------------------------------------------

def cbio_samples(study: str, sample_list: Optional[str] = None) -> List[str]:
    sl = sample_list or f"{study}_3way_complete"
    r = _get(f"{CBIO}/sample-lists/{sl}")
    return r.json()["sampleIds"]


def cbio_mutations(study: str, sample_id: str, force: bool = False) -> List[dict]:
    """DETAILED somatic mutation records for one sample (protein change, VAF counts)."""
    profile = f"{study}_mutations"

    def _go():
        return _post(
            f"{CBIO}/mutations/fetch?projection=DETAILED",
            {"sampleMolecularIdentifiers": [
                {"molecularProfileId": profile, "sampleId": sample_id}]},
        )

    return _cached_json(f"cbio_mut_{study}_{sample_id}.json", _go, force)


def cbio_expression(study: str, sample_id: str, entrez_ids: List[int],
                    profile_suffix: str = "rna_seq_v2_mrna",
                    force: bool = False) -> Dict[int, float]:
    """Per-gene tumor RNA value (RSEM by default) for one sample."""
    profile = f"{study}_{profile_suffix}"
    key = f"cbio_expr_{study}_{sample_id}_{profile_suffix}.json"

    def _go():
        out: Dict[str, float] = {}
        ids = sorted(set(int(i) for i in entrez_ids))
        for i in range(0, len(ids), 800):        # chunk to keep payloads sane
            chunk = ids[i:i + 800]
            data = _post(
                f"{CBIO}/molecular-profiles/{profile}/molecular-data/fetch",
                {"entrezGeneIds": chunk, "sampleIds": [sample_id]},
            )
            for rec in data:
                v = rec.get("value")
                if v is not None:
                    out[str(rec["entrezGeneId"])] = float(v)
        return out

    raw = _cached_json(key, _go, force)
    return {int(k): v for k, v in raw.items()}


def cbio_clinical(study: str, sample_id: str, force: bool = False) -> Dict[str, str]:
    def _go():
        r = _get(f"{CBIO}/studies/{study}/samples/{sample_id}/clinical-data")
        return {d["clinicalAttributeId"]: d["value"] for d in r.json()}

    return _cached_json(f"cbio_clin_{study}_{sample_id}.json", _go, force)


# --------------------------------------------------------------------------
# UniProt human reference proteome
# --------------------------------------------------------------------------

def human_proteome(force: bool = False) -> Dict[str, str]:
    """{gene_symbol: canonical protein sequence} for reviewed human entries."""
    path = _cache_path("uniprot_human_reviewed.fasta.gz")
    if not os.path.exists(path) or force:
        # The /stream endpoint truncates over long-haul proxies; the cursor-paged
        # /search endpoint is slower but finishes. 500 entries per page,
        # ~41 pages for the reviewed human proteome.
        url = (f"{UNIPROT}/uniprotkb/search?format=fasta&size=500&"
               "query=%28reviewed%3Atrue%29%20AND%20%28organism_id%3A9606%29")
        tmp, n_pages, n_entries = path + ".part", 0, 0
        with gzip.open(tmp, "wt") as out:
            while url:
                r = _get(url, tries=5, timeout=300)
                out.write(r.text)
                n_pages += 1
                n_entries += r.text.count(">")
                link = r.headers.get("link", "")
                url = None
                for part in link.split(","):
                    if 'rel="next"' in part:
                        url = part.split(";")[0].strip().strip("<>")
                        break
        if n_entries < 15000:
            os.remove(tmp)
            raise RuntimeError(f"UniProt returned only {n_entries} entries "
                               f"over {n_pages} pages; refusing a truncated proteome")
        os.replace(tmp, path)
    seqs: Dict[str, str] = {}
    name, buf = None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    seqs.setdefault(name, "".join(buf))
                gene = None
                for tok in line.split():
                    if tok.startswith("GN="):
                        gene = tok[3:]
                        break
                name = gene or line[1:].split("|")[1]
                buf = []
            elif line:
                buf.append(line)
    if name:
        seqs.setdefault(name, "".join(buf))
    return seqs


# NOTE: the self k-mer index lives in selfindex.SelfKmerIndex -- a plain Python
# set of 11.4M 9-mer strings costs ~1 GB, the encoded int64 array costs ~90 MB.


# --------------------------------------------------------------------------
# IEDB query-API: experimentally positive T-cell epitopes
# --------------------------------------------------------------------------

def _iedb_page(params: dict, timeout: int = 300):
    """PostgREST GET with proper encoding (the values contain spaces and parens)."""
    r = _get(f"{IEDB_QUERY}/tcell_search", timeout=timeout, params=params)
    return r.json()


def iedb_positive_epitopes(mhc_class: str = "I", lengths=(9, 10),
                           limit: int = 5000, force: bool = False) -> List[dict]:
    """Human-host, positive T-cell assay epitopes. Used as the TCR-recognition prior."""
    key = f"iedb_pos_{mhc_class}_{'_'.join(map(str, lengths))}_{limit}.json"

    def _go():
        rows: List[dict] = []
        for L in lengths:
            offset, page = 0, 1000
            while offset < limit:
                params = {
                    "qualitative_measure": "neq.Negative",
                    "host_organism_name": "eq.Homo sapiens (human)",
                    "mhc_class": f"eq.{mhc_class}",
                    "linear_sequence_length": f"eq.{L}",
                    "structure_type": "eq.Linear peptide",
                    "select": ("linear_sequence,mhc_allele_name,qualitative_measure,"
                               "disease_names,source_organism_name,"
                               "parent_source_antigen_name"),
                    # the API rejects offset without an explicit order
                    "order": "structure_id",
                    "limit": str(min(page, limit - offset)),
                    "offset": str(offset),
                }
                data = _iedb_page(params)
                if not data:
                    break
                rows.extend(data)
                offset += len(data)
                if len(data) < page:
                    break
        return rows

    return _cached_json(key, _go, force)


def iedb_human_tumor_epitopes(force: bool = False) -> List[dict]:
    """Positive human-source (self-antigen-derived) class-I epitopes -- the pool
    from which mutation-shaped, experimentally validated neoepitopes are mined."""

    def _go():
        rows, offset, page = [], 0, 1000
        while True:
            params = {
                "qualitative_measure": "neq.Negative",
                "host_organism_name": "eq.Homo sapiens (human)",
                "mhc_class": "eq.I",
                "source_organism_name": "eq.Homo sapiens (human)",
                "structure_type": "eq.Linear peptide",
                "select": ("linear_sequence,mhc_allele_name,qualitative_measure,"
                           "disease_names,parent_source_antigen_name,pubmed_id,"
                           "linear_sequence_length"),
                # the API rejects offset without an explicit order
                "order": "structure_id",
                "limit": str(page),
                "offset": str(offset),
            }
            data = _iedb_page(params)
            if not data:
                break
            rows.extend(data)
            offset += len(data)
            if len(data) < page:
                break
        return rows

    return _cached_json("iedb_human_source_classI.json", _go, force)
