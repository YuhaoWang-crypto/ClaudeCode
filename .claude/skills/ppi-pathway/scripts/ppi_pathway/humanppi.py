"""
humanppi.py — Cong-lab deep-learning human PPI predictions (Zhang, Humphreys,
Pei, ... Baker, Cong, *Science* 2025), from http://prodata.swmed.edu/humanPPI/ .

This is a structure-based (RoseTTAFold2-PPI / AlphaFold) screen of the human
proteome: ~29k confident predicted interactions with per-pair confidence scores
and rich metadata (subcellular locality, disease, GO process/function, PDB
templates). We download the *final predictions* table (not the multi-GB neural-
network weights — the predictions are the reusable product) and load it into a
networkx graph, so it can join STRING and BioGRID as a third, high-confidence,
structurally-supported evidence layer.

Two precision tiers are published:
    * final_predictions_90.tsv — expected precision 90% (~17.8k pairs)
    * final_predictions_80.tsv — expected precision 80% (~29.3k pairs, superset)

Key columns: Protein1/2 (UniProt AC), Name1/2 (gene symbol), RFprob (RF2-PPI
probability), AFprob/AFprob5/AFMprob (AlphaFold confidences), STRING (STRING
score), Source, Locality/Disease/Process/Function annotations.
"""
from __future__ import annotations

import io
import ssl
import tarfile
import urllib.request
from pathlib import Path

from .download import data_dir

HOST = "conglab.swmed.edu"
FINAL_PREDICTIONS_URL = f"https://{HOST}/humanPPI/downloads/final_predictions.tar.gz"
# Many .edu servers omit this intermediate from their TLS chain, which breaks
# strict verification. We fetch it once (plain HTTP, no chain needed) and add it
# to the trust store — verification stays ON, we just complete the chain.
INCOMMON_INTERMEDIATE_URL = "http://crt.usertrust.com/InCommonRSAServerCA2.crt"


def _ssl_context() -> ssl.SSLContext:
    """Default trust (system + proxy CA bundle) plus the InCommon intermediate
    that ``conglab.swmed.edu`` fails to send. Never disables verification."""
    ctx = ssl.create_default_context()
    try:
        cache = data_dir() / "InCommonRSAServerCA2.pem"
        if not cache.exists():
            with urllib.request.urlopen(INCOMMON_INTERMEDIATE_URL, timeout=60) as r:
                der = r.read()
            pem = ssl.DER_cert_to_PEM_cert(der)
            cache.write_text(pem)
        ctx.load_verify_locations(cafile=str(cache))
    except Exception as e:  # if it fails, the default context still applies
        print(f"  [warn] could not add InCommon intermediate: {e}")
    return ctx


def download(precision: int = 80, force: bool = False) -> Path:
    """Download + extract the final-predictions TSV for the given precision tier
    (80 or 90). Returns the path to the extracted ``final_predictions_<p>.tsv``."""
    if precision not in (80, 90):
        raise ValueError("precision must be 80 or 90")
    d = data_dir()
    out = d / f"final_predictions_{precision}.tsv"
    if out.exists() and not force:
        print(f"  [cached] {out.name}")
        return out

    print(f"  [get   ] {FINAL_PREDICTIONS_URL}")
    req = urllib.request.Request(FINAL_PREDICTIONS_URL,
                                 headers={"User-Agent": "ppi_pathway-skill/0.1"})
    with urllib.request.urlopen(req, timeout=180, context=_ssl_context()) as r:
        blob = r.read()
    with tarfile.open(fileobj=io.BytesIO(blob), mode="r:gz") as tf:
        for p in (80, 90):
            member = next((m for m in tf.getmembers()
                           if m.name.endswith(f"final_predictions_{p}.tsv")), None)
            if member is None:
                continue
            dest = d / f"final_predictions_{p}.tsv"
            with tf.extractfile(member) as src:
                dest.write_bytes(src.read())
            print(f"  [unpack] {dest.name}")
    return out


def load(precision: int = 80, min_rfprob: float = 0.0,
         path: str | Path | None = None):
    """Load humanPPI predictions into a networkx.Graph keyed by gene symbol.

    Edge attributes: ``rfprob`` (RF2-PPI probability), ``afprob`` (AlphaFold),
    ``string`` (STRING score if any), ``source`` (D/S/P evidence flags).
    ``min_rfprob`` filters to the most confident structural predictions.
    """
    import networkx as nx

    if path is None:
        path = data_dir() / f"final_predictions_{precision}.tsv"
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"{path} missing — run `python -m ppi_pathway.download --source humanppi`")

    g = nx.Graph()
    with open(path) as fh:
        header = None
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if header is None:
                header = cols
                idx = {name: i for i, name in enumerate(header)}
                continue
            a, b = cols[idx["Name1"]], cols[idx["Name2"]]
            if not a or not b or a == b or a == "none" or b == "none":
                continue
            try:
                rf = float(cols[idx["RFprob"]])
            except (ValueError, KeyError):
                rf = 0.0
            if rf < min_rfprob:
                continue

            def _f(key):
                try:
                    return float(cols[idx[key]])
                except Exception:
                    return None

            g.add_edge(a, b, rfprob=rf, afprob=_f("AFprob"),
                       string=_f("STRING"), source=cols[idx.get("Source", -1)]
                       if idx.get("Source") is not None else None,
                       up1=cols[idx["Protein1"]], up2=cols[idx["Protein2"]])
    g.graph["source"] = f"humanPPI (Cong lab, precision {precision}%)"
    g.graph["taxon"] = 9606
    return g
