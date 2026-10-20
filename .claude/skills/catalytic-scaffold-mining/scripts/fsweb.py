"""Stage 1 (remote): the Foldseek web server as a stand-in for a DALI sweep.

https://search.foldseek.com searches the whole AlphaFold Database server-side,
so a proteome-wide structural sweep needs no bulk download. The response
carries the target's full Calpha trace (`tCa`) and sequence (`tSeq`) inline,
which is enough to run the chemistry filter and the catalytic-shell geometry
check without ever fetching a structure file.

Two limitations to keep in view, both of which belong in any write-up:

* Hits are capped (afdb50 uses --max-accept 1000), and the cap fills with the
  query's closest relatives. Breadth across taxa therefore comes from running
  *several different queries* and *several databases*, not from one big search.
* The web results carry no TM-score and no per-residue pLDDT. TM-score is
  recomputed here from the coordinates; pLDDT is simply unavailable in remote
  mode and must be reported as missing rather than guessed.

HTTP goes through `curl` on purpose: it already has the session's proxy and CA
configuration, so there is nothing to re-derive in Python.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import time

import numpy as np

from structure import THREE_TO_ONE, Residue, Structure

API = "https://search.foldseek.com/api"
ONE_TO_THREE = {v: k for k, v in THREE_TO_ONE.items() if len(k) == 3}
ONE_TO_THREE.setdefault("X", "UNK")

# Databases worth combining for a cross-taxon sweep.
DEFAULT_DATABASES = ("afdb-proteome", "afdb-swissprot", "afdb50")


def _curl(args, timeout=600):
    if not shutil.which("curl"):
        raise RuntimeError("curl not found; it carries the proxy/CA config")
    out = subprocess.run(
        ["curl", "-sS", "--max-time", str(timeout)] + args,
        capture_output=True, text=True,
    )
    if out.returncode != 0:
        raise RuntimeError(f"curl failed: {out.stderr[:400]}")
    return out.stdout


def submit(pdb_path, databases=DEFAULT_DATABASES, mode="3diaa", timeout=300):
    args = [f"{API}/ticket", "-F", f"q=@{pdb_path}", "-F", f"mode={mode}"]
    for db in databases:
        args += ["-F", f"database[]={db}"]
    resp = json.loads(_curl(args, timeout=timeout))
    if "id" not in resp:
        raise RuntimeError(f"no ticket returned: {resp}")
    return resp["id"]


def wait(ticket, poll=15, limit=3600, verbose=True):
    waited = 0
    while waited < limit:
        st = json.loads(_curl([f"{API}/ticket/{ticket}"], timeout=120))
        status = st.get("status")
        if status == "COMPLETE":
            return True
        if status in ("ERROR", "RATELIMIT"):
            raise RuntimeError(f"foldseek job {ticket}: {status}")
        if verbose:
            print(f"    [{waited:4d}s] {status}", flush=True)
        time.sleep(poll)
        waited += poll
    raise TimeoutError(f"foldseek job {ticket} exceeded {limit}s")


def fetch(ticket, out_path, entry=0, timeout=900):
    _curl(["-o", str(out_path), f"{API}/result/{ticket}/{entry}"], timeout=timeout)
    with open(out_path) as fh:
        return json.load(fh)


def search(pdb_path, out_path, databases=DEFAULT_DATABASES, verbose=True):
    """Submit, wait, and download in one call."""
    ticket = submit(pdb_path, databases=databases)
    if verbose:
        print(f"  ticket {ticket} ({', '.join(databases)})", flush=True)
    wait(ticket, verbose=verbose)
    return fetch(ticket, out_path)


def iter_hits(raw):
    """Flatten the nested result structure into normalised hit dicts."""
    for block in raw.get("results", []):
        db = block.get("db")
        alignments = block.get("alignments") or []
        if alignments and isinstance(alignments[0], list):
            alignments = [a for grp in alignments for a in grp]
        for a in alignments:
            target = a.get("target", "")
            acc, desc = _parse_target(target)
            yield {
                "db": db,
                "target": target,
                "accession": acc,
                "description": desc,
                "tax_id": a.get("taxId"),
                "tax_name": a.get("taxName"),
                "seq_id": a.get("seqId"),
                "eval": a.get("eval"),
                "prob": a.get("prob"),
                "score": a.get("score"),
                "aln_length": a.get("alnLength"),
                "qstart": a.get("qStartPos"),
                "qend": a.get("qEndPos"),
                "tstart": a.get("dbStartPos"),
                "tend": a.get("dbEndPos"),
                "qlen": a.get("qLen"),
                "tlen": a.get("dbLen"),
                "qaln": a.get("qAln"),
                "taln": a.get("dbAln"),
                "tseq": a.get("tSeq"),
                "tca": a.get("tCa"),
            }


def _parse_target(target):
    """'AF-Q8N884-F1-model_v6 Cyclic GMP-AMP synthase' -> ('Q8N884', 'Cyclic...')"""
    head, _, desc = target.partition(" ")
    acc = head
    parts = head.split("-")
    if len(parts) >= 2 and parts[0] == "AF":
        acc = parts[1]
    return acc, desc.strip()


def target_structure(hit) -> Structure:
    """Rebuild a Calpha-only Structure from the inline tCa / tSeq fields.

    Residues are numbered 1..dbLen, which is the numbering Foldseek's
    tstart/tend use and -- for a single-fragment AlphaFold model -- is also the
    UniProt residue numbering. Multi-fragment entries (F2+, proteins over
    ~2700 aa) are offset; such hits are flagged rather than silently misnumbered.
    """
    seq = hit.get("tseq") or ""
    raw = hit.get("tca") or ""
    vals = [v for v in raw.split(",") if v]
    n = len(vals) // 3
    coords = np.array(vals[: n * 3], dtype=float).reshape(n, 3)
    residues = []
    for i in range(min(n, len(seq))):
        one = seq[i]
        res = Residue(
            resnum=i + 1,
            resname=ONE_TO_THREE.get(one, "UNK"),
            seq_index=i + 1,
        )
        res.atoms["CA"] = coords[i]
        residues.append(res)
    st = Structure(residues, name=hit.get("accession", "target"))
    st.fragment = _fragment_of(hit.get("target", ""))
    return st


def _fragment_of(target):
    for part in target.split("-"):
        if part.startswith("F") and part[1:].isdigit():
            return int(part[1:])
    return 1


def tm_score(query, target, mapping, norm_len=None):
    """TM-score over the aligned pairs, normalised by query length by default.

    Recomputed locally because the web API does not return it. Standard
    Zhang/Skolnick formulation; superposition is the one that optimises RMSD
    over the aligned set, which slightly under-estimates a true TM-align
    search over rotations -- fine for ranking, not a substitute for TM-align
    when a published number is needed.
    """
    pairs = [(query.by_index(qi), target.by_index(ti))
             for qi, ti in mapping.items()]
    pairs = [(q, t) for q, t in pairs
             if q is not None and t is not None
             and q.ca is not None and t.ca is not None]
    if len(pairs) < 3:
        return float("nan")
    L = norm_len or len(query)
    if L <= 21:
        d0 = 0.5
    else:
        d0 = 1.24 * (L - 15) ** (1.0 / 3.0) - 1.8
    d0 = max(d0, 0.5)
    P = np.array([t.ca for _, t in pairs])
    Q = np.array([q.ca for q, _ in pairs])
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, _, Wt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(V @ Wt))
    rot = V @ np.diag([1.0, 1.0, d]) @ Wt
    di = np.linalg.norm((Pc @ rot) - Qc, axis=1)
    return float(np.sum(1.0 / (1.0 + (di / d0) ** 2)) / L)
