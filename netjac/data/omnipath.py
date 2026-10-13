"""
Prior-knowledge signed networks (all open access):

  OmniPath interactions  https://omnipathdb.org/interactions  (datasets
      omnipath + kinaseextra + pathwayextra; directed, with a literature
      consensus sign: stimulation / inhibition)
  DoRothEA (levels A,B)  TF -> target, signed by mode of regulation
  OmniPath complexes     CORUM, ComplexPortal, hu.MAP, PDB, ...
  SIGNOR                 https://signor.uniroma2.it/getData.php?organism=9606
                         (used for the mechanism of each causal edge)

Cache: figures/_data_cache/omnipath/ (git-ignored).
"""
import os
import subprocess

import numpy as np

from netjac import DATA_CACHE

DIR = os.path.join(DATA_CACHE, "omnipath")
URLS = {
    "omnipath_signed.tsv": "https://omnipathdb.org/interactions?genesymbols=1"
                           "&datasets=omnipath,kinaseextra,pathwayextra&fields=sources",
    "dorothea_AB.tsv": "https://omnipathdb.org/interactions?genesymbols=1"
                       "&datasets=dorothea&dorothea_levels=A,B&fields=sources",
    "complexes.tsv": "https://omnipathdb.org/complexes",
    "signor_human.tsv": "https://signor.uniroma2.it/getData.php?organism=9606",
}


def _file(name):
    fn = os.path.join(DIR, name)
    if not os.path.exists(fn) or os.path.getsize(fn) < 10000:
        os.makedirs(DIR, exist_ok=True)
        subprocess.run(["curl", "-sS", "--max-time", "600", "-o", fn, URLS[name]],
                       check=True)
    return fn


def _read(name):
    with open(_file(name)) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return head, [ln.rstrip("\n").split("\t") for ln in fh]


def signed_edges(which="omnipath"):
    """Directed signed edges as a list of (source, target, sign +1/-1).
    Edges whose consensus is both or neither stimulation/inhibition are
    dropped; self-loops are dropped."""
    fn = {"omnipath": "omnipath_signed.tsv", "dorothea": "dorothea_AB.tsv"}[which]
    head, rows = _read(fn)
    h = {k: i for i, k in enumerate(head)}
    out = set()
    for r in rows:
        if len(r) < len(head) or r[h["is_directed"]] != "True":
            continue
        st = r[h["consensus_stimulation"]] == "True"
        inh = r[h["consensus_inhibition"]] == "True"
        if st == inh:
            continue
        a, b = r[h["source_genesymbol"]], r[h["target_genesymbol"]]
        if a == b or "_" in a or "_" in b:     # skip complexes as nodes
            continue
        out.add((a, b, 1 if st else -1))
    # a pair annotated with both signs is ambiguous -> drop
    sign = {}
    for a, b, s in out:
        sign.setdefault((a, b), set()).add(s)
    return [(a, b, next(iter(s))) for (a, b), s in sign.items() if len(s) == 1]


def complexes(min_size=2, max_size=30):
    head, rows = _read("complexes.tsv")
    h = {k: i for i, k in enumerate(head)}
    seen, out = set(), []
    for r in rows:
        if len(r) <= h["components_genesymbols"]:
            continue
        g = tuple(sorted(set(r[h["components_genesymbols"]].split("_"))))
        if min_size <= len(g) <= max_size and g not in seen:
            seen.add(g)
            out.append(g)
    return out


def signor_mechanisms():
    """{(source, target): mechanism} for protein->protein SIGNOR edges."""
    _, rows = _read("signor_human.tsv")
    out = {}
    for r in rows:
        if len(r) > 9 and r[1] == "protein" and r[5] == "protein":
            out[(r[0], r[4])] = r[9] or "unknown"
    return out


def feedback_pairs(edges):
    """2-node loops A->B, B->A from a signed edge list.
    Returns list of (A, B, s_AB, s_BA, kind) with kind 'positive' (s_AB*s_BA>0)
    or 'negative'."""
    d = {(a, b): s for a, b, s in edges}
    out = []
    for (a, b), s in d.items():
        if a < b and (b, a) in d:
            t = d[(b, a)]
            out.append((a, b, s, t, "positive" if s * t > 0 else "negative"))
    return out


def complex_pairs(cplx, genes):
    """Unique within-complex gene pairs restricted to `genes`."""
    gs = set(genes)
    out = set()
    for c in cplx:
        m = [g for g in c if g in gs]
        for i in range(len(m)):
            for j in range(i + 1, len(m)):
                out.add(tuple(sorted((m[i], m[j]))))
    return sorted(out)


def index_pairs(pairs, genes):
    ix = {g: i for i, g in enumerate(genes)}
    keep = [(ix[a], ix[b]) for a, b, *_ in pairs if a in ix and b in ix]
    return np.array(keep, dtype=int).reshape(-1, 2)
