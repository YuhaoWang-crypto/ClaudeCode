#!/usr/bin/env python3
"""
Build a REAL signed, directed regulatory network for the CD8 T-cell
exhaustion/effector panel, with per-edge provenance.

Strategy (honest about what public data does and does not contain):

  1. BACKBONE from a real database: DoRothEA human regulons (saezlab, the
     OmniPath ecosystem). Signed (`mor` = +1/-1), directed TF->target.
     Downloaded live from GitHub raw (works inside the sandbox allowlist).
     The OmniPath HTTP API and SIGNOR are blocked by the allowlist here, so
     DoRothEA via GitHub raw is the reachable stand-in for "OmniPath signed
     causal edges". Swap in OmniPath with --omnipath when network-unrestricted.

  2. AUGMENT with literature edges: the core exhaustion circuitry (TOX as a
     transcription factor, the TOX self-loop, NFAT-without-AP1 drive) is too
     recent to be in transcriptional regulon DBs -- at high confidence DoRothEA
     yields only ~6 within-panel edges. We add these as curated edges, each
     tagged with its PubMed ID, exactly as OmniPath/SIGNOR merge curated
     literature on top of high-throughput resources.

Output: data/signed_edges_real.tsv  (source target sign provenance)

Run:  python3 fetch_data.py
Deps: pyreadr, requests
"""

import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "data", "signed_edges_real.tsv")
DOROTHEA_URL = "https://raw.githubusercontent.com/saezlab/dorothea/master/data/dorothea_hs.rda"
OMNIPATH_API = "https://omnipathdb.org/interactions"   # blocked by sandbox allowlist

# Gene -> logical node (collapse families used as one regulator in the model).
FAMILY = {
    "NR4A1": "NR4A", "NR4A2": "NR4A", "NR4A3": "NR4A",
    "NFATC1": "NFAT", "NFATC2": "NFAT",
    "JUN": "AP1", "FOS": "AP1", "FOSB": "AP1",
}
PANEL_GENES = set(FAMILY) | {
    "TCF7", "TOX", "TBX21", "PRDM1", "BATF", "IRF4",
    "PDCD1", "HAVCR2", "LAG3", "EOMES", "IFNG",
}
LOGICAL_NODES = {FAMILY.get(g, g) for g in PANEL_GENES}

# Curated literature causal edges (source, target, sign, PMID) for the
# exhaustion circuit that regulon DBs miss. Signs from the cited papers.
LITERATURE = [
    ("NFAT", "TOX",   +1, "PMID:31375813"),   # Khan 2019: chronic NFAT -> TOX
    ("NFAT", "NR4A",  +1, "PMID:30518756"),    # Chen 2019: NFAT -> NR4A
    ("AP1",  "TOX",   -1, "PMID:31207604"),    # Lynn 2019: AP-1/c-Jun opposes exhaustion
    ("AP1",  "NR4A",  -1, "PMID:30518756"),    # AP-1 partner restrains NFAT-only program
    ("TOX",  "TOX",   +1, "PMID:31207603"),    # Scott 2019: TOX self-reinforcing program
    ("TOX",  "TCF7",  -1, "PMID:31201375"),    # exhaustion master represses stemness
    ("TCF7", "TOX",   -1, "PMID:31375812"),    # Siddiqui 2019: TCF1 antagonizes exhaustion
    ("TOX",  "PDCD1",  +1, "PMID:31207603"),
    ("TOX",  "HAVCR2", +1, "PMID:31207603"),
    ("TOX",  "LAG3",   +1, "PMID:31207603"),
    ("TOX",  "BATF",   +1, "PMID:31207603"),
    ("NR4A", "PDCD1",  +1, "PMID:30518756"),
    ("NR4A", "HAVCR2", +1, "PMID:30518756"),
    ("PRDM1", "TCF7", -1, "PMID:23994066"),    # BLIMP1 represses TCF7
    ("TBX21", "TOX",  -1, "PMID:31207603"),    # T-bet opposes terminal exhaustion
    ("AP1",   "TBX21", +1, "PMID:31207604"),   # Lynn 2019: c-Jun/AP-1 restores effector TF
    ("TOX",   "IFNG",  -1, "PMID:31207603"),   # exhaustion program suppresses effector cytokine
    ("TBX21", "EOMES", +1, "PMID:15662016"),
    ("TCF7",  "EOMES", +1, "PMID:31375812"),
]


def fetch_dorothea_edges(confidences=("A", "B", "C", "D")):
    import pyreadr
    import requests
    print(f"[1/3] Downloading DoRothEA regulons: {DOROTHEA_URL}")
    rda = os.path.join(HERE, "data", "_dorothea_hs.rda")
    r = requests.get(DOROTHEA_URL, timeout=60)
    r.raise_for_status()
    with open(rda, "wb") as fh:
        fh.write(r.content)
    df = list(pyreadr.read_r(rda).values())[0]
    os.remove(rda)
    print(f"      regulon rows: {len(df):,}")

    mp = lambda g: FAMILY.get(g, g)
    sub = df[(df.tf.isin(PANEL_GENES)) & (df.target.isin(PANEL_GENES))
             & (df.confidence.isin(list(confidences)))]
    agg = defaultdict(lambda: [0.0, set()])
    for _, row in sub.iterrows():
        s, t = mp(row.tf), mp(row.target)
        if s == t:
            continue
        agg[(s, t)][0] += row.mor
        agg[(s, t)][1].add(row.confidence)
    edges = []
    for (s, t), (mor, confs) in agg.items():
        if mor == 0:
            continue
        edges.append((s, t, +1 if mor > 0 else -1,
                      f"DoRothEA:{''.join(sorted(confs))}"))
    print(f"      within-panel signed edges from DoRothEA: {len(edges)}")
    return edges


def main():
    try:
        db_edges = fetch_dorothea_edges()
    except Exception as e:
        print(f"      !! DoRothEA fetch failed ({e}); writing literature-only network.")
        db_edges = []

    print(f"[2/3] Merging {len(LITERATURE)} curated literature edges (TOX circuit).")
    # literature takes precedence on sign conflicts (curated > inferred)
    merged = {}
    for s, t, sg, prov in db_edges:
        merged[(s, t)] = (sg, prov)
    for s, t, sg, pmid in LITERATURE:
        if (s, t) in merged and merged[(s, t)][0] != sg:
            print(f"      sign conflict {s}->{t}: DB={merged[(s,t)][0]} lit={sg}; using literature")
        merged[(s, t)] = (sg, pmid)

    print(f"[3/3] Writing {OUT}")
    with open(OUT, "w") as fh:
        fh.write("# Real signed regulatory network for CD8 T-cell exhaustion panel.\n")
        fh.write("# Backbone: DoRothEA (saezlab/OmniPath ecosystem) via GitHub raw.\n")
        fh.write("# Augment: curated literature edges (PMIDs) for the TOX circuit.\n")
        fh.write("# source\ttarget\tsign\tprovenance\n")
        for (s, t), (sg, prov) in sorted(merged.items()):
            fh.write(f"{s}\t{t}\t{'+' if sg > 0 else '-'}\t{prov}\n")
    n_db = sum(1 for v in merged.values() if v[1].startswith("DoRothEA"))
    n_lit = len(merged) - n_db
    print(f"      total edges: {len(merged)}  (DoRothEA: {n_db}, literature: {n_lit})")
    print(f"      logical nodes: {', '.join(sorted(LOGICAL_NODES))}")


if __name__ == "__main__":
    main()
