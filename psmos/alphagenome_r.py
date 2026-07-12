"""
alphagenome_r.py — R layer (regulatory grammar) via AlphaGenome, human/mouse only.

The R layer of the PSMOS framework — do the promoter/enhancer, TF-motif grammar,
and cell-type expression of a pathway gene look the same across species? — is the
one Evo2 explicitly *cannot* supply (Evo2 scores sequence naturalness, not
regulatory output). AlphaGenome predicts, from DNA sequence, the regulatory
tracks (RNA-seq expression, ATAC/DNase accessibility, TF ChIP, splicing) — but
**only for human and mouse**. So this adapter is, by construction, a human↔mouse
regulatory-conservation probe, and is labelled as such everywhere.

What this module does:
  1. LIVE now: resolve the human + mouse genomic interval of each pathway gate
     gene from Ensembl (real; runnable in this environment).
  2. GATED on ALPHAGENOME_API_KEY: call AlphaGenome for a TSS-centred window in
     each species, pull the predicted tracks, and define the R score as the
     human↔mouse concordance of the predicted regulatory profile.

AlphaGenome is an API-only model (weights not released), so — unlike Evo2 — it
cannot be self-hosted on Modal; it needs a DeepMind API key (free for
non-commercial). Without the key this writes the resolved intervals and a clear
"needs key" status, exactly the way Evo2 was "adapter-ready" before Modal.

    pip install alphagenome
    export ALPHAGENOME_API_KEY=...      # from https://deepmind.google.com/science/alphagenome
    python3 -m psmos.alphagenome_r Hippo
"""

import json
import os
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass, asdict

ENSEMBL = "https://rest.ensembl.org"
CACHE_DIR = os.path.join(os.path.dirname(__file__), "data", "alphagenome")

# AlphaGenome supports human and mouse only.
AG_SPECIES = {"human": "homo_sapiens", "mouse": "mus_musculus"}
# TSS-centred window; AlphaGenome accepts 2KB/16KB/100KB/500KB/1MB.
WINDOW = 100_000
# Representative ontology term for the readout (context-dependent; documented).
DEFAULT_ONTOLOGY = ["UBERON:0002107"]  # liver — Hippo/YAP-active regenerative tissue


@dataclass
class GeneInterval:
    family: str
    species_key: str
    gene: str
    chrom: str
    tss: int
    start: int
    end: int
    strand: int
    assembly: str
    found: bool
    note: str = ""


def _get(url, retries=4):
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=40) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            if e.code in (400, 404):
                return None
            last = e; time.sleep(2 ** i)
        except Exception as e:
            last = e; time.sleep(2 ** i)
    raise RuntimeError(f"Ensembl lookup failed: {last}")


def resolve_interval(family, species_key, symbols):
    esp = AG_SPECIES.get(species_key)
    if not esp:
        return GeneInterval(family, species_key, symbols[0] if symbols else "",
                            "", 0, 0, 0, 0, "", False, "not human/mouse (AlphaGenome scope)")
    for sym in symbols:
        d = _get(f"{ENSEMBL}/lookup/symbol/{esp}/{urllib.parse.quote(sym)}?content-type=application/json")
        if d and "seq_region_name" in d:
            strand = d.get("strand", 1)
            tss = d["start"] if strand == 1 else d["end"]
            half = WINDOW // 2
            return GeneInterval(
                family, species_key, sym, f"chr{d['seq_region_name']}",
                tss, max(1, tss - half), tss + half, strand,
                d.get("assembly_name", ""), True)
    return GeneInterval(family, species_key, symbols[0] if symbols else "",
                        "", 0, 0, 0, 0, "", False, "no Ensembl gene")


def resolve_pathway_intervals(pathway):
    os.makedirs(CACHE_DIR, exist_ok=True)
    out = []
    for fam_key, per_species in pathway.ortholog_seed.items():
        for sk in ("human", "mouse"):
            if sk in per_species:
                gi = resolve_interval(fam_key, sk, per_species[sk])
                out.append(gi)
                print(f"  [{fam_key:12s}] {sk:6s} {gi.gene:8s} -> "
                      f"{gi.chrom}:{gi.start}-{gi.end} ({gi.assembly}) "
                      f"{'ok' if gi.found else gi.note}")
    json.dump([asdict(g) for g in out],
              open(os.path.join(CACHE_DIR, f"{pathway.key}_intervals.json"), "w"), indent=1)
    return out


def _predict_tracks(model, gi: GeneInterval):
    """Call AlphaGenome for one interval; return a compact track summary."""
    from alphagenome.models import dna_client
    from alphagenome.data import genome

    interval = genome.Interval(chromosome=gi.chrom, start=gi.start, end=gi.end)
    out = model.predict_interval(
        interval=interval,
        requested_outputs=[
            dna_client.OutputType.RNA_SEQ,
            dna_client.OutputType.ATAC,
            dna_client.OutputType.CHIP_TF,
            dna_client.OutputType.SPLICE_SITES,
        ],
        ontology_terms=DEFAULT_ONTOLOGY,
    )
    import numpy as np
    def prof(track):
        try:
            return np.asarray(track.values).mean(axis=1)  # position profile
        except Exception:
            return None
    return {
        "rna_seq": prof(out.rna_seq),
        "atac": prof(out.atac),
        "chip_tf": prof(getattr(out, "chip_tf", None)),
    }


def compute_R(pathway):
    """Compute the R-layer human↔mouse regulatory concordance per family, if the
    AlphaGenome key is present. Returns (status, {family: R_score 0..1})."""
    key = os.environ.get("ALPHAGENOME_API_KEY")
    intervals = resolve_pathway_intervals(pathway)
    by = {}
    for gi in intervals:
        by.setdefault(gi.family, {})[gi.species_key] = gi

    if not key:
        return ("needs-key: set ALPHAGENOME_API_KEY (adapter ready; intervals resolved)", {})

    from alphagenome.models import dna_client
    import numpy as np
    model = dna_client.create(key)
    R = {}
    for fam, sp in by.items():
        if not (sp.get("human", GeneInterval(*[None]*10, False)).found and
                sp.get("mouse", GeneInterval(*[None]*10, False)).found):
            continue
        th = _predict_tracks(model, sp["human"])
        tm = _predict_tracks(model, sp["mouse"])
        # R = mean cross-species correlation of the predicted regulatory profiles
        # (TSS-anchored; mouse profile flipped for +/- strand alignment).
        cors = []
        for k in ("rna_seq", "atac", "chip_tf"):
            a, b = th.get(k), tm.get(k)
            if a is None or b is None:
                continue
            if sp["mouse"].strand != sp["human"].strand:
                b = b[::-1]
            n = min(len(a), len(b))
            if n > 8:
                c = np.corrcoef(a[:n], b[:n])[0, 1]
                if np.isfinite(c):
                    cors.append((c + 1) / 2)  # map [-1,1] -> [0,1]
        if cors:
            R[fam] = round(float(np.mean(cors)), 3)
    json.dump(R, open(os.path.join(CACHE_DIR, f"{pathway.key}_R.json"), "w"), indent=1)
    return ("live", R)


if __name__ == "__main__":
    import sys
    os.environ.setdefault("SSL_CERT_FILE", "/root/.ccr/ca-bundle.crt")
    from psmos.pathways import get_pathway
    pk = sys.argv[1] if len(sys.argv) > 1 else "Hippo"
    print(f"Resolving human/mouse regulatory intervals for {pk}:")
    status, R = compute_R(get_pathway(pk))
    print(f"\nAlphaGenome R-layer status: {status}")
    if R:
        for fam, r in R.items():
            print(f"  {fam:12s} R(human↔mouse) = {r}")
