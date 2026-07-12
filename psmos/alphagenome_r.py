"""
alphagenome_r.py — R layer (regulatory grammar) via AlphaGenome, human/mouse only.

The R layer of the PSMOS framework — do the promoter/enhancer, TF-motif grammar
and cell-type expression of a pathway gene look the same across species? — is the
one Evo2 explicitly *cannot* supply (Evo2 scores sequence naturalness, not
regulatory output). AlphaGenome predicts, from DNA sequence, the regulatory
tracks (RNA-seq, ATAC accessibility, TF ChIP, splicing) — but **only for human
and mouse**. So this adapter is, by construction, a human↔mouse regulatory-
conservation probe, and is labelled as such everywhere.

Pipeline (all live):
  1. Resolve each gate gene's genomic interval from Ensembl. Human is GRCh38
     (AlphaGenome's human reference). Mouse is fetched as GRCm39 and lifted over
     to GRCm38/mm10 (AlphaGenome's mouse reference) via Ensembl /map.
  2. Predict a TSS-centred window in each species (organism=HOMO_SAPIENS /
     MUS_MUSCULUS), pull RNA_SEQ / ATAC / CHIP_TF tracks.
  3. R(gene) = human↔mouse concordance of the TSS-anchored, strand-oriented
     regulatory profiles (mean Pearson across tracks, mapped to 0..1). R is a
     *reference-comparison*: human is the reference (R=1); mouse's R is the
     computed concordance.

AlphaGenome is an API-only model (weights not released), so — unlike Evo2 — it
cannot be self-hosted on Modal; it needs a DeepMind API key (free for
non-commercial). Without the key this writes the resolved intervals and a clear
"needs key" status.

    pip install alphagenome
    export ALPHAGENOME_API_KEY=...   export SSL_CERT_FILE=/root/.ccr/ca-bundle.crt
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

AG_SPECIES = {"human": "homo_sapiens", "mouse": "mus_musculus"}
# AlphaGenome reference assemblies. Human interval coords must be GRCh38, mouse
# must be GRCm38 (mm10). Ensembl current is GRCh38 for human but GRCm39 for
# mouse, so mouse coords are lifted over.
AG_ASSEMBLY = {"human": "GRCh38", "mouse": "GRCm38"}
WINDOW = 16_384  # AlphaGenome supported length (16KB), TSS-centred
DEFAULT_ONTOLOGY = ["UBERON:0002107"]  # liver — a Hippo/YAP-active regenerative tissue


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
    raise RuntimeError(f"Ensembl request failed: {last}")


def _liftover_tss(species_key, chrom, tss):
    """Lift a mouse TSS from GRCm39 to GRCm38 (AlphaGenome's mouse reference)."""
    if species_key != "mouse":
        return tss, AG_ASSEMBLY[species_key]
    url = f"{ENSEMBL}/map/mouse/GRCm39/{chrom}:{tss}..{tss}/GRCm38?content-type=application/json"
    d = _get(url)
    maps = (d or {}).get("mappings", [])
    if maps:
        return maps[0]["mapped"]["start"], "GRCm38"
    return tss, "GRCm39(no-liftover)"


def resolve_interval(family, species_key, symbols):
    esp = AG_SPECIES.get(species_key)
    if not esp:
        return GeneInterval(family, species_key, symbols[0] if symbols else "",
                            "", 0, 0, 0, 0, "", False, "not human/mouse (AlphaGenome scope)")
    for sym in symbols:
        d = _get(f"{ENSEMBL}/lookup/symbol/{esp}/{urllib.parse.quote(sym)}?content-type=application/json")
        if d and "seq_region_name" in d:
            strand = d.get("strand", 1)
            chrom = str(d["seq_region_name"])
            tss = d["start"] if strand == 1 else d["end"]
            tss, assembly = _liftover_tss(species_key, chrom, tss)
            half = WINDOW // 2
            return GeneInterval(family, species_key, sym, f"chr{chrom}", tss,
                                max(1, tss - half), tss - half + WINDOW, strand, assembly, True)
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


def _profiles(model, gi, dna_client):
    """Predict tracks for one interval; return TSS-anchored, 5'->3' oriented
    per-position profiles for RNA_SEQ / ATAC / CHIP_TF."""
    import numpy as np
    from alphagenome.data import genome

    organism = (dna_client.Organism.HOMO_SAPIENS if gi.species_key == "human"
                else dna_client.Organism.MUS_MUSCULUS)
    iv = genome.Interval(chromosome=gi.chrom, start=gi.start, end=gi.end)
    out = model.predict_interval(
        interval=iv, organism=organism,
        requested_outputs=[dna_client.OutputType.RNA_SEQ,
                           dna_client.OutputType.ATAC,
                           dna_client.OutputType.CHIP_TF],
        ontology_terms=DEFAULT_ONTOLOGY)

    def prof(track):
        if track is None:
            return None
        a = np.asarray(track.values)
        v = a.mean(axis=1) if a.ndim == 2 else a           # collapse tracks -> position profile
        if gi.strand == -1:                                # orient all profiles 5'->3'
            v = v[::-1]
        return v

    return {"rna_seq": prof(out.rna_seq),
            "atac": prof(getattr(out, "atac", None)),
            "chip_tf": prof(getattr(out, "chip_tf", None))}


def compute_R(pathway):
    """Compute per-family human↔mouse regulatory concordance R (0..1) if the
    AlphaGenome key is present. Returns (status, {family: R})."""
    key = os.environ.get("ALPHAGENOME_API_KEY")
    intervals = resolve_pathway_intervals(pathway)
    by = {}
    for gi in intervals:
        by.setdefault(gi.family, {})[gi.species_key] = gi

    if not key:
        return ("needs-key: set ALPHAGENOME_API_KEY (adapter ready; intervals resolved)", {})

    import numpy as np
    from alphagenome.models import dna_client
    model = dna_client.create(key)

    R, detail = {}, {}
    for fam, sp in by.items():
        h, m = sp.get("human"), sp.get("mouse")
        if not (h and h.found and m and m.found):
            continue
        ph, pm = _profiles(model, h, dna_client), _profiles(model, m, dna_client)
        cors = {}
        for k in ("rna_seq", "atac", "chip_tf"):
            a, b = ph.get(k), pm.get(k)
            if a is None or b is None:
                continue
            n = min(len(a), len(b))
            if n > 8 and np.std(a[:n]) > 0 and np.std(b[:n]) > 0:
                c = float(np.corrcoef(a[:n], b[:n])[0, 1])
                if np.isfinite(c):
                    cors[k] = round((c + 1) / 2, 3)     # [-1,1] -> [0,1]
        if cors:
            R[fam] = round(float(np.mean(list(cors.values()))), 3)
            detail[fam] = cors
            print(f"  {fam:12s} R(human↔mouse) = {R[fam]}   tracks={cors}")
    json.dump({"R": R, "detail": detail},
              open(os.path.join(CACHE_DIR, f"{pathway.key}_R.json"), "w"), indent=1)
    return ("live", R)


if __name__ == "__main__":
    import sys
    os.environ.setdefault("SSL_CERT_FILE", "/root/.ccr/ca-bundle.crt")
    os.environ.setdefault("GRPC_DEFAULT_SSL_ROOTS_FILE_PATH", "/root/.ccr/ca-bundle.crt")
    from psmos.pathways import get_pathway
    pk = sys.argv[1] if len(sys.argv) > 1 else "Hippo"
    print(f"Resolving human/mouse regulatory intervals for {pk}:")
    status, R = compute_R(get_pathway(pk))
    print(f"\nAlphaGenome R-layer status: {status}")
