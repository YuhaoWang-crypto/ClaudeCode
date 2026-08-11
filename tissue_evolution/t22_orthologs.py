"""
T22 - orthologue panel for the candidate proteins, projected onto human coordinates.

Everything downstream (contact-restricted coevolution, the subunit-swap
compatibility matrix) needs the same thing: for a given human residue position,
what amino acid does each other species have there? A pairwise human-vs-X
alignment answers that directly, so instead of building a real multiple
alignment I fetch Ensembl Compara's pairwise protein alignments and project
each one onto the ungapped human sequence.

This throws away insertions relative to human. That is the right trade here -
every position I will use is defined by a human PDB structure, so a column
that does not exist in human is a column I cannot place on a contact map
anyway. It also means all genes end up on a common coordinate system for free.

Only one2one orthologues are kept. one2many and many2many would silently mix
paralogues into what is supposed to be a single ortholog column, and for this
gene set (SUCLA2 vs SUCLG2, IDH1 vs IDH2, ACO1 vs ACO2) that is exactly the
error that would manufacture a coevolution signal out of nothing.
"""
import json
import os
import time
import urllib.error
import urllib.request

from .common import RAW, log

REST = "https://rest.ensembl.org"
CACHE = os.path.join(RAW, "ensembl_homology")

# The candidate set. The TCA catalytic core plus the four groups that were
# nominated for perturbation analysis, plus their structural partners - a
# subunit is useless here without the subunits it touches.
GENES = [
    # TCA catalytic core (T18 stratum 'catalytic')
    "CS", "ACO2", "IDH2", "IDH3A", "IDH3B", "IDH3G", "OGDH", "DLST", "DLD",
    "SUCLG1", "SUCLA2", "SUCLG2", "SDHA", "SDHB", "SDHC", "SDHD", "FH", "MDH2",
    # cytosolic paralogues (division-of-labour contrast)
    "ACO1", "IDH1", "MDH1", "GOT1", "GOT2", "OGDHL",
    # group 2: redox tolerance
    "SIRT3", "FXN", "NNT", "SOD2", "IDH3A",
    # group 3: SDH assembly + quinone pool
    "SDHAF1", "SDHAF2", "SDHAF3", "SDHAF4", "COQ2", "COQ7", "COQ9",
    # group 4: nuclear-encoded OXPHOS subunits that touch mtDNA-encoded ones
    "COX4I1", "COX5A", "COX5B", "COX6B1", "COX6C", "COX7A2", "COX7B", "COX8A",
    "UQCRC1", "UQCRC2", "UQCRFS1", "UQCRB", "UQCRQ", "CYC1",
    "ATP5F1A", "ATP5F1B", "ATP5PO", "ATP5MC1",
    "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS7", "NDUFS8", "NDUFV1", "NDUFA9",
    # group 4: the mtDNA-encoded partners themselves
    "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB", "MT-ATP6", "MT-ATP8",
    "MT-ND1", "MT-ND2", "MT-ND4", "MT-ND5",
]


def _get(url, tries=5):
    for k in range(tries):
        try:
            req = urllib.request.Request(
                url, headers={"Content-Type": "application/json",
                              "User-Agent": "tissue-evolution/1.0"})
            with urllib.request.urlopen(req, timeout=180) as fh:
                return json.load(fh)
        except urllib.error.HTTPError as exc:
            if exc.code in (400, 404):          # gene genuinely absent
                return None
            if exc.code == 429:                  # Compara rate limit
                time.sleep(float(exc.headers.get("Retry-After", 2)) + 1)
                continue
            if k == tries - 1:
                raise
        except Exception:                        # noqa: BLE001
            if k == tries - 1:
                raise
        time.sleep(2 ** k)
    return None


def fetch(symbol):
    os.makedirs(CACHE, exist_ok=True)
    p = os.path.join(CACHE, f"{symbol.replace('/', '_')}.json")
    if os.path.exists(p):
        with open(p) as fh:
            return json.load(fh)
    url = (f"{REST}/homology/symbol/homo_sapiens/{symbol}"
           "?content-type=application/json&type=orthologues"
           "&sequence=protein&aligned=1&format=full")
    d = _get(url)
    out = {}
    if d and d.get("data"):
        for h in d["data"][0].get("homologies", []):
            if h.get("type") != "ortholog_one2one":
                continue
            sp = h["target"]["species"]
            out[sp] = dict(src=h["source"]["align_seq"],
                           tgt=h["target"]["align_seq"],
                           pid=h["target"].get("perc_id"))
    with open(p, "w") as fh:
        json.dump(out, fh)
    return out


def project(pairs):
    """Pairwise alignments -> {species: string indexed by human residue}."""
    if not pairs:
        return {}, ""
    human = None
    cols = {}
    for sp, h in pairs.items():
        src, tgt = h["src"], h["tgt"]
        if len(src) != len(tgt):
            continue
        keep = [(s, t) for s, t in zip(src, tgt) if s != "-"]
        seq = "".join(s for s, _ in keep)
        if human is None:
            human = seq
        elif seq != human:
            # different Compara alignments can trim the human sequence
            # differently; only keep the ones on the majority coordinate frame
            continue
        cols[sp] = "".join(t for _, t in keep)
    if human is not None:
        cols["homo_sapiens"] = human
    return cols, human or ""


def build(genes=None, verbose=True):
    genes = genes or GENES
    out = {}
    for i, g in enumerate(dict.fromkeys(genes)):
        try:
            cols, human = project(fetch(g))
        except Exception as exc:                 # noqa: BLE001
            log(f"  {g}: {exc}")
            continue
        if human:
            out[g] = dict(human=human, cols=cols)
        if verbose and (i + 1) % 10 == 0:
            log(f"  {i + 1}/{len(genes)} genes, last {g}: "
                f"{len(cols)} species x {len(human)} aa")
    return out


def common_species(panel, min_frac=0.85):
    """Species present in at least min_frac of the genes."""
    from collections import Counter
    c = Counter()
    for d in panel.values():
        c.update(d["cols"].keys())      # count species, not sequence strings
    n = len(panel)
    return sorted(s for s, k in c.items() if k >= min_frac * n)


if __name__ == "__main__":
    panel = build()
    sp = common_species(panel)
    log(f"panel: {len(panel)} genes, {len(sp)} species shared by >=85% of genes")
    with open(os.path.join(CACHE, "_panel_species.json"), "w") as fh:
        json.dump(sp, fh, indent=1)
