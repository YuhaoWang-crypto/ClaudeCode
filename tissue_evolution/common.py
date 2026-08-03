"""
Shared paths, species ladder and small helpers for the tissue-evolution pipeline.

Everything downloaded lives under data/raw (never committed); everything the
pipeline computes lands in data/derived (committed when small) and
figures/tissue_evolution.
"""
import os
import sys
import time
import gzip
import urllib.request
import urllib.parse

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")
DERIVED = os.path.join(ROOT, "data", "derived")
FIGDIR = os.path.join(ROOT, "figures", "tissue_evolution")
for _d in (RAW, DERIVED, FIGDIR):
    os.makedirs(_d, exist_ok=True)


# --------------------------------------------------------------------------
# Species ladder.
#
# `divergence_mya` are TimeTree v5 median estimates for the human-<species>
# split; they are used only to ORDER the ladder and to convert dN into a rate
# per unit time, never as a calibrated molecular clock.
# --------------------------------------------------------------------------
SPECIES = [
    dict(key="macaque", biomart="mmulatta", ftp="macaca_mulatta",
         latin="Macaca mulatta", divergence_mya=28.8),
    dict(key="mouse", biomart="mmusculus", ftp="mus_musculus",
         latin="Mus musculus", divergence_mya=87.0),
    dict(key="opossum", biomart="mdomestica", ftp="monodelphis_domestica",
         latin="Monodelphis domestica", divergence_mya=159.0),
    dict(key="chicken", biomart="ggallus", ftp="gallus_gallus",
         latin="Gallus gallus", divergence_mya=319.0),
]
SPECIES_BY_KEY = {s["key"]: s for s in SPECIES}

HUMAN_CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]

# www.ensembl.org's BioMart returns HTTP 500 / an HTML status page under load
# often enough to kill a multi-hour job, so queries rotate over the mirrors.
BIOMART_HOSTS = ["https://www.ensembl.org", "https://useast.ensembl.org",
                 "https://asia.ensembl.org"]
ENSEMBL_FTP = "https://ftp.ensembl.org/pub/current_fasta"


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def http_get(url, timeout=600, retries=5, backoff=2.0):
    """GET with exponential backoff. Returns decoded text."""
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as fh:
                return fh.read().decode("utf-8", "replace")
        except Exception as exc:  # noqa: BLE001 - network layer, retry anything
            last = exc
            wait = backoff ** attempt
            log(f"  GET failed ({exc}); retry {attempt + 1}/{retries} in {wait:.0f}s")
            time.sleep(wait)
    raise RuntimeError(f"GET {url[:120]} failed after {retries} tries: {last}")


def download(url, dest, timeout=1800, retries=5, backoff=2.0):
    """Stream a (possibly large) file to disk, skipping if already present."""
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return dest
    tmp = dest + ".part"
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as fh, open(tmp, "wb") as out:
                while True:
                    chunk = fh.read(1 << 20)
                    if not chunk:
                        break
                    out.write(chunk)
            os.replace(tmp, dest)
            return dest
        except Exception as exc:  # noqa: BLE001
            last = exc
            wait = backoff ** attempt
            log(f"  download failed ({exc}); retry {attempt + 1}/{retries} in {wait:.0f}s")
            time.sleep(wait)
    raise RuntimeError(f"download {url[:120]} failed: {last}")


def biomart_query(dataset, attributes, filters=(), timeout=900, retries=12):
    """Run a BioMart TSV query, rotating over the Ensembl mirrors.

    BioMart signals overload two different ways - an HTTP 500, or a 200 whose
    body is an HTML status page - and signals a malformed query with a one-line
    `Query ERROR:`. All three are detected here rather than being parsed as
    data; only the last is fatal, since retrying a malformed query is pointless.
    """
    attr_xml = "".join(f'<Attribute name="{a}"/>' for a in attributes)
    filt_xml = "".join(f'<Filter name="{n}" value="{v}"/>' for n, v in filters)
    xml = (
        '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>'
        '<Query virtualSchemaName="default" formatter="TSV" header="1" '
        'uniqueRows="1" count="0" datasetConfigVersion="0.6">'
        f'<Dataset name="{dataset}" interface="default">{filt_xml}{attr_xml}</Dataset>'
        "</Query>"
    )
    qs = "/biomart/martservice?query=" + urllib.parse.quote(xml)
    for attempt in range(retries):
        host = BIOMART_HOSTS[attempt % len(BIOMART_HOSTS)]
        try:
            txt = http_get(host + qs, timeout=timeout, retries=1)
        except Exception as exc:                     # noqa: BLE001
            txt = f"<transport error {exc}>"
        if txt.startswith("Query ERROR"):
            raise RuntimeError(txt.strip().splitlines()[0])
        if not txt.lstrip().startswith("<"):
            return txt
        wait = min(60.0, 2.0 ** (attempt // len(BIOMART_HOSTS) + 1))
        log(f"  BioMart {host.split('//')[1].split('.')[0]} unavailable; "
            f"rotating, retry in {wait:.0f}s")
        time.sleep(wait)
    raise RuntimeError("BioMart unavailable on all mirrors after retries")


def to_tsv_gz(df, path, **kwargs):
    """Write a gzipped TSV deterministically.

    pandas stamps the current time into the gzip header, so re-running a module
    rewrites every .gz output with identical contents but a different checksum,
    and git reports a diff on files that did not change. Pinning mtime=0 makes
    a re-run a no-op in version control.
    """
    kwargs.setdefault("sep", "\t")
    return df.to_csv(path, compression={"method": "gzip", "mtime": 0}, **kwargs)


def read_fasta(path):
    """Yield (header, sequence) from a plain or gzipped FASTA."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        header, chunks = None, []
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:].strip(), []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, "".join(chunks)


def apply_figure_style():
    """Consistent, colour-blind-safe defaults for every figure in this package."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "figure.dpi": 130,
        "savefig.dpi": 160,
        "savefig.bbox": "tight",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linewidth": 0.6,
        "legend.frameon": False,
    })
    return plt


# Okabe-Ito, colour-blind safe.
PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7",
           "#E69F00", "#56B4E9", "#F0E442", "#666666"]


def die(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    raise SystemExit(1)
