"""
download.py — fetch the bulk STRING and BioGRID data files once and cache them.

Design goals
------------
* Idempotent: a file that is already fully downloaded (size matches the server's
  Content-Length) is skipped, so re-running is cheap.
* Ephemeral-container friendly: everything lands under a single cache dir
  (default ``./ppi_data``, overridable via the ``PPI_DATA_DIR`` env var) that is
  git-ignored. Nothing large is ever committed to the repo.
* No API key needed for STRING (public downloads). BioGRID's bulk *files* are
  also public; only the BioGRID *webservice* needs a key (not used here).

Usage
-----
    python -m ppi_pathway.download --source string          # human STRING
    python -m ppi_pathway.download --source biogrid          # human BioGRID
    python -m ppi_pathway.download --source all --taxon 9606

    from ppi_pathway import download
    paths = download.string_files(taxon=9606)
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from urllib.request import urlopen, Request

from . import STRING_VERSION, HUMAN_TAXON

STRING_BASE = "https://stringdb-downloads.org/download"
# BioGRID publishes per-organism tab3 files inside a single "by organism" zip
# for the current release. The tab3 format carries Entrez + gene symbols and
# rich metadata (experimental system, throughput, publication).
BIOGRID_ORGANISM_ZIP = (
    "https://downloads.thebiogrid.org/Download/BioGRID/"
    "Latest-Release/BIOGRID-ORGANISM-LATEST.tab3.zip"
)
# Common taxon -> BioGRID organism-file name fragment (used to pick one member
# out of the per-organism zip). Human is the default.
BIOGRID_ORGANISM_NAME = {
    9606: "Homo_sapiens",
    10090: "Mus_musculus",
    10116: "Rattus_norvegicus",
    7227: "Drosophila_melanogaster",
    6239: "Caenorhabditis_elegans",
    4932: "Saccharomyces_cerevisiae",
}


def data_dir() -> Path:
    d = Path(os.environ.get("PPI_DATA_DIR", "ppi_data"))
    d.mkdir(parents=True, exist_ok=True)
    return d


def _remote_size(url: str) -> int | None:
    try:
        req = Request(url, method="HEAD")
        with urlopen(req, timeout=60) as r:
            cl = r.headers.get("Content-Length")
            return int(cl) if cl else None
    except Exception:
        return None


def _download(url: str, dest: Path, force: bool = False) -> Path:
    """Stream ``url`` to ``dest`` unless a complete copy already exists."""
    remote = _remote_size(url)
    if dest.exists() and not force:
        if remote is None or dest.stat().st_size == remote:
            print(f"  [cached] {dest.name} ({dest.stat().st_size:,} bytes)")
            return dest
        print(f"  [stale ] {dest.name}: re-downloading")

    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"  [get   ] {url}")
    req = Request(url, headers={"User-Agent": "ppi_pathway-skill/0.1"})
    with urlopen(req, timeout=120) as r, open(tmp, "wb") as fh:
        total, chunk = 0, 1 << 20
        while True:
            buf = r.read(chunk)
            if not buf:
                break
            fh.write(buf)
            total += len(buf)
            if remote:
                pct = 100 * total / remote
                print(f"\r  [get   ] {dest.name}: {total:,}/{remote:,} "
                      f"({pct:4.1f}%)", end="", flush=True)
    print()
    tmp.replace(dest)
    return dest


def string_files(taxon: int = HUMAN_TAXON, version: str = STRING_VERSION,
                 force: bool = False) -> dict[str, Path]:
    """Download the three STRING files we need for offline analysis.

    Returns a dict with keys ``links`` (scored PPIs), ``info`` (preferred gene
    names + descriptions) and ``aliases`` (external id -> STRING id mapping).
    """
    d = data_dir()
    files = {
        "links": f"protein.links.v{version}/{taxon}.protein.links.v{version}.txt.gz",
        "info": f"protein.info.v{version}/{taxon}.protein.info.v{version}.txt.gz",
        "aliases": f"protein.aliases.v{version}/{taxon}.protein.aliases.v{version}.txt.gz",
    }
    out: dict[str, Path] = {}
    print(f"STRING v{version} for taxon {taxon} -> {d}")
    for key, rel in files.items():
        url = f"{STRING_BASE}/{rel}"
        dest = d / Path(rel).name
        out[key] = _download(url, dest, force=force)
    return out


def biogrid_file(taxon: int = HUMAN_TAXON, force: bool = False) -> Path:
    """Download the current BioGRID per-organism zip and extract the tab3 file
    for ``taxon``. Returns the path to the extracted ``.tab3.txt``."""
    import zipfile

    d = data_dir()
    org = BIOGRID_ORGANISM_NAME.get(taxon)
    if org is None:
        raise ValueError(
            f"No BioGRID organism-name mapping for taxon {taxon}; add it to "
            f"BIOGRID_ORGANISM_NAME in download.py")

    zip_path = d / "BIOGRID-ORGANISM-LATEST.tab3.zip"
    print(f"BioGRID (organism zip) -> {d}")
    _download(BIOGRID_ORGANISM_ZIP, zip_path, force=force)

    with zipfile.ZipFile(zip_path) as zf:
        member = next((n for n in zf.namelist() if org in n and n.endswith(".txt")), None)
        if member is None:
            raise RuntimeError(
                f"Could not find a '{org}' member in {zip_path.name}. "
                f"Members: {zf.namelist()[:5]}...")
        out = d / member
        if not out.exists() or force:
            print(f"  [unzip ] {member}")
            with zf.open(member) as src, open(out, "wb") as dst:
                dst.write(src.read())
        else:
            print(f"  [cached] {member}")
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Cache STRING/BioGRID bulk data.")
    ap.add_argument("--source", choices=["string", "biogrid", "all"],
                    default="all")
    ap.add_argument("--taxon", type=int, default=HUMAN_TAXON,
                    help="NCBI taxonomy id (9606 = human, the default)")
    ap.add_argument("--force", action="store_true",
                    help="re-download even if a cached copy exists")
    args = ap.parse_args(argv)

    if args.source in ("string", "all"):
        string_files(taxon=args.taxon, force=args.force)
    if args.source in ("biogrid", "all"):
        biogrid_file(taxon=args.taxon, force=args.force)
    print("done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
