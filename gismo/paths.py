"""Where GISMO raw inputs and reference data live, and how to fetch them.

Raw inputs (~150 MB Zenodo archive, ~100 MB gnomAD tables) are NOT committed.
They are cached under $GISMO_DATA (default: <repo>/data/gismo, git-ignored).
The small derived score table IS committed at gismo/data/gismo_scores.tsv.gz,
so lookups work offline without re-downloading anything.
"""
import os
import pathlib
import subprocess
import zipfile

REPO = pathlib.Path(__file__).resolve().parent.parent
PKG_DATA = pathlib.Path(__file__).resolve().parent / "data"
SCORES_TSV = PKG_DATA / "gismo_scores.tsv.gz"
FIG_DIR = REPO / "figures"

DATA = pathlib.Path(os.environ.get("GISMO_DATA", REPO / "data" / "gismo"))
ZENODO_ZIP = DATA / "GISMO-v1.0.0.zip"
ZENODO_DIR = DATA / "zenodo"
REF = DATA / "ref"

# Zenodo record 20547722 (concept) -> 20547723 (v1.0.0). Code + inputs, CC-BY 4.0.
ZENODO_URL = ("https://zenodo.org/api/records/20547723/files/"
              "cliao5/GISMO-GISMO.zip/content")

GL = "https://raw.githubusercontent.com/macarthur-lab/gene_lists/master/lists/"
KO = "https://raw.githubusercontent.com/macarthur-lab/gnomad_lof/master/R/ko_gene_lists/"
REF_URLS = {
    # human constraint
    "gnomad_v2.1.1_lof.txt.bgz": "https://storage.googleapis.com/gcp-public-data--gnomad/"
        "release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz",
    "gnomad_v4.1_constraint.tsv": "https://storage.googleapis.com/gcp-public-data--gnomad/"
        "release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv",
    # dosage sensitivity (Collins 2022)
    "phaplo_ptriplo.tsv.gz": "https://zenodo.org/api/records/6347673/files/"
        "Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz/content",
    # AlphaMissense gene-level (CC BY-NC-SA: use locally, never commit)
    "am_gene.tsv.gz": "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_gene_hg38.tsv.gz",
    # gnomAD benchmark gene sets (same ones the paper used)
    "mouse_het_lethal.tsv": KO + "list_mouse_het_lethal_genes.tsv",
    "mouse_lethal.tsv": KO + "list_mouse_lethal_genes.tsv",
    "CEGv2.tsv": GL + "CEGv2_subset_universe.tsv",
    "NEGv1.tsv": GL + "NEGv1_subset_universe.tsv",
    "olfactory_receptors.tsv": GL + "olfactory_receptors.tsv",
    "all_ar.tsv": GL + "all_ar.tsv",   # 1,183 OMIM autosomal recessive (= paper's set)
    "all_ad.tsv": GL + "all_ad.tsv",   # 709 OMIM autosomal dominant (= paper's set)
    "mgi_essential.tsv": GL + "mgi_essential.tsv",
}


def _curl(url, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    subprocess.run(["curl", "-sS", "-L", "-f", "--retry", "4", "-m", "900",
                    "-o", str(tmp), url], check=True)
    tmp.rename(dest)


def zenodo_root():
    """Return the unpacked author repo root (download + unzip on first use)."""
    if not ZENODO_ZIP.exists():
        print(f"[gismo] downloading Zenodo archive -> {ZENODO_ZIP}")
        _curl(ZENODO_URL, ZENODO_ZIP)
    if not ZENODO_DIR.exists():
        with zipfile.ZipFile(ZENODO_ZIP) as z:
            z.extractall(ZENODO_DIR)
    roots = [p for p in ZENODO_DIR.iterdir() if p.is_dir()]
    assert len(roots) == 1, roots
    return roots[0]


def ref(name):
    """Return a local path to a reference file, downloading it if needed."""
    p = REF / name
    if not p.exists():
        print(f"[gismo] downloading {name}")
        _curl(REF_URLS[name], p)
    return p


def fetch_all():
    zenodo_root()
    for k in REF_URLS:
        ref(k)
    return DATA
