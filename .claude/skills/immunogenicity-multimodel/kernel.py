"""Kernel sidecar for the immunogenicity-multimodel skill.

Verified backends (this session):
  - BioLib DTU/ImmunoGeNN  : MHC-II population immunogenicity (pIRS) + deimmunize
  - IEDB NetMHCpan (el)    : MHC-I binding %rank  [cloud REST]
  - standalone NetMHCpan   : MHC-I/II binding via mhctools [if binaries on PATH]
Every helper defers third-party imports to its body (starter kernel has no biolib).
"""

UNIFIED_COLS = ["id", "source", "model", "allele", "length", "peptide",
                "position", "affinity_nM", "percentile_rank", "score", "is_binder"]
RANK_STRONG_I = 0.5
RANK_WEAK_I = 2.0
RANK_STRONG_II = 2.0
RANK_WEAK_II = 10.0
IMMUNOGENN_THRESHOLD = 83.0   # pIRS_rank >= 83 -> immunogenic (per ImmunoGeNN paper)
IEDB_MHCI_URL = "http://tools-cluster-interface.iedb.org/tools_api/mhci/"


def read_fasta(path):
    """Parse a FASTA file -> {name: sequence}."""
    seqs, name, buf = {}, None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        elif line:
            buf.append(line)
    if name:
        seqs[name] = "".join(buf)
    return seqs


def classify_binder(rank, mhc2=False):
    """%rank -> strong/weak/none using DTU/IEDB thresholds."""
    import pandas as pd
    strong = RANK_STRONG_II if mhc2 else RANK_STRONG_I
    weak = RANK_WEAK_II if mhc2 else RANK_WEAK_I
    if rank is None or pd.isna(rank):
        return "none"
    if rank <= strong:
        return "strong"
    if rank <= weak:
        return "weak"
    return "none"


def run_immunogenn(fasta_path, mode="screen", timeout=1200):
    """BioLib DTU/ImmunoGeNN. mode='screen' -> pIRS per 15mer;
    mode='deimmunize' -> single-point variant scan of the first sequence.
    Returns (output_dir, unified_df_or_None). Anonymous run, no token."""
    import os, tempfile, biolib, pandas as pd
    os.environ.setdefault("BIOLIB_CACHE_DIR", "/tmp/biolib_cache")
    os.makedirs(os.environ["BIOLIB_CACHE_DIR"], exist_ok=True)
    app = biolib.load("DTU/ImmunoGeNN", suppress_version_warning=True)
    args = f"--fasta_file {os.path.basename(fasta_path)}"
    if mode == "deimmunize":
        args += " --mode deimmunize"
    td = tempfile.mkdtemp()
    job = app.cli(args=args, files=[fasta_path], blocking=True, timeout=timeout)
    if job.get_status() != "completed":
        return td, None
    job.save_files(td)
    if mode == "deimmunize":
        return td, None  # caller reads SAVs.csv / scores.csv / deimmunized_variants.fasta
    pirs = pd.read_csv(os.path.join(td, "pIRS.csv"))
    rows = []
    for _, r in pirs.iterrows():
        rank = r["pIRS_rank"]
        rows.append({"id": r["id"], "source": "biolib", "model": "immunogenn",
                     "allele": "DRB1(global)", "length": len(str(r["peptide_seq"])),
                     "peptide": r["peptide_seq"], "position": int(r["peptide_pos"]) - 1,
                     "affinity_nM": None, "percentile_rank": rank, "score": r["pIRS"],
                     "is_binder": "immunogenic" if rank >= IMMUNOGENN_THRESHOLD else "none"})
    return td, pd.DataFrame(rows)


def run_iedb_mhci(sequences, alleles, lengths=(9,), method="netmhcpan_el", timeout=180):
    """IEDB cloud NetMHCpan (verified method: netmhcpan_el). Returns unified DataFrame."""
    import io, requests, pandas as pd
    pairs = [(a, l) for a in alleles for l in lengths]
    allele_str = ",".join(a for a, _ in pairs)
    length_str = ",".join(str(l) for _, l in pairs)
    frames = []
    for name, seq in sequences.items():
        data = {"method": method, "sequence_text": seq,
                "allele": allele_str, "length": length_str}
        r = requests.post(IEDB_MHCI_URL, data=data, timeout=timeout)
        r.raise_for_status()
        df = pd.read_csv(io.StringIO(r.text), sep="\t")
        df["__id"] = name
        frames.append(df)
    raw = pd.concat(frames, ignore_index=True)
    rows = []
    for _, r in raw.iterrows():
        rank = float(r["percentile_rank"]) if pd.notna(r.get("percentile_rank")) else None
        rows.append({"id": r.get("__id"), "source": "iedb", "model": method,
                     "allele": r.get("allele"),
                     "length": int(r["length"]) if pd.notna(r.get("length")) else len(str(r.get("peptide", ""))),
                     "peptide": r.get("peptide"),
                     "position": int(r["start"]) - 1 if pd.notna(r.get("start")) else None,
                     "affinity_nM": None,
                     "percentile_rank": rank,
                     "score": float(r["score"]) if pd.notna(r.get("score")) else None,
                     "is_binder": classify_binder(rank, mhc2=False)})
    return pd.DataFrame(rows)


def run_multimodel(fasta_path, alleles_mhci=("HLA-A*02:01",), mhci_lengths=(9,),
                   use=("immunogenn", "iedb")):
    """One protein -> run selected backends -> stacked unified table."""
    import pandas as pd
    seqs = read_fasta(fasta_path)
    parts = []
    if "immunogenn" in use:
        _, df = run_immunogenn(fasta_path, mode="screen")
        if df is not None:
            parts.append(df)
    if "iedb" in use:
        parts.append(run_iedb_mhci(seqs, list(alleles_mhci), tuple(mhci_lengths)))
    merged = pd.concat(parts, ignore_index=True)
    return merged[UNIFIED_COLS]


def per_position_matrix(merged, seq_len):
    """Project peptide-level calls onto residues. MHC binding -> 100-rank;
    immunogenn -> pIRS_rank. Rows=models, intensity 0-100 (higher=stronger)."""
    import numpy as np, pandas as pd
    real = merged[~merged["is_binder"].astype(str).str.startswith("SKIPPED")]
    models = sorted(real["model"].dropna().unique())
    M = np.zeros((len(models), seq_len))
    for i, model in enumerate(models):
        sub = real[real["model"] == model]
        for _, r in sub.iterrows():
            if pd.isna(r["position"]) or pd.isna(r["percentile_rank"]):
                continue
            s, L = int(r["position"]), int(r["length"])
            sig = float(r["percentile_rank"]) if model == "immunogenn" else 100.0 - float(r["percentile_rank"])
            for p in range(s, min(s + L, seq_len)):
                M[i, p] = max(M[i, p], sig)
    return models, M
