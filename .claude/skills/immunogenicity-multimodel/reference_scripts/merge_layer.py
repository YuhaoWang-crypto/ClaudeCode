#!/usr/bin/env python3
"""
Merge layer for the DTU multi-model immunogenicity orchestrator.

One protein in -> runs every available backend in parallel -> one unified
long-format table (one row per peptide x model x allele) + a per-position
summary + a combined heatmap.

Backends (auto-detected; a backend that can't run is skipped with a note,
never a silent pass):
  - BioLib  : DTU/ImmunoGeNN         -> MHC-II population immunogenicity (pIRS)
  - IEDB    : NetMHCpan (netmhcpan_el) -> MHC-I binding %rank        [cloud]
  - standalone: NetMHCpan/IIpan       -> MHC-I/II binding            [if on PATH]

Unified per-peptide schema:
  id, source, model, allele, length, peptide, position,
  affinity_nM, percentile_rank, score, is_binder
"""
import os, json, tempfile
import pandas as pd

os.environ.setdefault("BIOLIB_CACHE_DIR", "/tmp/biolib_cache")

UNIFIED_COLS = ["id","source","model","allele","length","peptide","position",
                "affinity_nM","percentile_rank","score","is_binder"]


def _read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name: seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        elif line:
            buf.append(line)
    if name: seqs[name] = "".join(buf)
    return seqs


# ---- backend: ImmunoGeNN via BioLib --------------------------------------
def backend_immunogenn(fasta_path):
    """MHC-II population immunogenicity (pIRS), per 15mer peptide."""
    import biolib
    app = biolib.load("DTU/ImmunoGeNN", suppress_version_warning=True)
    with tempfile.TemporaryDirectory() as td:
        job = app.cli(args=f"--fasta_file {os.path.basename(fasta_path)}",
                      files=[fasta_path], blocking=True, timeout=1200)
        if job.get_status() != "completed":
            return _skip_row("immunogenn", "biolib", "job did not complete")
        job.save_files(td)
        pirs = pd.read_csv(os.path.join(td, "pIRS.csv"))
    rows = []
    for _, r in pirs.iterrows():
        rank = r["pIRS_rank"]                      # 0-100 percentile; >=83 = immunogenic
        rows.append({
            "id": r["id"], "source": "biolib", "model": "immunogenn",
            "allele": "DRB1(global)", "length": len(str(r["peptide_seq"])),
            "peptide": r["peptide_seq"], "position": int(r["peptide_pos"]) - 1,
            "affinity_nM": None, "percentile_rank": rank, "score": r["pIRS"],
            "is_binder": "immunogenic" if rank >= 83 else "none",
        })
    return pd.DataFrame(rows)


# ---- backend: NetMHCpan via IEDB -----------------------------------------
def backend_iedb_mhci(sequences, alleles, lengths=(9,)):
    import mhc_iedb_backend as IB
    try:
        return IB.run_iedb_mhci(sequences, alleles=alleles, lengths=lengths,
                                method="netmhcpan_el")
    except Exception as e:
        return _skip_row("netmhcpan_el", "iedb", str(e)[:120])


# ---- backend: standalone NetMHCpan ---------------------------------------
def backend_standalone_mhci(sequences, alleles, lengths=(8,9,10,11)):
    import mhc_standalone_backend as SB
    if SB.which("netmhcpan") is None:
        return _skip_row("netmhcpan", "standalone", "binary not on PATH")
    return SB.run_netmhcpan(sequences, alleles=alleles, lengths=lengths)


def _skip_row(model, source, why):
    return pd.DataFrame([{**{c: None for c in UNIFIED_COLS},
                          "model": model, "source": source,
                          "is_binder": f"SKIPPED:{why}"}])


# ---- orchestration --------------------------------------------------------
def run_merged(fasta_path, alleles_mhci=("HLA-A*02:01",),
               mhci_lengths=(9,), use=("immunogenn","iedb")):
    """Run selected backends and stack their unified outputs."""
    seqs = _read_fasta(fasta_path)
    parts = []
    if "immunogenn" in use:
        parts.append(backend_immunogenn(fasta_path))
    if "iedb" in use:
        parts.append(backend_iedb_mhci(seqs, list(alleles_mhci), mhci_lengths))
    if "standalone" in use:
        parts.append(backend_standalone_mhci(seqs, list(alleles_mhci)))
    merged = pd.concat([p for p in parts if p is not None], ignore_index=True)
    return merged[UNIFIED_COLS] if set(UNIFIED_COLS).issubset(merged.columns) else merged


def per_position_summary(merged, seq_len, seq_id):
    """Collapse peptide-level calls onto residue positions for a combined view.
    A residue inherits the strongest signal from any peptide window covering it."""
    real = merged[~merged["is_binder"].astype(str).str.startswith("SKIPPED")].copy()
    models = sorted(real["model"].unique())
    mat = pd.DataFrame(0.0, index=range(seq_len), columns=models)
    for _, r in real.iterrows():
        if pd.isna(r["position"]) or pd.isna(r["percentile_rank"]):
            continue
        start = int(r["position"]); L = int(r["length"])
        # signal: for MHC binding, lower %rank = stronger -> use (100-rank) as intensity;
        # for immunogenn pIRS_rank, higher = more immunogenic -> use rank directly
        if r["model"] == "immunogenn":
            sig = float(r["percentile_rank"])
        else:
            sig = 100.0 - float(r["percentile_rank"])
        for pos in range(start, min(start + L, seq_len)):
            mat.at[pos, r["model"]] = max(mat.at[pos, r["model"]], sig)
    mat.index.name = "position"
    return mat


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--alleles", nargs="+", default=["HLA-A*02:01"])
    ap.add_argument("--use", nargs="+", default=["immunogenn","iedb"])
    ap.add_argument("--out", default="merged_results.csv")
    a = ap.parse_args()
    m = run_merged(a.fasta, alleles_mhci=a.alleles, use=a.use)
    m.to_csv(a.out, index=False)
    print(m.groupby(["source","model"]).size().to_string())
    print("\nsaved:", a.out)
