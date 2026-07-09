#!/usr/bin/env python3
"""
Standalone MHC backend for the DTU multi-model orchestrator.

Wraps the license-gated DTU command-line tools (NetMHCpan / NetMHCIIpan /
NetChop / NetMHCstabpan) via the `mhctools` library, and normalizes their
output to the SAME schema the BioLib backend emits, so both plug into one
merged results table.

  columns: id, source, model, allele, length, peptide, position,
           affinity_nM, percentile_rank, score, is_binder

WHY standalone (not BioLib): NetMHCpan-family predictors are NOT published on
BioLib. They are downloaded from
https://services.healthtech.dtu.dk/ under an academic license (accept the
license; the download link expires ~4 h). Install once, put the binaries on
PATH (or point env vars at them), and this backend drives them.

Binaries expected on PATH (override via env):
  NETMHCPAN_BIN     (default: netMHCpan)      -> MHC-I binding
  NETMHCIIPAN_BIN   (default: netMHCIIpan)    -> MHC-II binding
  NETCHOP_BIN       (default: netchop)        -> proteasomal cleavage
  NETMHCSTABPAN_BIN (default: netMHCstabpan)  -> pMHC-I stability

If a binary is absent, that model is skipped with an explicit note in the
manifest — never a silent pass.
"""
import os, shutil
import pandas as pd

# ---- binding thresholds (DTU/IEDB conventions) ----------------------------
RANK_STRONG_I  = 0.5    # %rank <= 0.5  -> strong binder (MHC-I)
RANK_WEAK_I    = 2.0    # %rank <= 2.0  -> weak binder   (MHC-I)
RANK_STRONG_II = 2.0    # MHC-II strong
RANK_WEAK_II   = 10.0   # MHC-II weak

BINARIES = {
    "netmhcpan":     os.environ.get("NETMHCPAN_BIN", "netMHCpan"),
    "netmhciipan":   os.environ.get("NETMHCIIPAN_BIN", "netMHCIIpan"),
    "netchop":       os.environ.get("NETCHOP_BIN", "netchop"),
    "netmhcstabpan": os.environ.get("NETMHCSTABPAN_BIN", "netMHCstabpan"),
}


def which(model):
    """Return resolved binary path, or None if not installed."""
    return shutil.which(BINARIES[model])


def _to_rows(binding_predictions, model, source):
    """mhctools BindingPrediction list -> normalized rows.

    Verified field names (mhctools 3.15 BindingPrediction._fields):
      source_sequence_name, offset, peptide, allele, score, affinity,
      percentile_rank, prediction_method_name
    Note: there is NO `length` field — derive it from the peptide; and the
    id field is `source_sequence_name` (not `source_sequence_key`).
    """
    is_mhc2 = "ii" in model
    strong = RANK_STRONG_II if is_mhc2 else RANK_STRONG_I
    weak   = RANK_WEAK_II if is_mhc2 else RANK_WEAK_I
    rows = []
    for p in binding_predictions:
        d = p.to_dict() if hasattr(p, "to_dict") else dict(p._asdict())
        rank = d.get("percentile_rank")
        binder = "none"
        if rank is not None and pd.notna(rank):
            if rank <= strong: binder = "strong"
            elif rank <= weak: binder = "weak"
        pep = d.get("peptide") or ""
        rows.append({
            "id": d.get("source_sequence_name"),
            "source": source, "model": model,
            "allele": d.get("allele"),
            "length": len(pep),
            "peptide": pep,
            "position": d.get("offset"),
            "affinity_nM": d.get("affinity"),
            "percentile_rank": rank,
            "score": d.get("score"),
            "is_binder": binder,
        })
    return rows


def run_netmhcpan(sequences, alleles, lengths=(8, 9, 10, 11)):
    """MHC-I binding. `sequences`: {name: protein_seq}. Returns normalized DataFrame."""
    from mhctools import NetMHCpan41
    if which("netmhcpan") is None:
        return _skipped("netmhcpan")
    pred = NetMHCpan41(alleles=list(alleles),
                       default_peptide_lengths=list(lengths),
                       program_name=BINARIES["netmhcpan"])
    res = pred.predict_subsequences(sequences)
    return pd.DataFrame(_to_rows(res, "netmhcpan", "standalone"))


def run_netmhciipan(sequences, alleles, lengths=(15,)):
    """MHC-II binding."""
    from mhctools import NetMHCIIpan43
    if which("netmhciipan") is None:
        return _skipped("netmhciipan")
    pred = NetMHCIIpan43(alleles=list(alleles),
                         default_peptide_lengths=list(lengths),
                         program_name=BINARIES["netmhciipan"])
    res = pred.predict_subsequences(sequences)
    return pd.DataFrame(_to_rows(res, "netmhciipan", "standalone"))


def run_netmhcstabpan(sequences, alleles, lengths=(9,)):
    """pMHC-I stability (half-life)."""
    from mhctools import NetMHCstabpan
    if which("netmhcstabpan") is None:
        return _skipped("netmhcstabpan")
    pred = NetMHCstabpan(alleles=list(alleles),
                         default_peptide_lengths=list(lengths),
                         program_name=BINARIES["netmhcstabpan"])
    res = pred.predict_subsequences(sequences)
    return pd.DataFrame(_to_rows(res, "netmhcstabpan", "standalone"))


def _skipped(model):
    return pd.DataFrame([{
        "id": None, "source": "standalone", "model": model, "allele": None,
        "length": None, "peptide": None, "position": None, "affinity_nM": None,
        "percentile_rank": None, "score": None,
        "is_binder": f"SKIPPED:{BINARIES[model]}_not_on_PATH",
    }])


def status_report():
    """Which standalone binaries are installed here."""
    return {m: (which(m) or "NOT FOUND") for m in BINARIES}


if __name__ == "__main__":
    import json
    print("Standalone MHC binary status:")
    print(json.dumps(status_report(), indent=2))
    print("\nInstall via https://services.healthtech.dtu.dk/ (academic license).")
    print("Then export e.g.:  export NETMHCPAN_BIN=/path/to/netMHCpan")
