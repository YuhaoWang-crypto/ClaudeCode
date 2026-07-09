#!/usr/bin/env python3
"""
IEDB cloud REST backend for the DTU multi-model orchestrator.

No local install and no license needed: IEDB hosts NetMHCpan / NetMHCIIpan /
NetChop as a public REST API that stays in sync with the DTU tools. Output is
normalized to the SAME schema as the BioLib and standalone backends:

  id, source, model, allele, length, peptide, position,
  affinity_nM, percentile_rank, score, is_binder

Endpoints:
  MHC-I  : http://tools-cluster-interface.iedb.org/tools_api/mhci/
  MHC-II : http://tools-cluster-interface.iedb.org/tools_api/mhcii/

Verified in-session: netmhcpan_el (MHC-I) — returns a real 200 with the
tab-delimited schema below. Other IEDB methods (netmhcpan_ba, ann, smm, and the
MHC-II methods) use the same request shape but were NOT exercised here; note
that netmhcpan_ba returned HTTP 500 from this endpoint at test time, so confirm
the exact method name against the live IEDB method list before relying on it.
%rank thresholds follow DTU/IEDB conventions (MHC-I strong<=0.5, weak<=2.0).
"""
import io, requests
import pandas as pd

MHCI_URL  = "http://tools-cluster-interface.iedb.org/tools_api/mhci/"
MHCII_URL = "http://tools-cluster-interface.iedb.org/tools_api/mhcii/"
RANK_STRONG_I, RANK_WEAK_I = 0.5, 2.0
RANK_STRONG_II, RANK_WEAK_II = 2.0, 10.0


def _classify(rank, mhc2):
    strong = RANK_STRONG_II if mhc2 else RANK_STRONG_I
    weak   = RANK_WEAK_II if mhc2 else RANK_WEAK_I
    if rank is None or pd.isna(rank): return "none"
    if rank <= strong: return "strong"
    if rank <= weak:   return "weak"
    return "none"


def _post(url, data, timeout=180):
    r = requests.post(url, data=data, timeout=timeout)
    r.raise_for_status()
    return pd.read_csv(io.StringIO(r.text), sep="\t")


def run_iedb_mhci(sequences, alleles, lengths=(9,), method="netmhcpan_el", seq_id_map=None):
    """MHC-I binding via IEDB. `sequences`: {name: seq}. Returns normalized DataFrame."""
    # IEDB pairs each allele with a length: allele="A,B" needs length="9,9".
    # Expand the allele x length grid so every combination is requested.
    pairs = [(a, l) for a in alleles for l in lengths]
    allele_str = ",".join(a for a, _ in pairs)
    length_str = ",".join(str(l) for _, l in pairs)
    frames = []
    for name, seq in sequences.items():
        data = {"method": method, "sequence_text": seq,
                "allele": allele_str, "length": length_str}
        df = _post(MHCI_URL, data)
        df["__id"] = name
        frames.append(df)
    raw = pd.concat(frames, ignore_index=True)
    return _normalize(raw, model=method, mhc2=False)


def run_iedb_mhcii(sequences, alleles, method="netmhciipan_el", lengths=(15,), seq_id_map=None):
    """MHC-II binding via IEDB."""
    frames = []
    for name, seq in sequences.items():
        data = {"method": method, "sequence_text": seq,
                "allele": ",".join(alleles),
                "length": ",".join(str(l) for l in lengths)}
        df = _post(MHCII_URL, data)
        df["__id"] = name
        frames.append(df)
    raw = pd.concat(frames, ignore_index=True)
    return _normalize(raw, model=method, mhc2=True)


def _normalize(raw, model, mhc2):
    """IEDB tsv -> unified schema. IEDB columns:
       allele, seq_num, start, end, length, peptide, core, icore,
       score, percentile_rank  (+ ic50 when a BA method is used)."""
    rows = []
    for _, r in raw.iterrows():
        rank = r.get("percentile_rank")
        rank = float(rank) if pd.notna(rank) else None
        ic50 = None
        for c in ("ic50", "affinity", "aff(nm)"):
            if c in raw.columns and pd.notna(r.get(c)):
                ic50 = float(r[c]); break
        rows.append({
            "id": r.get("__id"),
            "source": "iedb", "model": model,
            "allele": r.get("allele"),
            "length": int(r["length"]) if pd.notna(r.get("length")) else len(str(r.get("peptide",""))),
            "peptide": r.get("peptide"),
            "position": int(r["start"]) - 1 if pd.notna(r.get("start")) else None,  # 0-based
            "affinity_nM": ic50,
            "percentile_rank": rank,
            "score": float(r["score"]) if pd.notna(r.get("score")) else None,
            "is_binder": _classify(rank, mhc2),
        })
    return pd.DataFrame(rows)


if __name__ == "__main__":
    seqs = {"seq1": "SIINFEHLLKAPWDVKELS"}
    df = run_iedb_mhci(seqs, alleles=["HLA-A*02:01"], lengths=[9], method="netmhcpan_el")
    print(df.to_string())
    print("\nstrong binders:", int((df.is_binder=="strong").sum()),
          "| weak:", int((df.is_binder=="weak").sum()))
