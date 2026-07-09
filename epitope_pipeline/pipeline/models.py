"""Registry of epitope / immunogenicity models and their CLI adapters.

Each model is described by a `ModelSpec`. An adapter knows how to:
  1. report whether the tool is installed (`resolve_binary`),
  2. build the shell command for a given FASTA + allele set (`build_command`),
  3. parse the tool's stdout into normalized `EpitopeRecord`s (`parse`).

The command templates below match the documented CLI of the current
stand-alone packages. Versions drift, so every template is overridable from
config.yaml / the CLI without touching code.
"""

from __future__ import annotations

import re
import shutil
from dataclasses import dataclass, field
from typing import Callable, List, Optional


# --------------------------------------------------------------------------- #
# Normalized output record shared by every model
# --------------------------------------------------------------------------- #
@dataclass
class EpitopeRecord:
    model: str            # e.g. "NetMHCpan-4.1"
    category: str         # mhc_i | mhc_ii | bcell_linear | bcell_conf | aux
    seq_id: str           # FASTA identifier the peptide came from
    peptide: str          # the epitope / peptide sequence (may be "" for per-residue)
    position: Optional[int] = None
    allele: Optional[str] = None
    core: Optional[str] = None
    score: Optional[float] = None       # raw model score (higher = better binder)
    rank: Optional[float] = None        # %Rank (lower = stronger); the key metric
    affinity_nM: Optional[float] = None
    bind_level: Optional[str] = None    # SB / WB / "" as reported
    extra: dict = field(default_factory=dict)

    def as_row(self) -> dict:
        return {
            "model": self.model,
            "category": self.category,
            "seq_id": self.seq_id,
            "position": self.position,
            "allele": self.allele,
            "peptide": self.peptide,
            "core": self.core,
            "score": self.score,
            "rank_pct": self.rank,
            "affinity_nM": self.affinity_nM,
            "bind_level": self.bind_level,
        }


# --------------------------------------------------------------------------- #
# Model specification
# --------------------------------------------------------------------------- #
@dataclass
class ModelSpec:
    name: str
    category: str
    binary: str                       # executable name expected on PATH
    # build_command(spec, fasta, alleles, lengths, workdir) -> list[str]
    build_command: Callable
    # parse(spec, stdout, fasta_ids) -> list[EpitopeRecord]
    parse: Callable
    needs_alleles: bool = True
    needs_structure: bool = False     # DiscoTope etc. need a PDB, not FASTA
    default_lengths: tuple = ()
    homepage: str = ""
    install_hint: str = ""

    def resolve_binary(self, override: Optional[str] = None) -> Optional[str]:
        cand = override or self.binary
        return shutil.which(cand)


# --------------------------------------------------------------------------- #
# Command builders
# --------------------------------------------------------------------------- #
def _cmd_netmhcpan(spec, fasta, alleles, lengths, workdir):
    lens = ",".join(str(x) for x in (lengths or spec.default_lengths))
    return [spec.binary, "-f", fasta, "-a", ",".join(alleles),
            "-l", lens, "-BA"]


def _cmd_netmhciipan(spec, fasta, alleles, lengths, workdir):
    lens = ",".join(str(x) for x in (lengths or spec.default_lengths))
    return [spec.binary, "-f", fasta, "-a", ",".join(alleles),
            "-length", lens, "-BA"]


def _cmd_netctlpan(spec, fasta, alleles, lengths, workdir):
    # NetCTLpan takes one allele at a time; the runner iterates alleles for it.
    return [spec.binary, "-f", fasta, "-a", ",".join(alleles)]


def _cmd_bepipred(spec, fasta, alleles, lengths, workdir):
    # BepiPred-3.0 CLI: bepipred3_CLI.py -i in.fasta -o outdir -pred mjv_pred
    return [spec.binary, "-i", fasta, "-o", workdir, "-pred", "mjv_pred"]


def _cmd_signalp(spec, fasta, alleles, lengths, workdir):
    return [spec.binary, "-fasta", fasta, "-org", "euk", "-format", "short"]


# --------------------------------------------------------------------------- #
# Parsers
# --------------------------------------------------------------------------- #
_NUM = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"


def _to_float(tok: str) -> Optional[float]:
    try:
        return float(tok)
    except (TypeError, ValueError):
        return None


def _parse_netmhcpan(spec, stdout, fasta_ids):
    """Parse the fixed-width table emitted by NetMHCpan / NetMHCpan-BA.

    Columns (4.1): Pos MHC Peptide Core Of Gp Gl Ip Il Icore Identity
                   Score_EL %Rank_EL Score_BA %Rank_BA Aff(nM) BindLevel
    We locate columns by header when present, else fall back to positional.
    """
    records: List[EpitopeRecord] = []
    for line in stdout.splitlines():
        s = line.strip()
        if not s or s.startswith("#") or s.startswith("-"):
            continue
        if s.startswith("Pos") or s.lower().startswith("protein"):
            continue
        toks = s.split()
        # A data row starts with an integer position and has >= 12 tokens.
        if len(toks) < 12 or not toks[0].isdigit():
            continue
        try:
            pos = int(toks[0])
            allele = toks[1]
            peptide = toks[2]
            core = toks[3]
            identity = toks[10]
            rank_el = _to_float(toks[12]) if len(toks) > 12 else None
            aff = None
            rank_ba = None
            # BA columns appear only with -BA; grab trailing numeric fields.
            nums = [t for t in toks[11:] if re.fullmatch(_NUM, t)]
            if len(nums) >= 5:
                # Score_EL %Rank_EL Score_BA %Rank_BA Aff(nM)
                rank_el = _to_float(nums[1])
                rank_ba = _to_float(nums[3])
                aff = _to_float(nums[4])
            bind = ""
            if "<= SB" in s:
                bind = "SB"
            elif "<= WB" in s:
                bind = "WB"
            records.append(EpitopeRecord(
                model=spec.name, category=spec.category, seq_id=identity,
                peptide=peptide, position=pos, allele=allele, core=core,
                score=rank_el, rank=rank_ba if rank_ba is not None else rank_el,
                affinity_nM=aff, bind_level=bind,
                extra={"rank_el": rank_el},
            ))
        except (ValueError, IndexError):
            continue
    return records


def _parse_netmhciipan(spec, stdout, fasta_ids):
    """NetMHCIIpan table: Pos MHC Peptide Of Core Core_Rel Identity
    Score_EL %Rank_EL Exp_Bind Score_BA Affinity %Rank_BA BindLevel."""
    records: List[EpitopeRecord] = []
    for line in stdout.splitlines():
        s = line.strip()
        if not s or s.startswith("#") or s.startswith("-"):
            continue
        toks = s.split()
        if len(toks) < 9 or not toks[0].isdigit():
            continue
        try:
            pos = int(toks[0])
            allele = toks[1]
            peptide = toks[2]
            core = toks[4]
            identity = toks[6]
            nums = [t for t in toks[6:] if re.fullmatch(_NUM, t)]
            rank_el = _to_float(nums[1]) if len(nums) > 1 else None
            rank_ba = _to_float(nums[-1]) if len(nums) >= 5 else None
            aff = _to_float(nums[-2]) if len(nums) >= 5 else None
            bind = "SB" if "<=SB" in s.replace(" ", "") else (
                "WB" if "<=WB" in s.replace(" ", "") else "")
            records.append(EpitopeRecord(
                model=spec.name, category=spec.category, seq_id=identity,
                peptide=peptide, position=pos, allele=allele, core=core,
                score=rank_el, rank=rank_ba if rank_ba is not None else rank_el,
                affinity_nM=aff, bind_level=bind,
            ))
        except (ValueError, IndexError):
            continue
    return records


def _parse_netctlpan(spec, stdout, fasta_ids):
    """NetCTLpan combined score table.
    Columns: N pep_id allele peptide MHC TAP Cle Comb %Rank."""
    records: List[EpitopeRecord] = []
    for line in stdout.splitlines():
        s = line.strip()
        if not s or s.startswith("#") or s.startswith("-") or s.startswith("N "):
            continue
        toks = s.split()
        if len(toks) < 9 or not toks[0].isdigit():
            continue
        try:
            records.append(EpitopeRecord(
                model=spec.name, category=spec.category, seq_id=toks[1],
                peptide=toks[3], position=int(toks[0]), allele=toks[2],
                score=_to_float(toks[7]), rank=_to_float(toks[8]),
                bind_level="SB" if "<-E" in s else "",
            ))
        except (ValueError, IndexError):
            continue
    return records


def _parse_bepipred(spec, stdout, fasta_ids):
    """BepiPred writes per-residue scores; we surface residues above threshold
    as candidate linear B-cell epitope positions. Output is read from the
    raw_output.csv the runner captures, but we also tolerate stdout CSV."""
    records: List[EpitopeRecord] = []
    for line in stdout.splitlines():
        parts = line.strip().split(",")
        if len(parts) < 3:
            continue
        seq_id, residue, score = parts[0], parts[1], _to_float(parts[-1])
        if score is None:
            continue
        records.append(EpitopeRecord(
            model=spec.name, category=spec.category, seq_id=seq_id,
            peptide=residue, score=score, rank=None,
            bind_level="epitope" if score >= 0.1512 else "",
        ))
    return records


def _parse_signalp(spec, stdout, fasta_ids):
    records: List[EpitopeRecord] = []
    for line in stdout.splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        toks = s.split()
        if len(toks) < 3:
            continue
        records.append(EpitopeRecord(
            model=spec.name, category=spec.category, seq_id=toks[0],
            peptide=toks[1], score=None,
            extra={"prediction": " ".join(toks[1:])},
        ))
    return records


# --------------------------------------------------------------------------- #
# The registry
# --------------------------------------------------------------------------- #
REGISTRY = {
    "netmhcpan": ModelSpec(
        name="NetMHCpan-4.1", category="mhc_i", binary="netMHCpan",
        build_command=_cmd_netmhcpan, parse=_parse_netmhcpan,
        needs_alleles=True, default_lengths=(8, 9, 10, 11),
        homepage="https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/",
        install_hint="Download stand-alone tarball (academic), add netMHCpan to PATH.",
    ),
    "netmhciipan": ModelSpec(
        name="NetMHCIIpan-4.3", category="mhc_ii", binary="netMHCIIpan",
        build_command=_cmd_netmhciipan, parse=_parse_netmhciipan,
        needs_alleles=True, default_lengths=(15,),
        homepage="https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/",
        install_hint="Download stand-alone tarball (academic), add netMHCIIpan to PATH.",
    ),
    "netctlpan": ModelSpec(
        name="NetCTLpan-1.1", category="mhc_i", binary="netCTLpan",
        build_command=_cmd_netctlpan, parse=_parse_netctlpan,
        needs_alleles=True, default_lengths=(9,),
        homepage="https://services.healthtech.dtu.dk/services/NetCTLpan-1.1/",
        install_hint="Download stand-alone tarball (academic), add netCTLpan to PATH.",
    ),
    "bepipred": ModelSpec(
        name="BepiPred-3.0", category="bcell_linear", binary="bepipred3_CLI.py",
        build_command=_cmd_bepipred, parse=_parse_bepipred,
        needs_alleles=False,
        homepage="https://services.healthtech.dtu.dk/services/BepiPred-3.0/",
        install_hint="pip install the BepiPred-3.0 package; exposes bepipred3_CLI.py.",
    ),
    "signalp": ModelSpec(
        name="SignalP-6.0", category="aux", binary="signalp6",
        build_command=_cmd_signalp, parse=_parse_signalp,
        needs_alleles=False,
        homepage="https://services.healthtech.dtu.dk/services/SignalP-6.0/",
        install_hint="pip install the SignalP-6.0 package; exposes signalp6.",
    ),
}


def get_models(keys: List[str]) -> List[ModelSpec]:
    out = []
    for k in keys:
        k = k.lower()
        if k not in REGISTRY:
            raise KeyError(f"Unknown model '{k}'. Known: {', '.join(REGISTRY)}")
        out.append(REGISTRY[k])
    return out
