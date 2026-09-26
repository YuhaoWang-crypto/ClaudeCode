"""Humanness and T-cell-epitope-risk scoring against the real human germline set.

WHAT THIS MEASURES, PRECISELY
-----------------------------
`germline_9mer_coverage(seq)` returns the fraction of overlapping 9-mers in a
variable-domain sequence that also occur somewhere in the human germline
IGHV/IGKV/IGLV V-region repertoire (143 reviewed UniProt entries).

The rationale is the standard one behind peptide-level humanness scores: a
9-mer that already exists in the human germline antibody repertoire is a 9-mer
the human T-cell compartment is tolerised to, so it is unlikely to be presented
as a foreign MHC-II core.

This is NOT:
  * an OASis score (that is computed against the Observed Antibody Space
    repertoire database of ~10^9 sequences, not the germline set),
  * an MHC-II binding prediction,
  * an anti-drug-antibody probability.

A high coverage means "few foreign-looking 9-mers", not "not immunogenic".
Fully human antibodies still raise ADA. Actual MHC-II risk needs
NetMHCIIpan/IEDB; `mhc2_peptides()` emits the exact peptide set to submit.
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA

# A V domain starts at the mature N-terminus, right after the ~19-20 aa signal
# peptide. UniProt germline entries include that signal peptide, so trim it.
V_START = re.compile(r"[EQDA][VIL][QKMES][LMV]")
KMER = 9


def _read_multi_fasta(path):
    out, name, buf = {}, None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name:
                out[name] = "".join(buf)
            name, buf = line[1:].strip(), []
        else:
            buf.append(line.strip())
    if name:
        out[name] = "".join(buf)
    return out


def _trim_signal_peptide(seq):
    """Return (trimmed_sequence, was_trimmed)."""
    m = V_START.search(seq, 10)
    if m and m.start() <= 30:
        return seq[m.start():], True
    return seq, False


def load_germline_kmers(k=KMER):
    """Build the human germline k-mer set. Returns (kmer_set, stats)."""
    kmers, n_seq, n_trim = set(), 0, 0
    for locus in ("IGHV", "IGKV", "IGLV"):
        path = DATA / f"germline_{locus}.fasta"
        if not path.exists():
            continue
        for _, seq in _read_multi_fasta(path).items():
            seq = seq.upper().replace("X", "")
            if len(seq) < k:
                continue
            trimmed, ok = _trim_signal_peptide(seq)
            n_seq += 1
            n_trim += int(ok)
            for i in range(len(trimmed) - k + 1):
                kmers.add(trimmed[i:i + k])
    return kmers, {"germline_sequences": n_seq, "signal_peptides_trimmed": n_trim,
                   "distinct_kmers": len(kmers), "k": k}


def germline_9mer_coverage(seq, kmers, k=KMER):
    """Fraction of the sequence's 9-mers present in the human germline set."""
    if len(seq) < k:
        return 1.0, []
    windows = [seq[i:i + k] for i in range(len(seq) - k + 1)]
    foreign = [(i + 1, w) for i, w in enumerate(windows) if w not in kmers]
    return (len(windows) - len(foreign)) / len(windows), foreign


def _aligner():
    from Bio import Align
    from Bio.Align import substitution_matrices
    al = Align.PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score = -11
    al.extend_gap_score = -1
    # local: the germline V gene covers only the V region of the query, and the
    # query may carry extra N-terminal or J-derived residues.
    al.mode = "local"
    return al


def nearest_germline(seq, locus_files=("IGHV", "IGKV", "IGLV")):
    """Closest human germline V gene by local alignment identity.

    Uses BLOSUM62 local alignment, so it tolerates the offsets and indels that
    make a naive index-by-index comparison assign a lambda chain to an IGHV
    gene. Still a reporting aid, not a substitute for IMGT/ANARCI numbering.
    """
    al = _aligner()
    best = {"identity": 0.0, "gene": None, "locus": None, "aligned_len": 0}
    for locus in locus_files:
        path = DATA / f"germline_{locus}.fasta"
        if not path.exists():
            continue
        for header, g in _read_multi_fasta(path).items():
            g, _ = _trim_signal_peptide(g.upper().replace("X", ""))
            if len(g) < 60:
                continue
            aln = al.align(seq, g)[0]
            n_ident = n_tot = 0
            for (qs, qe), (gs, ge) in zip(*aln.aligned):
                for off in range(qe - qs):
                    n_tot += 1
                    n_ident += int(seq[qs + off] == g[gs + off])
            if n_tot < 50:
                continue
            ident = n_ident / n_tot
            if ident > best["identity"]:
                gene = (header.split("GN=")[-1].split()[0]
                        if "GN=" in header else header[:24])
                best = {"identity": round(ident, 4), "gene": gene,
                        "locus": locus, "aligned_len": int(n_tot)}
    return best


def mhc2_peptides(seq, length=15, step=4, region=None):
    """15-mers (step 4) covering `region` (1-based inclusive) or the whole domain.

    This is the input format for NetMHCIIpan / IEDB MHC-II binding prediction.
    Submitting these to a third party discloses the sequences, so it needs the
    user's explicit authorisation - hence this only writes the file.
    """
    lo, hi = region if region else (1, len(seq))
    out = []
    for start in range(1, len(seq) - length + 2, step):
        end = start + length - 1
        if end < lo or start > hi:
            continue
        out.append({"start": start, "end": end, "peptide": seq[start - 1:end]})
    return out


def chemical_liabilities(seq):
    """Sequence-intrinsic chemical instability motifs in a variable domain."""
    pats = {
        "deamidation_NG_NS": r"N[GS]",
        "isomerisation_DG_DP": r"D[GP]",
        "fragmentation_DP": r"DP",
        "Nglyc_sequon": r"N[^P][ST]",
        "unpaired_Cys": r"C",
        "oxidation_M": r"M",
        "oxidation_W": r"W",
        "hydrophobic_run4": r"[ILVFWMY]{4,}",
    }
    found = {}
    for name, p in pats.items():
        hits = [m.start() + 1 for m in re.finditer(p, seq)]
        if hits:
            found[name] = hits
    return found
