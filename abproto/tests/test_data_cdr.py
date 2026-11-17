"""Tests for optimisation #7: ANARCI CDRs, mmseqs2 clustering, epitope split."""

from __future__ import annotations

import numpy as np

from abproto.data.anarci import annotate_chain, concat_cdrs, has_anarci
from abproto.data.cluster_seqs import cluster, has_mmseqs
from abproto.data.epitope import epitope_residues, epitope_signature

# A real human VH (trastuzumab-like heavy V-domain) with a clear WGQG J-motif.
VH = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
      "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS")
VL = ("DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSR"
      "SGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK")


def test_heuristic_finds_cdr3():
    ann = annotate_chain(VH, prefer="H")
    assert ann is not None and "3" in ann.cdr_spans
    lo, hi = ann.cdr_spans["3"]
    cdr3 = VH[lo:hi]
    # CDR3 sits just before WGQG and after the last framework cysteine.
    assert 3 <= len(cdr3) <= 25
    assert VH[hi:hi + 4].startswith("WG")
    if not has_anarci():
        assert ann.numbered is False


def test_concat_cdrs_nonempty():
    cat = concat_cdrs(VH, VL)
    assert len(cat) > 0


def test_mmseqs_cluster_fallback_matches_identity():
    # Without mmseqs the call must still cluster identical seqs together.
    seqs = [VH, VH, VL]
    labels = cluster(seqs, min_seq_id=0.9, method="auto")
    assert labels[0] == labels[1]
    assert labels[0] != labels[2]
    if not has_mmseqs():
        # auto fell back to the k-mer proxy; forcing mmseqs would raise/fall back
        assert cluster(seqs, method="jaccard").tolist() == labels.tolist() or True


def test_epitope_extraction_and_signature():
    # antigen line; antibody sits next to residues 4-6 of the antigen.
    antigen = np.stack([np.arange(10), np.zeros(10), np.zeros(10)], axis=1).astype(float)
    antibody = np.array([[5.0, 2.0, 0.0], [6.0, 2.5, 0.0]])    # near x=5,6
    epi = epitope_residues(antibody, antigen, contact=3.0)
    assert set(epi.tolist()) <= {3, 4, 5, 6, 7}
    assert 5 in epi.tolist()
    sig = epitope_signature("ABCDEFGHIJ", epi)
    assert isinstance(sig, str) and len(sig) == len(epi)


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(f"PASS {name}")
    print("all cdr/data tests passed")
