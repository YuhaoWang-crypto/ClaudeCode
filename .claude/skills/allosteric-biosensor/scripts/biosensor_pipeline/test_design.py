"""
test_design.py -- reproducibility + correctness tests for the design engine.

Run:  python3 -m biosensor_pipeline.test_design
Pure, deterministic, offline.  No network, no Boltz.
"""

from __future__ import annotations
from .design import circular_permute, build_chimera, verify_chimera
from .screen import build_library
from .systems import SYSTEMS, TEM1


def test_circular_permutation_conserves_residues():
    seq = "ACDEFGHIKL"
    site = 4  # remove 'F'
    cp = circular_permute(seq, site, "GS")
    assert cp == "GHIKL" + "GS" + "ACDE", cp
    body = cp.replace("GS", "", 1)
    assert sorted(body) == sorted(seq[:site] + seq[site + 1:])
    print("[OK] circular permutation removes exactly the site residue, conserves the rest")


def test_insertion_preserves_reporter():
    reporter = "MREPORTER" * 3
    receptor = "ABCDEFGHIK"
    ch = build_chimera(
        name="t", reporter_name="R", reporter_seq=reporter,
        receptor_name="X", receptor_seq=receptor,
        insertion_index=5, permutation_site=4, gs_linker="GGS", flank_linker="G",
    )
    assert ch.sequence.startswith(reporter[:5])
    assert ch.sequence.endswith(reporter[5:])
    chk = verify_chimera(ch, reporter, receptor)
    assert chk["all_ok"], chk
    print("[OK] insertion preserves reporter head/tail and conserves receptor residues")


def test_all_library_constructions_valid():
    for key, sysm in SYSTEMS.items():
        lib = build_library(sysm)
        assert 1 <= len(lib) < 10, f"{key}: library not <10 ({len(lib)})"
        for c in lib:
            chk = verify_chimera(c, sysm.reporter.seq, sysm.receptor.seq)
            assert chk["all_ok"], (key, c.name, chk)
        print(f"[OK] {key}: {len(lib)} valid variants (<10, as in the paper)")


def test_determinism():
    """Same inputs -> identical sequence, always."""
    sysm = SYSTEMS["dig"]
    a = build_library(sysm)
    b = build_library(sysm)
    assert [c.sequence for c in a] == [c.sequence for c in b]
    print("[OK] construction is deterministic (identical output across runs)")


def test_catalytic_residues_present():
    """The TEM-1 catalytic serines are where the motif says they are."""
    assert TEM1.seq[TEM1.catalytic["S70"]] == "S"
    assert TEM1.seq[TEM1.catalytic["S130"]] == "S"
    assert TEM1.seq[TEM1.catalytic["E166"]] == "E"
    print("[OK] TEM-1 catalytic residues verified by position")


if __name__ == "__main__":
    test_circular_permutation_conserves_residues()
    test_insertion_preserves_reporter()
    test_all_library_constructions_valid()
    test_determinism()
    test_catalytic_residues_present()
    print("\nALL TESTS PASSED ✅")
