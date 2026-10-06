"""Tests for the AP Novo package audit.

The load-bearing ones:

* `test_bare_integers_are_read_as_designable_lengths` -- pins AP Novo's own
  grammar rule. Everything the audit concludes follows from it.
* `test_broken_indexed_job_is_unsatisfiable` / `test_broken_unindexed_job_...`
  -- reproduce the two shipped failures, with the same numbers the upstream
  sampler produced.
* `test_fixed_manifest_jobs_are_all_satisfiable` -- the replacement package
  stays runnable.
* `test_adenylate_split_puts_n346_c348_on_the_ncaa_side` -- the cross-source
  consistency check between this motif and the directed-evolution literature.
"""

from __future__ import annotations

import json
import pathlib
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from trnaplat.apnovo import audit_apnovo_package as audit  # noqa: E402

FIXED = ROOT / "trnaplat/apnovo/fixed"

#: The shipped manifest's indexed job, verbatim.
BROKEN_INDEXED = (
    "A300,302,305,306,330,332,338,339,342,344,346,348,384,396,399,401,417,426,"
    "437,60-120,40-80,40-80,40-80/B1"
)


# ------------------------------------------------------------- the grammar

def test_bare_integers_are_read_as_designable_lengths():
    assert audit.is_designable("302")
    assert audit.is_designable("60-120")
    assert not audit.is_designable("A300")
    assert not audit.is_designable("B1")


def test_broken_indexed_job_is_unsatisfiable():
    check = audit.check_motif_str(BROKEN_INDEXED, "260-360", declared_residues=19)
    assert check["n_fixed"] == 1, "only A300 carries a chain letter"
    assert check["silently_dropped"] == 18
    assert sum(check["exact_length_segments"]) == 6492
    assert (check["designable_min"], check["designable_max"]) == (6672, 6852)
    assert not check["ok"]
    assert "Cannot sample values in range [260, 360]" in check["error"]


def test_broken_unindexed_job_is_unsatisfiable():
    check = audit.check_motif_str("60-250/B1", "260-360")
    assert not check["ok"]
    assert "[60, 250]" in check["error"]


def test_a_correctly_interleaved_motif_is_satisfiable():
    motif = "20-100,A56,20-70,A100,20-70,A125,20-100/A2001"   # the repo's example
    check = audit.check_motif_str(motif, None)
    assert check["ok"]
    assert check["fixed_residues"] == ["A56", "A100", "A125"]
    assert not check["exact_length_segments"]


def test_placeholder_form_counts_its_prefix_motifs():
    motif = "A56,A100,A125|20-100,{},20-70,{},20-70,{},20-100/A2001"
    check = audit.check_motif_str(motif, None)
    assert check["ok"]
    assert check["n_fixed"] == 3


def test_more_than_one_designable_chain_is_rejected():
    check = audit.check_motif_str("50-100/60-120/B1", None)
    assert not check["ok"]
    assert "exactly 1 designable chain" in check["error"]


def test_seq_length_bounds_only_the_designable_sum():
    """The subtlety behind defect 2: fixed residues never enter the sum."""
    check = audit.check_motif_str("100-100,A5,100-100/B1", "200-200")
    assert check["ok"], "200 designable + 1 fixed must satisfy seq_length 200"
    assert (check["designable_min"], check["designable_max"]) == (200, 200)


# ---------------------------------------------------------------- the CIF

@pytest.fixture(scope="module")
def full_motif():
    path = FIXED / "motif_active_site_full.cif"
    if not path.exists():
        pytest.skip("run `python3 -m trnaplat.apnovo.build_fixed_package` first")
    return audit.read_motif_cif(path)


def test_motif_residues_match_wild_type_mm_pylrs(full_motif):
    rows = audit.check_numbering(full_motif["residues"])
    assert len(rows) == 19
    assert all(r["match"] for r in rows), [r for r in rows if not r["match"]]


def test_ligand_is_the_adenylate_not_the_free_ncaa(full_motif):
    comps = {a["comp"] for a in full_motif["ligand_atoms"]}
    assert comps == {"YLY"}
    # CCD YLY is C22 H35 N8 O9 P; heavy atoms only here.
    counts = {}
    for atom in full_motif["ligand_atoms"]:
        counts[atom["element"]] = counts.get(atom["element"], 0) + 1
    assert counts == {"C": 22, "N": 8, "O": 9, "P": 1}
    assert len(full_motif["ligand_atoms"]) == 40, "free pyrrolysine has 18"


def test_adenylate_split_puts_n346_c348_on_the_ncaa_side(full_motif):
    rows = audit.partition_by_moiety(full_motif)
    side = {r["residue"]: r["side"] for r in rows}
    assert side["ASN346"] == "ncAA"
    assert side["CYS348"] == "ncAA"
    n_amp = sum(1 for r in rows if r["side"] == "AMP")
    assert n_amp == 10, "half the 'ncAA pocket' motif is the ATP site"


def test_ncaa_pocket_cif_keeps_only_the_amino_acid_side():
    path = FIXED / "motif_ncaa_pocket.cif"
    if not path.exists():
        pytest.skip("run `python3 -m trnaplat.apnovo.build_fixed_package` first")
    info = audit.read_motif_cif(path)
    ids = sorted(r for _, r in info["residues"])
    assert ids == [300, 302, 305, 306, 346, 348, 384, 401, 417]
    # the ligand must stay whole, or CCD YLY no longer matches it
    assert len(info["ligand_atoms"]) == 40


# ------------------------------------------------------- the fixed package

@pytest.mark.parametrize("name", ["manifest_fixed.json",
                                  "manifest_partial_diffusion.json"])
def test_fixed_manifest_jobs_are_all_satisfiable(name):
    path = FIXED / name
    if not path.exists():
        pytest.skip("run `python3 -m trnaplat.apnovo.build_fixed_package` first")
    spec = json.loads(path.read_text())
    defaults = spec.get("defaults", {})
    assert spec["designs"], "a manifest with no jobs proves nothing"
    for job in spec["designs"]:
        seq_length = job.get("seq_length", defaults.get("seq_length"))
        check = audit.check_motif_str(job["motif_str"], seq_length)
        assert check["ok"], f"{job['name']}: {check['error']}"
        assert not check["exact_length_segments"], (
            f"{job['name']} still has a bare integer that would be read as a "
            "designable length"
        )


@pytest.mark.parametrize("name", ["manifest_fixed.json",
                                  "manifest_partial_diffusion.json"])
def test_fixed_manifests_have_no_unknown_top_level_keys(name):
    """`Manifest` rejects unknown top-level keys -- caught by the upstream check."""
    path = FIXED / name
    if not path.exists():
        pytest.skip("run `python3 -m trnaplat.apnovo.build_fixed_package` first")
    allowed = {"defaults", "designs", "evaluation", "folding", "resequence",
               "settings"}
    extra = set(json.loads(path.read_text())) - allowed
    assert not extra, f"{name} has keys the schema will reject: {sorted(extra)}"


def test_fixed_manifest_declares_the_ligand_state_by_ccd_code():
    path = FIXED / "manifest_fixed.json"
    if not path.exists():
        pytest.skip("run `python3 -m trnaplat.apnovo.build_fixed_package` first")
    spec = json.loads(path.read_text())
    states = {s["name"]: s for s in spec["folding"]["states"]}
    assert states["apo_monomer"]["ligands"] == []
    assert states["adenylate_complex"]["ligands"][0]["ccd_code"] == "YLY"


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
