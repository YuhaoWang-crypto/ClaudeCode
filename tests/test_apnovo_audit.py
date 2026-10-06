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

import collections
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


def test_cif_header_check_rejects_a_space_after_data():
    """The defect AlphaFold 3 found: `data <name>` instead of `data_<name>`."""
    problem = audit.check_cif_header(["data pylrs_pyl_motif", "_entry.id x"])
    assert problem is not None
    assert "data_<name>" in problem


def test_cif_header_check_accepts_a_legal_header():
    assert audit.check_cif_header(["# a comment", "", "data_motif",
                                   "_entry.id motif"]) is None


def test_cif_header_check_rejects_a_bare_data_token():
    assert audit.check_cif_header(["data_"]) is not None


def test_read_motif_cif_refuses_a_malformed_header_by_default(tmp_path):
    path = tmp_path / "bad.cif"
    path.write_text("data bad\nloop_\n_atom_site.group_PDB\n_atom_site.type_symbol\n"
                    "_atom_site.label_atom_id\n_atom_site.label_comp_id\n"
                    "_atom_site.label_asym_id\n_atom_site.label_seq_id\n"
                    "_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n"
                    "_atom_site.auth_seq_id\n_atom_site.auth_asym_id\n"
                    "ATOM N N MET A 1 0.0 0.0 0.0 1 A\n")
    with pytest.raises(audit.CifHeaderError):
        audit.read_motif_cif(path)
    assert audit.read_motif_cif(path, strict=False)["n_atoms"] == 1


def test_fixed_cifs_have_a_legal_data_block_header():
    for name in ("motif_active_site_full.cif", "motif_ncaa_pocket.cif"):
        path = FIXED / name
        if not path.exists():
            pytest.skip("run `python3 -m trnaplat.apnovo.build_fixed_package` first")
        assert audit.check_cif_header(path.read_text().splitlines()) is None, name


def test_fixed_cifs_carry_every_column_af3_reads():
    """Enumerated from alphafold3/structure/parsing.py, not guessed."""
    from trnaplat.apnovo.build_fixed_package import ATOM_SITE_COLUMNS
    af3_reads = {
        "group_PDB", "type_symbol", "label_atom_id", "label_comp_id",
        "label_asym_id", "label_entity_id", "label_seq_id",
        "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z", "occupancy",
        "B_iso_or_equiv", "auth_seq_id", "pdbx_PDB_model_num",
    }
    assert af3_reads <= set(ATOM_SITE_COLUMNS)
    for name in ("motif_active_site_full.cif", "motif_ncaa_pocket.cif"):
        path = FIXED / name
        if not path.exists():
            pytest.skip("run `python3 -m trnaplat.apnovo.build_fixed_package` first")
        present = {l.strip().split(".", 1)[1]
                   for l in path.read_text().splitlines()
                   if l.strip().startswith("_atom_site.")}
        assert af3_reads <= present, f"{name} missing {sorted(af3_reads - present)}"


def test_fixed_cifs_keep_author_numbering():
    """auth_seq_id holds the motif's identity, which the manifest refers to."""
    path = FIXED / "motif_ncaa_pocket.cif"
    if not path.exists():
        pytest.skip("run `python3 -m trnaplat.apnovo.build_fixed_package` first")
    lines = path.read_text().splitlines()
    cols = [l.strip().split(".", 1)[1] for l in lines
            if l.strip().startswith("_atom_site.")]
    idx = {name: i for i, name in enumerate(cols)}
    auth, label, ligand_label = set(), set(), set()
    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        row = line.split()
        if row[idx["label_comp_id"]] in {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN",
                                         "GLU", "GLY", "HIS", "ILE", "LEU", "LYS",
                                         "MET", "PHE", "PRO", "SER", "THR", "TRP",
                                         "TYR", "VAL"}:
            auth.add(int(row[idx["auth_seq_id"]]))
            label.add(row[idx["label_seq_id"]])
        else:
            ligand_label.add(row[idx["label_seq_id"]])
    assert auth == {300, 302, 305, 306, 346, 348, 384, 401, 417}
    assert len(label) == 9, "each residue needs its own label_seq_id"
    assert all(v.isdigit() for v in label)
    assert ligand_label <= {".", "?"}, "a non-polymer atom has no label_seq_id"


def test_normalise_header_repairs_the_shipped_defect():
    from trnaplat.apnovo.build_fixed_package import normalise_header
    out = normalise_header(["data pylrs_pyl_motif", "_entry.id pylrs_pyl_motif"])
    assert out[0] == "data_pylrs_pyl_motif"
    assert out[1] == "_entry.id pylrs_pyl_motif"
    # already legal headers are left alone
    assert normalise_header(["data_ok", "x"])[0] == "data_ok"


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
    counts = collections.Counter(r["side"] for r in rows)
    # ✅ 9 unambiguously the amino-acid pocket, 6 within 0.5 A of both halves,
    # 4 the nucleotide site. Reproduced identically from the audited CIF and
    # from an independent gemmi extraction of 2Q7H.
    assert counts["ncAA"] == 9
    assert counts["both"] == 6
    assert counts["AMP"] == 4
    assert counts["ncAA"] + counts["both"] + counts["AMP"] == 19


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


# --------------------------------------------------- the Modal licence gate

def _modal_app_source() -> str:
    return (ROOT / "trnaplat/apnovo/modal_app.py").read_text()


def test_generate_is_gated_on_an_explicit_attestation():
    """The weighted step must refuse by default, not merely warn."""
    source = _modal_app_source()
    assert "attest_non_commercial: bool = False" in source, \
        "the gate must default to refusing"
    gate = source.split("def generate(")[1].split("def ")[0]
    assert "if not attest_non_commercial:" in gate
    assert "raise RuntimeError" in gate
    # the refusal must come before any weight download
    refuse_at = gate.index("raise RuntimeError")
    assert gate.index("wget") > refuse_at, \
        "weights must not be fetched before the attestation is checked"


def test_generate_keeps_folding_inputs_consistent_with_resequence():
    """The pipeline refuses folding.inputs=["resequenced"] with resequence off."""
    gate = _modal_app_source().split("def generate(")[1].split("\ndef ")[0]
    assert '"enabled": False' in gate
    assert '["generated"]' in gate, (
        "disabling resequence without switching folding.inputs to 'generated' "
        "makes run_pipeline.py refuse the manifest"
    )


def test_weight_urls_point_at_the_official_buckets():
    source = _modal_app_source()
    assert "storage.googleapis.com/alphaprotein_novo/generator.bin.zst" in source
    assert "storage.googleapis.com/alphafold3/af3_leaving_atom.bin.zst" in source


def test_licence_free_functions_declare_no_weight_volume():
    """featurize and gpu_probe must not even mount the weights volume."""
    source = _modal_app_source()
    for name in ("featurize", "gpu_probe"):
        block = source.split(f"def {name}(")[0].rsplit("@app.function", 1)[1]
        assert "volumes=" not in block, f"{name} should not mount the weights volume"


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
