"""Tests for the TCR triage pipeline.

Unit tests cover the pure logic and run offline. Integration tests need the
cached PDB entries and the proteome FASTA, and skip cleanly when those are
absent so the suite still passes without network access.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from tcr_design.constants import AA20, HLA_A0201_ANCHORS, MAX_SASA, VOLUME
from tcr_design.liability import assess_liabilities, mutate_sequence
from tcr_design.offtarget import _similarity
from tcr_design.scan import score_mutation
from tcr_design.structure import Residue, annotate_loops
from tcr_design.validate import _rank, _spearman, locate_in_scan

CACHE = Path.home() / ".cache" / "tcr_design"
PDB_CACHE = CACHE / "pdb"
PROTEOME = CACHE / "proteome" / "human_reviewed.fasta.gz"

# Published CDR sequences for the A6 TCR (PDB 1AO7), used as ground truth for
# the loop annotation.
A6_ALPHA = (
    "KEVEQNSGPLSVPEGAIASLNCTYSDRGSQSFFWYRQYSGKSPELIMSIYSNGDKEDGRFTAQLNKASQY"
    "VSLLIRDSQPSDSATYLCAVTTDSWGKLQFGAGTQVVVTPDIQNP"
)
A6_BETA = (
    "GVTQTPKFQVLKTGQSMTLQCAQDMNHEYMSWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRSTTE"
    "DFPLRLLSAAPSQTSVYFCASRPGLAGGRPEQYFGPGTRLTVTEDLKNVFPPEVAVF"
)


# --------------------------------------------------------------------------
# Reference data integrity
# --------------------------------------------------------------------------

def test_constant_tables_cover_all_twenty_residues():
    for table in (VOLUME, MAX_SASA):
        assert set(table) == set(AA20)


def test_glycine_is_the_smallest_and_tryptophan_the_largest():
    assert min(VOLUME, key=VOLUME.get) == "G"
    assert max(VOLUME, key=VOLUME.get) == "W"


# --------------------------------------------------------------------------
# Loop annotation
# --------------------------------------------------------------------------

def test_cdr3_alpha_matches_the_published_a6_junction():
    loops = annotate_loops(A6_ALPHA, "D")
    chain, start, end = loops["CDR3"]
    assert chain == "D"
    assert A6_ALPHA[start:end] == "CAVTTDSWGKLQF"


def test_cdr3_beta_matches_the_published_a6_junction():
    loops = annotate_loops(A6_BETA, "E")
    _chain, start, end = loops["CDR3"]
    assert A6_BETA[start:end] == "CASRPGLAGGRPEQYF"


def test_cdr1_matches_the_published_a6_loops():
    alpha = annotate_loops(A6_ALPHA, "D")
    _c, start, end = alpha["CDR1"]
    assert A6_ALPHA[start:end] == "DRGSQS"

    beta = annotate_loops(A6_BETA, "E")
    _c, start, end = beta["CDR1"]
    assert A6_BETA[start:end] == "MNHEY"


def test_annotate_loops_returns_nothing_for_a_non_tcr_sequence():
    assert annotate_loops("MKVLAAGIVGLNLGGK", "A") == {}


# --------------------------------------------------------------------------
# Scoring behaviour
# --------------------------------------------------------------------------

def _residue(**kwargs) -> Residue:
    defaults = dict(
        chain="D", number=100, insertion=" ", aa="S", index=99,
        loop="A-CDR3", sasa_complex=10.0, sasa_free=60.0,
        contacts_peptide=10, contacts_mhc=2, void_slots=3,
        growth_reaches_partner=True, partner_charge=0.0,
        partner_polar_atoms=3, backbone_nh_bonded=False, phi=-60.0,
    )
    defaults.update(kwargs)
    return Residue(**defaults)


def test_growth_into_open_space_beats_growth_into_a_crowded_pocket():
    roomy = score_mutation(_residue(void_slots=3), "W")
    crowded = score_mutation(_residue(void_slots=0), "W")
    assert roomy.packing > crowded.packing
    assert crowded.packing < 0


def test_growth_pointing_at_solvent_is_discounted():
    toward_partner = score_mutation(_residue(growth_reaches_partner=True), "F")
    toward_solvent = score_mutation(_residue(growth_reaches_partner=False), "F")
    assert toward_partner.packing > toward_solvent.packing


def test_apolar_burial_only_pays_off_at_a_buried_position():
    buried = score_mutation(_residue(sasa_complex=2.0, sasa_free=60.0), "I")
    exposed = score_mutation(_residue(sasa_complex=58.0, sasa_free=60.0), "I")
    assert buried.apolar_burial > exposed.apolar_burial


def test_charge_complementarity_is_rewarded_and_repulsion_penalised():
    near_negative = _residue(partner_charge=-0.8)
    assert score_mutation(near_negative, "R").electrostatic > 0
    assert score_mutation(near_negative, "D").electrostatic < 0


def test_burying_a_charge_with_no_counter_charge_is_penalised():
    lonely = score_mutation(
        _residue(partner_charge=0.0, sasa_complex=1.0, sasa_free=60.0), "R"
    )
    assert lonely.electrostatic < 0


def test_proline_is_penalised_and_more_so_when_the_amide_is_bonded():
    free = score_mutation(_residue(backbone_nh_bonded=False), "P")
    bonded = score_mutation(_residue(backbone_nh_bonded=True), "P")
    assert free.conformational < 0
    assert bonded.conformational < free.conformational


def test_mutating_away_from_left_handed_glycine_is_penalised():
    left_handed = score_mutation(_residue(aa="G", phi=70.0), "L")
    normal = score_mutation(_residue(aa="G", phi=-70.0), "L")
    assert left_handed.conformational < normal.conformational


def test_specificity_flag_follows_the_contact_split():
    peptide = score_mutation(_residue(contacts_peptide=10, contacts_mhc=1), "A")
    leaning = score_mutation(_residue(contacts_peptide=4, contacts_mhc=6), "A")
    mhc_only = score_mutation(_residue(contacts_peptide=0, contacts_mhc=10), "A")
    assert peptide.specificity_flag == "peptide-focused"
    assert leaning.specificity_flag == "mhc-leaning"
    assert mhc_only.specificity_flag == "mhc-dominated"


def test_a_mutation_keeps_its_identity_in_its_name():
    mutation = score_mutation(_residue(aa="S", number=100, chain="D"), "W")
    assert mutation.name == "D:S100W"


# --------------------------------------------------------------------------
# Liabilities
# --------------------------------------------------------------------------

def test_mutate_sequence_replaces_exactly_one_position():
    assert mutate_sequence("ACDEF", 2, "W") == "ACWEF"


class _FakeComplex:
    """Minimal stand-in exposing the attributes the liability check needs."""

    def __init__(self, sequence: str, chain: str = "D"):
        self.sequences = {chain: sequence}
        self.loops = {"A-CDR3": (chain, 0, len(sequence))}
        self._chain = chain

    def residue_by_label(self, chain, number):
        return Residue(
            chain=chain, number=number, insertion=" ",
            aa=self.sequences[chain][number], index=number,
        )


def test_a_new_glycosylation_sequon_is_detected():
    # Introducing N at index 2 creates N-A-T.
    complex_obj = _FakeComplex("GGSATGG")
    mutation = score_mutation(_residue(aa="S", number=2, chain="D"), "N")
    assert "new N-glycosylation sequon" in assess_liabilities(complex_obj, mutation)


def test_a_pre_existing_sequon_is_not_blamed_on_the_mutation():
    complex_obj = _FakeComplex("NATGGGG")
    mutation = score_mutation(_residue(aa="G", number=4, chain="D"), "A")
    assert "new N-glycosylation sequon" not in assess_liabilities(complex_obj, mutation)


def test_introducing_cysteine_is_always_flagged():
    complex_obj = _FakeComplex("GGSGGGG")
    mutation = score_mutation(_residue(aa="S", number=2, chain="D"), "C")
    assert "new free cysteine" in assess_liabilities(complex_obj, mutation)


def test_proline_after_aspartate_is_flagged_as_isomerisation_risk():
    complex_obj = _FakeComplex("GGDSGGG")
    mutation = score_mutation(_residue(aa="S", number=3, chain="D"), "P")
    assert any("isomerisation" in item for item in assess_liabilities(complex_obj, mutation))


# --------------------------------------------------------------------------
# Off-target similarity
# --------------------------------------------------------------------------

def test_identical_peptide_scores_one():
    import numpy as np

    weights = np.array([1.0, 0, 5.0, 0, 2.0])
    similarity, identical, pattern = _similarity("ACDEF", "ACDEF", weights)
    assert similarity == pytest.approx(1.0)
    assert identical == 3
    assert pattern == "A.D.F"


def test_unweighted_positions_do_not_affect_similarity():
    import numpy as np

    weights = np.array([1.0, 0.0, 1.0])
    same = _similarity("ACD", "ACD", weights)[0]
    changed_only_at_unweighted = _similarity("ACD", "AWD", weights)[0]
    assert same == pytest.approx(changed_only_at_unweighted)


def test_conservative_substitution_scores_above_a_non_conservative_one():
    import numpy as np

    weights = np.array([1.0])
    conservative = _similarity("I", "L", weights)[0]
    radical = _similarity("I", "D", weights)[0]
    assert conservative > radical


def test_anchor_convention_covers_both_ends_of_the_peptide():
    assert 2 in HLA_A0201_ANCHORS
    assert -1 in HLA_A0201_ANCHORS
    # P2 and the C-terminal anchor both prefer small hydrophobics.
    assert "L" in HLA_A0201_ANCHORS[2]
    assert "V" in HLA_A0201_ANCHORS[-1]


# --------------------------------------------------------------------------
# Validation helpers
# --------------------------------------------------------------------------

def test_rank_averages_ties():
    assert _rank([1.0, 2.0, 2.0, 3.0]) == [1.0, 2.5, 2.5, 4.0]


def test_spearman_recovers_perfect_monotonic_relationships():
    assert _spearman([1, 2, 3, 4], [10, 20, 30, 40]) == pytest.approx(1.0)
    assert _spearman([1, 2, 3, 4], [40, 30, 20, 10]) == pytest.approx(-1.0)


def test_spearman_needs_at_least_three_points():
    assert _spearman([1, 2], [1, 2]) is None


def test_locate_in_scan_finds_a_mutation_and_reports_its_rank():
    mutations = [
        score_mutation(_residue(aa="S", number=100), aa) for aa in "WAIG"
    ]
    found = locate_in_scan(mutations, "D", 100, "W")
    assert found is not None
    assert found.name == "D:S100W"
    assert 1 <= found.rank <= len(mutations)
    assert found.position_total == len(mutations)


def test_locate_in_scan_returns_none_for_an_absent_mutation():
    mutations = [score_mutation(_residue(aa="S", number=100), "W")]
    assert locate_in_scan(mutations, "D", 999, "W") is None


# --------------------------------------------------------------------------
# Integration: these need cached structures
# --------------------------------------------------------------------------

requires_pdb = pytest.mark.skipif(
    not (PDB_CACHE / "1ao7.pdb").exists(),
    reason="1AO7 not cached; run the pipeline once with network access first",
)


@requires_pdb
def test_1ao7_chain_typing_and_interface_geometry():
    from tcr_design.structure import load_complex

    complex_obj = load_complex("1AO7", PDB_CACHE)
    assert complex_obj.chains["peptide"] == "C"
    assert complex_obj.peptide == "LLFGYPVYV"
    assert complex_obj.chains["tcr_alpha"] == "D"
    assert complex_obj.chains["tcr_beta"] == "E"

    paratope = complex_obj.paratope()
    assert 15 <= len(paratope) <= 40

    buried = sum(r.sasa_free - r.sasa_complex for r in paratope)
    # TCR-pMHC interfaces bury roughly 600-1000 A^2 per side.
    assert 500 < buried < 1500


@requires_pdb
def test_1ao7_cdr3_loops_are_recovered_from_the_structure():
    from tcr_design.structure import load_complex

    complex_obj = load_complex("1AO7", PDB_CACHE)
    sequences = complex_obj.sequences
    chain, start, end = complex_obj.loops["A-CDR3"]
    assert sequences[chain][start:end] == "CAVTTDSWGKLQF"
    chain, start, end = complex_obj.loops["B-CDR3"]
    assert sequences[chain][start:end] == "CASRPGLAGGRPEQYF"


@requires_pdb
def test_tax_peptide_tyrosine_five_is_the_dominant_tcr_contact():
    from tcr_design.offtarget import peptide_contact_profile
    from tcr_design.structure import load_complex

    complex_obj = load_complex("1AO7", PDB_CACHE)
    profile = peptide_contact_profile(complex_obj)
    assert len(profile) == 9
    # Tyr5 is the well-described critical TCR contact of the Tax epitope, and
    # the C-terminal valine is an MHC anchor the TCR should barely touch.
    assert profile.index(max(profile)) == 4
    assert profile[8] < profile[4]


@requires_pdb
def test_scan_is_deterministic_and_gates_every_mutation():
    from tcr_design.liability import apply_gates
    from tcr_design.scan import scan_complex
    from tcr_design.structure import load_complex

    complex_obj = load_complex("1AO7", PDB_CACHE)
    first = apply_gates(complex_obj, scan_complex(complex_obj))
    second = apply_gates(complex_obj, scan_complex(complex_obj))

    assert [m.name for m in first] == [m.name for m in second]
    assert [m.triage_score for m in first] == [m.triage_score for m in second]
    assert all(m.verdict != "pending" for m in first)
    assert all(not math.isnan(m.triage_score) for m in first)


@requires_pdb
def test_disulfide_cysteines_are_never_scanned():
    from tcr_design.scan import scan_complex
    from tcr_design.structure import load_complex

    complex_obj = load_complex("1AO7", PDB_CACHE)
    mutations = scan_complex(complex_obj, min_contacts=1)
    disulfides = {
        (r.chain, r.number) for r in complex_obj.residues if r.in_disulfide
    }
    assert not disulfides & {(m.chain, m.number) for m in mutations}


@pytest.mark.skipif(
    not (PDB_CACHE / "3qdg.pdb").exists() or not (PDB_CACHE / "4l3e.pdb").exists(),
    reason="DMF5 structure pair not cached",
)
def test_dmf5_structure_pair_recovers_the_two_designed_mutations():
    from tcr_design.structure import load_complex
    from tcr_design.validate import diff_structures

    wild_type = load_complex("3QDG", PDB_CACHE)
    variant = load_complex("4L3E", PDB_CACHE)
    differences = diff_structures(wild_type, variant)
    tcr = {(wt, number, var) for role, _c, number, wt, var, _l in differences
           if role.startswith("tcr")}
    # The designed high-affinity DMF5 variant differs by alpha D26Y and
    # beta L98W.
    assert ("D", 26, "Y") in tcr
    assert ("L", 98, "W") in tcr


@pytest.mark.skipif(not PROTEOME.exists(), reason="human proteome not cached")
@requires_pdb
def test_proteome_screen_requires_identity_at_the_hottest_position():
    from tcr_design.offtarget import scan_proteome
    from tcr_design.structure import load_complex

    complex_obj = load_complex("1AO7", PDB_CACHE)
    hits, motif = scan_proteome(complex_obj, PROTEOME, top_n=10)

    assert motif["identity_required_at"] == ["P5Y"]
    # Every reported hit must carry the required tyrosine and satisfy the
    # anchor filter, and none may be the target peptide itself.
    for hit in hits:
        assert hit.peptide[4] == "Y"
        assert hit.peptide != complex_obj.peptide
        assert hit.peptide[1] in HLA_A0201_ANCHORS[2]
        assert hit.peptide[-1] in HLA_A0201_ANCHORS[-1]


# --------------------------------------------------------------------------
# Benchmark / calibration
# --------------------------------------------------------------------------

def test_ddg_sign_convention_is_positive_for_weaker_binding():
    from tcr_design.benchmark import _ddg_from_affinities

    # A mutant with a ten-fold larger Kd binds ten-fold more weakly.
    weaker = _ddg_from_affinities("1e-5", "1e-6", 298.0)
    tighter = _ddg_from_affinities("1e-7", "1e-6", 298.0)
    assert weaker > 0
    assert tighter < 0
    assert weaker == pytest.approx(-tighter, rel=1e-6)
    # RT ln(10) at 298 K is about 1.36 kcal/mol.
    assert weaker == pytest.approx(1.364, abs=0.01)


def test_ddg_rejects_nonsense_affinities():
    from tcr_design.benchmark import _ddg_from_affinities

    assert _ddg_from_affinities("", "1e-6", 298.0) is None
    assert _ddg_from_affinities("0", "1e-6", 298.0) is None
    assert _ddg_from_affinities("abc", "1e-6", 298.0) is None


def test_atlas_measured_values_tolerate_error_bars():
    from tcr_design.benchmark import _parse_measured_number

    assert _parse_measured_number("0.19 ± 0.11") == pytest.approx(0.19)
    assert _parse_measured_number("1.2 +/- 0.3") == pytest.approx(1.2)
    assert _parse_measured_number("-0.04") == pytest.approx(-0.04)
    assert _parse_measured_number("N/A") is None
    assert _parse_measured_number("n.d.") is None


def test_mutation_strings_parse_in_both_database_formats():
    from tcr_design.benchmark import ATLAS_MUTATION, MUTATION_PATTERN

    # SKEMPI carries the chain; ATLAS keeps it in a separate column.
    wt, chain, number, insertion, mut = MUTATION_PATTERN.match("AE99S").groups()
    assert (wt, chain, number, insertion, mut) == ("A", "E", "99", "", "S")
    wt, number, insertion, mut = ATLAS_MUTATION.match("D26W").groups()
    assert (wt, number, insertion, mut) == ("D", "26", "", "W")


def test_repeat_measurements_are_averaged_not_double_counted():
    from tcr_design.benchmark import BenchmarkEntry, _average_duplicates

    def entry(ddg):
        return BenchmarkEntry(
            pdb_id="1AO7", chain="D", number=26, insertion="", wild_type="D",
            mutant="W", ddg=ddg, side="tcr", location="COR", source="test",
        )

    merged = _average_duplicates([entry(1.0), entry(2.0), entry(3.0)])
    assert len(merged) == 1
    assert merged[0].ddg == pytest.approx(2.0)
    assert "3 measurements" in merged[0].note


def test_rank_data_averages_ties():
    import numpy as np

    from tcr_design.benchmark import _rankdata

    assert list(_rankdata(np.array([1.0, 2.0, 2.0, 3.0]))) == [1.0, 2.5, 2.5, 4.0]


def test_roc_auc_is_one_for_perfect_separation_and_half_for_none():
    from tcr_design.benchmark import roc_auc

    assert roc_auc([3, 4, 1, 2], [True, True, False, False]) == pytest.approx(1.0)
    assert roc_auc([1, 2, 3, 4], [True, True, False, False]) == pytest.approx(0.0)
    assert roc_auc([1, 1, 1, 1], [True, True, False, False]) == pytest.approx(0.5)


def test_roc_auc_needs_both_classes():
    from tcr_design.benchmark import roc_auc

    assert roc_auc([1, 2, 3], [True, True, True]) is None


def test_ridge_recovers_a_known_linear_relationship():
    import numpy as np

    from tcr_design.benchmark import _fit_ridge

    rng = np.random.default_rng(0)
    X = rng.normal(size=(200, 3))
    truth = np.array([1.5, -2.0, 0.5])
    y = X @ truth + 0.25

    coefficients, intercept = _fit_ridge(X, y, ridge_lambda=1e-6)
    assert coefficients == pytest.approx(truth, abs=0.01)
    assert intercept == pytest.approx(0.25, abs=0.01)


def test_stronger_ridge_shrinks_coefficients():
    import numpy as np

    from tcr_design.benchmark import _fit_ridge

    rng = np.random.default_rng(1)
    X = rng.normal(size=(100, 3))
    y = X @ np.array([2.0, 2.0, 2.0])
    weak, _ = _fit_ridge(X, y, 1e-6)
    strong, _ = _fit_ridge(X, y, 1000.0)
    assert np.abs(strong).sum() < np.abs(weak).sum()


def test_spearman_and_pearson_agree_on_a_monotonic_relationship():
    from tcr_design.benchmark import pearson, spearman

    x = [1, 2, 3, 4, 5]
    y = [1, 4, 9, 16, 25]
    assert spearman(x, y) == pytest.approx(1.0)
    assert pearson(x, y) < 1.0  # monotonic but not linear


SKEMPI_CACHE = CACHE / "benchmark" / "skempi_v2.csv"
ATLAS_CACHE = CACHE / "benchmark" / "atlas_Mutants.csv"


@pytest.mark.skipif(not SKEMPI_CACHE.exists(), reason="SKEMPI not cached")
def test_skempi_loads_only_tcr_side_single_mutations():
    from tcr_design.benchmark import load_skempi

    entries = load_skempi(SKEMPI_CACHE)
    assert len(entries) > 100
    assert all(e.side == "tcr" for e in entries)
    assert all(e.source == "SKEMPI2" for e in entries)
    assert all(e.wild_type in AA20 and e.mutant in AA20 for e in entries)


@pytest.mark.skipif(not ATLAS_CACHE.exists(), reason="ATLAS not cached")
def test_atlas_loads_tcr_mutants_with_wild_type_peptide_and_mhc():
    from tcr_design.benchmark import load_atlas

    entries = load_atlas(ATLAS_CACHE)
    assert len(entries) > 100
    assert all(e.side == "tcr" for e in entries)
    assert all(e.source == "ATLAS" for e in entries)
    # The A6 CDR1-alpha aromatic substitution is in ATLAS and tightens binding.
    a6 = [
        e for e in entries
        if e.pdb_id == "1AO7" and e.number == 26 and e.mutant == "W"
    ]
    assert a6, "expected the A6 D26W mutation in ATLAS"
    assert a6[0].ddg < -0.5

    # Multi-point mutants are separated by '|' and must not leak in as singles.
    assert all(
        "|" not in e.wild_type and "|" not in e.mutant for e in entries
    )


@pytest.mark.skipif(
    not (SKEMPI_CACHE.exists() and ATLAS_CACHE.exists()),
    reason="both benchmark sources needed",
)
def test_pooling_deduplicates_overlapping_mutations():
    from tcr_design.benchmark import load_combined

    entries, provenance = load_combined(SKEMPI_CACHE, ATLAS_CACHE)
    assert provenance["overlapping_mutations"] > 0
    assert provenance["pooled_unique"] == len(entries)
    assert provenance["pooled_unique"] < (
        provenance["atlas_entries"] + provenance["skempi_entries"]
    )
    keys = [e.key for e in entries]
    assert len(keys) == len(set(keys))


# --------------------------------------------------------------------------
# Calibrated model
# --------------------------------------------------------------------------

MODEL_PATH = Path(__file__).resolve().parents[1] / "tcr_design" / "data" / "calibrated_ddg_model.json"


def test_ranking_key_prefers_predicted_ddg_when_available():
    from tcr_design.scan import ranking_key

    tighter = score_mutation(_residue(), "W")
    weaker = score_mutation(_residue(), "A")
    tighter.predicted_ddg, weaker.predicted_ddg = -1.5, +0.8
    # Sorting ascending by the key must put the tighter binder first.
    assert ranking_key(tighter) < ranking_key(weaker)


def test_ranking_key_falls_back_to_the_z_score_without_a_model():
    from tcr_design.scan import ranking_key

    good = score_mutation(_residue(), "W")
    bad = score_mutation(_residue(), "A")
    good.peptide_z, bad.peptide_z = 2.0, -1.0
    assert good.predicted_ddg is None
    assert ranking_key(good) < ranking_key(bad)


def test_gates_use_predicted_ddg_when_the_model_is_present():
    from tcr_design.liability import apply_gates

    complex_obj = _FakeComplex("GGSGGGG")

    def gated(predicted, z):
        mutation = score_mutation(_residue(aa="S", number=2, chain="D"), "W")
        mutation.predicted_ddg = predicted
        mutation.peptide_z = z
        return apply_gates(complex_obj, [mutation], ddg_threshold=-0.25)[0].verdict

    # A strongly negative predicted ddG advances even with a poor z-score,
    # and a positive one is held back even with a good z-score.
    assert gated(-1.0, -5.0) == "advance"
    assert gated(+1.0, +5.0) == "below-threshold"


@pytest.mark.skipif(not MODEL_PATH.exists(), reason="calibrated model not generated")
def test_shipped_model_loads_and_declares_its_performance():
    from tcr_design.calibrated import load_model

    model = load_model(MODEL_PATH)
    assert model is not None
    assert model.feature_set == "extended"
    held_out = model.performance["held_out"]
    # The model is weak but must beat predicting the mean, or it should not ship.
    assert held_out["rmse"] < held_out["rmse_predicting_the_mean"]
    assert held_out["spearman"] > 0.2
    assert "Spearman" in model.summary()


@pytest.mark.skipif(not MODEL_PATH.exists(), reason="calibrated model not generated")
@requires_pdb
def test_shipped_model_predicts_physically_sensible_directions():
    from tcr_design.calibrated import load_model
    from tcr_design.structure import load_complex

    model = load_model(MODEL_PATH)
    complex_obj = load_complex("1AO7", PDB_CACHE)
    # A large, well-buried, peptide-contacting residue: deleting it to alanine
    # should be predicted to weaken binding more than a conservative swap.
    target = max(
        complex_obj.paratope(),
        key=lambda r: r.contacts_peptide * r.buried_fraction,
    )
    to_alanine = score_mutation(target, "A")
    to_similar = score_mutation(target, target.aa if target.aa != "A" else "S")
    assert model.predict(target, to_alanine) > model.predict(target, to_similar)


@pytest.mark.skipif(not MODEL_PATH.exists(), reason="calibrated model not generated")
@requires_pdb
def test_scan_attaches_predicted_ddg_only_when_a_model_is_passed():
    from tcr_design.calibrated import load_model
    from tcr_design.scan import scan_complex
    from tcr_design.structure import load_complex

    complex_obj = load_complex("1AO7", PDB_CACHE)
    without = scan_complex(complex_obj, min_contacts=10)
    with_model = scan_complex(
        complex_obj, min_contacts=10, ddg_model=load_model(MODEL_PATH)
    )
    assert without and with_model
    assert all(m.predicted_ddg is None for m in without)
    assert all(m.predicted_ddg is not None for m in with_model)
    assert all(-15.0 < m.predicted_ddg < 15.0 for m in with_model)


# --------------------------------------------------------------------------
# Epistasis / additivity
# --------------------------------------------------------------------------

def test_mutation_set_splits_on_both_separators():
    from tcr_design.epistasis import _mutation_set

    assert _mutation_set("D26W | L98W") == {"D26W", "L98W"}
    assert _mutation_set("A99M.G100S") == {"A99M", "G100S"}
    assert _mutation_set("D26W") == {"D26W"}


def test_decomposition_error_sign_and_labels():
    from tcr_design.epistasis import Decomposition

    def made(sum_of_parts, measured):
        return Decomposition(
            tcr="X", peptide="P", study="1", whole=frozenset({"a", "b"}),
            parts=(frozenset({"a"}), frozenset({"b"})),
            sum_of_parts=sum_of_parts, measured=measured,
        )

    # Measured minus predicted: positive means the stack underdelivered.
    assert made(-2.0, -1.0).error == pytest.approx(1.0)
    assert made(-2.0, -1.0).behaviour == "antagonistic"
    assert made(-2.0, -3.0).behaviour == "synergistic"
    assert made(-2.0, -2.1).behaviour == "additive"


def test_decompositions_must_partition_the_whole_exactly():
    from tcr_design.epistasis import find_decompositions

    # {a,b} and {b,c} overlap, so {a,b,c} has no valid two-part decomposition
    # from them; {a} + {b,c} does partition it.
    index = {
        ("T", "P"): {
            frozenset({"a", "b", "c"}): -3.0,
            frozenset({"a", "b"}): -2.0,
            frozenset({"b", "c"}): -2.0,
            frozenset({"a"}): -1.0,
            frozenset({"b", "c"}): -2.0,
        }
    }
    found = find_decompositions(index)
    for item in found:
        union = set()
        for part in item.parts:
            assert not (union & part), "parts must be disjoint"
            union |= part
        assert frozenset(union) == item.whole


def test_decomposition_sums_the_parts_it_names():
    from tcr_design.epistasis import find_decompositions

    index = {
        ("T", "P"): {
            frozenset({"a", "b"}): -3.0,
            frozenset({"a"}): -1.0,
            frozenset({"b"}): -1.5,
        }
    }
    found = find_decompositions(index)
    assert len(found) == 1
    assert found[0].sum_of_parts == pytest.approx(-2.5)
    assert found[0].measured == pytest.approx(-3.0)
    assert found[0].error == pytest.approx(-0.5)


@pytest.mark.skipif(not ATLAS_CACHE.exists(), reason="ATLAS not cached")
def test_atlas_index_keeps_different_peptides_apart():
    from tcr_design.epistasis import index_atlas

    index = index_atlas(ATLAS_CACHE, by_study=False)
    # The same receptor against a different epitope must be a separate group,
    # or changing the target would masquerade as epistasis.
    dmf5 = [key for key in index if key[0] == "DMF5"]
    assert len({key[1] for key in dmf5}) > 1


@pytest.mark.skipif(not ATLAS_CACHE.exists(), reason="ATLAS not cached")
def test_additivity_holds_for_pairs_and_fails_for_deep_stacks():
    from tcr_design.epistasis import run_additivity_analysis

    results = run_additivity_analysis(ATLAS_CACHE)
    within = results["within_study"]["summary"]["by_stack_size"]
    assert within[2]["n"] >= 3
    # Pairs add up.
    assert abs(within[2]["mean_error"]) < 0.5
    assert within[2]["within_half_kcal"] == within[2]["n"]
    # Deep stacks systematically underdeliver against the sum of their parts.
    deep = [s for size, s in within.items() if size >= 8]
    assert deep, "expected some deeply stacked variants in ATLAS"
    assert max(s["mean_error"] for s in deep) > 1.0
