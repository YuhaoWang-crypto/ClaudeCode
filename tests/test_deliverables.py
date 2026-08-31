"""Companion files, the scoreboard, the validation gate, and the end-to-end run."""

import json

import pandas as pd
import pytest

from binder_campaign.companions import (
    assert_companion_coverage,
    build_instrument_realization,
    build_per_seed_metrics,
    recompute_final_scores_from_companion,
)
from binder_campaign.demo import run_demo
from binder_campaign.gates import (
    InstrumentGateRow,
    ScoringConstruct,
    ValidationGate,
    gate_is_pass,
    write_gate,
)
from binder_campaign.ledger import LedgerTotals
from binder_campaign.schema import default_method_vocab, default_sheet_schema
from binder_campaign.scoreboard import (
    SCOREBOARD_COLUMNS,
    build_scoreboard,
    scoreboard_gaps,
)
from binder_campaign.scoring import (
    DEFAULT_ARMS,
    DesignScore,
    InstrumentMask,
    SeedRecord,
    aggregate_arm,
)

MASK = InstrumentMask()


# --------------------------------------------------------------------------- #
# validation gate
# --------------------------------------------------------------------------- #


def _rows(ca_rmsd=1.0, **kw):
    base = dict(ca_rmsd=ca_rmsd, ca_rmsd_threshold=2.0,
                control_name="natural ligand", control_score=0.6,
                negative_control_scores=(0.2, 0.15))
    base.update(kw)
    return [InstrumentGateRow(arm=a, **base) for a in DEFAULT_ARMS]


def _gate(**kw):
    base = dict(target="PD-L1", instruments=_rows(),
                construct=ScoringConstruct(chains=("A",), residue_range="18-134"))
    base.update(kw)
    return ValidationGate(**base)


def test_gate_passes_on_fold_recapitulation_plus_control_separation():
    g = _gate()
    assert g.status == "PASS"
    assert g.instruments[0].separation_value == pytest.approx(0.4)


def test_gate_fails_when_the_fold_is_not_recapitulated():
    g = _gate(instruments=_rows(ca_rmsd=5.0))
    assert g.status == "FAIL"
    assert "fold recapitulation" in g.blockers()[0]


def test_separation_failure_bars_self_authorisation_until_remedies_are_tried():
    failing = _rows(control_score=0.1, negative_control_scores=(0.5,))
    g = _gate(instruments=failing)
    assert g.status == "FAIL"
    assert "remedies not yet tried" in g.blockers()[0]

    tried = _rows(control_score=0.1, negative_control_scores=(0.5,),
                  remedies_attempted=(
                      "rerun_control_with_control_ligand_msa",
                      "alternative_cofolder",
                      "af_unmasked_template_injection"))
    assert "after all three remedies" in _gate(instruments=tried).blockers()[0]


def test_no_literature_control_drops_condition_b():
    g = _gate(instruments=_rows(control_name=None, control_score=None,
                                negative_control_scores=()),
              no_literature_control=True)
    assert g.status == "PASS"


def test_an_antibody_only_control_that_fails_still_passes_on_fold_alone():
    g = _gate(instruments=_rows(control_score=0.1, negative_control_scores=(0.5,),
                                control_is_antibody=True))
    assert g.status == "PASS"


def test_gdf8_must_record_a_gdf11_arm():
    g = _gate(target="GDF-8", counter_target="GDF-11")
    assert g.status == "FAIL"
    assert "must record that arm" in g.blockers()[0]

    g.counter_target_instruments = _rows()
    assert g.status == "PASS"


def test_gate_file_round_trips_and_is_what_submit_gate_reads(tmp_path):
    write_gate(str(tmp_path), _gate())
    assert gate_is_pass(str(tmp_path), "PD-L1")
    assert not gate_is_pass(str(tmp_path), "TNFa")   # absent file is not a pass
    write_gate(str(tmp_path), _gate(target="TNFa", instruments=_rows(ca_rmsd=9.0)))
    assert not gate_is_pass(str(tmp_path), "TNFa")   # FAIL is not a pass


def test_scoring_construct_asserts_chain_count_and_cofactor_atoms():
    c = ScoringConstruct(chains=("A", "B"), residue_range="1-100",
                         cofactors=("ZN",), n_target_chains=2)
    assert c.matches(n_chains_folded=3, cofactor_atoms=3)      # 2 target + binder
    assert not c.matches(n_chains_folded=2, cofactor_atoms=3)  # a monomer crop
    assert not c.matches(n_chains_folded=3, cofactor_atoms=0)  # cofactor-blind


# --------------------------------------------------------------------------- #
# companions
# --------------------------------------------------------------------------- #


def _design(design_id, target="PD-L1", ipsae=0.5, dockq=0.5, n_seeds=5, off=None):
    on = {
        a: aggregate_arm(a, [
            SeedRecord(seed=100 * i + k, ipsae_ab=ipsae, ipsae_ba=ipsae, dockq=dockq)
            for k in range(n_seeds)
        ])
        for i, a in enumerate(DEFAULT_ARMS)
    }
    off_book = {}
    if off is not None:
        off_book = {
            a: aggregate_arm(a, [
                SeedRecord(seed=100 * i + k, ipsae_ab=off, ipsae_ba=off, dockq=0.3)
                for k in range(n_seeds)
            ])
            for i, a in enumerate(DEFAULT_ARMS)
        }
    return DesignScore(design_id, target, on, off_book)


def _sheet(designs, mask=MASK):
    return pd.DataFrame([{
        "design_id": d.design_id,
        "target": d.target,
        "final_score": d.final_score(mask),
        "rank": i + 1,
    } for i, d in enumerate(designs)])


def test_companion_covers_every_ranked_design_id():
    designs = [_design(f"d{i}") for i in range(4)]
    comp = build_per_seed_metrics(designs, MASK)
    assert_companion_coverage(_sheet(designs), comp, MASK)


def test_incomplete_companion_coverage_is_a_deliverable_defect():
    designs = [_design(f"d{i}") for i in range(4)]
    comp = build_per_seed_metrics(designs[:2], MASK)
    with pytest.raises(AssertionError, match="deliverable defect"):
        assert_companion_coverage(_sheet(designs), comp, MASK)


def test_every_sheet_score_reproduces_from_the_companion_alone():
    designs = [_design(f"d{i}", ipsae=0.3 + 0.1 * i, dockq=0.4 + 0.05 * i)
               for i in range(5)]
    frame = recompute_final_scores_from_companion(
        _sheet(designs), build_per_seed_metrics(designs, MASK), MASK
    )
    assert (frame["abs_diff"] < 1e-4).all()


def test_a_corrupted_sheet_score_fails_the_companion_recompute():
    designs = [_design(f"d{i}") for i in range(3)]
    sheet = _sheet(designs)
    sheet.loc[0, "final_score"] = 0.99
    with pytest.raises(AssertionError, match="do not reproduce from the companion"):
        recompute_final_scores_from_companion(
            sheet, build_per_seed_metrics(designs, MASK), MASK
        )


def test_counter_screened_companion_must_be_seed_matched():
    mask = InstrumentMask(counter_screened=True)
    designs = [_design(f"d{i}", target="GDF-8", off=0.2) for i in range(3)]
    comp = build_per_seed_metrics(designs, mask)
    assert_companion_coverage(_sheet(designs, mask), comp, mask)

    # drop the off-target arm entirely
    comp_on_only = comp[comp["side"] == "on_target"]
    with pytest.raises(AssertionError, match="no off-target arm"):
        assert_companion_coverage(_sheet(designs, mask), comp_on_only, mask)


def test_instrument_realization_is_derived_and_recounted(tmp_path):
    designs = [_design(f"d{i}") for i in range(4)]
    sheet = _sheet(designs)
    comp = build_per_seed_metrics(designs, MASK)
    gate = _gate()
    gates = {"PD-L1": json.loads(json.dumps(gate.to_json()))}

    frame = build_instrument_realization(sheet, comp, gates, {"PD-L1": MASK})
    assert set(frame["arm_name"]) == set(DEFAULT_ARMS)
    assert (frame["gate_status"] == "PASS").all()
    assert (frame["n_ranked_with_score"] == 4).all()
    assert (frame["used_in_final_score"]).all()
    assert frame["control_separation_value"].notna().all()


def test_a_dropped_arm_must_carry_a_drop_reason():
    designs = [_design(f"d{i}") for i in range(3)]
    reduced = InstrumentMask(arms=("ef2full", "ptxv2"))
    comp = build_per_seed_metrics(designs, MASK)  # all three arms were run
    gates = {"PD-L1": json.loads(json.dumps(_gate().to_json()))}

    frame = build_instrument_realization(_sheet(designs), comp, gates,
                                         {"PD-L1": reduced})
    dropped = frame[~frame["used_in_final_score"]]
    assert len(dropped) == 1
    assert dropped.iloc[0]["arm_name"] == "ef2fast"
    assert dropped.iloc[0]["drop_reason"]


def test_no_literature_control_writes_NA_and_names_the_path():
    designs = [_design(f"d{i}") for i in range(3)]
    gate = _gate(instruments=_rows(control_name=None, control_score=None,
                                   negative_control_scores=()),
                 no_literature_control=True)
    gates = {"PD-L1": json.loads(json.dumps(gate.to_json()))}
    frame = build_instrument_realization(
        _sheet(designs), build_per_seed_metrics(designs, MASK), gates,
        {"PD-L1": MASK}
    )
    assert (frame["control_separation_value"] == "NA").all()
    assert frame["drop_reason"].str.contains("no_literature_control").all()


# --------------------------------------------------------------------------- #
# scoreboard
# --------------------------------------------------------------------------- #


def _ranked(n=5):
    return [{
        "design_id": f"d{i}", "sequence": "ACDEFGHIKLMNPQRSTVWY" * (2 + i % 2),
        "final_score": 0.5 + 0.01 * i, "pose_PASS": i % 2 == 0,
        "root_backbone_id": f"bb{i % 3}", "tm90_cluster_id": f"c{i % 2}",
        "structure_method": "rfdiffusion" if i % 2 else "boltzgen",
        "fold_class": "not_all_alpha" if i == 0 else "all_alpha",
    } for i in range(n)]


def test_scoreboard_columns_are_exactly_the_frozen_set():
    totals = LedgerTotals(0, 0, 0, 0.0,
                          {"PD-L1": {"designs_generated": 100,
                                     "designs_screened": 90,
                                     "designs_ranked": 30}}, {})
    board = build_scoreboard({"PD-L1": _ranked()}, totals,
                             {"PD-L1": 1200.0}, {"PD-L1": 6.0})
    assert list(board.columns) == list(SCOREBOARD_COLUMNS)


def test_every_non_target_scoreboard_column_is_a_real_numeric():
    totals = LedgerTotals(0, 0, 0, 0.0, {"PD-L1": {}}, {})
    board = build_scoreboard({"PD-L1": _ranked()}, totals,
                             {"PD-L1": 1200.0}, {"PD-L1": 6.0})
    row = board.iloc[0]
    for col in SCOREBOARD_COLUMNS[1:]:
        assert pd.api.types.is_number(row[col]) and pd.notna(row[col])
    assert scoreboard_gaps(board, "f1").empty


def test_diversity_columns_come_from_the_ranked_set():
    totals = LedgerTotals(0, 0, 0, 0.0, {"PD-L1": {}}, {})
    board = build_scoreboard({"PD-L1": _ranked(6)}, totals,
                             {"PD-L1": 0.0}, {"PD-L1": 0.0})
    row = board.iloc[0]
    assert row["n_distinct_root_backbones"] == 3
    assert row["max_seqs_per_root_backbone"] == 2
    assert row["n_tm90_clusters"] == 2
    assert row["n_structure_methods"] == 2
    assert row["top_method_share"] == pytest.approx(0.5)
    assert row["n_non_all_alpha"] == 1


# --------------------------------------------------------------------------- #
# schema
# --------------------------------------------------------------------------- #


def test_method_vocab_canonicalises_aliases_and_rejects_unknown_tokens():
    v = default_method_vocab()
    assert v.canonical_structure_method("RFdiffusion") == "rfdiffusion"
    assert v.canonical_seq_method("SolubleMPNN") == "solublempnn"
    with pytest.raises(ValueError, match="not in the frozen enum"):
        v.canonical_structure_method("secret_model_v9")


def _schema_row(schema):
    row = {c: "x" for c in schema.columns}
    row["structure_method"] = "rfdiffusion"
    row["seq_method"] = "solublempnn"
    row["fold_class"] = "all_alpha"
    return row


def test_schema_refuses_renamed_missing_or_extra_columns():
    schema = default_sheet_schema()
    row = _schema_row(schema)
    schema.validate_row(row, ranked=False)
    with pytest.raises(ValueError, match="unknown sheet columns"):
        schema.validate_row({**row, "surprise": 1}, ranked=False)
    with pytest.raises(ValueError, match="not in frozen enum"):
        schema.validate_row({**row, "structure_method": "secret_v9"}, ranked=False)
    del row["rank"]
    with pytest.raises(ValueError, match="missing sheet columns"):
        schema.validate_row(row, ranked=False)


def test_mandatory_nonnull_is_fail_loud_on_a_ranked_row():
    schema = default_sheet_schema()
    row = _schema_row(schema)
    row["pose_dockq"] = None
    with pytest.raises(ValueError, match="mandatory_nonnull column 'pose_dockq'"):
        schema.validate_row(row, ranked=True)


def test_schema_has_one_predicted_structure_path_per_arm():
    schema = default_sheet_schema()
    for arm in DEFAULT_ARMS:
        assert f"predicted_structure_path_{arm}" in schema.columns
        assert f"ipsae_{arm}" in schema.columns
        assert f"sc_DockQ_{arm}" in schema.columns


# --------------------------------------------------------------------------- #
# end to end
# --------------------------------------------------------------------------- #


def test_demo_produces_every_deliverable(tmp_path):
    summary = run_demo(str(tmp_path))

    for name in ("design_sheet.csv", "per_seed_metrics.parquet",
                 "instrument_realization.csv", "scoreboard.csv",
                 "sheet_schema.json", "method_vocab.json", "deviations.jsonl"):
        assert (tmp_path / name).exists(), name
    assert (tmp_path / "state" / "governor.json").exists()
    assert (tmp_path / "state" / "gates" / "PD-L1.json").exists()
    assert (tmp_path / "state" / "ledger_agg.json").exists()

    sheet = pd.read_csv(tmp_path / "design_sheet.csv")
    # exactly 30 ranked rows per target, ranks contiguous
    for target, grp in sheet.groupby("target"):
        assert len(grp) == 30
        assert sorted(grp["rank"]) == list(range(1, 31))
        assert grp["pose_dockq"].notna().all()
        assert grp["final_score"].notna().all()
        assert grp["sequence"].nunique() == 30

    assert summary["even_pace_ceiling"] == 237
    assert summary["n_scoreboard_gaps"] == 0

    board = pd.read_csv(tmp_path / "scoreboard.csv")
    assert list(board.columns) == list(SCOREBOARD_COLUMNS)
    assert len(board) == 2


def test_demo_is_deterministic(tmp_path):
    a = run_demo(str(tmp_path / "a"), seed=3)
    b = run_demo(str(tmp_path / "b"), seed=3)
    assert a["per_target"] == b["per_target"]
    assert (pd.read_csv(tmp_path / "a" / "design_sheet.csv")["design_id"].tolist()
            == pd.read_csv(tmp_path / "b" / "design_sheet.csv")["design_id"].tolist())
