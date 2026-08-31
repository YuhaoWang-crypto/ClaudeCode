"""The design-count / job-metadata / deviations ledger tree."""

from datetime import datetime, timedelta, timezone

import pytest

from binder_campaign.ledger import (
    METHOD_FLOOR,
    Deviation,
    DesignCountLedger,
    JobMetadataRow,
    LedgerRow,
)


def row(**kw):
    base = dict(job_id="j1", target="PD-L1", structure_method="rfdiffusion",
                stage="gen_screen", n_generated=10, n_scored=8, gpu_seconds=100.0,
                writer_frame_id="f1")
    base.update(kw)
    return LedgerRow(**base)


def test_stage_must_come_from_the_fixed_vocabulary():
    with pytest.raises(ValueError, match="unknown stage"):
        row(stage="scoring")


def test_writes_are_idempotent_on_job_id_and_stage(tmp_path):
    led = DesignCountLedger(str(tmp_path), "f1")
    assert led.append(row()) is True
    assert led.append(row()) is False           # same (job_id, stage)
    assert led.append(row(stage="final")) is True  # different stage
    assert led.totals().designs_generated == 20


def test_a_reopened_handle_still_refuses_a_duplicate(tmp_path):
    DesignCountLedger(str(tmp_path), "f1").append(row())
    assert DesignCountLedger(str(tmp_path), "f1").append(row()) is False


def test_one_file_per_writing_frame(tmp_path):
    DesignCountLedger(str(tmp_path), "f1").append(row(job_id="a"))
    DesignCountLedger(str(tmp_path), "f2").append(row(job_id="b"))
    names = sorted(p.name for p in tmp_path.glob("*.jsonl"))
    assert names == ["f1.jsonl", "f2.jsonl"]


def test_totals_use_the_kickoff_definitions(tmp_path):
    led = DesignCountLedger(str(tmp_path), "f1")
    led.append(row(job_id="g", stage="generate", n_generated=100, n_scored=0))
    led.append(row(job_id="s", stage="screen", n_generated=0, n_scored=90))
    led.append(row(job_id="i", stage="intermediate", n_generated=0, n_scored=40))
    led.append(row(job_id="f", stage="final", n_generated=0, n_scored=30))

    t = led.totals()
    assert t.designs_generated == 100
    assert t.designs_screened == 90    # intermediate does NOT count
    assert t.designs_ranked == 30


def test_gen_screen_counts_toward_both_halves(tmp_path):
    led = DesignCountLedger(str(tmp_path), "f1")
    led.append(row(job_id="gs", stage="gen_screen", n_generated=50, n_scored=45))
    t = led.totals()
    assert t.designs_generated == 50
    assert t.designs_screened == 45


def test_ledger_agg_is_what_submit_gate_step_d2_reads(tmp_path):
    led = DesignCountLedger(str(tmp_path), "f1")
    led.append(row(job_id="a", target="TNFa", n_generated=300, n_scored=20))
    agg = led.totals().ledger_agg()
    assert agg["TNFa"] == {"gen": 300, "scr": 20}


def test_floor_matrix_flags_every_cell_below_fifty(tmp_path):
    led = DesignCountLedger(str(tmp_path), "f1")
    led.append(row(job_id="a", structure_method="rfdiffusion", n_generated=60))
    led.append(row(job_id="b", structure_method="boltzgen", n_generated=10))

    m = led.totals().floor_matrix(["PD-L1"], ["rfdiffusion", "boltzgen", "genie3"])
    assert m["floor"] == METHOD_FLOOR
    assert m["matrix"]["PD-L1"]["rfdiffusion"] == 60
    open_cells = {(c["structure_method"], c["shortfall"]) for c in m["open_obligations"]}
    assert ("boltzgen", 40) in open_cells
    assert ("genie3", 50) in open_cells   # a method that contributed zero
    assert not any(c["structure_method"] == "rfdiffusion"
                   for c in m["open_obligations"])


def test_ratecard_is_computed_from_the_job_metadata_ledger(tmp_path):
    led = DesignCountLedger(str(tmp_path), "f1")
    t0 = datetime(2026, 6, 1, tzinfo=timezone.utc)
    led.append_job_metadata(JobMetadataRow(
        job_id="j1", campaign="c", target="PD-L1", launching_frame_id="f1",
        gpu_type="H100", est_hourly_usd=4.39,
        submitted_at_utc=t0.isoformat(),
        terminal_at_utc=(t0 + timedelta(hours=2)).isoformat(),
    ))
    assert led.ratecard_usd(now=t0 + timedelta(hours=5)) == pytest.approx(8.78)


def test_a_still_running_job_bills_to_now(tmp_path):
    led = DesignCountLedger(str(tmp_path), "f1")
    t0 = datetime(2026, 6, 1, tzinfo=timezone.utc)
    led.append_job_metadata(JobMetadataRow(
        job_id="j1", campaign="c", target="PD-L1", launching_frame_id="f1",
        gpu_type="H100", est_hourly_usd=4.0, submitted_at_utc=t0.isoformat(),
    ))
    assert led.ratecard_usd(now=t0 + timedelta(hours=3)) == pytest.approx(12.0)


def test_deviations_are_append_only_and_deduplicated(tmp_path):
    led = DesignCountLedger(str(tmp_path), "f1")
    dev = Deviation(kind="relaxation", target="PD-L1", what="caps relaxed",
                    action="disclosed", at_utc="2026-06-01T00:00:00+00:00")
    led.append_deviation(dev)
    led.append_deviation(dev)
    rows = led.deviations()
    assert len(rows) == 1
    assert set(rows[0]) == {"at_utc", "kind", "target", "what", "action"}


def test_a_read_only_handle_refuses_to_write(tmp_path):
    led = DesignCountLedger(str(tmp_path))
    with pytest.raises(RuntimeError, match="read-only"):
        led.append(row())
    assert led.totals().designs_generated == 0
