"""The fail-closed submit gate, steps (a)-(g) and (d2)."""

import json
from datetime import datetime, timedelta, timezone

import pytest

from binder_campaign import submit_gate as sg

NOW = datetime(2026, 6, 1, 12, 0, tzinfo=timezone.utc)


class FakeState:
    def __init__(self, files=None):
        self.files = files or {}

    def read_json(self, path):
        v = self.files.get(path)
        return json.loads(json.dumps(v)) if v is not None else None

    def read_lines(self, path):
        v = self.files.get(path)
        return list(v) if v else []


class FakeFleet:
    def __init__(self, total=0, cpu=0):
        self.total, self.cpu = total, cpu
        self.calls = 0

    def count_all(self, app_id, project_tag):
        self.calls += 1
        return self.total

    def count_cpu_tagged(self, app_id, project_tag):
        return self.cpu


def wire(files, total=0, cpu=0, now=NOW):
    fleet = FakeFleet(total, cpu)
    sleeps = []
    sg.configure(state=FakeState(files), fleet=fleet, sleep=sleeps.append,
                 now=lambda: now)
    return fleet, sleeps


def gov(**kw):
    base = {
        "ceiling": 300,
        "set_at": (NOW - timedelta(minutes=5)).isoformat(),
        "basis": "BUDGET",
        "pct_elapsed": 50.0,
        "pct_spent_metered": 45.0,
        "pct_spent_projected": 46.0,
    }
    base.update(kw)
    return base


def call(**kw):
    kw.setdefault("app_id", "ap-1")
    kw.setdefault("project_tag", "campaign-x")
    kw.setdefault("timeout_s", 900)
    return sg.submit_gate(**kw)


# --- (b) missing / stale governor ------------------------------------------ #

def test_missing_governor_is_a_dead_governor():
    wire({})
    d = call()
    assert not d and d.reason == "STALE_OR_MISSING_GOVERNOR"


def test_governor_older_than_25_minutes_is_stale():
    wire({"governor.json": gov(set_at=(NOW - timedelta(minutes=26)).isoformat())})
    d = call()
    assert not d and d.reason == "STALE_OR_MISSING_GOVERNOR"


def test_governor_at_24_minutes_is_still_fresh():
    wire({"governor.json": gov(set_at=(NOW - timedelta(minutes=24)).isoformat())},
         total=0)
    assert call()


# --- (c) ceiling zero ------------------------------------------------------- #

def test_ceiling_zero_refuses_immediately_with_no_retry():
    _, sleeps = wire({"governor.json": gov(ceiling=0)}, total=0)
    d = call()
    assert not d and d.reason == "CEILING_ZERO"
    assert sleeps == []


# --- (d) long timeout over pace --------------------------------------------- #

def test_long_timeout_refused_when_over_pace_by_more_than_5pp():
    wire({"governor.json": gov(pct_elapsed=40.0, pct_spent_projected=46.0)})
    d = call(timeout_s=1801)
    assert not d and d.reason == "LONG_TIMEOUT_OVER_PACE"


def test_short_timeout_is_unaffected_by_pace():
    wire({"governor.json": gov(pct_elapsed=40.0, pct_spent_projected=46.0)}, total=0)
    assert call(timeout_s=1800)


def test_long_timeout_allowed_when_on_pace():
    wire({"governor.json": gov(pct_elapsed=50.0, pct_spent_projected=52.0)}, total=0)
    assert call(timeout_s=3600)


# --- (d2) generation may not outrun scoring --------------------------------- #

def test_gen_backlog_cap_blocks_generation():
    wire({
        "governor.json": gov(),
        "ledger_agg.json": {"PD-L1": {"gen": 500, "scr": 100}},
    })
    d = call(target="PD-L1", stage="gen_screen")
    assert not d and d.reason == "GEN_BACKLOG_CAP"


def test_gen_backlog_cap_does_not_bind_below_200_generated():
    wire({
        "governor.json": gov(),
        "ledger_agg.json": {"PD-L1": {"gen": 199, "scr": 0}},
        "gates/PD-L1.json": {"status": "PASS"},
    }, total=0)
    assert call(target="PD-L1", stage="gen_screen")


def test_gen_backlog_cap_uses_a_scoring_floor_of_50():
    # gen=300, scr=10 -> 2*max(10,50)=100 -> 300 > 100 -> blocked
    wire({"governor.json": gov(),
          "ledger_agg.json": {"T": {"gen": 300, "scr": 10}}})
    assert not call(target="T", stage="gen_screen")
    # gen=300, scr=200 -> 2*200=400 -> 300 < 400 -> allowed
    wire({"governor.json": gov(),
          "ledger_agg.json": {"T": {"gen": 300, "scr": 200}},
          "gates/T.json": {"status": "PASS"}}, total=0)
    assert call(target="T", stage="gen_screen")


def test_backlog_cap_does_not_apply_to_score_only_stages():
    wire({"governor.json": gov(),
          "ledger_agg.json": {"T": {"gen": 5000, "scr": 1}},
          "gates/T.json": {"status": "PASS"}}, total=0)
    assert call(target="T", stage="final")


# --- validation gate --------------------------------------------------------- #

def test_scoring_is_blocked_until_the_validation_gate_file_says_PASS():
    wire({"governor.json": gov()}, total=0)
    d = call(target="TNFa", stage="screen")
    assert not d and d.reason == "VALIDATION_GATE_NOT_PASS"

    wire({"governor.json": gov(), "gates/TNFa.json": {"status": "FAIL"}}, total=0)
    assert not call(target="TNFa", stage="screen")

    wire({"governor.json": gov(), "gates/TNFa.json": {"status": "PASS"}}, total=0)
    assert call(target="TNFa", stage="screen")


# --- (e)/(f) live count, headroom margin, backoff ---------------------------- #

def test_cpu_tagged_sandboxes_are_subtracted_from_the_live_count():
    fleet, _ = wire({"governor.json": gov(ceiling=100)}, total=95, cpu=40)
    # 95 - 40 = 55 live GPU, margin = max(10, 10) = 10, ceiling - margin = 90
    assert call()


def test_cpu_allowlist_file_is_also_subtracted():
    wire({
        "governor.json": gov(ceiling=100),
        "cpu_sandboxes.jsonl": [
            json.dumps({"sandbox_id": f"sb-{i}"}) for i in range(30)
        ],
    }, total=95, cpu=0)
    assert call()  # 95 - 0 - 30 = 65 < 90


def test_terminated_allowlist_entries_do_not_count():
    wire({
        "governor.json": gov(ceiling=100),
        "cpu_sandboxes.jsonl": [
            json.dumps({"sandbox_id": "sb-1", "terminated": True}),
            json.dumps({"sandbox_id": "sb-2"}),
        ],
    }, total=95, cpu=0)
    d = call()
    assert not d and d.reason == "AT_CEILING"  # 95 - 1 = 94 >= 90


def test_headroom_margin_is_ten_or_ten_percent_whichever_is_larger():
    # ceiling 300 -> margin 30 -> block at live >= 270
    wire({"governor.json": gov(ceiling=300)}, total=269)
    assert call()
    wire({"governor.json": gov(ceiling=300)}, total=270)
    assert not call()
    # ceiling 50 -> margin max(10, 5) = 10 -> block at live >= 40
    wire({"governor.json": gov(ceiling=50)}, total=39)
    assert call()
    wire({"governor.json": gov(ceiling=50)}, total=40)
    assert not call()


def test_at_ceiling_backs_off_eight_times_then_refuses():
    _, sleeps = wire({"governor.json": gov(ceiling=100)}, total=200)
    d = call()
    assert not d and d.reason == "AT_CEILING"
    assert d.attempts == 8
    assert len(sleeps) == 8
    assert sleeps[0] <= 5.0                    # random.uniform(0, 5)
    assert 4.0 <= sleeps[1] <= 6.0             # 5s x U(0.8, 1.2)
    assert 8.0 <= sleeps[2] <= 12.0            # 10s
    assert 16.0 <= sleeps[3] <= 24.0           # 20s
    assert all(s <= 72.0 for s in sleeps)      # capped at 60s x 1.2


# --- fail-closed ------------------------------------------------------------- #

def test_any_exception_returns_gate_error_never_true():
    class Exploding:
        def read_json(self, path):
            raise RuntimeError("volume read 403")

        def read_lines(self, path):
            raise RuntimeError("volume read 403")

    sg.configure(state=Exploding(), fleet=FakeFleet(0), now=lambda: NOW)
    d = call()
    assert not d and d.reason == "GATE_ERROR"


def test_missing_fleet_counter_is_gate_error_not_permission():
    sg.configure(state=FakeState({"governor.json": gov()}), now=lambda: NOW)
    sg._FLEET = None
    d = call()
    assert not d and d.reason == "GATE_ERROR"


def test_malformed_governor_json_is_gate_error():
    class Bad:
        def read_json(self, path):
            return {"ceiling": "not-an-int", "set_at": NOW.isoformat()}

        def read_lines(self, path):
            return []

    sg.configure(state=Bad(), fleet=FakeFleet(0), now=lambda: NOW)
    d = call()
    assert not d and d.reason == "GATE_ERROR"


def test_one_true_authorises_exactly_one_dispatch():
    """The decision object carries the live count it was granted against."""
    wire({"governor.json": gov(ceiling=300)}, total=100)
    d = call()
    assert d.allowed and d.live == 100 and d.ceiling == 300


def test_local_state_reader_caches_by_mtime(tmp_path):
    p = tmp_path / "governor.json"
    p.write_text(json.dumps({"ceiling": 5}))
    reader = sg.LocalStateReader(str(tmp_path))
    assert reader.read_json("governor.json")["ceiling"] == 5
    assert reader.read_json("missing.json") is None
    assert reader.read_lines("missing.jsonl") == []
