"""BUDGET's cycle: burn rates, the pace band, calibration, the fail-safes."""

from datetime import datetime, timedelta, timezone

import pytest

from binder_campaign.governor import (
    BOOTSTRAP_CEILING,
    BudgetReading,
    GovernorState,
    SandboxSpec,
    budget_cycle,
    campaign_hourly_target,
    daemon_limit,
    derive_max_instances,
    is_calibrated,
    orchestrator_resume_clear,
    sandbox_hourly_rate,
    seed_governor,
    watchdog_dead_orchestrator,
    watchdog_runaway,
    weighted_per_sandbox_rate,
)

T0 = datetime(2026, 6, 1, tzinfo=timezone.utc)
END = T0 + timedelta(hours=48)
BUDGET = 50_000.0


def reading(**kw):
    base = dict(
        now=T0 + timedelta(hours=10),
        t0=T0,
        campaign_end=END,
        budget_usd=BUDGET,
        metered_usd=None,
        metered_asof_utc=None,
        ratecard_usd=0.0,
        live_sb=0,
    )
    base.update(kw)
    return BudgetReading(**base)


# --- burn-rate arithmetic ---------------------------------------------------- #

def test_typical_h100_box_is_about_4_39_per_hour():
    """gpu_rate + cores x $0.047/h + mem_gib x $0.008/h, per the prompt."""
    rate = sandbox_hourly_rate(SandboxSpec("H100", cores=4, mem_gib=32))
    assert rate == pytest.approx(4.394, abs=0.005)


def test_even_pace_ceiling_is_about_237_h100_equivalents():
    hourly = campaign_hourly_target(BUDGET, END - T0)
    assert hourly == pytest.approx(1041.67, abs=0.5)
    assert derive_max_instances(hourly, [SandboxSpec()]) == 237


def test_bootstrap_cap_runs_about_1_37x_even_pace():
    hourly = campaign_hourly_target(BUDGET, END - T0)
    even = derive_max_instances(hourly, [SandboxSpec()])
    assert BOOTSTRAP_CEILING / even == pytest.approx(1.37, abs=0.01)


def test_weighted_rate_follows_the_live_fleet_mix():
    mixed = [SandboxSpec("H100", count=1), SandboxSpec("A100-40GB", count=1)]
    r = weighted_per_sandbox_rate(mixed)
    assert sandbox_hourly_rate(SandboxSpec("A100-40GB")) < r
    assert r < sandbox_hourly_rate(SandboxSpec("H100"))


def test_unknown_gpu_type_is_refused_rather_than_guessed():
    with pytest.raises(KeyError, match="modal.com/pricing"):
        sandbox_hourly_rate(SandboxSpec("B200"))


# --- bootstrap stepping ------------------------------------------------------ #

def test_setup_seeds_325_with_basis_bootstrap():
    g = seed_governor(T0)
    assert g.ceiling == 325 and g.basis == "BOOTSTRAP"


def test_uncalibrated_bootstrap_steps_by_at_most_50_per_cycle():
    prev = GovernorState(ceiling=100, set_at=T0, basis="BUDGET", metered_asof_utc=T0)
    res = budget_cycle(reading(now=T0 + timedelta(minutes=30), ratecard_usd=50.0), prev)
    assert res.state.ceiling == 150
    assert res.state.calibrated is False


def test_bootstrap_never_exceeds_325():
    prev = GovernorState(ceiling=300, set_at=T0, basis="BUDGET", metered_asof_utc=T0)
    res = budget_cycle(reading(now=T0 + timedelta(minutes=30)), prev)
    assert res.state.ceiling == 325


def test_raises_stop_after_three_hours_without_a_calibrated_reading():
    prev = GovernorState(ceiling=100, set_at=T0, basis="BUDGET", metered_asof_utc=T0)
    res = budget_cycle(reading(now=T0 + timedelta(hours=3, minutes=1)), prev)
    assert res.state.ceiling == 100  # holds, lowering remains allowed
    assert res.raise_blocked_reason == "BILLING_DARK_OVER_3H"
    assert any("BILLING_DARK" in a for a in res.alerts)


def test_billing_dark_warn_at_t0_plus_90_minutes():
    prev = GovernorState(ceiling=100, set_at=T0, basis="BUDGET", metered_asof_utc=T0)
    res = budget_cycle(reading(now=T0 + timedelta(minutes=95)), prev)
    assert any("BILLING_DARK_WARN" in a for a in res.alerts)


# --- calibration ------------------------------------------------------------- #

def _prev_metered(usd, asof):
    return GovernorState(ceiling=200, set_at=T0, basis="BUDGET",
                         metered_usd=usd, metered_asof_utc=asof, calibrated=True)


def test_calibration_requires_all_four_conditions():
    prev = _prev_metered(0.0, T0)
    asof = T0 + timedelta(hours=2)

    # canary absent
    ok, why = is_calibrated(reading(metered_usd=100.0, metered_asof_utc=asof,
                                    canary_job_in_metered=False), prev)
    assert not ok and why == "CANARY_NOT_IN_METERED"

    # below the $20 floor and not covering 90% of ratecard
    ok, why = is_calibrated(reading(metered_usd=10.0, metered_asof_utc=asof,
                                    ratecard_usd=1000.0,
                                    canary_job_in_metered=True), prev)
    assert not ok and why == "METERED_BELOW_CALIBRATION_FLOOR"

    # unchanged from the previous reading
    ok, why = is_calibrated(reading(metered_usd=100.0, metered_asof_utc=asof,
                                    canary_job_in_metered=True),
                            _prev_metered(98.0, T0))
    assert not ok and why == "METERED_UNCHANGED"

    # as-of timestamp did not advance
    ok, why = is_calibrated(reading(metered_usd=100.0, metered_asof_utc=T0,
                                    canary_job_in_metered=True), prev)
    assert not ok and why == "METERED_ASOF_NOT_ADVANCED"

    # all four satisfied
    ok, why = is_calibrated(reading(metered_usd=100.0, metered_asof_utc=asof,
                                    canary_job_in_metered=True), prev)
    assert ok and why is None


def test_no_raise_while_ratecard_materially_trails_metered():
    prev = _prev_metered(100.0, T0)
    res = budget_cycle(
        reading(now=T0 + timedelta(hours=4), metered_usd=200.0,
                metered_asof_utc=T0 + timedelta(hours=3), ratecard_usd=900.0,
                canary_job_in_metered=True),
        prev,
    )
    assert res.raise_blocked_reason == "RATECARD_TRAILS_METERED"
    assert res.state.ceiling == prev.ceiling


def test_no_raise_while_two_calibrated_readings_disagree_by_more_than_20_percent():
    prev = _prev_metered(1000.0, T0)
    res = budget_cycle(
        reading(now=T0 + timedelta(hours=4), metered_usd=1500.0,
                metered_asof_utc=T0 + timedelta(hours=3), ratecard_usd=1500.0,
                canary_job_in_metered=True),
        prev,
    )
    assert res.raise_blocked_reason == "METERED_READINGS_DISAGREE"


def test_no_raise_while_live_exceeds_the_current_ceiling():
    prev = GovernorState(ceiling=100, set_at=T0, basis="BUDGET", metered_asof_utc=T0)
    res = budget_cycle(reading(now=T0 + timedelta(minutes=30), live_sb=150), prev)
    assert res.raise_blocked_reason == "LIVE_OVER_CEILING"
    assert res.state.ceiling == 100


# --- the pace band ----------------------------------------------------------- #

def test_over_band_by_15pp_writes_ceiling_zero_for_one_cycle():
    prev = GovernorState(ceiling=300, set_at=T0, basis="BUDGET", metered_asof_utc=T0,
                         calibrated=True)
    # 10h of 48h elapsed = 20.8%; spend 40% -> over by ~19pp
    res = budget_cycle(
        reading(metered_usd=0.40 * BUDGET, metered_asof_utc=T0 + timedelta(hours=10),
                canary_job_in_metered=True),
        prev,
    )
    assert res.state.ceiling == 0
    assert res.lower_daemon_limit_to == 4  # max(2, standing singletons)
    assert any("OVER band" in a for a in res.alerts)


def test_ceiling_zero_binds_on_the_projection_while_billing_is_dark():
    """metered_usd=0 but the rate-card reconstruction is already over the band."""
    prev = GovernorState(ceiling=300, set_at=T0, basis="BUDGET", metered_asof_utc=T0)
    res = budget_cycle(reading(metered_usd=None, ratecard_usd=0.45 * BUDGET), prev)
    assert res.state.ceiling == 0
    assert res.state.pct_spent_projected == pytest.approx(45.0)


def test_live_over_ceiling_for_two_cycles_lowers_to_80_percent():
    prev = GovernorState(ceiling=300, set_at=T0, basis="BUDGET", metered_asof_utc=T0)
    res = budget_cycle(
        reading(now=T0 + timedelta(hours=10), live_sb=400,
                consecutive_over_ceiling=1),
        prev,
    )
    assert res.state.ceiling == 240  # floor(0.8 * 300)


def test_under_pace_for_three_cycles_demands_reallocation():
    prev = GovernorState(ceiling=300, set_at=T0, basis="BUDGET", metered_asof_utc=T0,
                         calibrated=True)
    # 50% elapsed, 10% spent -> 40pp under
    res = budget_cycle(
        reading(now=T0 + timedelta(hours=24), metered_usd=0.10 * BUDGET,
                metered_asof_utc=T0 + timedelta(hours=24),
                canary_job_in_metered=True, live_sb=280,
                consecutive_under_pace=2),
        prev,
    )
    assert res.wake_scale_up is True
    assert res.consecutive_under_pace == 3


def test_under_pace_streak_resets_when_back_in_band():
    prev = GovernorState(ceiling=300, set_at=T0, basis="BUDGET", metered_asof_utc=T0,
                         calibrated=True)
    res = budget_cycle(
        reading(now=T0 + timedelta(hours=24), metered_usd=0.49 * BUDGET,
                metered_asof_utc=T0 + timedelta(hours=24),
                canary_job_in_metered=True, live_sb=290,
                consecutive_under_pace=2),
        prev,
    )
    assert res.wake_scale_up is False
    assert res.consecutive_under_pace == 0


# --- fail-safes -------------------------------------------------------------- #

def test_watchdog_runaway_halts_and_records_the_pre_halt_ceiling():
    prev = GovernorState(ceiling=200, set_at=T0, basis="BUDGET", metered_asof_utc=T0)
    assert watchdog_runaway(prev, live_gpu=350, now=T0) is None  # 350 <= max(400, 300)
    halted = watchdog_runaway(prev, live_gpu=450, now=T0)
    assert halted.ceiling == 0 and halted.basis == "WATCHDOG_RUNAWAY"
    assert halted.pre_halt_ceiling == 200
    # never re-fires, never overwrites pre_halt_ceiling
    assert watchdog_runaway(halted, live_gpu=900, now=T0) is None


def test_budget_preserves_a_runaway_halt_until_the_fleet_drains():
    halted = GovernorState(ceiling=0, set_at=T0, basis="WATCHDOG_RUNAWAY",
                           metered_asof_utc=T0, pre_halt_ceiling=200)
    still = budget_cycle(reading(live_sb=400, consecutive_over_ceiling=1), halted)
    assert still.state.ceiling == 0 and still.state.basis == "WATCHDOG_RUNAWAY"

    cleared = budget_cycle(reading(live_sb=150), halted)
    assert cleared.state.basis == "RUNAWAY_CLEARED"
    assert cleared.state.ceiling == 200


def test_budget_preserves_a_dead_orchestrator_halt_while_the_heartbeat_is_stale():
    dead = watchdog_dead_orchestrator(
        GovernorState(ceiling=300, set_at=T0, basis="BUDGET", metered_asof_utc=T0), T0
    )
    assert dead.ceiling == 0 and dead.basis == "DEAD_ORCHESTRATOR"
    held = budget_cycle(reading(heartbeat_age_s=20 * 60), dead)
    assert held.state.ceiling == 0 and held.state.basis == "DEAD_ORCHESTRATOR"


def test_only_the_orchestrator_clears_dead_orchestrator():
    dead = GovernorState(ceiling=0, set_at=T0, basis="DEAD_ORCHESTRATOR",
                         metered_asof_utc=T0, pct_elapsed=50.0,
                         pct_spent_metered=45.0, pct_spent_projected=45.0)
    resumed = orchestrator_resume_clear(dead, T0 + timedelta(hours=1))
    assert resumed.basis == "ORCHESTRATOR_RESUMED"
    assert resumed.ceiling > 0


def test_daemon_limit_never_exceeds_ceiling_plus_100_or_450():
    assert daemon_limit(300, 10) == 310
    assert daemon_limit(300, 500) == 400   # N capped at 100
    assert daemon_limit(400, 100) == 450   # hard cap
    assert daemon_limit(0, 0) == 1         # max(1, ...)


def test_published_governor_json_carries_every_required_field():
    prev = seed_governor(T0)
    res = budget_cycle(reading(), prev)
    payload = res.state.to_json()
    for field in ("ceiling", "set_at", "basis", "metered_usd", "metered_asof_utc",
                  "ratecard_usd", "live_sb", "pct_elapsed", "pct_spent_metered",
                  "pct_spent_projected", "calibrated"):
        assert field in payload
