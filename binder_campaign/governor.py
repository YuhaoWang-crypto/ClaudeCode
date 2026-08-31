"""The concurrency governor and the BUDGET sub-agent's 15-minute cycle.

The governor is a SINGLE campaign-wide integer with exactly four authorised
writers (SETUP's T0 seed, BUDGET, the WATCHDOG fail-safe, and the orchestrator's
resume-clear).  This module implements BUDGET's arithmetic and the fail-safe
transitions; :mod:`binder_campaign.submit_gate` implements the enforcement side.

Everything here is pure: :func:`budget_cycle` takes a reading and the previous
published state and returns the next :class:`GovernorState` plus the actions the
BUDGET agent must take (WAKE files, Slack alerts, daemon-limit lowering).  That
keeps the pace band, the calibration ladder and the runaway brake unit-testable
without a Modal account.

Prices are the prompt's published per-second Modal rates; the prompt requires
verifying them at kickoff, so :func:`gpu_hourly_rate` accepts an override table.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from datetime import datetime, timedelta, timezone
from typing import Literal, Mapping

__all__ = [
    "GPU_RATES_PER_SECOND",
    "CPU_RATE_PER_CORE_HOUR",
    "MEM_RATE_PER_GIB_HOUR",
    "BOOTSTRAP_CEILING",
    "SandboxSpec",
    "GovernorState",
    "BudgetReading",
    "BudgetCycleResult",
    "sandbox_hourly_rate",
    "weighted_per_sandbox_rate",
    "derive_max_instances",
    "is_calibrated",
    "budget_cycle",
    "watchdog_runaway",
    "watchdog_dead_orchestrator",
    "orchestrator_resume_clear",
]

# --- published Modal rates (modal.com/pricing; verify at kickoff) ----------- #
GPU_RATES_PER_SECOND: dict[str, float] = {
    "H100": 0.001097,
    "A100-80GB": 0.000694,
    "A100-40GB": 0.000583,
    "L40S": 0.000542,
    "A10G": 0.000306,
    "L4": 0.000222,
    "T4": 0.000164,
    "CPU": 0.0,
}
CPU_RATE_PER_CORE_HOUR = 0.047  # $0.0000131/core/s; physical core = 2 vCPU
MEM_RATE_PER_GIB_HOUR = 0.008  # $0.00000222/GiB/s
MIN_BILLED_CORES = 0.125

BOOTSTRAP_CEILING = 325
BOOTSTRAP_STEP = 50
DAEMON_LIMIT_HARD_CAP = 450

PACE_BAND_OVER_PP = 15.0  # >= pct_elapsed + 15 for one cycle -> ceiling 0
PACE_BAND_UNDER_PP = 20.0  # < pct_elapsed - 20 for 3 cycles -> scale-up WAKE
CALIBRATION_MIN_USD = 20.0
CALIBRATION_MIN_DELTA_USD = 5.0
CALIBRATION_COVERAGE = 0.90
RATECARD_DIVERGENCE_FRAC = 0.50
RATECARD_DIVERGENCE_ABS = 200.0
METERED_DISAGREEMENT_FRAC = 0.20

Basis = Literal[
    "BOOTSTRAP",
    "BUDGET",
    "DEAD_ORCHESTRATOR",
    "WATCHDOG_RUNAWAY",
    "RUNAWAY_CLEARED",
    "ORCHESTRATOR_RESUMED",
]


# --------------------------------------------------------------------------- #
# per-sandbox burn rate
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class SandboxSpec:
    """A container spec, for the burn-rate arithmetic."""

    gpu_type: str = "H100"
    cores: float = 4.0
    mem_gib: float = 32.0
    count: int = 1


def gpu_hourly_rate(gpu_type: str, rates: Mapping[str, float] | None = None) -> float:
    table = rates or GPU_RATES_PER_SECOND
    if gpu_type not in table:
        raise KeyError(f"unknown GPU type {gpu_type!r}; verify against modal.com/pricing")
    return table[gpu_type] * 3600.0


def sandbox_hourly_rate(
    spec: SandboxSpec, rates: Mapping[str, float] | None = None
) -> float:
    """``gpu_rate + cores x $0.047/h + mem_gib x $0.008/h``.

    A 4-core / 32-GiB H100 box is roughly $4.39/h all-in, as the prompt states.
    """
    cores = max(spec.cores, MIN_BILLED_CORES)
    return (
        gpu_hourly_rate(spec.gpu_type, rates)
        + cores * CPU_RATE_PER_CORE_HOUR
        + spec.mem_gib * MEM_RATE_PER_GIB_HOUR
    )


def weighted_per_sandbox_rate(
    fleet: list[SandboxSpec], rates: Mapping[str, float] | None = None
) -> float:
    """Live-fleet-weighted mean of ``(gpu_rate + cpu_adder + mem_adder)``.

    Must be re-derived whenever the GPU mix or container spec changes.
    """
    total_count = sum(s.count for s in fleet)
    if total_count == 0:
        return sandbox_hourly_rate(SandboxSpec(), rates)
    return sum(sandbox_hourly_rate(s, rates) * s.count for s in fleet) / total_count


def derive_max_instances(
    target_hourly_rate: float,
    fleet: list[SandboxSpec],
    rates: Mapping[str, float] | None = None,
) -> int:
    """``floor(target_hourly_rate / weighted_per_sandbox_rate)``.

    The prompt's arithmetic: $50,000 / 48 h = $1,042/h against a $4.39/h H100
    box gives ~237 H100-equivalents (the even-pace ceiling), which is why the
    325-instance bootstrap cap runs about 1.37x even pace.
    """
    rate = weighted_per_sandbox_rate(fleet, rates)
    if rate <= 0:
        raise ValueError("weighted per-sandbox rate must be positive")
    return max(0, math.floor(target_hourly_rate / rate))


# --------------------------------------------------------------------------- #
# governor.json
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class GovernorState:
    """The published ``/state/governor.json``."""

    ceiling: int
    set_at: datetime
    basis: str
    metered_usd: float = 0.0
    metered_asof_utc: datetime | None = None
    ratecard_usd: float = 0.0
    live_sb: int = 0
    pct_elapsed: float = 0.0
    pct_spent_metered: float = 0.0
    pct_spent_projected: float = 0.0
    calibrated: bool = False
    pre_halt_ceiling: int | None = None

    def to_json(self) -> dict:
        d = {
            "ceiling": self.ceiling,
            "set_at": self.set_at.astimezone(timezone.utc).isoformat(),
            "basis": self.basis,
            "metered_usd": self.metered_usd,
            "metered_asof_utc": (
                self.metered_asof_utc.astimezone(timezone.utc).isoformat()
                if self.metered_asof_utc
                else None
            ),
            "ratecard_usd": self.ratecard_usd,
            "live_sb": self.live_sb,
            "pct_elapsed": self.pct_elapsed,
            "pct_spent_metered": self.pct_spent_metered,
            "pct_spent_projected": self.pct_spent_projected,
            "calibrated": self.calibrated,
        }
        if self.pre_halt_ceiling is not None:
            d["pre_halt_ceiling"] = self.pre_halt_ceiling
        return d


def seed_governor(t0: datetime) -> GovernorState:
    """SETUP's one-time T0 seed: ``{ceiling: 325, basis: BOOTSTRAP}``."""
    return GovernorState(
        ceiling=BOOTSTRAP_CEILING,
        set_at=t0,
        basis="BOOTSTRAP",
        metered_asof_utc=t0,
    )


# --------------------------------------------------------------------------- #
# BUDGET's cycle
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class BudgetReading:
    """What BUDGET gathers on one 15-minute cycle."""

    now: datetime
    t0: datetime
    campaign_end: datetime
    budget_usd: float
    #: metered spend summed over rows where object_id == this campaign's app_id
    #: ONLY (never the workspace total).  None while billing is dark.
    metered_usd: float | None
    metered_asof_utc: datetime | None
    #: rate card x elapsed lifetime, from /ledger/job_metadata/*.jsonl only.
    ratecard_usd: float
    live_sb: int
    fleet: list[SandboxSpec] = field(default_factory=lambda: [SandboxSpec()])
    canary_job_in_metered: bool = False
    heartbeat_age_s: float = 0.0
    #: consecutive prior cycles where live_sb exceeded the ceiling
    consecutive_over_ceiling: int = 0
    #: consecutive prior cycles that met the -20pp / idle scale-up condition
    consecutive_under_pace: int = 0


@dataclass(frozen=True)
class BudgetCycleResult:
    state: GovernorState
    #: True when BUDGET must write /state/WAKE with a reallocation demand
    wake_scale_up: bool = False
    #: one-line Slack alerts this cycle must post
    alerts: tuple[str, ...] = ()
    #: when set, BUDGET calls host.compute.set_concurrency_limit(this)
    lower_daemon_limit_to: int | None = None
    #: why the ceiling did not go up, for the status block
    raise_blocked_reason: str | None = None
    consecutive_over_ceiling: int = 0
    consecutive_under_pace: int = 0


def is_calibrated(
    reading: BudgetReading, previous: GovernorState
) -> tuple[bool, str | None]:
    """The four-part CALIBRATED test.  Returns ``(calibrated, blocking_reason)``."""
    if reading.metered_usd is None:
        return False, "NO_METERED_READING"
    # (i) the canary job's cost is present in the app-scoped report
    if not reading.canary_job_in_metered:
        return False, "CANARY_NOT_IN_METERED"
    # (ii) value > $20 or covers >= 90% of ledger-derived ratecard_usd
    covers = (
        reading.ratecard_usd > 0
        and reading.metered_usd >= CALIBRATION_COVERAGE * reading.ratecard_usd
    )
    if not (reading.metered_usd > CALIBRATION_MIN_USD or covers):
        return False, "METERED_BELOW_CALIBRATION_FLOOR"
    # (iii) differs from the previous reading by > $5
    if abs(reading.metered_usd - previous.metered_usd) <= CALIBRATION_MIN_DELTA_USD:
        return False, "METERED_UNCHANGED"
    # (iv) its as-of timestamp has advanced
    if (
        reading.metered_asof_utc is None
        or previous.metered_asof_utc is None
        or reading.metered_asof_utc <= previous.metered_asof_utc
    ):
        return False, "METERED_ASOF_NOT_ADVANCED"
    return True, None


def _pace(reading: BudgetReading, rate: float) -> tuple[float, float, float]:
    """``(pct_elapsed, pct_spent_metered, pct_spent_projected)``."""
    window = (reading.campaign_end - reading.t0).total_seconds()
    pct_elapsed = 100.0 * (reading.now - reading.t0).total_seconds() / window
    pct_elapsed = min(max(pct_elapsed, 0.0), 100.0)

    metered = reading.metered_usd or 0.0
    asof = reading.metered_asof_utc or reading.t0
    pct_metered = 100.0 * metered / reading.budget_usd

    hours_since = max((reading.now - asof).total_seconds() / 3600.0, 0.0)
    projected_usd = metered + reading.live_sb * rate * hours_since
    pct_projected = 100.0 * projected_usd / reading.budget_usd

    if reading.metered_usd is None:
        # while billing is dark the band binds on dispatch-ledger data
        pct_projected = max(
            pct_projected, 100.0 * reading.ratecard_usd / reading.budget_usd
        )
    return pct_elapsed, pct_metered, pct_projected


def budget_cycle(
    reading: BudgetReading,
    previous: GovernorState,
    standing_singletons: int = 4,
) -> BudgetCycleResult:
    """One BUDGET cycle: publish the governor, decide raise / hold / lower / halt."""
    rate = weighted_per_sandbox_rate(reading.fleet)
    pct_elapsed, pct_metered, pct_projected = _pace(reading, rate)
    calibrated, calib_reason = is_calibrated(reading, previous)

    base = replace(
        previous,
        set_at=reading.now,
        metered_usd=reading.metered_usd if reading.metered_usd is not None else 0.0,
        metered_asof_utc=reading.metered_asof_utc or reading.t0,
        ratecard_usd=reading.ratecard_usd,
        live_sb=reading.live_sb,
        pct_elapsed=pct_elapsed,
        pct_spent_metered=pct_metered,
        pct_spent_projected=pct_projected,
        calibrated=calibrated,
        basis="BUDGET",
    )

    alerts: list[str] = []
    hours_since_t0 = (reading.now - reading.t0).total_seconds() / 3600.0
    if not calibrated:
        if hours_since_t0 >= 1.5 and not reading.canary_job_in_metered:
            alerts.append("BILLING_DARK_WARN: canary not in metered by T0+90m")
        if hours_since_t0 >= 3.0:
            alerts.append(
                "BILLING_DARK: no calibrated metered reading by T0+3h; "
                "governor raises are suspended"
            )

    # --- fail-safe bases are preserved, not overwritten --------------------- #
    if previous.basis == "DEAD_ORCHESTRATOR" and reading.heartbeat_age_s > 15 * 60:
        return BudgetCycleResult(
            state=replace(base, ceiling=0, basis="DEAD_ORCHESTRATOR"),
            alerts=tuple(alerts),
            raise_blocked_reason="DEAD_ORCHESTRATOR",
            consecutive_over_ceiling=reading.consecutive_over_ceiling,
            consecutive_under_pace=reading.consecutive_under_pace,
        )

    if previous.basis == "WATCHDOG_RUNAWAY":
        pre_halt = previous.pre_halt_ceiling or BOOTSTRAP_CEILING
        cleared = (
            reading.live_sb <= pre_halt and reading.consecutive_over_ceiling == 0
        )
        if not cleared:
            return BudgetCycleResult(
                state=replace(base, ceiling=0, basis="WATCHDOG_RUNAWAY",
                              pre_halt_ceiling=pre_halt),
                alerts=tuple(alerts),
                raise_blocked_reason="WATCHDOG_RUNAWAY",
                consecutive_over_ceiling=reading.consecutive_over_ceiling,
                consecutive_under_pace=reading.consecutive_under_pace,
            )
        return BudgetCycleResult(
            state=replace(base, ceiling=pre_halt, basis="RUNAWAY_CLEARED",
                          pre_halt_ceiling=None),
            alerts=tuple(alerts),
            consecutive_over_ceiling=0,
            consecutive_under_pace=reading.consecutive_under_pace,
        )

    band_value = max(pct_metered, pct_projected)

    # --- OVER the band: ceiling 0, pause NEW dispatch only ------------------ #
    if band_value >= pct_elapsed + PACE_BAND_OVER_PP:
        return BudgetCycleResult(
            state=replace(base, ceiling=0),
            alerts=tuple(alerts + [
                f"OVER band: spent {band_value:.1f}% vs elapsed {pct_elapsed:.1f}%; "
                f"ceiling=0 until the next reading is back inside the band"
            ]),
            lower_daemon_limit_to=max(2, standing_singletons),
            raise_blocked_reason="OVER_PACE_BAND",
            consecutive_over_ceiling=reading.consecutive_over_ceiling,
            consecutive_under_pace=0,
        )

    # --- live over the ceiling for two consecutive cycles: LOWER ------------ #
    over_ceiling = reading.live_sb > previous.ceiling
    n_over = reading.consecutive_over_ceiling + 1 if over_ceiling else 0
    if n_over >= 2:
        return BudgetCycleResult(
            state=replace(base, ceiling=math.floor(0.8 * previous.ceiling)),
            alerts=tuple(alerts),
            raise_blocked_reason="LIVE_OVER_CEILING",
            consecutive_over_ceiling=n_over,
            consecutive_under_pace=reading.consecutive_under_pace,
        )

    # --- may we RAISE? ------------------------------------------------------ #
    blocked = _raise_blockers(reading, previous, calibrated, calib_reason,
                              hours_since_t0, over_ceiling)

    ceiling = previous.ceiling
    if blocked is None and ceiling < BOOTSTRAP_CEILING:
        step = BOOTSTRAP_STEP if not calibrated else BOOTSTRAP_CEILING
        ceiling = min(BOOTSTRAP_CEILING, previous.ceiling + step)

    # --- UNDER the band for three cycles: demand reallocation --------------- #
    under = pct_metered < pct_elapsed - PACE_BAND_UNDER_PP or (
        reading.live_sb < 0.5 * previous.ceiling and pct_elapsed > 10.0
    )
    n_under = reading.consecutive_under_pace + 1 if under else 0

    return BudgetCycleResult(
        state=replace(base, ceiling=ceiling),
        wake_scale_up=n_under >= 3,
        alerts=tuple(alerts),
        raise_blocked_reason=blocked,
        consecutive_over_ceiling=n_over,
        consecutive_under_pace=n_under,
    )


def _raise_blockers(
    reading: BudgetReading,
    previous: GovernorState,
    calibrated: bool,
    calib_reason: str | None,
    hours_since_t0: float,
    over_ceiling: bool,
) -> str | None:
    """Every rule that forbids RAISING the governor.  ``None`` means raise is OK."""
    # never raise while the live count exceeds the current ceiling
    if over_ceiling:
        return "LIVE_OVER_CEILING"

    # never raise while ratecard materially exceeds metered (a stale figure)
    metered = reading.metered_usd
    if metered is not None:
        if (
            reading.ratecard_usd > metered * (1 + RATECARD_DIVERGENCE_FRAC)
            and reading.ratecard_usd - metered > RATECARD_DIVERGENCE_ABS
        ):
            return "RATECARD_TRAILS_METERED"
        # never raise while two consecutive metered readings disagree by > 20%
        if previous.metered_usd > 0:
            rel = abs(metered - previous.metered_usd) / previous.metered_usd
            if previous.calibrated and calibrated and rel > METERED_DISAGREEMENT_FRAC:
                return "METERED_READINGS_DISAGREE"

    if not calibrated:
        # BOOTSTRAP EXCEPTION: step up on UNCALIBRATED rate-card readings, but
        # only until T0+3h; after that the ceiling holds where it is.
        if hours_since_t0 >= 3.0:
            return "BILLING_DARK_OVER_3H"
        return None  # bootstrap step permitted, disclosed as UNCALIBRATED
    return None


# --------------------------------------------------------------------------- #
# WATCHDOG fail-safes and the orchestrator's resume-clear
# --------------------------------------------------------------------------- #


def watchdog_runaway(
    previous: GovernorState, live_gpu: int, now: datetime
) -> GovernorState | None:
    """Runaway brake: fires on ``live_gpu > max(2*ceiling, ceiling+100)``.

    The caller applies the "two consecutive cycles" condition.  Never re-fires
    while the basis is already ``WATCHDOG_RUNAWAY`` and never overwrites
    ``pre_halt_ceiling``; a ``pre_halt_ceiling`` of 0 is never written.
    """
    if previous.basis == "WATCHDOG_RUNAWAY":
        return None
    if live_gpu <= max(2 * previous.ceiling, previous.ceiling + 100):
        return None
    pre_halt = previous.ceiling if previous.ceiling > 0 else BOOTSTRAP_CEILING
    return replace(
        previous,
        ceiling=0,
        basis="WATCHDOG_RUNAWAY",
        set_at=now,
        pre_halt_ceiling=pre_halt,
        live_sb=live_gpu,
    )


def watchdog_dead_orchestrator(
    previous: GovernorState, now: datetime
) -> GovernorState:
    """Heartbeat stale by > 30 minutes: ceiling 0, basis DEAD_ORCHESTRATOR."""
    return replace(previous, ceiling=0, basis="DEAD_ORCHESTRATOR", set_at=now)


def orchestrator_resume_clear(
    previous: GovernorState, now: datetime
) -> GovernorState:
    """The orchestrator (and only it) clears DEAD_ORCHESTRATOR on resume.

    The restored ceiling is what the linear pace band yields at the most recent
    calibrated reading.
    """
    headroom = max(0.0, previous.pct_elapsed + PACE_BAND_OVER_PP -
                   max(previous.pct_spent_metered, previous.pct_spent_projected))
    fraction = min(1.0, headroom / PACE_BAND_OVER_PP) if PACE_BAND_OVER_PP else 1.0
    ceiling = max(1, int(BOOTSTRAP_CEILING * fraction))
    return replace(
        previous, ceiling=ceiling, basis="ORCHESTRATOR_RESUMED", set_at=now
    )


def daemon_limit(ceiling: int, n_live_subagents: int) -> int:
    """``min(450, max(1, ceiling + N))`` with N capped at 100, per the prompt.

    The orchestrator is the only process permitted to RAISE the daemon limit,
    and only to ``governor.ceiling + N`` where N counts its running sub-agents
    (one CPU state-reader each).
    """
    n = min(max(n_live_subagents, 0), 100)
    return min(DAEMON_LIMIT_HARD_CAP, min(ceiling + 100, max(1, ceiling + n)))


def campaign_hourly_target(budget_usd: float, window: timedelta) -> float:
    return budget_usd / (window.total_seconds() / 3600.0)
