"""End-to-end dry run of the campaign machinery on mock scores.

Runs the whole close-out path for two targets — one monomeric, one
counter-screened — with synthetic designs and synthetic co-folding scores, and
writes the real deliverables:

* ``design_sheet.csv`` (exactly 30 ranked rows per target, or the actual N with
  a deviation row if the ladder cannot reach 30),
* ``per_seed_metrics.parquet``,
* ``instrument_realization.csv``,
* ``scoreboard.csv``,
* ``deviations.jsonl``,
* ``sheet_schema.json`` / ``method_vocab.json``,
* ``gates/<target>.json``, ``governor.json``, ``ledger_agg.json``.

No GPU, no Modal, no network.  The point is to exercise every invariant the
prompt makes mechanical, so a real campaign can plug real scorers into the same
interfaces.  Run with ``python -m binder_campaign.demo``.
"""

from __future__ import annotations

import json
import os
import random
from datetime import datetime, timedelta, timezone

import pandas as pd

from .companions import (
    assert_companion_coverage,
    build_instrument_realization,
    build_per_seed_metrics,
    recompute_final_scores_from_companion,
)
from .filters import GateThresholds, run_prescoring_gates
from .gates import (
    InstrumentGateRow,
    ScoringConstruct,
    ValidationGate,
    read_gate,
    write_gate,
)
from .governor import (
    BudgetReading,
    SandboxSpec,
    budget_cycle,
    campaign_hourly_target,
    derive_max_instances,
    seed_governor,
)
from .lcp import lcp_score
from .ledger import DesignCountLedger, LedgerRow
from .schema import default_method_vocab, default_sheet_schema
from .scoreboard import build_scoreboard, scoreboard_gaps
from .scoring import (
    DEFAULT_ARMS,
    DesignScore,
    InstrumentMask,
    SeedRecord,
    aggregate_arm,
    score_pool,
)
from .sheet_writer import ROWS_PER_TARGET, select_and_rank

AA = "ACDEFGHIKLMNPQRSTVWY"
STRUCT_METHODS = ("rfdiffusion", "rfdiffusion3", "freebindcraft", "boltzgen",
                  "pxdesign", "genie3", "proteina_complexa")
SEQ_METHODS = ("solublempnn", "solublecaliby", "solublecaliby_ensemble")


def _random_sequence(rng: random.Random, n: int) -> str:
    # biased toward a native-like composition so LCP does not fire on everything
    weights = [8, 2, 6, 7, 4, 7, 2, 6, 6, 9, 2, 4, 5, 4, 5, 6, 5, 7, 1, 3]
    return "".join(rng.choices(AA, weights=weights, k=n))


def _make_designs(
    rng: random.Random, target: str, n: int, counter_screened: bool
) -> tuple[list[DesignScore], list[dict]]:
    """Synthetic designs plus their sheet-side metadata."""
    scores: list[DesignScore] = []
    meta: list[dict] = []
    n_roots = max(n // 3, 8)

    for i in range(n):
        root = f"{target}_bb{i % n_roots:03d}"
        method = STRUCT_METHODS[i % len(STRUCT_METHODS)]
        seq_method = SEQ_METHODS[i % len(SEQ_METHODS)]
        design_id = f"{target}_{method}_{i:04d}"
        seq = _random_sequence(rng, rng.randint(55, 110))

        quality = rng.betavariate(2.0, 3.0)  # the design's latent goodness
        on: dict[str, object] = {}
        off: dict[str, object] = {}
        for arm in DEFAULT_ARMS:
            n_seeds = 5 if i < int(0.8 * n) else 3
            seeds = rng.sample(range(1000), n_seeds)
            recs = []
            for s in seeds:
                base = quality * 0.6 + rng.gauss(0, 0.05)
                recs.append(SeedRecord(
                    seed=s,
                    ipsae_ab=max(0.0, base + rng.gauss(0, 0.02)),
                    ipsae_ba=max(0.0, base + rng.gauss(0, 0.02)),
                    dockq=min(1.0, max(0.0, quality * 0.8 + rng.gauss(0, 0.1))),
                    structure_path=(
                        f"ordered_structures/{target}/{design_id}/"
                        f"predicted_{arm}_seedbest.cif"
                    ),
                ))
            on[arm] = aggregate_arm(arm, recs)
            if counter_screened:
                off_recs = [
                    SeedRecord(
                        seed=s,
                        ipsae_ab=max(0.0, quality * 0.35 + rng.gauss(0, 0.05)),
                        ipsae_ba=max(0.0, quality * 0.35 + rng.gauss(0, 0.05)),
                        dockq=min(1.0, max(0.0, quality * 0.5 + rng.gauss(0, 0.1))),
                    )
                    for s in seeds
                ]
                off[arm] = aggregate_arm(arm, off_recs)

        scores.append(DesignScore(design_id, target, on, off))  # type: ignore[arg-type]

        helical_fraction = rng.uniform(0.4, 1.0)
        meta.append({
            "design_id": design_id,
            "target": target,
            "sequence": seq,
            "binder_len": len(seq),
            "structure_method": method,
            "seq_method": seq_method,
            "opt_round": i % 3,
            "root_backbone_id": root,
            "parent_design_id": "" if i % 3 == 0 else f"{target}_bb{i % n_roots:03d}_r0",
            "tm90_cluster_id": f"{target}_c{i % max(n // 4, 6):03d}",
            "fold_class": "all_alpha" if helical_fraction >= 0.7 else "not_all_alpha",
            "designed_structure_path":
                f"ordered_structures/{target}/{design_id}/designed.pdb",
            "novelty_verdict_path": f"novelty/{design_id}.json",
            "monomer_plddt": rng.uniform(65, 95),
            "lcp_score": lcp_score(seq),
            "esmc_ll": rng.gauss(-1.2, 0.3),
            "construct_status": "OK",
            "relaxation_step": "",
            "binder_binder_clashes_NN": 0,
        })
    return scores, meta


def run_demo(out_dir: str, seed: int = 7) -> dict:
    rng = random.Random(seed)
    os.makedirs(out_dir, exist_ok=True)
    state_root = os.path.join(out_dir, "state")
    ledger_root = os.path.join(out_dir, "ledger")
    os.makedirs(state_root, exist_ok=True)

    thresholds = GateThresholds()
    vocab = default_method_vocab()
    schema = default_sheet_schema(vocab=vocab)

    with open(os.path.join(out_dir, "method_vocab.json"), "w") as fh:
        json.dump(vocab.to_json(), fh, indent=2, sort_keys=True)
    with open(os.path.join(out_dir, "sheet_schema.json"), "w") as fh:
        fh.write(schema.dumps())

    targets = {
        "PD-L1": InstrumentMask(),
        "GDF-8": InstrumentMask(counter_screened=True),
    }

    # ---- governor: T0 seed and one BUDGET cycle ---------------------------- #
    t0 = datetime(2026, 6, 1, tzinfo=timezone.utc)
    end = t0 + timedelta(hours=48)
    gov = seed_governor(t0)
    hourly = campaign_hourly_target(50_000.0, end - t0)
    even_pace = derive_max_instances(hourly, [SandboxSpec()])
    cycle = budget_cycle(
        BudgetReading(
            now=t0 + timedelta(minutes=15), t0=t0, campaign_end=end,
            budget_usd=50_000.0, metered_usd=None, metered_asof_utc=None,
            ratecard_usd=180.0, live_sb=120,
        ),
        gov,
    )
    with open(os.path.join(state_root, "governor.json"), "w") as fh:
        json.dump(cycle.state.to_json(), fh, indent=2, sort_keys=True)

    # ---- ledger ------------------------------------------------------------ #
    led = DesignCountLedger(ledger_root, writer_frame_id="demo_frame_0001")

    results: dict[str, object] = {}
    ranked_by_target: dict[str, list[dict]] = {}
    all_sheet_rows: list[dict] = []
    all_companion: list[pd.DataFrame] = []
    all_designs: dict[str, list[DesignScore]] = {}
    gates_json: dict[str, dict] = {}
    all_deviations: list[dict] = []

    for target, mask in targets.items():
        # ---- validation gate ----------------------------------------------- #
        gate = ValidationGate(
            target=target,
            construct=ScoringConstruct(
                chains=("A",) if target == "PD-L1" else ("A", "B"),
                residue_range="18-134" if target == "PD-L1" else "268-375",
                cofactors=(),
                n_target_chains=1 if target == "PD-L1" else 2,
                native_oligomer_n=1 if target == "PD-L1" else 2,
                reference_structure=None,  # a real campaign fills this from a lookup
            ),
            mask=mask,
            instruments=[
                InstrumentGateRow(
                    arm=arm, ca_rmsd=1.1, ca_rmsd_threshold=2.0,
                    control_name="PD-1 ectodomain" if target == "PD-L1" else "ActRIIB",
                    control_score=0.62, negative_control_scores=(0.18, 0.21),
                    n_target_chains_folded=1 if target == "PD-L1" else 2,
                )
                for arm in mask.arms
            ],
            counter_target="GDF-11" if mask.counter_screened else None,
            counter_target_instruments=[
                InstrumentGateRow(arm=arm, ca_rmsd=1.3, ca_rmsd_threshold=2.0)
                for arm in mask.arms
            ] if mask.counter_screened else [],
        )
        write_gate(state_root, gate)
        gates_json[target] = read_gate(state_root, target) or {}
        assert gates_json[target]["status"] == "PASS", gates_json[target]["blockers"]

        # ---- generate + pre-scoring gates ---------------------------------- #
        n_pool = 240
        designs, meta = _make_designs(rng, target, n_pool, mask.counter_screened)
        all_designs[target] = designs

        verdicts: dict[str, dict] = {}
        for m in meta:
            verdicts[m["design_id"]] = run_prescoring_gates(
                m["design_id"], m["sequence"],
                mean_plddt=m["monomer_plddt"],
                thresholds=thresholds,
            )
        rejected = [d for d, v in verdicts.items() if v["verdict"] == "REJECT"]

        led.append(LedgerRow(
            job_id=f"{target}_gen_0001", target=target,
            structure_method="rfdiffusion", stage="gen_screen",
            n_generated=n_pool, n_scored=n_pool - len(rejected),
            gpu_seconds=3600.0, writer_frame_id="demo_frame_0001",
        ))
        for m_name in STRUCT_METHODS[1:]:
            led.append(LedgerRow(
                job_id=f"{target}_gen_{m_name}", target=target,
                structure_method=m_name, stage="gen_screen",
                n_generated=60, n_scored=55, gpu_seconds=1800.0,
                writer_frame_id="demo_frame_0001",
            ))
        led.append(LedgerRow(
            job_id=f"{target}_final_0001", target=target,
            structure_method="rfdiffusion", stage="final",
            n_generated=0, n_scored=ROWS_PER_TARGET, gpu_seconds=900.0,
            writer_frame_id="demo_frame_0001",
        ))

        # ---- score the pool (transductive z over exactly this pool) -------- #
        survivors = [d for d in designs if verdicts[d.design_id]["verdict"] == "PASS"]
        scored = score_pool(survivors, mask)
        by_id = {m["design_id"]: m for m in meta}
        rows = []
        for s in scored:
            row = {**by_id[s["design_id"]], **s}
            v = verdicts[s["design_id"]]
            row["novelty_verdict"] = v["gates"]["novelty"]["verdict"]
            row["liability_verdict"] = v["gates"]["liability"]["verdict"]
            row["monomer_foldability_verdict"] = \
                v["gates"]["monomer_foldability"]["verdict"]
            row["structural_plausibility_verdict"] = \
                v["gates"]["structural_plausibility"]["verdict"]
            for arm in mask.arms:
                row.setdefault(f"ipsae_NN_{arm}", float("nan"))
                row.setdefault(f"sc_DockQ_NN_{arm}", float("nan"))
                if not mask.counter_screened:
                    row.setdefault(f"selectivity_delta_{arm}", float("nan"))
                    row.setdefault(f"ipsae_offtarget_{arm}", float("nan"))
            rows.append(row)

        # ---- sheet writer --------------------------------------------------- #
        result = select_and_rank(
            rows, target=target, mask=mask, schema=schema, vocab=vocab,
            thresholds=thresholds,
        )
        ranked_by_target[target] = result.ranked
        all_sheet_rows.extend(result.ranked)
        all_deviations.extend(result.deviations)
        for dev in result.deviations:
            from .ledger import Deviation
            led.append_deviation(Deviation(**dev))
        results[target] = result.diagnostics

        # ---- companion ------------------------------------------------------ #
        ranked_ids = {r["design_id"] for r in result.ranked}
        all_companion.append(build_per_seed_metrics(
            [d for d in designs if d.design_id in ranked_ids], mask
        ))

    # ---- deliverables ------------------------------------------------------- #
    sheet = pd.DataFrame(all_sheet_rows)
    companion = pd.concat(all_companion, ignore_index=True)
    mask_for = targets

    for target, mask in targets.items():
        # each target is checked under its OWN frozen mask, so the
        # counter-screened target's seed-matched off-target arm is verified too
        sheet_t = sheet[sheet["target"] == target]
        comp_t = companion[companion["target"] == target]
        assert_companion_coverage(sheet_t, comp_t, mask)
        recompute_final_scores_from_companion(sheet_t, comp_t, mask)

    realization = build_instrument_realization(sheet, companion, gates_json, mask_for)

    totals = led.totals()
    with open(os.path.join(state_root, "ledger_agg.json"), "w") as fh:
        json.dump(totals.ledger_agg(), fh, indent=2, sort_keys=True)
    with open(os.path.join(state_root, "floor_matrix.json"), "w") as fh:
        json.dump(totals.floor_matrix(targets, STRUCT_METHODS), fh, indent=2)

    board = build_scoreboard(
        ranked_by_target,
        totals,
        gpu_dollars={t: 1200.0 + 300.0 * i for i, t in enumerate(targets)},
        active_compute_hours={t: 6.0 + i for i, t in enumerate(targets)},
    )
    gaps = scoreboard_gaps(board, owning_frame_id="demo_frame_0001")

    sheet_cols = [c for c in schema.columns if c in sheet.columns]
    sheet[sheet_cols].to_csv(os.path.join(out_dir, "design_sheet.csv"), index=False)
    companion.to_parquet(os.path.join(out_dir, "per_seed_metrics.parquet"),
                         index=False)
    realization.to_csv(os.path.join(out_dir, "instrument_realization.csv"),
                       index=False)
    board.to_csv(os.path.join(out_dir, "scoreboard.csv"), index=False)
    gaps.to_csv(os.path.join(out_dir, "scoreboard_gaps.csv"), index=False)
    with open(os.path.join(out_dir, "deviations.jsonl"), "w") as fh:
        for dev in led.deviations():
            fh.write(json.dumps(dev, sort_keys=True) + "\n")

    return {
        "out_dir": out_dir,
        "even_pace_ceiling": even_pace,
        "hourly_target_usd": hourly,
        "governor": cycle.state.to_json(),
        "raise_blocked_reason": cycle.raise_blocked_reason,
        "totals": {
            "designs_generated": totals.designs_generated,
            "designs_screened": totals.designs_screened,
            "designs_ranked": totals.designs_ranked,
        },
        "per_target": results,
        "n_sheet_rows": len(sheet),
        "n_companion_rows": len(companion),
        "n_scoreboard_gaps": len(gaps),
        "n_deviations": len(all_deviations),
    }


def main() -> None:  # pragma: no cover - CLI
    out = os.environ.get("BINDER_DEMO_OUT", "demo_out")
    summary = run_demo(out)
    print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":  # pragma: no cover
    main()
