"""A reference implementation of the Anthropic protein-binder-design campaign prompt.

The prompt (``prompts/multi_target_binder_design_prompt.md``, CC BY 4.0, from
the HuggingFace dataset ``Anthropic/claude-protein-binder-design``) specifies a
48-hour, $50,000, >100-agent autonomous campaign.  Most of it is *biology* that
needs GPUs, Modal, and a specific agent harness.  A large part of it, though, is
**mechanically specified**: a scoring instrument with exact aggregation rules, a
concurrency governor with exact arithmetic, a fail-closed submit gate, a
sequence restraint defined by an equation, a ledger with fixed totals, and a
sheet writer with a fixed cap and relaxation ladder.

This package implements that mechanical core, exactly as written, with the
GPU-bound steps behind injected interfaces:

===========================  =============================================
:mod:`~binder_campaign.lcp`             Local Composition Perplexity (Figure 1)
:mod:`~binder_campaign.scoring`         the three-arm instrument, ipSAE/DockQ,
                                        ``final_score`` and ``rank_zscore``
:mod:`~binder_campaign.submit_gate`     the fail-closed dispatch gate (a)-(g)
:mod:`~binder_campaign.governor`        BUDGET's cycle, pace band, calibration
:mod:`~binder_campaign.ledger`          design-count / job-metadata / deviations
:mod:`~binder_campaign.filters`         the four pre-scoring gates
:mod:`~binder_campaign.gates`           ``/state/gates/{target}.json``
:mod:`~binder_campaign.schema`          frozen vocabularies and sheet schema
:mod:`~binder_campaign.sheet_writer`    caps, rank, relaxation ladder
:mod:`~binder_campaign.companions`      per-seed metrics, instrument realization
:mod:`~binder_campaign.scoreboard`      the scoreboard CSV
:mod:`~binder_campaign.demo`            an end-to-end dry run on mock scores
===========================  =============================================
"""

from __future__ import annotations

__version__ = "0.1.0"

from . import (  # noqa: F401
    companions,
    filters,
    gates,
    governor,
    lcp,
    ledger,
    schema,
    scoreboard,
    scoring,
    sheet_writer,
    submit_gate,
)

__all__ = [
    "companions",
    "filters",
    "gates",
    "governor",
    "lcp",
    "ledger",
    "schema",
    "scoreboard",
    "scoring",
    "sheet_writer",
    "submit_gate",
]
