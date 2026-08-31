# `binder_campaign` — a reference implementation of Anthropic's protein-binder-design campaign prompt

Anthropic released the prompts from its *de novo* miniprotein binder design
campaign as the HuggingFace dataset
[`Anthropic/claude-protein-binder-design`](https://huggingface.co/datasets/Anthropic/claude-protein-binder-design)
(CC BY 4.0). This directory implements the mechanically-specified core of that
prompt as tested Python.

## What was downloaded

Everything the prompt release bundle contains except the 1.16 GB external
reference corpus, which is a convenience copy of third-party sources (the
`source_url` in its index is authoritative anyway):

```
prompts/
  multi_target_binder_design_prompt.md   the 14-target, 48 h / $50k campaign prompt (114 KB, 288 paragraphs)
  single_target/<TARGET>.md              16 single-target, 24 h / $10k prompts
  kickoff/{multi,single}_target_kickoff.md   the T0 user messages
  Figure 1.jpg, Figure 2.jpg             the LCP definition and LCP results figures
  README.md, REDISTRIBUTION_NOTES.csv, MANIFEST.sha256
```

The corpus (`prompts/protein_binder_design_prompts_release.zip`, 1.16 GB, 317
files) was **not** downloaded. Fetch it with:

```bash
curl -L -O https://huggingface.co/datasets/Anthropic/claude-protein-binder-design/resolve/main/prompts/protein_binder_design_prompts_release.zip
```

## What "implement" can and cannot mean here

The prompt describes a 48-hour, $50,000, >100-concurrent-agent campaign that
generates 10^5–10^6 designs on ~325 concurrent H100s, orders 30 designs per
target from two CROs, and posts to Slack and Google Drive from a bespoke
"Claude Science" agent harness (`host.delegate`, `host.compute`,
`wait_for_notification`, …). None of that can be executed here: there is no
Modal account, no GPU, no harness, no CRO, no 48-hour window.

But a large fraction of the prompt is not biology at all — it is **exactly
specified mechanism**, and that part is implementable and testable:

| Prompt requirement | Module | Tests |
|---|---|---|
| LCP restraint, "implement the per-position penalty exactly as defined in Figure 1" | `lcp.py` | `test_lcp.py` |
| The three-arm instrument: `ipSAE_min`, `sc_DockQ`, `final_score`, `rank_zscore`, reduced masks, selectivity | `scoring.py` | `test_scoring.py` |
| `submit_gate()` steps (a)–(g) + (d2), fail-closed | `submit_gate.py` | `test_submit_gate.py` |
| Governor arithmetic, pace band, calibration ladder, WATCHDOG fail-safes | `governor.py` | `test_governor.py` |
| Design-count / job-metadata / deviations ledger tree | `ledger.py` | `test_ledger.py` |
| The four pre-scoring gates (novelty, liability, foldability, plausibility) | `filters.py` | `test_filters.py` |
| `/state/gates/{target}.json` validation gate | `gates.py` | `test_deliverables.py` |
| Frozen `method_vocab.json` / `sheet_schema.json` | `schema.py` | `test_deliverables.py` |
| Selection caps (a)–(f), rank order, relaxation ladder, write-time gate recompute | `sheet_writer.py` | `test_sheet_writer.py` |
| `per_seed_metrics` / `instrument_realization` companions | `companions.py` | `test_deliverables.py` |
| The 17-column scoreboard | `scoreboard.py` | `test_deliverables.py` |

Everything that needs a GPU or a cloud account is behind an injected
interface — `HomologySearcher` (MMseqs2), `StateReader`, `FleetCounter` (Modal
SDK), and the `SeedRecord` stream that a real ESMFold2 / Protenix run would
produce. Swap in real implementations and the same invariants hold.

## Run it

```bash
pip install numpy pandas pyarrow pytest
python3 -m pytest tests/ -q          # 153 tests
python3 -m binder_campaign.demo      # end-to-end dry run -> demo_out/
```

The demo runs two targets (PD-L1, and GDF-8 with its GDF-11 counter-screen)
through the full close-out path — validation gate, pre-scoring filters, the
three-arm instrument, the sheet writer, both companions, the scoreboard — and
writes real `design_sheet.csv` (30 ranked rows per target),
`per_seed_metrics.parquet`, `instrument_realization.csv`, `scoreboard.csv`,
`deviations.jsonl`, and the frozen schema/vocab/gate/governor JSON.

## Notes on faithfulness

Three places where the prompt needed interpretation, resolved and documented in
the code:

1. **LCP is defined per *window*, the prompt asks for a per-*position*
   penalty.** `lcp_score` is Figure 1's `C₃` verbatim;
   `lcp_per_position_penalty` spreads each window's penalty over the residues
   it covers, and sums back to `C₃`.

2. **Order of the eligibility routing and the write-time recompute.** The
   prompt says a mismatch between carried and recomputed values "halts the
   writer", *and* that rows with missing or zero scores are "automatically"
   routed to an unranked section. Routing runs first: a *missing* measurement
   goes to the unranked set, a *present but wrong* value halts. Recomputing
   first would turn every missing value into a campaign halt.

3. **Novelty relaxation is not automated.** The ladder implements steps (i)
   diversity caps and (ii) liability flags. `NOVELTY_LAST_RESORT` is reachable
   only by the caller explicitly supplying novelty-flagged rows, because the
   prompt makes it conditional on regeneration having been tried first and
   keeps the target-mimic and natural-sequence-copy bans absolute.

Two prompt behaviours are deliberately **not** implemented, because they are
measurements rather than rules: TM-score and DockQ themselves (both take
structures and belong to DockQ/TM-align), and the MMseqs2 UniRef90 search. The
*decision rules* on top of them are implemented and tested, including a real
affine-gap Smith-Waterman for the "≥30 % gapped local identity over ≥40 aligned
residues" rule the prompt calls out for catching Ubiquitin with terminal
extensions.

## Licence

The prompt files under `prompts/` are Anthropic's, CC BY 4.0. This
implementation is original work.
