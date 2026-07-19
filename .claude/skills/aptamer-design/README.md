# aptamer-design — Claude Code Skill

Structure-informed DNA/RNA aptamer design + Boltz-2.1 in-silico ranking, packaged as a
reusable Claude Code skill.

## Install

**Personal (all your projects):**
```
unzip aptamer-design.zip -d ~/.claude/skills/
# -> ~/.claude/skills/aptamer-design/SKILL.md
```

**Project (share via repo):**
```
unzip aptamer-design.zip -d .claude/skills/
# -> .claude/skills/aptamer-design/SKILL.md
```

Restart / reload Claude Code. It auto-discovers the skill from `SKILL.md`'s frontmatter
and invokes it when you ask to design or rank aptamers.

## Use
Just ask, e.g.:
> "Design DNA/RNA aptamers for <protein>, it's for a diagnostic probe."

Or run the ranker standalone on Boltz metrics:
```
python3 scripts/rank_candidates.py examples/gfra1_candidates.json --use-case diagnostic
```

## Requirements
- **Boltz_API MCP server** connected (for co-folding; free estimates, ~$0.025/paid run).
- Web access for UniProt / RCSB. Optional: PubMed / bioRxiv MCP.
- Python 3 (standard library only) for the ranking script.

## Contents
```
aptamer-design/
  SKILL.md                         # entry point (frontmatter + workflow)
  README.md                        # this file
  references/
    design-principles.md           # scaffolds, SELEX library, chemistry
    boltz-workflow.md              # Boltz-2.1 co-folding how-to
    metric-interpretation.md       # ipTM/pLDDT meaning + mandatory caveats
    t-selex-integration.md         # UPGRADED: T-SELEX × Boltz consensus in-silico SELEX
  templates/
    report_template.md             # deliverable skeleton
  scripts/
    rank_candidates.py             # use-case-weighted composite ranker (quick mode)
    consensus_rank.py              # HDOCK ⊕ Boltz Borda consensus + specificity gate
    tselex_boltz_pipeline.py       # orchestrator scaffold for the full SELEX loop
  model/                           # data integration + baseline LM (see model/README.md)
  examples/
    gfra1_candidates.json          # worked example (GFRα1, 16 candidates)
    egfr_boltz_results.json        # EGFR (counter-screen killed all — no specific hit)
    pdl1_consensus.json            # PD-L1 consensus demo (RNA-HP is PD-L2-selective lead)
```

## Scope & honesty
This skill produces a **computationally prioritized starting pool for SELEX / binding
assays**, not validated binders. ipTM/pLDDT are relative confidence signals, not
affinities. See `references/metric-interpretation.md`.
