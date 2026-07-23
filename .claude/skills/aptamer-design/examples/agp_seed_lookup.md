# Seed lookup: can we start AGP / α1-acid glycoprotein (PDB 3KQ0) from a database?

**Question:** the CTLA-4 campaign showed the *real* aptamer `aptamerd6` is promiscuous.
Rather than cold-start, is there a database seed we can recognize for AGP (ORM1, P02763, 3KQ0)?

## Run

```bash
python3 model/seed_from_db.py --target "alpha-1-acid glycoprotein" \
        --uniprot P02763 --family "lipocalin,orosomucoid,ORM1,ORM2,acute-phase glycoprotein"
```

```
VERDICT: NO_SEED
No literature aptamer for this target OR its family in the corpus. -> fall back to de novo.
```

## Why NO_SEED is the correct, honest answer

A literature check (2026) finds **no published DNA/RNA aptamer against AGP / orosomucoid**.
SELEX aptamers exist for many plasma proteins (thrombin, VEGF, IgE, PDGF, nucleolin…) and
for lookalike acute-phase / lipocalin targets, but **not for AGP**. So the seed corpus has
nothing to return — and the tool says so instead of fabricating a hit. That absence is itself
information: AGP has no validated nucleic-acid binder to lean on.

This contrasts cleanly with the two other verdict classes the mode produces:

| target | verdict | what the DB gives us |
|---|---|---|
| **Thrombin** (P00734) | `DIRECT_SEED` | TBA `GGTTGGTGTGGTTGG` + HD22 — real seeds, plus a **PDB anchor** (1HAO) for calibration |
| **CTLA-4** (P16410) | `DIRECT_METADATA_ONLY` | `aptamerd6` lineage exists but is tagged **NEGATIVE-ANCHOR** (this skill cross-validated it as promiscuous) → calibrate, don't seed |
| **AGP** (P02763) | `NO_SEED` | nothing → genuine cold-start |

## Consequence for AGP

`seed_from_db` sends AGP straight to **de novo** design — which is exactly the campaign already
run in `agp_campaign_report.md`, and which **failed the decoy gate** (best on-target ipTM 0.845
< scramble decoy 0.862; AGP's glycan-decorated, sticky β-barrel is an electrostatic sink). So:

- **There is no database shortcut for AGP.** The warm-start path is empty; the honest route is
  wet-lab SELEX with counter-selection vs ORM2 + serum albumin and attention to the glycan shield.
- The `NO_SEED` result **explains** the de novo failure rather than contradicting it: no prior
  binder exists because AGP is a genuinely hard, under-explored aptamer target — not because the
  design was under-powered.

## Takeaway

The database-seed mode helps **only when a seed exists**. Its value on a hard target like AGP is
negative-space honesty: it tells you the cupboard is bare and routes you to de novo + the gate,
instead of dressing up a random pool as "database-backed."
