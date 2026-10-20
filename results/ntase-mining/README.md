# Results — structure-first mining of the de novo RNA synthase scaffold

Output of `.claude/skills/catalytic-scaffold-mining/` applied to the catalytic
architecture of OsNTP4 / OsNTP10 / AtURT1 (Li *et al.*, *Vita* 2026, 10.0078).

Full write-up and honesty labelling: [`../../NTASE_REPORT.md`](../../NTASE_REPORT.md).

**Every row in every file is a ⚠️ nominated candidate.** None of these proteins
has been shown to synthesise RNA. Structural similarity plus conservation of an
acidic carboxylate configuration nominates a protein for functional analysis; it
says nothing about its substrates, products or biological role.

## Files

| File | Query | Database | Rows |
|---|---|---|---|
| `urt1_proteome.annotated.tsv` | AtURT1 (PDB 6L8K, catalytic 491/493/547) | `afdb-proteome` + `afdb-swissprot` | 1004 hits / 816 unique proteins |
| `cgas_proteome.candidates.tsv` | human cGAS (PDB 4KM5, catalytic 225/227/319) — the paper's own query | `afdb-proteome` + `afdb-swissprot` | 190 |
| `ARATH_local.candidates.tsv` | AtURT1 | *A. thaliana* proteome (27,402 models), local Foldseek — the paper's own species | 128 |
| `yeast_local.candidates.tsv` | AtURT1 | *S. cerevisiae* proteome, local Foldseek | 18 |
| `SCHPO_local.candidates.tsv` | AtURT1 | *S. pombe* proteome, local Foldseek | 19 |
| `PLAF7_local.candidates.tsv` | AtURT1 | *P. falciparum* proteome, local Foldseek | 13 |
| `metja_local.candidates.tsv` | AtURT1 | *M. jannaschii* proteome, local Foldseek | 83 |

Searches ran with `--window 2`. The four local files carry pLDDT and
side-chain geometry; the two web files cannot (see `reference/data-access.md`).

Two runs demonstrate that the filter is both sensitive and specific.

*A. thaliana* is the end-to-end reproduction, on one of the paper's own three
species: of **27,402** predicted structures, **20** pass the triad + geometry
filter (0.07%), every registration offset zero, catalytic-shell pLDDT 85–97.
Those 20 are the plant NTP family the paper describes — URT1 (the validated
enzyme), HESO1, TRL, AtNTP2/4/6/7, PAPS1/PAPS3, MEE44 and DNA polymerase
lambda. Four of them carry Glu rather than Asp at the third position, exactly
the Glu/Asp interchange the paper allows, and the chemical-class filter keeps
them rather than discarding them.

*S. pombe* is the cleanest specificity check: of 5,196 predicted structures,
exactly 9 pass, and those 9 are precisely the organism's known terminal
nucleotidyltransferases — the whole Cid family (cid1, cid11, cid12, cid13,
cid14, cid16), pla1, and both CCA-adding enzymes. Nothing missing, nothing
spurious.

### Stage 5 outputs

| File | Contents |
|---|---|
| `urt1_clades.clades.tsv` | 72 representative catalytic domains assigned to subfamilies; `with_reference` marks the clade containing the validated enzymes |
| `urt1_clades.tree.nwk` | the dendrogram, Newick |

**This is distance-based clustering on the catalytic domain, not a
maximum-likelihood phylogeny** — no branch support, and a weaker method than
the source paper's MEGA tree. It reproduces the structure (the three validated
enzymes form one tight, plant-restricted clade; pairwise distance 0.31–0.37
against an all-pairs median of 0.766) and ranks the nearest non-plant
relatives, with *Dictyostelium* `Q54Y43` first at 0.586 — ahead of human TUT4
(0.663), *S. pombe* Cid1 (0.668) and human TENT4B (0.736).

## Columns

| Column | Meaning |
|---|---|
| `core_rmsd`, `n_core`, `core_frac` | **rank by this.** Outlier-trimmed RMSD over the catalytic shell, how many residue pairs survived trimming, and what fraction of the shell that is. A tight RMSD over few residues is weak evidence — always read `n_core` alongside. |
| `local_rmsd`, `n_shell` | Untrimmed shell RMSD. Inflated by flexible loops next to the active site; kept for transparency. |
| `n_conserved` / `n_conserved_exact` / `n_catalytic` | Catalytic positions conserved (windowed / at offset 0 only) out of the total. |
| `max_offset` | Largest registration offset used. `0` is strongest; non-zero means the structural alignment was shifted and the counterpart was found nearby. |
| `catalytic_residues` | Explicit per-position mapping, e.g. `D491->D100;D493->D102;A547->D157`, with `[+2]` marking an offset. |
| `site_passed` | All catalytic positions conserved. |
| `tm_score_qnorm` | TM-score over aligned pairs, normalised by query length. Penalises targets much shorter than the query — informative but not the ranking metric. |
| `shell_plddt` | Mean AlphaFold pLDDT **over the catalytic shell only**. Local mode only; empty in web results. Never substitute whole-chain pLDDT. |
| `carboxylate_geometry` | Distance between catalytic carboxylate carbons in the target, with the deviation from the query in parentheses. Local mode only. |
| `lineage_group`, `lineage`, `known_class` | Added by `annotate.py`. `known_class = unknown` means *the annotation is uninformative*, **not** that the protein is novel. |
| `label` | Always `hypothesis`. |

## Reproduce

Commands are listed at the end of [`../../NTASE_REPORT.md`](../../NTASE_REPORT.md).

## Positive controls embedded in these files

`urt1_proteome.annotated.tsv` recovers all three experimentally validated
enzymes in the top 6 by `core_rmsd`, and independently derives OsNTP4's
catalytic triad as **D100 / D102 / D157** — the exact residues the paper reports
and mutates. That reproduction is the reason to trust the rest of the table's
*geometry*; it is not a reason to trust any claim about *activity*.
