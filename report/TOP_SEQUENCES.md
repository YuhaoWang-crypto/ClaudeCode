# Tirzepatide binders — TOP recommended sequences

Target: Tirzepatide (Aib2/Aib13 + K20 C20-diacid-γGlu-(AEEA)₂). Epitope: F22/V23/L26/I27 (+E3, PPPS tail).

| ID | src | AA | topology | key score | sequence | note |
|---|---|---|---|---|---|---|
| **TZP-S1.1** | S1o_short | 66 | 3-helix bundle (opt) | bind(screen)=0.75, fold=0.87 | `SLEELSKLLEEISKLIEEFSKLGSPGEELLKKLEEALKKLEEHLKKLGNSGSEEIRKLAEEIKKLG` | Stable — best of optimization; shortest confident bundle |
| **TZP-S1** | LS3_bundle2 | 75 | 3-helix bundle | protein_iptm=0.96, pTM=0.96 | `SLEELSKLLEEISKLIEEWSKLGSPGSEELLKKLEEALKKLEEHLKKLGSDGTSEEIRKLAEEIKKLAEELKKLG` | Stable — atomistic-confirmed |
| **TZP-S2** | LS5_c | 49 | mixed alpha/beta | protein_iptm=0.97, pTM=0.96 | `KTVTLTVEGKTYTVTVSDGSKLEGSPSLEELSKLLEEISKLIEEFSKLG` | Stable — has beta-sheet(33%)+helix |
| **TZP-S3** | LS4_hairpinF | 47 | alpha-hairpin | protein_iptm=0.97, pTM=0.96 | `SLEELSKLLEEISKLIEEFSKLGPNGKELKKLAEELKKLAEELKKLA` | Stable — compact two-helix |
| **TZP-S4** | LS2_b | 34 | disulfide-stapled helix | protein_iptm=0.96, pTM=0.97 | `SAEELLKKLCELSKLLEEICKLIEEWLKKLEELG` | Stable — stapled (C10-C18), protease-resistant |
| **TZP-P3** | R2_035 | 34 | amphipathic helix | bind(un-forced)=0.41 | `SLSTLENELSTLENEISTIENEWSTGLENEISTG` | Peptide — highest peptide signal |
| **TZP-P1** | R2_008 | 26 | amphipathic helix | bind(un-forced)=0.37 | `SLSTLENELSTLENEISTIENEFSTG` | Peptide — best-balanced 26mer, F-anchor |
| **TZP-P2** | R2_002 | 26 | amphipathic helix | bind(un-forced)=0.37 | `SLSTLENELSTLENEISTIENEWSTG` | Peptide — cleanest, W-anchor, 0 liabilities |
| **TZP-P4** | R2_010 | 26 | amphipathic helix | bind(un-forced)=0.35 | `SLSTLENELSTLENEISTIENEYSTG` | Peptide — Tyr-anchor |
| **TZP-P5** | R2_007 | 26 | amphipathic helix | bind(un-forced)=0.34 | `SLSTWENELSTLENEISTIENEWSTG` | Peptide — di-Trp |
| **TZP-P6** | R2_034 | 21 | amphipathic helix | bind(un-forced)=0.33 | `SLSTLENELSTLENEISTIEG` | Peptide — minimal 21mer, cheapest SPPS |
| **TZP-B1** | design_spec_13 | 96 | mixed alpha/beta mini-protein | bind(un-forced)=0.20 | `MVSKTFTVTLPTGATLTVTVTYDYETNVLTVKVTLPQTGKSDEATFKVAAGVTVETILPGVGLVVKAHYDAEKNTITITLECPALGATFTFTLNLS` | Existing — best legacy binder, orthogonal scaffold |

**Fc-fusion:** append `-GGGGSGGGGSGGGGS-[human IgG1 hinge-CH2-CH3]` to any binder (full sequences in `results/wetlab_constructs.fasta`; IgG4-S228P variant also provided).

**Interpretation:** `protein_iptm` (atomistic binder↔peptide interface) is the cleanest metric; `binding_confidence (un-forced)` is the honest peptide metric; `fold/pTM` = autonomous-fold confidence. Scores are triage signals, confirm by SPR/BLI + CD.
