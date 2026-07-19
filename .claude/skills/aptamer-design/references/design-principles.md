# Aptamer design principles (scaffolds, library, chemistry)

## Length & composition
- Typical aptamer length: **26–45 nt**. Shorter = cheaper/more stable; longer = larger
  interface potential.
- GC ~40–62%. Avoid homopolymer runs > 4 nt (except intended G-tracts).

## Scaffold families (design several; diversity matters)

### G-quadruplex (GQ) — top choice for diagnostics
- Four G-tracts of ≥3 G separated by short loops: `G≥3 N1-3 G≥3 N1-3 G≥3 N1-3 G≥3`.
- Compact, rigid, **intrinsically nuclease-resistant and thermostable**; presents a flat
  tetrad face + diversifiable loops.
- Example core: `GGGAGGGTGGGAGGG` + a variable tail/loop.

### Hairpin / stem-loop
- 6–8 bp stem + 12–18 nt apical loop. The **loop presents the recognition surface**;
  randomize the loop in SELEX. Verify stem base-pairing (reverse-complement) explicitly.

### Dual-hairpin / three-way junction
- Two stem-loops off a junction → larger, more shape-specific interface. Good when a
  single loop is too small (small interfaces score noisily).

### Pseudoknot
- Interleaved base pairing → rigid 3D surface; highest fold confidence in practice but
  harder/costlier to synthesize (better for therapeutic RNA than cheap diagnostics).

## Epitope-directed biasing
- **Polyanion complementarity**: aptamers are polyanionic (phosphate backbone), like
  heparan sulfate. A basic (Arg/Lys) or heparin/GAG groove on the target is a natural
  docking site — bias compact folds (GQ) toward it.
- **Ligand-competitive**: to block a protein–protein interaction, target the partner's
  binding footprint (from a complex structure). *For diagnostics this is usually NOT
  wanted* — endogenous ligand competes with the probe, and electrostatic grooves are
  cross-reactive.

## SELEX-biased starting library
Full random N40 converges slowly. Bias it:
```
5'-[fwd primer 18nt]-NNN-[structured core, e.g. GGGAGGGN]-N12-[complementary arm]-NNN-[rev primer 18nt]-3'
```
- Embed a pre-folding structural motif; keep an inner fully-random region for diversity.
- Library complexity ~4^15–4^20.
- Run 8–12 SELEX rounds with **counter-SELEX** against paralogs/off-targets for specificity.

## Chemistry & stabilization
| Modification | Purpose |
|---|---|
| 2'-F / 2'-OMe pyrimidines | nuclease resistance (essential for RNA in vivo) |
| 3'-inverted dT | blocks 3' exonuclease |
| 5'-PEG / cholesterol / 40 kDa PEG | half-life extension |
| phosphorothioate (partial) | backbone stability |
| G-quadruplex scaffold | intrinsic nuclease resistance, minimal mods needed |

## Common pitfalls
- Small predicted interfaces (< ~10 residues) with "100% on-site" overlap are **small-denominator
  artifacts** — down-weight them.
- A single MFE structure ignores conformational ensembles and Mg²⁺-dependent tertiary folds.
- Electrostatic/heparin-groove binders can be promiscuous — always counter-screen.
