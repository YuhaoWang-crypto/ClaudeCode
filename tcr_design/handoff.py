"""Hand the shortlist to the tools that can actually rank it.

The scan in this package is a filter. It exists to decide what deserves real
compute. This module writes the input files those real tools expect, so the
expensive stage is a copy-paste away rather than a reimplementation:

* **Boltz-2 / AlphaFold3** - co-fold each mutant complex and read interface
  confidence (ipTM, interface PAE/pLDDT). Needs a GPU.
* **Rosetta flex ddG** - the best-validated physics-based interface ddG
  protocol for point mutations. Resfiles are written here; the protocol
  capture itself comes from the Kortemme lab repository.
* **FoldX BuildModel** - fast empirical ddG, useful as an orthogonal check.
* **ProteinMPNN / LigandMPNN** - inverse folding over the paratope when you
  want loop redesign rather than point mutations.

Nothing here runs those tools. It writes their inputs and records the exact
commands, so the provenance of every downstream number is traceable.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from .liability import mutate_sequence
from .scan import Mutation
from .structure import Complex

# Chain order used in the co-folding inputs. Keeping it stable across mutants
# makes the resulting confidence metrics comparable.
CHAIN_ROLE_ORDER = ("mhc", "b2m", "peptide", "tcr_alpha", "tcr_beta")


@dataclass
class Handoff:
    directory: Path
    written: list[Path]
    commands: dict[str, str]


def _mutant_sequences(complex_obj: Complex, mutation: Mutation) -> dict[str, str]:
    """Chain sequences for a complex carrying one substitution."""
    residue = complex_obj.residue_by_label(mutation.chain, mutation.number)
    if residue is None:
        raise ValueError(f"{mutation.name}: residue not found in structure")
    sequences = dict(complex_obj.sequences)
    sequences[mutation.chain] = mutate_sequence(
        sequences[mutation.chain], residue.index, mutation.mutant
    )
    return sequences


def _ordered_chains(complex_obj: Complex) -> list[tuple[str, str]]:
    ordered = []
    for role in CHAIN_ROLE_ORDER:
        chain_id = complex_obj.chains.get(role)
        if chain_id:
            ordered.append((chain_id, role))
    return ordered


def write_boltz_yaml(
    complex_obj: Complex, mutation: Mutation | None, path: Path
) -> Path:
    """Write a Boltz YAML input for the wild type or one mutant."""
    sequences = (
        complex_obj.sequences if mutation is None else _mutant_sequences(complex_obj, mutation)
    )
    lines = ["version: 1", "sequences:"]
    for chain_id, role in _ordered_chains(complex_obj):
        lines.append("  - protein:")
        lines.append(f"      id: {chain_id}")
        lines.append(f"      sequence: {sequences[chain_id]}")
        lines.append(f"      # role: {role}")
    path.write_text("\n".join(lines) + "\n")
    return path


def write_alphafold3_json(
    complex_obj: Complex, mutation: Mutation | None, path: Path, seeds: tuple[int, ...] = (1,)
) -> Path:
    """Write an AlphaFold3 JSON input for the wild type or one mutant."""
    sequences = (
        complex_obj.sequences if mutation is None else _mutant_sequences(complex_obj, mutation)
    )
    name = complex_obj.pdb_id.lower() if mutation is None else mutation.name.replace(":", "_")
    payload = {
        "name": name,
        "modelSeeds": list(seeds),
        "sequences": [
            {"protein": {"id": chain_id, "sequence": sequences[chain_id]}}
            for chain_id, _ in _ordered_chains(complex_obj)
        ],
        "dialect": "alphafold3",
        "version": 1,
    }
    path.write_text(json.dumps(payload, indent=2) + "\n")
    return path


def write_rosetta_resfile(mutation: Mutation, path: Path) -> Path:
    """Write a Rosetta resfile restricting design to the single mutated position.

    Resfile grammar: a ``NATAA`` default, a ``start`` line, then one record per
    designed position as ``<resnum> <chain> PIKAA <aa>``.
    """
    lines = [
        f"# {mutation.name}: {mutation.wild_type}{mutation.number} -> {mutation.mutant}",
        "NATAA",
        "start",
        f"{mutation.number} {mutation.chain} PIKAA {mutation.mutant}",
    ]
    path.write_text("\n".join(lines) + "\n")
    return path


def write_foldx_individual_list(mutations: list[Mutation], path: Path) -> Path:
    """Write a FoldX BuildModel mutation list.

    One mutation per line as ``<wt><chain><resnum><mut>;`` which is the format
    FoldX ``individual_list.txt`` expects.
    """
    lines = [
        f"{m.wild_type}{m.chain}{m.number}{m.mutant};" for m in mutations
    ]
    path.write_text("\n".join(lines) + "\n")
    return path


def write_proteinmpnn_spec(
    complex_obj: Complex, mutations: list[Mutation], directory: Path
) -> list[Path]:
    """Write ProteinMPNN chain and fixed-position inputs for paratope redesign.

    ProteinMPNN designs the positions that are *not* fixed, so every position
    except the shortlisted paratope ones is pinned. Indices are 1-based within
    each chain, which is what the ``--fixed_positions_jsonl`` loader expects.
    """
    key = complex_obj.pdb_id.lower()
    designable: dict[str, list[int]] = {}
    for mutation in mutations:
        residue = complex_obj.residue_by_label(mutation.chain, mutation.number)
        if residue is None:
            continue
        designable.setdefault(mutation.chain, [])
        if residue.index + 1 not in designable[mutation.chain]:
            designable[mutation.chain].append(residue.index + 1)
    for chain_id in designable:
        designable[chain_id].sort()

    tcr_chains = complex_obj.tcr_chain_ids()
    fixed: dict[str, list[int]] = {}
    for chain_id in tcr_chains:
        length = len(complex_obj.sequences[chain_id])
        free = set(designable.get(chain_id, []))
        fixed[chain_id] = [i for i in range(1, length + 1) if i not in free]

    written = []
    fixed_path = directory / "proteinmpnn_fixed_positions.jsonl"
    fixed_path.write_text(json.dumps({key: fixed}) + "\n")
    written.append(fixed_path)

    chains_path = directory / "proteinmpnn_chains.jsonl"
    designed = [c for c in tcr_chains]
    context = [c for c in complex_obj.partner_chain_ids()]
    chains_path.write_text(json.dumps({key: [designed, context]}) + "\n")
    written.append(chains_path)

    readable_path = directory / "designable_positions.json"
    readable_path.write_text(
        json.dumps(
            {
                "designable_positions_1based_per_chain": designable,
                "designed_chains": designed,
                "fixed_context_chains": context,
                "note": (
                    "Designed chains are redesigned in the context of the fixed "
                    "chains. Keep the pMHC chains as context so the model sees "
                    "the epitope."
                ),
            },
            indent=2,
        )
        + "\n"
    )
    written.append(readable_path)
    return written


TIER2_README = """# Tier 2: running the shortlist through real models

The candidates in this directory were selected by an **uncalibrated structural
triage heuristic**. Nothing here is an affinity prediction yet. These are the
inputs for the tools that can produce one.

Run at least two orthogonal methods and treat disagreement as the signal it is.

## 1. Co-folding (needs a GPU)

AlphaFold3 is the better choice for TCR-pMHC and should be run first if you
have access to the weights:

```bash
python run_alphafold.py --json_path=alphafold3/{first_af3} \\
  --model_dir=$AF3_WEIGHTS --output_dir=af3_out
```

Boltz inputs are provided as an open-weights fallback:

```bash
boltz predict boltz/ --use_msa_server --out_dir boltz_out
```

**Read this before trusting either.** Two independent 2025 benchmarks on
TCR-pMHC complexes, both including structures deposited after the models'
training cutoffs, found AlphaFold3 clearly ahead, with roughly 90% of models
at acceptable quality or better. The same benchmark reported that Boltz-2 and
tFold-TCR produced **zero** medium- or high-quality models across 20
complexes: they recover the fold (TCR interface RMSD under 2 A) but place the
receptor in the wrong docking geometry, which is precisely the thing interface
engineering depends on. Treat Boltz output here as a structural sanity check,
not as a ranking signal.

Second caveat: **co-folding confidence is not affinity.** ipTM measures
modelling self-consistency, not Kd. A Boltz-2 affinity-head adaptation scored
about r = 0.15 against TCR3d affinity data, and a 2026 analysis argues that the
apparent predictive power of AlphaFold-derived scores tracks how the TCR is
positioned over the pMHC rather than the chemistry of the contacts. Use
co-folding to discard mutants that break the interface, not to rank survivors.

## 2. Rosetta flex ddG (CPU, best-validated physics option)

Resfiles are in `rosetta/`. The protocol capture is not vendored here; get it
from the Kortemme lab repository (`flex_ddG_tutorial`) and point it at the
structure plus each resfile. Expect roughly 35 backrub models per mutation and
a few CPU-hours per mutation for converged numbers.

flex ddG was benchmarked on 1,240 mutations of general protein-protein
interfaces, with its largest gains on small-to-large substitutions and on
antibody-antigen interfaces. **No TCR-specific flex ddG benchmark appears to
exist**, so its accuracy on these interfaces is an extrapolation. TCR-pMHC
interfaces are unusually flat and low-affinity, which is the regime where ddG
methods tend to do worst.

## 3. MM/PBSA (the best TCR-specific energetics number available)

If you can run MD, this is the method with real TCR validation behind it.
Benchmarked against measured mutant affinities in the ATLAS database it reached
Pearson r = 0.9 across 7 mutants with solved structures and r = 0.74 across 35
modelled mutants, with balanced accuracy 0.79 for calling ddG above
0.5 kcal/mol. Those sample sizes are small, and the protocol omitted explicit
interface waters and entropy corrections.

## 4. FoldX BuildModel (fast, orthogonal, less accurate)

```bash
foldx --command=BuildModel --pdb={pdb_file} \\
  --mutant-file=foldx/individual_list.txt --numberOfRuns=5
```

Repair the structure first (`--command=RepairPDB`). Average over runs. FoldX is
quicker than flex ddG and noticeably less accurate on charged and buried-polar
substitutions.

## 5. Inverse folding for loop redesign

For redesigning whole CDR loops rather than single positions:

```bash
python protein_mpnn_run.py --pdb_path {pdb_file} \\
  --pdb_path_chains "{designed_chains}" \\
  --chain_id_jsonl proteinmpnn_chains.jsonl \\
  --fixed_positions_jsonl proteinmpnn_fixed_positions.jsonl \\
  --num_seq_per_target 64 --sampling_temp 0.2 --out_folder mpnn_out
```

On TCR-pMHC interfaces specifically, ProteinMPNN and ESM-IF1 recover 44% and
50% of native CDR3-interface residues respectively, and removing the pMHC from
the input costs about 20% of that, so always keep the pMHC chains as context.
That benchmark had **no wet-lab validation**, and it found that randomly chosen
dissimilar sequences sometimes still scored high model confidence. Confidence
alone is not a filter. ProteinMPNN also optimises sequence for a fixed
backbone, so it will propose sequences that need a different loop conformation
than the one it was given; re-fold anything it suggests.

## 6. The step that is not optional

Whatever survives, screen it for cross-reactivity properly. The structural
screen in this package is a cheap first pass on one crystal structure. The
established cascade is:

1. Determine the recognition motif of the *actual engineered receptor*
   experimentally: a glycine/alanine or combinatorial X-scan, or randomised-CDR3
   yeast selection. Residue-importance rankings from such selections correlate
   with a Rosetta computational alanine scan, so run that too.
2. Search the motif broadly against the human proteome. CrossDome's
   TCR-centred mode, which weights positions by TCR contact hotspots in the
   same spirit as the screen in this package, raised top-50 enrichment of
   validated cross-reactive peptides from 63% to 82% across 16 known cases, and
   moved the titin peptide from rank 27 to rank 6 of 36,000 in the MAGE-A3
   case. Extend the search to common SNP variants and layer tissue expression.
3. Run an unbiased experimental screen. A yeast-display pHLA library of around
   10^8 peptides is the only method that has prospectively caught a peptide
   that killed patients.
4. Confirm function in physiologically realistic cells from the at-risk
   tissue. In the MAGE-A3 case titin was undetectable in cardiac-derived cell
   lines and absent from mouse models; it showed up only in beating
   iPSC-derived cardiomyocytes.

**Do not substitute a sequence-based TCR-epitope predictor for any of this.**
On the task that off-target screening actually requires, scoring the effect of
single amino-acid changes in a peptide, these models perform at chance: across
eight predictors tested on 21 TCRs against 172 single-residue mutants of one
epitope, none achieved a significant correlation with measured binding.
"""


def write_handoff(
    complex_obj: Complex,
    shortlist: list[Mutation],
    directory: Path,
    include_wild_type: bool = True,
) -> Handoff:
    """Write every tier-2 input for a shortlist of mutations."""
    directory.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    boltz_dir = directory / "boltz"
    af3_dir = directory / "alphafold3"
    rosetta_dir = directory / "rosetta"
    foldx_dir = directory / "foldx"
    for sub in (boltz_dir, af3_dir, rosetta_dir, foldx_dir):
        sub.mkdir(parents=True, exist_ok=True)

    if include_wild_type:
        written.append(write_boltz_yaml(complex_obj, None, boltz_dir / "wild_type.yaml"))
        written.append(
            write_alphafold3_json(complex_obj, None, af3_dir / "wild_type.json")
        )

    first_af3 = "wild_type.json"
    for mutation in shortlist:
        stem = mutation.name.replace(":", "_")
        written.append(write_boltz_yaml(complex_obj, mutation, boltz_dir / f"{stem}.yaml"))
        written.append(
            write_alphafold3_json(complex_obj, mutation, af3_dir / f"{stem}.json")
        )
        written.append(write_rosetta_resfile(mutation, rosetta_dir / f"{stem}.resfile"))
    if shortlist:
        first_af3 = shortlist[0].name.replace(":", "_") + ".json"
        written.append(
            write_foldx_individual_list(shortlist, foldx_dir / "individual_list.txt")
        )
    written += write_proteinmpnn_spec(complex_obj, shortlist, directory)

    designed_chains = " ".join(complex_obj.tcr_chain_ids())
    readme = directory / "RUN_TIER2.md"
    readme.write_text(
        TIER2_README.format(
            first_af3=first_af3,
            pdb_file=complex_obj.path.name,
            designed_chains=designed_chains,
        )
    )
    written.append(readme)

    commands = {
        "boltz": "boltz predict boltz/ --use_msa_server --out_dir boltz_out",
        "alphafold3": (
            f"python run_alphafold.py --json_path=alphafold3/{first_af3} "
            "--model_dir=$AF3_WEIGHTS --output_dir=af3_out"
        ),
        "foldx": (
            f"foldx --command=BuildModel --pdb={complex_obj.path.name} "
            "--mutant-file=foldx/individual_list.txt --numberOfRuns=5"
        ),
        "proteinmpnn": (
            f"python protein_mpnn_run.py --pdb_path {complex_obj.path.name} "
            f'--pdb_path_chains "{designed_chains}" '
            "--chain_id_jsonl proteinmpnn_chains.jsonl "
            "--fixed_positions_jsonl proteinmpnn_fixed_positions.jsonl "
            "--num_seq_per_target 64 --sampling_temp 0.2 --out_folder mpnn_out"
        ),
    }
    return Handoff(directory=directory, written=written, commands=commands)
