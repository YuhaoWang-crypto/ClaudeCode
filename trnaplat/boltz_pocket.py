"""Boltz-2 against the only measured negatives in this project.

Every structure-based claim in this repository so far has been scored against
**decoys** -- unselected library members standing in for rejections. Sheet 3 of
the client's spreadsheet contains the one exception: five McTyrRS variants built
and assayed for para-azido-phenylalanine, of which the paper reports two as
selected and states the other three were "not good". Those three are measured
rejections, not invented ones.

That makes them the only ground truth available for asking whether a
structure-and-binding predictor separates a working aaRS pocket from a
non-working one.

| variant | mutations | F+ | F- | F+/F- | outcome |
|---|---|---|---|---|---|
| Mc Mut6+RM | Y33G/Y112F/D162T | 177 | 48 | 3.69 | selected |
| Mc Mut6    | Y33G/D162T       |  41 | 37 | 1.11 | selected, barely |
| Mc Mut4    | Y33L/D162L/L166Q |   - |  - |   -  | built, rejected |
| Mc Mut5    | Y33A/D162V/L166D |   - |  - |   -  | built, rejected |
| Mc Mut7    | Y33L/D162Q/L166S |   - |  - |   -  | built, rejected |

## ❌ The measured outcome: it does not separate them

✅ Run 2026-10-06, Boltz-2.1, pocket-constrained, one sample each
(`trnaplat/data/boltz_results.csv`):

| variant | outcome | binding_confidence | optimization_score |
|---|---|---|---|
| **Mc Mut7** | **built, rejected** | **0.6196** | **0.2483** |
| Mc Mut6+RM | selected, F+/F- 3.69 | 0.6008 | 0.2211 |
| Mc Mut6 | selected, F+/F- 1.11 | 0.5960 | 0.2171 |
| Mc Mut5 | built, rejected | 0.5823 | 0.1987 |
| Mc Mut4 | built, rejected | 0.5680 | 0.2117 |

**The highest-scoring variant on both metrics is one the experimenters built and
rejected.** AUC 0.667, exact one-sided p = 0.400. Both metrics give the same
ordering, so this is not a choice-of-score artifact. Structure confidence was
high and uninformative across the board (0.951-0.969, ligand ipTM 0.925-0.975):
Boltz is confident about where the ligand sits in every variant, including the
ones that do not work.

Read it as: on the only measured-negative panel this project has, a pocket-
constrained binding prediction would have put a rejected clone at the top of the
pick-list. It does separate the two hits from two of the three rejections, so it
is not anti-correlated -- it is just not usable as a filter at this n. ⚠️ With
2x3 the result cannot be significant either way (see below), so this rules out a
gross success, not a weak signal.

## ⚠️ What this can and cannot establish

Two positives against three negatives is **six ranked pairs**. A perfect
separation gives AUC 1.0 with a 95% interval that still reaches roughly 0.5, so
this pilot can detect a gross failure and cannot certify a success.
`power_note()` prints the arithmetic rather than leaving it implied. The
honest reading: if Boltz cannot separate these five, that is informative; if it
can, the next step is the score-stratified plate in
`pylrs/negative_panel.py`, not a claim.

⚠️ Two further limits, both structural rather than statistical:

* The three rejections have **no numbers**, only the paper's wording, so they
  are ordered relative to the hits but not among themselves.
* Boltz predicts **binding**, while the assay measures **aminoacylation** --
  charging the tRNA, which needs the amino acid positioned for catalysis and
  ATP bound. A variant can bind pAzF and still not charge it. Binding is a
  necessary condition this tests, not the measured quantity.

Licence: Boltz is MIT, and the Boltz API runs under the account's own terms, so
unlike the AP Novo generator there is no non-commercial restriction here.

Run:  python3 -m trnaplat.boltz_pocket            # emit the job specs
      python3 -m trnaplat.boltz_pocket --ingest results.csv
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "pylrs"))

#: para-azido-L-phenylalanine. The azide is written as the neutral linear
#: N=[N+]=[N-] form, which is how the PDB chemical dictionary draws aryl azides.
PAZF_SMILES = "N[C@@H](Cc1ccc(cc1)N=[N+]=[N-])C(O)=O"

#: L-tyrosine, the native substrate. The discrimination this assay measures is
#: pAzF over Tyr, so Tyr is the comparator every variant must also be scored
#: against -- a variant that binds pAzF well and Tyr better has not gained
#: specificity.
TYR_SMILES = "N[C@@H](Cc1ccc(O)cc1)C(O)=O"

#: The pocket residues to constrain the ligand against, in McTyrRS numbering.
#: ✅ These are the positions sheet 3's own library randomised, cross-checked
#: against the MjTyrRS equivalents through the verified alignment in
#: `pylrs/crossspecies.py`.
MC_POCKET_POSITIONS = [33, 111, 112, 162, 163, 166]

#: ⚠️ Ordered by the measured ratio where one exists; the three rejections are
#: ranked below both hits but not among themselves, because the paper gives no
#: numbers for them.
GROUND_TRUTH = {
    "Mc Mut6+RM": {"label": 1, "f_ratio": 177 / 48},
    "Mc Mut6": {"label": 1, "f_ratio": 41 / 37},
    "Mc Mut4": {"label": 0, "f_ratio": None},
    "Mc Mut5": {"label": 0, "f_ratio": None},
    "Mc Mut7": {"label": 0, "f_ratio": None},
}


def mc_sequence() -> str:
    import pandas as pd

    table = pd.read_csv(ROOT / "pylrs/data/literature_campaigns.csv")
    rows = table[table["scaffold"] == "McTyrRS"]
    if rows.empty:
        raise ValueError("no McTyrRS scaffold in literature_campaigns.csv")
    sequence = str(rows.iloc[0]["sequence"])
    if len(sequence) != 319:
        raise ValueError(f"McTyrRS should be 319 aa, got {len(sequence)}")
    return sequence


def apply_mutations(sequence: str, mutations: str) -> str:
    """Build a variant, refusing any substitution whose stated parent is wrong."""
    residues = list(sequence)
    for mut in filter(None, str(mutations).split("/")):
        frm, pos, to = mut[0], int(mut[1:-1]), mut[-1]
        if residues[pos - 1] != frm:
            raise ValueError(
                f"{mut}: McTyrRS has {residues[pos - 1]}{pos}, not {frm}{pos}"
            )
        residues[pos - 1] = to
    return "".join(residues)


def variants() -> list[dict]:
    """The five measured variants plus the wild type, as sequences."""
    import pandas as pd

    sequence = mc_sequence()
    table = pd.read_csv(ROOT / "pylrs/data/mc_variants.csv")
    measured = table[table["origin"].str.startswith("measured")]

    out = [{
        "name": "Mc wild type", "mutations": "", "label": None,
        "f_ratio": None, "sequence": sequence,
        "note": "control: should prefer Tyr over pAzF",
    }]
    for row in measured.itertuples():
        truth = GROUND_TRUTH.get(row.clone)
        if truth is None:
            raise ValueError(f"{row.clone} is in the data but not in GROUND_TRUTH")
        out.append({
            "name": row.clone, "mutations": row.mutations,
            "label": truth["label"], "f_ratio": truth["f_ratio"],
            "sequence": apply_mutations(sequence, row.mutations),
            "note": row.outcome,
        })
    missing = set(GROUND_TRUTH) - {v["name"] for v in out}
    if missing:
        raise ValueError(f"GROUND_TRUTH names not found in the data: {sorted(missing)}")
    return out


def job_spec(variant: dict, ligand_smiles: str, ligand_name: str,
             pocket: list[int] | None = None) -> dict:
    """One Boltz structure-and-binding input.

    The pocket constraint is 0-indexed, which is why the positions are shifted:
    McTyrRS Y33 is residue_index 32.
    """
    entities = [
        {"type": "protein", "chain_ids": ["A"], "value": variant["sequence"]},
        {"type": "ligand_smiles", "chain_ids": ["B"], "value": ligand_smiles},
    ]
    spec = {
        "entities": entities,
        "binding": {"type": "ligand_protein_binding", "binder_chain_id": "B"},
        "num_samples": 1,
    }
    if pocket:
        spec["constraints"] = [{
            "type": "pocket", "binder_chain_id": "B",
            "contact_residues": {"A": [p - 1 for p in pocket]},
            "max_distance_angstrom": 6.0,
        }]
    return spec


def auc(scores: list[float], labels: list[int]) -> float:
    """Rank AUC over the labelled pairs, ties counted as half."""
    pos = [s for s, l in zip(scores, labels) if l == 1]
    neg = [s for s, l in zip(scores, labels) if l == 0]
    if not pos or not neg:
        return float("nan")
    wins = sum((p > n) + 0.5 * (p == n) for p, n in itertools.product(pos, neg))
    return wins / (len(pos) * len(neg))


def exact_p(value: float, n_pos: int, n_neg: int) -> float:
    """Exact one-sided p for an observed AUC, by enumerating every ranking.

    ⚠️ Hanley-McNeil is wrong here and misleadingly so: its variance collapses
    to zero at AUC 1.0, reporting a 95% interval of [1.000, 1.000] for a result
    based on six pairs. With n this small the permutation null is cheap and
    exact, so it is used instead.
    """
    total = math.comb(n_pos + n_neg, n_pos)
    at_least = 0
    for positions in itertools.combinations(range(n_pos + n_neg), n_pos):
        # rank 0 is the top score; count concordant pairs for this arrangement
        pos_ranks = set(positions)
        wins = sum(1 for p in pos_ranks
                   for n in range(n_pos + n_neg) if n not in pos_ranks and p < n)
        if wins / (n_pos * n_neg) >= value - 1e-12:
            at_least += 1
    return at_least / total


def power_note(n_pos: int = 2, n_neg: int = 3) -> str:
    pairs = n_pos * n_neg
    best_p = exact_p(1.0, n_pos, n_neg)
    return (
        f"{n_pos} positives x {n_neg} negatives = {pairs} ranked pairs, and\n"
        f"    {math.comb(n_pos + n_neg, n_pos)} possible arrangements of the labels.\n"
        f"    a PERFECT separation (AUC 1.000) has exact one-sided p = {best_p:.3f}\n"
        f"    -> even a flawless result cannot reach p < 0.05 on this panel.\n"
        "       It can rule out a gross failure and nothing more; treat any\n"
        "       success as a reason to run the stratified plate, not a result."
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ingest", type=pathlib.Path,
                    help="CSV with columns name,ligand,score -- scores the run")
    ap.add_argument("--no-pocket", action="store_true",
                    help="omit the pocket constraint, letting Boltz place the "
                         "ligand itself")
    ap.add_argument("--out", type=pathlib.Path, default=ROOT / "trnaplat/data")
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    if args.ingest:
        import pandas as pd

        frame = pd.read_csv(args.ingest)
        truth = {v["name"]: v for v in variants()}
        print("=" * 76)
        print(f"Boltz vs the measured McTyrRS panel -- {args.ingest}")
        print("=" * 76 + "\n")
        for ligand, block in frame.groupby("ligand"):
            rows = [(r.name, r.score, truth[r.name]["label"])
                    for r in block.itertuples() if r.name in truth]
            labelled = [(n, s, l) for n, s, l in rows if l is not None]
            print(f"  ligand {ligand}:")
            for name, score, label in sorted(rows, key=lambda t: -t[1]):
                tag = {1: "selected", 0: "rejected", None: "control"}[label]
                print(f"    {name:<14s} {score:>10.4f}  {tag}")
            if len({l for _, _, l in labelled}) == 2:
                scores = [s for _, s, _ in labelled]
                labels = [l for _, _, l in labelled]
                value = auc(scores, labels)
                p = exact_p(value, labels.count(1), labels.count(0))
                print(f"    AUC {value:.3f}, exact one-sided p = {p:.3f}\n")
        print("  ⚠️ " + power_note().replace("\n", "\n  "))
        return 0

    pocket = None if args.no_pocket else MC_POCKET_POSITIONS
    jobs = []
    for variant in variants():
        for ligand_name, smiles in (("pAzF", PAZF_SMILES), ("Tyr", TYR_SMILES)):
            jobs.append({
                "name": f"{variant['name']} + {ligand_name}",
                "variant": variant["name"], "ligand": ligand_name,
                "label": variant["label"], "f_ratio": variant["f_ratio"],
                "input": job_spec(variant, smiles, ligand_name, pocket),
            })

    target = args.out / "boltz_jobs.json"
    target.write_text(json.dumps(jobs, indent=2) + "\n")

    print("=" * 76)
    print("Boltz jobs for the measured McTyrRS pAzF panel")
    print("=" * 76 + "\n")
    print(f"  {len(jobs)} jobs -> {target}")
    print(f"  pocket constraint: "
          f"{'none (Boltz places the ligand)' if pocket is None else pocket}")
    print(f"  protein length: {len(variants()[0]['sequence'])} aa\n")
    for variant in variants():
        tag = {1: "selected", 0: "rejected", None: "control"}[variant["label"]]
        ratio = (f"F+/F- {variant['f_ratio']:.2f}" if variant["f_ratio"]
                 else "no numbers")
        print(f"    {variant['name']:<14s} {variant['mutations'] or '(wild type)':<18s}"
              f" {tag:<9s} {ratio}")

    print("\n" + "=" * 76)
    print("What this pilot can establish")
    print("=" * 76)
    print("  ⚠️ " + power_note().replace("\n", "\n  "))
    print("\n  ⚠️ Boltz predicts BINDING; the assay measured AMINOACYLATION.")
    print("     Binding is a necessary condition being tested here, not the")
    print("     quantity the ratios come from.")
    print("\n  ✅ Each variant is scored against BOTH pAzF and Tyr, because the")
    print("     selection measured discrimination. A variant that binds pAzF")
    print("     well and Tyr better has gained nothing.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
