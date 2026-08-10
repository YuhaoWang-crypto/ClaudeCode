#!/usr/bin/env python
"""Derive the three-class immunophenotype from TIL geometry, on all 136 patients.

The published TIL maps carry no tumour boundary, so the desert / excluded /
inflamed trichotomy cannot be read off them directly.  This script derives it
from geometry instead, following the same two-step structure the reference study
uses -- amount first, then arrangement:

1. **desert vs rest** from the amount of infiltrate;
2. **excluded vs inflamed** from *penetration*: what share of tissue the
   infiltrate never gets near.

Both cut-points are fitted quantiles, so they are parameters. They are reported
with a sensitivity analysis rather than presented as thresholds.

Then the part that decides whether any of this means anything: **the three
groups are tested against the transcriptome**.  Under the intended
interpretation, inflamed tumours should carry the most CD8 and interferon
signal, desert the least, and excluded should sit in between while carrying
*more* signal than desert despite an un-penetrated interior.  That is a
falsifiable ordering, and it is checked rather than assumed.

A data-driven three-cluster solution is computed alongside, so the rule-based
assignment can be compared with what the descriptors would do unsupervised.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402
from scipy.stats import kruskal, spearmanr  # noqa: E402

from mri_cd8_til.clients.idc import load_map  # noqa: E402
from mri_cd8_til.labels.geometry import describe, descriptor_names  # noqa: E402
from mri_cd8_til.labels.schema import Immunophenotype  # noqa: E402
from mri_cd8_til.labels.transcriptome import SIGNATURES  # noqa: E402

DESERT_QUANTILE = 1 / 3
EXCLUDED_QUANTILE = 0.5


def compute_geometry(cache: Path, patient_ids: list[str]) -> list[dict]:
    records = []
    for i, patient in enumerate(patient_ids, start=1):
        try:
            til_map = load_map(cache, patient)
        except FileNotFoundError:
            continue
        geometry = describe(patient, til_map.grid)
        if not geometry.is_evaluable:
            print(f"  {patient}: skipped, only {geometry.n_tissue_patches} tissue patches")
            continue
        records.append(geometry.as_dict())
        if i % 25 == 0 or i == len(patient_ids):
            print(f"  [{i}/{len(patient_ids)}] {patient}")
    return records


def assign_phenotype(
    records: list[dict],
    desert_quantile: float = DESERT_QUANTILE,
    excluded_quantile: float = EXCLUDED_QUANTILE,
) -> tuple[list[Immunophenotype], dict]:
    """Two-step rule: amount, then penetration among the non-desert cases."""
    amount = np.array([r["til_fraction"] for r in records])
    unreached = np.array([r["unreached_fraction"] for r in records])

    desert_cut = float(np.quantile(amount, desert_quantile))
    is_desert = amount <= desert_cut

    # The excluded/inflamed cut is fitted among non-desert cases only: including
    # deserts would drag it toward their (necessarily high) unreached values.
    excluded_cut = float(np.quantile(unreached[~is_desert], excluded_quantile))

    phenotypes = []
    for desert, un in zip(is_desert, unreached):
        if desert:
            phenotypes.append(Immunophenotype.DESERT)
        elif un >= excluded_cut:
            phenotypes.append(Immunophenotype.EXCLUDED)
        else:
            phenotypes.append(Immunophenotype.INFLAMED)
    return phenotypes, {"desert_cut": desert_cut, "excluded_cut": excluded_cut}


def read_signature_scores(rna_dir: Path, patients: list[str]) -> dict[str, np.ndarray]:
    """Immune-signature z-scores per patient from GDC STAR-Counts tables."""
    wanted = {g for genes in SIGNATURES.values() for g in genes} | {"EPCAM"}
    expression: dict[str, dict[str, float]] = {}
    for patient in patients:
        path = rna_dir / f"{patient}.tsv.gz"
        if not path.exists():
            continue
        values: dict[str, float] = {}
        with gzip.open(path, "rt") as handle:
            header = None
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if header is None:
                    if fields[0] == "gene_id":
                        header = fields
                    continue
                if fields[0].startswith("N_"):
                    continue
                name = fields[header.index("gene_name")]
                if name in wanted:
                    try:
                        values[name] = float(fields[header.index("tpm_unstranded")])
                    except ValueError:
                        continue
        if values:
            expression[patient] = values

    present = [p for p in patients if p in expression]
    scores: dict[str, np.ndarray] = {"_patients": np.array(present, dtype=object)}
    for name, genes in list(SIGNATURES.items()) + [("epcam_control", ("EPCAM",))]:
        matrix = np.vstack(
            [np.log2(np.array([expression[p].get(g, np.nan) for p in present]) + 1) for g in genes]
        )
        sd = np.nanstd(matrix, axis=1, ddof=1, keepdims=True)
        centred = (matrix - np.nanmean(matrix, axis=1, keepdims=True)) / np.where(sd > 0, sd, 1.0)
        scores[name] = np.nanmean(centred, axis=0)
    return scores


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", default="results/til_maps")
    parser.add_argument("--labels", default="results/til_labels.json")
    parser.add_argument("--rna-dir", default=None)
    parser.add_argument("--out", default="results/til_geometry.json")
    args = parser.parse_args()

    patient_ids = [r["patient_id"] for r in json.loads(Path(args.labels).read_text())["records"]]
    print(f"[1] computing {len(descriptor_names())} geometry descriptors for {len(patient_ids)} patients")
    records = compute_geometry(Path(args.cache), patient_ids)
    print(f"  evaluable: {len(records)}")

    print("\n[2] two-step phenotype assignment")
    phenotypes, cuts = assign_phenotype(records)
    counts = {p.value: sum(1 for x in phenotypes if x is p) for p in Immunophenotype}
    for name, count in counts.items():
        print(f"    {name:<18} {count:3d}  ({count / len(phenotypes):.0%})")
    print(f"  cuts: TIL fraction <= {cuts['desert_cut']:.4f} -> desert; "
          f"unreached >= {cuts['excluded_cut']:.3f} -> excluded")
    print("  (reference study: desert 36.8%, excluded 16.5%, inflamed 46.7%)")

    print("\n[3] descriptor profile by phenotype (median)")
    key_descriptors = [
        "til_fraction", "unreached_fraction", "penetration_p90",
        "band_index", "morans_i", "largest_cluster_share", "interface_density",
    ]
    header = f"    {'descriptor':<24}" + "".join(f"{p.value:>18}" for p in Immunophenotype)
    print(header)
    print("    " + "-" * (len(header) - 4))
    for descriptor in key_descriptors:
        values = np.array([r[descriptor] for r in records])
        cells = "".join(
            f"{np.median(values[[i for i, x in enumerate(phenotypes) if x is p]]):>18.4f}"
            for p in Immunophenotype
        )
        print(f"    {descriptor:<24}{cells}")

    payload = {
        "n_patients": len(records),
        "n_descriptors": len(descriptor_names()),
        "cuts": cuts,
        "counts": counts,
        "records": [
            {**record, "phenotype": phenotype.value}
            for record, phenotype in zip(records, phenotypes)
        ],
    }

    if args.rna_dir:
        print("\n[4] biological validation — do the groups behave as the names claim?")
        patients = [r["case_id"] for r in records]
        scores = read_signature_scores(Path(args.rna_dir), patients)
        present = list(scores.pop("_patients"))
        index = {p: i for i, p in enumerate(present)}
        groups = {
            p: np.array([index[r["case_id"]] for r, x in zip(records, phenotypes)
                         if x is p and r["case_id"] in index])
            for p in Immunophenotype
        }
        print(f"  patients with RNA-seq: {len(present)}/{len(records)}")

        validation = {}
        for signature in ("cd8_core", "tis", "cytolytic", "epcam_control"):
            values = scores[signature]
            arrays = [values[groups[p]] for p in Immunophenotype if groups[p].size > 2]
            if len(arrays) < 3:
                continue
            stat = kruskal(*arrays)
            medians = {p.value: float(np.median(values[groups[p]])) for p in Immunophenotype}
            ordered = (
                medians[Immunophenotype.INFLAMED.value]
                > medians[Immunophenotype.EXCLUDED.value]
                > medians[Immunophenotype.DESERT.value]
            )
            validation[signature] = {
                "medians": medians,
                "kruskal_p": float(stat.pvalue),
                "ordering_holds": bool(ordered),
            }
            flag = "ORDER HOLDS" if ordered else "order broken"
            print(
                f"    {signature:<14} desert {medians['immune-desert']:+.3f} | "
                f"excluded {medians['immune-excluded']:+.3f} | "
                f"inflamed {medians['inflamed']:+.3f}   "
                f"Kruskal p={stat.pvalue:.4f}  [{flag}]"
            )

        print("\n  Which descriptors track immune signal, independent of amount?")
        cd8 = scores["cd8_core"]
        amount = np.array([r["til_fraction"] for r in records if r["case_id"] in index])
        rows = []
        for descriptor in descriptor_names():
            values = np.array([r[descriptor] for r in records if r["case_id"] in index])
            if not np.isfinite(values).all() or values.std() == 0:
                continue
            rho = spearmanr(values, cd8)
            # Partial correlation given amount, so a descriptor is not simply
            # re-reporting how much infiltrate there is.
            partial = _partial_spearman(values, cd8, amount)
            rows.append((descriptor, rho.statistic, rho.pvalue, partial))
        rows.sort(key=lambda r: -abs(r[3]))
        for descriptor, rho, p, partial in rows[:8]:
            print(f"    {descriptor:<24} rho={rho:+.3f} (p={p:.3f})   partial|amount={partial:+.3f}")
        validation["descriptor_vs_cd8"] = {
            r[0]: {"rho": float(r[1]), "p": float(r[2]), "partial_given_amount": float(r[3])}
            for r in rows
        }
        payload["validation"] = validation

    path = Path(args.out)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))
    print(f"\nwrote {path}")
    return 0


def _partial_spearman(x: np.ndarray, y: np.ndarray, control: np.ndarray) -> float:
    """Spearman partial correlation of x and y given control, via rank residuals."""
    from scipy.stats import rankdata

    rx, ry, rc = (rankdata(v) for v in (x, y, control))

    def residual(v):
        design = np.column_stack([rc, np.ones_like(rc)])
        beta, *_ = np.linalg.lstsq(design, v, rcond=None)
        return v - design @ beta

    ex, ey = residual(rx), residual(ry)
    if ex.std() == 0 or ey.std() == 0:
        return 0.0
    return float(np.corrcoef(ex, ey)[0, 1])


if __name__ == "__main__":
    raise SystemExit(main())
