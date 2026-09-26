"""Stage 1 - curate one training set per panel target from ChEMBL.

Design carried over from the 3CLpro virtual-assay work, plus the fixes its own
limitations section asked for.

What is enforced here:

* **Targets are pinned by UniProt accession**, and the resolved ChEMBL record
  must come back as SINGLE PROTEIN. The 3CLpro run had to keyword-filter assay
  descriptions after the fact because CHEMBL4523582 is the whole pp1ab
  polyprotein and mixed Mpro, PLpro and RdRp activities - 19% of its labels
  belonged to a different protein. Resolving from an accession makes that class
  of contamination structurally impossible instead of something to catch later.

* **A confidence-score floor.** ChEMBL's assay confidence_score says how
  securely an assay is assigned to a molecular target; anything below 8 is a
  homologous-protein or subunit-level assignment. Those are dropped, which is
  the systematic version of the keyword filtering the earlier work did by hand.

  confidence_score lives on the **assay** record, not the activity record - the
  /activity endpoint does not carry it. A first version read it off the activity
  rows, got None for every one, and silently filtered all 80,000 activities in
  the panel to zero. The distinct assay IDs are therefore fetched from /assay in
  batches and joined. Anything that reads a quality field off the wrong endpoint
  fails this way: not with an error, but with an empty dataset that looks like a
  real answer.

* **A data-reality gate.** Fewer than MIN_COMPOUNDS usable compounds and the
  target is reported as failed rather than modelled on a set too small to
  scaffold-split.

* **Censored values are dropped.** Only standard_relation "=" is kept. A ">"
  IC50 means the compound was inactive up to the top concentration tested,
  which is a different measurement, not a weak number.

* **Replicate aggregation is by median**, and an honest noise floor is measured
  from it. The spread has to be computed across *distinct assays*, not across
  raw rows: ChEMBL deposits the same measurement more than once, so a naive
  row-level spread came out at 0.01 log units - an assay reproducibility no
  biochemical assay has. Only compounds measured in two or more different assay
  records contribute, and each assay contributes its own median first.

  Even then a median is the wrong summary, because the distribution is bimodal.
  On PARP1, half of the compounds measured in two assays agree to within 0.05
  log, and those are concentrated in the pairs that come from two *different
  documents* with a range of exactly 0.000 - one primary measurement re-reported,
  not two independent assays. The remaining third disagree by more than 0.5 log,
  with a 90th percentile of 2.0. A median over both modes reads 0.04 and
  understates real assay noise by a factor of ~20.

  So the floor is reported as the full quantile set plus a "when two assays
  actually disagree" median (non-zero ranges only), and the report quotes that
  rather than a single flattering number. This matters because the noise floor is
  the ceiling model RMSE is judged against.

* **Assay-type mixing is tracked, not hidden.** IC50/Ki/Kd are pooled into one
  pAffinity, but the functional-vs-binding group of each record is kept so the
  per-group offset can be reported.

ChEMBL is CC BY-SA 3.0: derived datasets keep the ChEMBL IDs, cite the version
recorded in results/panel_provenance.json, and must be shared alike.
"""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
API = "https://www.ebi.ac.uk/chembl/api/data"

FORCE_REFETCH = False  # raw activities are cached; set True to re-pull from ChEMBL
MIN_COMPOUNDS = 300  # below this a scaffold-split benchmark is not meaningful
MIN_CONFIDENCE = 8  # ChEMBL assay confidence: 8/9 = direct single-protein
MAX_MW = 650.0
KEEP_TYPES = ("IC50", "Ki", "Kd")
PAGE = 1000


def api(path: str, **params) -> dict:
    url = f"{API}/{path}.json?" + urllib.parse.urlencode(params)
    last = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=120) as fh:
                return json.load(fh)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            last = exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"ChEMBL request failed after retries: {path} {params}: {last}")


def resolve_target(uniprot: str) -> dict | None:
    """UniProt accession -> the single-protein ChEMBL target record."""
    res = api(
        "target",
        target_components__accession=uniprot,
        target_type="SINGLE PROTEIN",
        limit=10,
    )
    for t in res.get("targets", []):
        # Require the accession to actually be a component of this target, so a
        # fuzzy backend match cannot slip a different protein through.
        accs = {
            c.get("accession")
            for comp in t.get("target_components", [])
            for c in [comp]
        }
        if uniprot in accs and t.get("organism") == "Homo sapiens":
            return t
    # Fall back to any organism but say so.
    for t in res.get("targets", []):
        if uniprot in {c.get("accession") for c in t.get("target_components", [])}:
            return t
    return None


def fetch_activities(target_chembl_id: str) -> pd.DataFrame:
    rows, offset = [], 0
    while True:
        res = api(
            "activity",
            target_chembl_id=target_chembl_id,
            standard_type__in=",".join(KEEP_TYPES),
            standard_relation="=",
            limit=PAGE,
            offset=offset,
        )
        batch = res.get("activities", [])
        rows.extend(batch)
        meta = res.get("page_meta", {})
        offset += PAGE
        if not batch or offset >= meta.get("total_count", 0):
            break
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(
        [
            {
                "molecule_chembl_id": r.get("molecule_chembl_id"),
                "canonical_smiles": r.get("canonical_smiles"),
                "standard_type": r.get("standard_type"),
                "standard_value": r.get("standard_value"),
                "standard_units": r.get("standard_units"),
                "pchembl_value": r.get("pchembl_value"),
                "potential_duplicate": r.get("potential_duplicate"),
                "document_chembl_id": r.get("document_chembl_id"),
                "assay_chembl_id": r.get("assay_chembl_id"),
                "assay_type": r.get("assay_type"),
                "bao_label": r.get("bao_label"),
                "assay_description": r.get("assay_description"),
                "data_validity_comment": r.get("data_validity_comment"),
            }
            for r in rows
        ]
    )


def fetch_assay_confidence(assay_ids: list[str]) -> dict[str, int]:
    """assay_chembl_id -> confidence_score, batched over the /assay endpoint."""
    conf: dict[str, int] = {}
    ids = sorted(set(a for a in assay_ids if isinstance(a, str)))
    for i in range(0, len(ids), 50):
        chunk = ids[i : i + 50]
        res = api("assay", assay_chembl_id__in=",".join(chunk), limit=len(chunk))
        for a in res.get("assays", []):
            cs = a.get("confidence_score")
            if cs is not None:
                conf[a["assay_chembl_id"]] = int(cs)
    return conf


_normalizer = rdMolStandardize.Normalizer()
_chooser = rdMolStandardize.LargestFragmentChooser()
_uncharger = rdMolStandardize.Uncharger()


def standardize(smiles: str) -> tuple[str | None, str | None, float | None]:
    """Canonical parent SMILES, Murcko scaffold, MW - or (None, None, None)."""
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if mol is None:
        return None, None, None
    try:
        mol = _uncharger.uncharge(_chooser.choose(_normalizer.normalize(mol)))
        Chem.RemoveStereochemistry(mol)  # stereo is inconsistently reported
        smi = Chem.MolToSmiles(mol)
        scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=mol) or smi
        return smi, scaf, float(Descriptors.MolWt(mol))
    except Exception:
        return None, None, None


def curate_one(target: dict) -> dict:
    name, uniprot = target["name"], target["uniprot"]
    rec = resolve_target(uniprot)
    if rec is None:
        return {"name": name, "status": "no_single_protein_target", "n_compounds": 0}
    tid = rec["target_chembl_id"]
    print(
        f"\n{name} ({uniprot}) -> {tid}  "
        f"{rec.get('pref_name', '?')[:52]}  [{rec.get('target_type')}, {rec.get('organism')}]"
    )

    cache = DATA / f"raw_{name}.csv.gz"
    if cache.exists() and not FORCE_REFETCH:
        raw = pd.read_csv(cache)
        print(f"  (reusing cached raw activities: {len(raw)} rows)")
    else:
        raw = fetch_activities(tid)
    if raw.empty:
        return {"name": name, "status": "no_activities", "n_compounds": 0}
    n_raw = len(raw)

    df = raw.copy()
    audit = {"n_raw_activities": n_raw}

    # --- assay-level quality filters -----------------------------------------
    conf = fetch_assay_confidence(df["assay_chembl_id"].tolist())
    df["confidence_score"] = df["assay_chembl_id"].map(conf)
    audit["n_distinct_assays"] = int(df["assay_chembl_id"].nunique())
    audit["n_assays_without_confidence"] = int(
        df["assay_chembl_id"].nunique() - len(conf)
    )
    # A missing confidence_score is treated as failing the floor, but it is
    # counted separately so an API change cannot quietly empty the dataset again.
    low_conf = df["confidence_score"].fillna(-1) < MIN_CONFIDENCE
    audit["dropped_low_confidence"] = int(low_conf.sum())
    audit["confidence_score_distribution"] = (
        df["confidence_score"].value_counts(dropna=False).sort_index().to_dict()
    )
    df = df[~low_conf]
    if df.empty:
        return {
            "name": name, "uniprot": uniprot, "target_chembl_id": tid,
            "bps_family": target["bps_family"], "assay_readout": target["assay_readout"],
            "status": "all_filtered_at_confidence", "n_compounds": 0, "audit": audit,
        }

    bad = df["data_validity_comment"].notna() & (
        df["data_validity_comment"].astype(str).str.strip() != ""
    )
    audit["dropped_data_validity_flagged"] = int(bad.sum())
    audit["data_validity_reasons"] = (
        df.loc[bad, "data_validity_comment"].value_counts().head(6).to_dict()
    )
    df = df[~bad]

    df["pchembl_value"] = pd.to_numeric(df["pchembl_value"], errors="coerce")
    no_p = df["pchembl_value"].isna()
    audit["dropped_no_pchembl"] = int(no_p.sum())
    df = df[~no_p]

    if df.empty:
        return {
            "name": name, "uniprot": uniprot, "target_chembl_id": tid,
            "bps_family": target["bps_family"], "assay_readout": target["assay_readout"],
            "status": "all_filtered", "n_compounds": 0, "audit": audit,
        }

    # --- structure standardisation ------------------------------------------
    std = df["canonical_smiles"].map(standardize)
    df["std_smiles"] = [s[0] for s in std]
    df["scaffold"] = [s[1] for s in std]
    df["mw"] = [s[2] for s in std]
    unparsed = df["std_smiles"].isna()
    audit["dropped_unparseable"] = int(unparsed.sum())
    df = df[~unparsed]
    heavy = df["mw"] > MAX_MW
    audit["dropped_mw_over_limit"] = int(heavy.sum())
    df = df[~heavy]

    df["assay_group"] = np.where(df["assay_type"].eq("B"), "binding", "functional")

    # --- aggregate replicates, and measure an inter-assay noise floor --------
    df.to_csv(DATA / f"raw_{name}.csv.gz", index=False)  # cache: re-analysis needs no refetch

    # Per (compound, assay) median first, so repeated deposition of a single
    # measurement cannot masquerade as agreement between independent assays.
    per_assay = (
        df.groupby(["molecule_chembl_id", "assay_chembl_id"])["pchembl_value"]
        .median()
        .reset_index()
    )
    spread = per_assay.groupby("molecule_chembl_id")["pchembl_value"].agg(
        n_assays="size",
        inter_assay_range=lambda v: float(v.max() - v.min()),
    )
    multi = spread.loc[spread["n_assays"] > 1, "inter_assay_range"]
    discordant = multi[multi > 0.0]
    q = lambda v, x: round(float(v.quantile(x)), 3) if len(v) else None
    noise = {
        "definition": (
            "range of per-assay medians over compounds measured in >1 distinct "
            "ChEMBL assay record. The distribution is bimodal: exact ties are "
            "re-deposition of one measurement, not independent replication, so "
            "noise_floor_log is the median over non-zero ranges."
        ),
        "n_compounds_in_multiple_assays": int(len(multi)),
        "frac_exact_ties": round(float((multi == 0).mean()), 3) if len(multi) else None,
        "quantiles_all": {str(x): q(multi, x) for x in (0.25, 0.5, 0.75, 0.9, 0.95)},
        "n_discordant": int(len(discordant)),
        "noise_floor_log": q(discordant, 0.5),
        "noise_floor_q75_log": q(discordant, 0.75),
        "n_rows_flagged_potential_duplicate": int(
            pd.to_numeric(df["potential_duplicate"], errors="coerce").fillna(0).sum()
        ),
    }

    grp = df.groupby("molecule_chembl_id")
    agg = grp.agg(
        pAffinity=("pchembl_value", "median"),
        n_meas=("pchembl_value", "size"),
        std_smiles=("std_smiles", "first"),
        scaffold=("scaffold", "first"),
        mw=("mw", "first"),
        assay_group=("assay_group", lambda v: v.mode().iat[0]),
        types=("standard_type", lambda v: "|".join(sorted(set(v)))),
    ).reset_index()
    agg = agg.merge(spread.reset_index(), on="molecule_chembl_id", how="left")

    by_group = agg.groupby("assay_group")["pAffinity"].agg(["size", "mean"]).to_dict("index")
    group_offset = None
    if {"binding", "functional"} <= set(by_group):
        group_offset = round(
            abs(by_group["binding"]["mean"] - by_group["functional"]["mean"]), 3
        )

    agg = agg.sort_values("molecule_chembl_id").reset_index(drop=True)
    status = "ok" if len(agg) >= MIN_COMPOUNDS else "below_data_gate"
    out = {
        "name": name,
        "uniprot": uniprot,
        "bps_family": target["bps_family"],
        "assay_readout": target["assay_readout"],
        "target_chembl_id": tid,
        "target_pref_name": rec.get("pref_name"),
        "target_type": rec.get("target_type"),
        "organism": rec.get("organism"),
        "status": status,
        "n_compounds": int(len(agg)),
        "n_scaffolds": int(agg["scaffold"].nunique()),
        "pAffinity_mean": round(float(agg["pAffinity"].mean()), 3),
        "pAffinity_sd": round(float(agg["pAffinity"].std()), 3),
        "pAffinity_min": round(float(agg["pAffinity"].min()), 3),
        "pAffinity_max": round(float(agg["pAffinity"].max()), 3),
        "assay_group_counts": {k: int(v["size"]) for k, v in by_group.items()},
        "binding_vs_functional_mean_offset_log": group_offset,
        "noise_floor": noise,
        "audit": audit,
    }
    if status == "ok":
        agg.to_csv(DATA / f"train_{name}.csv", index=False)
    print(
        f"  {n_raw} raw -> {len(agg)} compounds / {out['n_scaffolds']} scaffolds  "
        f"pAff {out['pAffinity_mean']:.2f}+-{out['pAffinity_sd']:.2f}  "
        f"noise floor {noise['noise_floor_log']} log "
        f"(n={noise['n_discordant']} discordant of {noise['n_compounds_in_multiple_assays']}, "
        f"{noise['frac_exact_ties']:.0%} exact ties)  [{status}]"
    )
    return out


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    panel = json.loads((ROOT / "configs" / "panel.json").read_text())

    status = api("status")
    print(
        f"ChEMBL release {status.get('chembl_release_version', '?')} "
        f"({status.get('chembl_db_version', '?')})"
    )

    summaries = [curate_one(t) for t in panel["targets"]]

    ok = [s for s in summaries if s["status"] == "ok"]
    print(f"\n{'=' * 78}\n{len(ok)}/{len(summaries)} targets passed the data gate "
          f"(>= {MIN_COMPOUNDS} compounds)")
    print(f"{'target':<10s} {'n':>6s} {'scaf':>5s} {'pAff':>12s} {'noise':>6s}  family")
    for s in sorted(summaries, key=lambda r: -r["n_compounds"]):
        if s["status"] != "ok":
            print(f"{s['name']:<10s} {s['n_compounds']:>6d} {'-':>5s} {'-':>12s} "
                  f"{'-':>6s}  FAILED: {s['status']}")
            continue
        nf = s["noise_floor"]["noise_floor_log"]
        print(
            f"{s['name']:<10s} {s['n_compounds']:>6d} {s['n_scaffolds']:>5d} "
            f"{s['pAffinity_mean']:>6.2f}+-{s['pAffinity_sd']:<5.2f} "
            f"{nf if nf is not None else '-':>6}  {s['bps_family']}"
        )

    (RESULTS / "panel_provenance.json").write_text(
        json.dumps(
            {
                "chembl_release": status.get("chembl_release_version"),
                "chembl_db_version": status.get("chembl_db_version"),
                "licence": "ChEMBL data is CC BY-SA 3.0 - attribute and share alike",
                "min_compounds_gate": MIN_COMPOUNDS,
                "min_assay_confidence_score": MIN_CONFIDENCE,
                "max_mw": MAX_MW,
                "activity_types": list(KEEP_TYPES),
                "censored_values": "dropped (standard_relation '=' only)",
                "targets": summaries,
            },
            indent=2,
        )
    )
    print(f"\nwrote {RESULTS / 'panel_provenance.json'}")


if __name__ == "__main__":
    main()
