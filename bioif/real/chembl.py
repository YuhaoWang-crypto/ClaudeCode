"""
bioif.real.chembl -- a REAL data source behind the affinity adapter.

This replaces the `affinity-model@0.1-stub` source in the demo with measured
bioactivity from the EMBL-EBI ChEMBL REST API. Nothing here is a stub.

Two things matter for the interface argument:

  * ChEMBL keeps assay type, units, standardised activity and the source
    document *on every activity record*. That is not bureaucracy -- it is
    the minimum needed to know whether two numbers may be compared. The
    contract in bioif.core asks for the same fields for the same reason.

  * A target name does not resolve to a target. `resolve_target("KRAS")`
    returns several ChEMBL targets of different `target_type`. Picking the
    first hit is the identity bug that no amount of model accuracy fixes, so
    `resolve_target` refuses to guess.

Network use is cached. A compact snapshot of what the demo needs is
committed under `_snapshot/`, so `python3 -m bioif.selftest` is deterministic
and runs offline; set BIOIF_REFRESH=1 to re-fetch from the live API.
"""
from __future__ import annotations

import csv
import json
import os
import subprocess
import urllib.parse
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

BASE = "https://www.ebi.ac.uk/chembl/api/data"
SNAPSHOT = Path(__file__).parent / "_snapshot"
CACHE = Path(__file__).parent / "_cache"


class OfflineError(RuntimeError):
    """Raised when live data is needed but neither cache nor network has it."""


def _fetch(path: str, **params) -> dict:
    """GET from ChEMBL, with an on-disk cache keyed by the full URL."""
    url = f"{BASE}/{path}?" + urllib.parse.urlencode(params)
    CACHE.mkdir(exist_ok=True)
    key = CACHE / (urllib.parse.quote(url, safe="")[-180:] + ".json")
    if key.exists() and not os.environ.get("BIOIF_REFRESH"):
        return json.loads(key.read_text())
    if os.environ.get("BIOIF_OFFLINE"):
        raise OfflineError(f"BIOIF_OFFLINE set and not cached: {url}")
    proc = subprocess.run(["curl", "-sS", "--max-time", "90", url],
                          capture_output=True, text=True)
    if proc.returncode != 0 or not proc.stdout:
        raise OfflineError(f"ChEMBL unreachable: {url} ({proc.stderr[:200]})")
    data = json.loads(proc.stdout)
    key.write_text(json.dumps(data))
    return data


# --------------------------------------------------------------------------
# Identity resolution -- the part that must be allowed to refuse
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class TargetHit:
    chembl_id: str
    pref_name: str
    target_type: str
    organism: str


@dataclass(frozen=True)
class Resolution:
    """Either a single resolved target, or an explicit refusal."""
    resolved: TargetHit | None
    candidates: tuple[TargetHit, ...]
    reason: str = ""

    @property
    def ok(self) -> bool:
        return self.resolved is not None


def resolve_target(query: str, organism: str = "Homo sapiens",
                   target_type: str = "SINGLE PROTEIN") -> Resolution:
    """
    Resolve a free-text target name against ChEMBL.

    Refuses unless exactly one candidate survives the organism and
    target_type constraints. A gene symbol is not an identifier: "KRAS"
    matches single proteins, protein complexes, protein-protein interactions
    and a protein family, and silently taking the first row is how a
    compound measured against a PPI ends up scored as a direct inhibitor.
    """
    try:
        rows = _fetch("target/search.json", q=query, limit=25).get("targets", [])
    except OfflineError:
        # fall back to the committed snapshot for the query it covers, so the
        # test suite is deterministic and needs no network
        if query.upper() != "KRAS" or not SNAP_TARGETS.exists():
            raise
        rows = json.loads(SNAP_TARGETS.read_text())
    hits = tuple(TargetHit(t["target_chembl_id"], t.get("pref_name") or "",
                           t.get("target_type") or "", t.get("organism") or "")
                 for t in rows)
    if not hits:
        return Resolution(None, (), f"no ChEMBL target matches {query!r}")
    keep = tuple(h for h in hits
                 if (not organism or h.organism == organism)
                 and (not target_type or h.target_type == target_type))
    if len(keep) == 1:
        return Resolution(keep[0], hits)
    if not keep:
        return Resolution(None, hits,
                          f"{len(hits)} candidates for {query!r}, none is a "
                          f"{organism} {target_type}".rstrip())
    kinds = sorted({h.target_type for h in keep})
    return Resolution(None, hits,
                      f"{query!r} is ambiguous: {len(keep)} candidates "
                      f"spanning {len(kinds)} target types ({', '.join(kinds)})")


# --------------------------------------------------------------------------
# Measured activities
# --------------------------------------------------------------------------

#: The fields that decide whether two activity records may be compared. This
#: list IS the interface contract, expressed in ChEMBL's vocabulary.
ACTIVITY_FIELDS = ("molecule_chembl_id", "assay_chembl_id", "target_chembl_id",
                   "standard_type", "standard_relation", "standard_units",
                   "pchembl_value", "assay_type", "assay_description",
                   "document_chembl_id", "target_organism")


def fetch_activities(target_chembl_id: str, molecule_chembl_id: str | None = None,
                     max_records: int = 5000) -> list[dict]:
    """All pChEMBL-bearing activities for a target (optionally one compound)."""
    out: list[dict] = []
    for offset in range(0, max_records, 1000):
        params = dict(target_chembl_id=target_chembl_id,
                      pchembl_value__isnull="false", limit=1000, offset=offset)
        if molecule_chembl_id:
            params["molecule_chembl_id"] = molecule_chembl_id
        page = _fetch("activity.json", **params)
        rows = page.get("activities", [])
        out += [{k: a.get(k) for k in ACTIVITY_FIELDS} for a in rows]
        if len(rows) < 1000:
            break
    return out


# --------------------------------------------------------------------------
# Committed snapshot -- keeps the tests deterministic and offline-capable
# --------------------------------------------------------------------------

SNAP_ACTIVITIES = SNAPSHOT / "kras_activities.csv"
SNAP_TARGETS = SNAPSHOT / "kras_target_search.json"
SNAP_META = SNAPSHOT / "PROVENANCE.json"

DEMO_TARGET = "CHEMBL2189121"          # GTPase KRas, Homo sapiens, SINGLE PROTEIN
DEMO_MOLECULE = "CHEMBL4072295"        # research compound, 16 records / 16 assays


def build_snapshot() -> dict:
    """Re-fetch from the live API and rewrite the committed snapshot."""
    SNAPSHOT.mkdir(exist_ok=True)
    acts = fetch_activities(DEMO_TARGET)
    with SNAP_ACTIVITIES.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=ACTIVITY_FIELDS, extrasaction="ignore")
        w.writeheader()
        for a in acts:
            a = dict(a)
            if a.get("assay_description"):
                a["assay_description"] = a["assay_description"][:160]
            w.writerow(a)
    # Keep only the fields that decide identity; the full records carry
    # megabytes of component cross-references the contract does not use.
    tgts = _fetch("target/search.json", q="KRAS", limit=25)["targets"]
    SNAP_TARGETS.write_text(json.dumps(
        [{k: t.get(k) for k in ("target_chembl_id", "pref_name",
                                "target_type", "organism")} for t in tgts],
        indent=1))
    meta = {
        "source": "EMBL-EBI ChEMBL REST API",
        "base_url": BASE,
        "target": DEMO_TARGET,
        "filter": "pchembl_value__isnull=false",
        "n_activities": len(acts),
        "fetched_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "note": ("Measured bioactivity, not model output. Fields retained are "
                 "exactly those needed to decide comparability."),
    }
    SNAP_META.write_text(json.dumps(meta, indent=1))
    return meta


def load_snapshot() -> list[dict]:
    if not SNAP_ACTIVITIES.exists():
        raise OfflineError("no snapshot; run python3 -m bioif.real.chembl")
    with SNAP_ACTIVITIES.open() as fh:
        return list(csv.DictReader(fh))


def snapshot_provenance() -> dict:
    return json.loads(SNAP_META.read_text()) if SNAP_META.exists() else {}


if __name__ == "__main__":
    print(json.dumps(build_snapshot(), indent=1))
