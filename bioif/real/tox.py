"""
bioif.real.tox -- real data for three more pairings of the interface map.

Covers the blueprint's C->F (compound -> in-vitro cellular readout), C->G
(compound -> mutagenicity) and, most usefully, the F->G edge that the
blueprint's own Demo B marks `blocked`: does a p53 reporter result tell you
anything about mutagenicity?

Two public datasets, both with measured labels:

  Tox21 (MoleculeNet release) -- 7,831 compounds x 12 binary assay endpoints.
    SR-p53    : p53 response-element reporter (the C->F endpoint of interest)
    SR-MMP    : mitochondrial membrane potential -- the general-cytotoxicity
                confound the blueprint warns about, kept as a control
    SR-ATAD5  : DNA-damage response reporter
  Hansen Ames benchmark (N=6512) -- bacterial reverse mutation, the C->G
    endpoint.

The two sets are joined on the InChIKey **skeleton** (first block), which
ignores salt form, charge and stereochemistry. That is a deliberate,
declared choice: joining on full InChIKey loses most of the overlap to salt
and stereo differences, and joining on name or CAS loses more. It is exactly
the identity decision the blueprint's bottleneck #1 is about, so it is
recorded in provenance rather than buried.

Network use is cached and a compact snapshot is committed, so everything
runs offline; BIOIF_REFRESH=1 re-fetches.
"""
from __future__ import annotations

import csv
import gzip
import json
import os
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

SNAPSHOT = Path(__file__).parent / "_snapshot"
CACHE = Path(__file__).parent / "_cache"

TOX21_URL = ("https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/"
             "tox21.csv.gz")
AMES_URL = ("https://doc.ml.tu-berlin.de/toxbenchmark/"
            "Mutagenicity_N6512.csv")

SNAP_TOX = SNAPSHOT / "tox21_keyed.csv"
SNAP_AMES = SNAPSHOT / "ames_keyed.csv"
SNAP_META = SNAPSHOT / "TOX_PROVENANCE.json"

TOX21_ENDPOINTS = ("NR-AR", "NR-AR-LBD", "NR-AhR", "NR-Aromatase", "NR-ER",
                   "NR-ER-LBD", "NR-PPAR-gamma", "SR-ARE", "SR-ATAD5",
                   "SR-HSE", "SR-MMP", "SR-p53")

#: The endpoint this work is about, and the two controls that decide whether
#: a positive means "genotoxic stress" or merely "the cells were dying".
P53 = "SR-p53"
CYTOTOX_CONTROL = "SR-MMP"
DDR_CONTROL = "SR-ATAD5"


class OfflineError(RuntimeError):
    pass


def _download(url: str, dest: Path) -> Path:
    CACHE.mkdir(exist_ok=True)
    if dest.exists() and not os.environ.get("BIOIF_REFRESH"):
        return dest
    if os.environ.get("BIOIF_OFFLINE"):
        raise OfflineError(f"BIOIF_OFFLINE set and not cached: {url}")
    p = subprocess.run(["curl", "-sSL", "--max-time", "180", "-o", str(dest),
                        url], capture_output=True, text=True)
    if p.returncode != 0 or not dest.exists() or dest.stat().st_size == 0:
        raise OfflineError(f"could not fetch {url}: {p.stderr[:200]}")
    return dest


# --------------------------------------------------------------------------
# Identity: the join key, and the fact that it is a choice
# --------------------------------------------------------------------------

def skeleton_key(smiles: str) -> str | None:
    """
    InChIKey skeleton (first block) of the salt-stripped, neutralised parent.

    Returns None when RDKit cannot parse or kekulise the input -- a refusal,
    not a guess. The counts of such refusals are reported in provenance.
    """
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem.MolStandardize import rdMolStandardize
        RDLogger.DisableLog("rdApp.*")
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return None
        try:
            m = rdMolStandardize.FragmentParent(m)
            m = rdMolStandardize.Uncharger().uncharge(m)
        except Exception:
            pass
        ik = Chem.MolToInchiKey(m)
        return ik.split("-")[0] if ik else None
    except Exception:
        return None


# --------------------------------------------------------------------------
# Snapshot
# --------------------------------------------------------------------------

def build_snapshot() -> dict:
    SNAPSHOT.mkdir(exist_ok=True)
    tox_gz = _download(TOX21_URL, CACHE / "tox21.csv.gz")
    ames_raw = _download(AMES_URL, CACHE / "ames.csv")

    n_tox_bad = n_ames_bad = 0
    with gzip.open(tox_gz, "rt") as fh, SNAP_TOX.open("w", newline="") as out:
        w = csv.writer(out)
        w.writerow(["skeleton", "mol_id", "smiles", *TOX21_ENDPOINTS])
        seen = set()
        for row in csv.DictReader(fh):
            k = skeleton_key(row["smiles"])
            if not k:
                n_tox_bad += 1
                continue
            if k in seen:
                continue
            seen.add(k)
            w.writerow([k, row["mol_id"], row["smiles"],
                        *[row[e] for e in TOX21_ENDPOINTS]])
    n_tox = len(seen)

    with open(ames_raw) as fh, SNAP_AMES.open("w", newline="") as out:
        w = csv.writer(out)
        w.writerow(["skeleton", "cas", "smiles", "ames"])
        seen_a = set()
        for row in csv.DictReader(fh):
            k = skeleton_key(row["Canonical_Smiles"])
            if not k:
                n_ames_bad += 1
                continue
            if k in seen_a:
                continue
            seen_a.add(k)
            w.writerow([k, row.get("CAS_NO", ""), row["Canonical_Smiles"],
                        row["Activity"]])
    n_ames = len(seen_a)

    meta = {
        "tox21": {"source": "MoleculeNet Tox21 release", "url": TOX21_URL,
                  "unique_skeletons": n_tox, "unparsed": n_tox_bad,
                  "endpoints": list(TOX21_ENDPOINTS)},
        "ames": {"source": "Hansen mutagenicity benchmark N=6512",
                 "url": AMES_URL, "unique_skeletons": n_ames,
                 "unparsed": n_ames_bad},
        "join_key": ("InChIKey skeleton (first block) of the salt-stripped, "
                     "neutralised parent; ignores salt form, charge and "
                     "stereochemistry"),
        "fetched_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "note": ("Measured assay labels, not model output. Tox21 labels are "
                 "qHTS actives/inactives; Ames is the bacterial reverse "
                 "mutation assay."),
    }
    SNAP_META.write_text(json.dumps(meta, indent=1))
    return meta


def load_tox21() -> list[dict]:
    if not SNAP_TOX.exists():
        raise OfflineError("no Tox21 snapshot; run python3 -m bioif.real.tox")
    with SNAP_TOX.open() as fh:
        return list(csv.DictReader(fh))


def load_ames() -> list[dict]:
    if not SNAP_AMES.exists():
        raise OfflineError("no Ames snapshot; run python3 -m bioif.real.tox")
    with SNAP_AMES.open() as fh:
        return list(csv.DictReader(fh))


def provenance() -> dict:
    return json.loads(SNAP_META.read_text()) if SNAP_META.exists() else {}


# --------------------------------------------------------------------------
# Layer 1 of the blueprint's evidence ladder: measured, no model at all
# --------------------------------------------------------------------------

@dataclass
class TwoByTwo:
    endpoint: str
    n: int
    a: int            # endpoint+ , Ames+
    b: int            # endpoint+ , Ames-
    c: int            # endpoint- , Ames+
    d: int            # endpoint- , Ames-

    @property
    def ppv(self) -> float:
        return self.a / max(self.a + self.b, 1)

    @property
    def base(self) -> float:
        return self.c / max(self.c + self.d, 1)

    @property
    def risk_ratio(self) -> float:
        return self.ppv / self.base if self.base else float("inf")

    @property
    def odds_ratio(self) -> float:
        return (self.a * self.d) / (self.b * self.c) if self.b and self.c \
            else float("inf")

    @property
    def sensitivity(self) -> float:
        return self.a / max(self.a + self.c, 1)

    def fisher_p(self) -> float:
        """Two-sided Fisher exact test, computed exactly (no scipy needed)."""
        from math import comb
        a, b, c, d = self.a, self.b, self.c, self.d
        r1, r2 = a + b, c + d
        c1, n = a + c, a + b + c + d
        def prob(x):
            return (comb(r1, x) * comb(r2, c1 - x)) / comb(n, c1)
        obs = prob(a)
        lo = max(0, c1 - r2)
        hi = min(r1, c1)
        return min(1.0, sum(prob(x) for x in range(lo, hi + 1)
                            if prob(x) <= obs * (1 + 1e-9)))


def overlap(tox=None, ames=None) -> dict[str, dict]:
    """Compounds measured in both datasets, keyed by skeleton."""
    tox = {r["skeleton"]: r for r in (tox or load_tox21())}
    ames = {r["skeleton"]: r for r in (ames or load_ames())}
    return {k: {"tox": tox[k], "ames": ames[k]} for k in set(tox) & set(ames)}


def two_by_two(endpoint: str = P53, ov: dict | None = None) -> TwoByTwo:
    ov = ov if ov is not None else overlap()
    a = b = c = d = 0
    for rec in ov.values():
        v = rec["tox"][endpoint]
        if v in ("", "NA"):
            continue
        pos = int(float(v)) == 1
        mut = int(float(rec["ames"]["ames"])) == 1
        if pos and mut:
            a += 1
        elif pos:
            b += 1
        elif mut:
            c += 1
        else:
            d += 1
    return TwoByTwo(endpoint, a + b + c + d, a, b, c, d)


def report() -> str:
    ov = overlap()
    p = provenance()
    L = ["Measured co-occurrence: Tox21 reporter endpoints vs Ames",
         f"  tox21 : {p.get('tox21', {}).get('unique_skeletons', '?')} unique "
         f"skeletons   ames: {p.get('ames', {}).get('unique_skeletons', '?')}",
         f"  joined on {p.get('join_key', '?')}",
         f"  OVERLAP: {len(ov)} compounds measured in both",
         ""]
    L.append(f"  {'endpoint':<12}{'n':>6}{'ep+ %Ames+':>12}{'ep- %Ames+':>12}"
             f"{'RR':>7}{'OR':>7}{'sens':>7}{'Fisher p':>11}")
    for e in (P53, DDR_CONTROL, CYTOTOX_CONTROL):
        t = two_by_two(e, ov)
        L.append(f"  {e:<12}{t.n:>6}{t.ppv * 100:>11.1f}%"
                 f"{t.base * 100:>11.1f}%{t.risk_ratio:>7.2f}"
                 f"{t.odds_ratio:>7.2f}{t.sensitivity:>7.2f}"
                 f"{t.fisher_p():>11.2e}")
    t = two_by_two(P53, ov)
    L += ["",
          "Reading (this is Layer 1: measured, no model involved):",
          f"  A p53-reporter positive raises the probability of being Ames-positive",
          f"  from {t.base:.1%} to {t.ppv:.1%} -- a real {t.risk_ratio:.2f}x enrichment.",
          f"  But {t.base:.1%} of p53-NEGATIVE compounds are still Ames-positive, and",
          f"  the reporter catches only {t.sensitivity:.0%} of the mutagens. So the edge is",
          "  usable for PRIORITISATION and is nowhere near a replacement for the",
          "  assay: as a screen it would miss most of what it is screening for."]
    return "\n".join(L)


if __name__ == "__main__":
    if not SNAP_TOX.exists() or os.environ.get("BIOIF_REFRESH"):
        print(json.dumps(build_snapshot(), indent=1))
    print(report())
