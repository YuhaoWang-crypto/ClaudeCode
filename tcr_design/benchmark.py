"""Calibration against measured TCR-pMHC binding free energies.

Where the data comes from
-------------------------
The intended source was ATLAS (Borrman et al., Proteins 2017), which curates
wild-type and mutant TCR-pMHC affinities alongside structures. Its host
(atlas.ibbr.umd.edu) currently fails TLS certificate validation, so the primary
source used here is **SKEMPI 2.0**, which tags 751 of its entries as
``TCR/pMHC`` and is the standard benchmark for protein-protein interface ddG:

    Jankauskaite J, Jimenez-Garcia B, Dapkunas J, Fernandez-Recio J,
    Moal IH. "SKEMPI 2.0: an updated benchmark of changes in protein-protein
    binding energy, kinetics and thermodynamics upon mutation."
    Bioinformatics 2019;35(3):462-469. https://life.bsc.es/pid/skempi2/

``load_atlas_csv`` accepts an ATLAS export in the same shape, so the two can be
pooled once ATLAS is reachable again.

How ddG is computed
-------------------
SKEMPI records dissociation constants for the wild type and the mutant, so

    ddG = R * T * ln(Kd_mut / Kd_wt)

in kcal/mol, with the usual convention that **positive means weaker binding**.
Note that the triage score is built so that *higher means tighter*, so a score
that works produces a **negative** correlation with ddG.

What is and is not honest to claim here
---------------------------------------
Fitting six weights on a few hundred mutations drawn from a dozen structures is
easy to overfit, and mutations from the same structure are not independent. So
every reported number is **leave-one-structure-out**: the model never sees any
mutation from the structure it is scored on. In-sample numbers are computed too,
and the gap between them is reported, because that gap is the overfitting.
"""

from __future__ import annotations

import csv
import math
import re
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .scan import WEIGHTS, Mutation, score_mutation
from .structure import Complex, load_complex

SKEMPI_URL = "https://life.bsc.es/pid/skempi2/database/download/skempi_v2.csv"
GAS_CONSTANT = 1.987204258e-3  # kcal / (mol K)
DEFAULT_TEMPERATURE = 298.0

# Fraction of a chain's benchmark mutations whose wild-type residue must match
# the structure before the chain's numbering is trusted at all.
NUMBERING_AGREEMENT_FLOOR = 0.6

# "AE99S" -> wild-type A, chain E, residue 99, mutant S. The residue number may
# carry an insertion code.
MUTATION_PATTERN = re.compile(r"^([A-Z])([A-Za-z0-9])(-?\d+)([A-Za-z]?)([A-Z])$")

# Order of the score components in the design matrix.
COMPONENTS = (
    "packing",
    "apolar_burial",
    "electrostatic",
    "hbond",
    "conformational",
    "substitution_prior",
)


@dataclass
class BenchmarkEntry:
    """One measured mutation with somewhere to hang the computed features."""

    pdb_id: str
    chain: str
    number: int
    insertion: str
    wild_type: str
    mutant: str
    ddg: float
    side: str  # "tcr" or "pmhc"
    location: str  # SKEMPI interface class: COR / SUP / RIM / INT / SUR
    source: str
    partner_chains: str = ""
    tcr_chains: str = ""

    scored: Mutation | None = None
    residue: object | None = None
    features: tuple[float, ...] = ()
    triage_score: float | None = None
    note: str = ""

    @property
    def name(self) -> str:
        return f"{self.chain}:{self.wild_type}{self.number}{self.insertion}{self.mutant}"

    @property
    def key(self) -> tuple:
        return (self.pdb_id, self.chain, self.number, self.insertion, self.mutant)

    def as_row(self) -> dict:
        row = {
            "pdb_id": self.pdb_id,
            "mutation": self.name,
            "chain": self.chain,
            "number": self.number,
            "wild_type": self.wild_type,
            "mutant": self.mutant,
            "ddg_measured": round(self.ddg, 3),
            "side": self.side,
            "interface_location": self.location,
            "source": self.source,
        }
        if self.triage_score is not None:
            row["triage_score"] = round(self.triage_score, 4)
        for name, value in zip(COMPONENTS, self.features):
            row[name] = round(value, 4)
        if self.scored is not None:
            row["contacts_peptide"] = self.scored.contacts_peptide
            row["contacts_mhc"] = self.scored.contacts_mhc
            row["buried_fraction"] = round(self.scored.buried_fraction, 3)
        if self.note:
            row["note"] = self.note
        return row


def download_skempi(cache_dir: Path) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / "skempi_v2.csv"
    if path.exists() and path.stat().st_size > 500_000:
        return path
    with urllib.request.urlopen(SKEMPI_URL, timeout=300) as response:
        path.write_bytes(response.read())
    return path


def _parse_temperature(value: str) -> float:
    cleaned = re.sub(r"\(.*?\)", "", value or "").strip()
    try:
        return float(cleaned)
    except ValueError:
        return DEFAULT_TEMPERATURE


def _ddg_from_affinities(mutant_kd: str, wild_kd: str, temperature: float) -> float | None:
    try:
        km, kw = float(mutant_kd), float(wild_kd)
    except (TypeError, ValueError):
        return None
    if km <= 0 or kw <= 0:
        return None
    return GAS_CONSTANT * temperature * math.log(km / kw)


def load_skempi(
    path: Path,
    max_mutations: int = 1,
    side: str | None = "tcr",
) -> list[BenchmarkEntry]:
    """Load the TCR/pMHC subset of SKEMPI 2.0 as benchmark entries.

    Only rows with both a wild-type and a mutant affinity are usable. Rows with
    more mutations than ``max_mutations`` are skipped, since the scoring
    function is additive over single substitutions and multi-mutant data would
    test a different claim.

    SKEMPI's ``#Pdb`` field looks like ``1AO7_ABC_DE``: the PDB id, then the
    chains of each binding partner. That second group is the TCR side, which is
    what tells us which mutations are on the receptor.
    """
    entries: list[BenchmarkEntry] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for row in csv.DictReader(handle, delimiter=";"):
            if row.get("Hold_out_type") != "TCR/pMHC":
                continue
            mutations = [m for m in row["Mutation(s)_PDB"].split(",") if m]
            if len(mutations) != max_mutations:
                continue
            ddg = _ddg_from_affinities(
                row["Affinity_mut_parsed"],
                row["Affinity_wt_parsed"],
                _parse_temperature(row.get("Temperature", "")),
            )
            if ddg is None:
                continue

            parts = row["#Pdb"].split("_")
            if len(parts) != 3:
                continue
            pdb_id, partner_chains, tcr_chains = parts

            match = MUTATION_PATTERN.match(mutations[0])
            if not match:
                continue
            wild_type, chain, number, insertion, mutant = match.groups()
            entry_side = (
                "tcr" if chain in tcr_chains
                else "pmhc" if chain in partner_chains
                else "unknown"
            )
            if side and entry_side != side:
                continue
            entries.append(
                BenchmarkEntry(
                    pdb_id=pdb_id.upper(),
                    chain=chain,
                    number=int(number),
                    insertion=insertion,
                    wild_type=wild_type,
                    mutant=mutant,
                    ddg=ddg,
                    side=entry_side,
                    location=row.get("iMutation_Location(s)", ""),
                    source="SKEMPI2",
                    partner_chains=partner_chains,
                    tcr_chains=tcr_chains,
                )
            )
    return _average_duplicates(entries)


# ATLAS serves a valid certificate but omits its InCommon intermediate, so a
# default trust store cannot build the chain. The intermediate is published at
# the address named in the leaf certificate's Authority Information Access
# extension; fetching it and appending it to the bundle fixes verification
# properly, with no weakening of TLS.
ATLAS_BASE = "https://atlas.ibbr.umd.edu/web/"
ATLAS_FILES = ("Mutants.csv", "TCRs.csv", "MHCs.csv")
ATLAS_INTERMEDIATE_URL = "https://repository.emsign.com/certs/EEEMIncommonOVG2C.crt"
ATLAS_INTERMEDIATE_SUBJECT = b"InCommon"

ATLAS_MUTATION = re.compile(r"^([A-Z])(-?\d+)([A-Za-z]?)([A-Z])$")
_MISSING = {"N/A", "n.d.", "", "WT"}


def build_atlas_ca_bundle(cache_dir: Path, system_bundle: Path | None = None) -> Path:
    """Build a CA bundle that can verify ATLAS, by adding the missing intermediate.

    This does not disable or weaken certificate validation: it supplies the
    intermediate the server should have sent, so the chain verifies normally.
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    combined = cache_dir / "atlas_ca_bundle.pem"
    if combined.exists() and combined.stat().st_size > 1000:
        return combined

    if system_bundle is None:
        for candidate in (
            Path("/root/.ccr/ca-bundle.crt"),
            Path("/etc/ssl/certs/ca-certificates.crt"),
            Path("/etc/pki/tls/certs/ca-bundle.crt"),
        ):
            if candidate.exists():
                system_bundle = candidate
                break
    if system_bundle is None:
        raise RuntimeError("no system CA bundle found to extend")

    with urllib.request.urlopen(ATLAS_INTERMEDIATE_URL, timeout=60) as response:
        raw = response.read()
    if ATLAS_INTERMEDIATE_SUBJECT not in raw:
        # DER-encoded certificates still carry the subject as readable bytes.
        raise RuntimeError(
            f"{ATLAS_INTERMEDIATE_URL} did not return the expected InCommon "
            "intermediate; refusing to add an unrecognised certificate"
        )
    if b"-----BEGIN CERTIFICATE-----" in raw:
        pem = raw.decode("ascii", errors="strict")
    else:
        import base64
        import textwrap

        body = "\n".join(textwrap.wrap(base64.b64encode(raw).decode("ascii"), 64))
        pem = f"-----BEGIN CERTIFICATE-----\n{body}\n-----END CERTIFICATE-----\n"

    combined.write_text(system_bundle.read_text() + "\n" + pem)
    return combined


def download_atlas(cache_dir: Path) -> Path:
    """Download the ATLAS mutant table, repairing the certificate chain first."""
    import ssl

    cache_dir.mkdir(parents=True, exist_ok=True)
    target = cache_dir / "atlas_Mutants.csv"
    if target.exists() and target.stat().st_size > 100_000:
        return target

    context = ssl.create_default_context(cafile=str(build_atlas_ca_bundle(cache_dir)))
    for name in ATLAS_FILES:
        destination = cache_dir / f"atlas_{name}"
        if destination.exists() and destination.stat().st_size > 500:
            continue
        with urllib.request.urlopen(
            ATLAS_BASE + name, timeout=300, context=context
        ) as response:
            destination.write_bytes(response.read())
    return target


def _parse_measured_number(value: str) -> float | None:
    """Pull a number out of an ATLAS cell such as ``0.19 +/- 0.11`` or ``1.2``."""
    if value is None:
        return None
    cleaned = value.replace("±", "+/-").split("+/-")[0].strip()
    if cleaned in _MISSING:
        return None
    try:
        return float(cleaned)
    except ValueError:
        return None


def load_atlas(path: Path, max_mutations: int = 1) -> list[BenchmarkEntry]:
    """Load TCR-side mutants from the ATLAS ``Mutants.csv`` export.

    Only rows whose MHC and peptide are wild type are kept, so the measured ddG
    is attributable to the TCR substitution alone. ``template_PDB`` names the
    wild-type structure the mutant was built from, which is the structure the
    score needs; ``true_PDB`` is used when the mutant itself was solved.
    """
    entries: list[BenchmarkEntry] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for row in csv.DictReader(handle):
            tcr_mutation = (row.get("TCR_mut") or "").strip()
            if tcr_mutation in _MISSING:
                continue
            if (row.get("MHC_mut") or "").strip() not in ("WT", "N/A"):
                continue
            if (row.get("PEP_mut") or "").strip() not in ("WT", "N/A"):
                continue

            # ATLAS writes multi-point mutants as "D26Y | L98W", and older rows
            # use a dot. Split on both so the count is right either way.
            parts = [p.strip() for p in re.split(r"[|.]", tcr_mutation) if p.strip()]
            if len(parts) != max_mutations:
                continue
            ddg = _parse_measured_number(row.get("Delta_DeltaG_kcal_per_mol", ""))
            if ddg is None:
                continue

            template = (row.get("template_PDB") or "").strip()
            true_pdb = (row.get("true_PDB") or "").strip()
            pdb_id = template if template not in _MISSING else true_pdb
            if pdb_id in _MISSING or len(pdb_id) != 4:
                continue

            # The chain column can list one chain per mutation; for a single
            # substitution the first token is the one that matters.
            chain_field = (row.get("TCR_PDB_chain") or "").strip()
            chain = chain_field.split("|")[0].strip()
            if len(chain) != 1:
                continue

            match = ATLAS_MUTATION.match(parts[0])
            if not match:
                continue
            wild_type, number, insertion, mutant = match.groups()
            entries.append(
                BenchmarkEntry(
                    pdb_id=pdb_id.upper(),
                    chain=chain,
                    number=int(number),
                    insertion=insertion,
                    wild_type=wild_type,
                    mutant=mutant,
                    ddg=ddg,
                    side="tcr",
                    location=(row.get("CDR") or "").strip(),
                    source="ATLAS",
                )
            )
    return _average_duplicates(entries)


def load_combined(
    skempi_path: Path | None = None,
    atlas_path: Path | None = None,
) -> tuple[list[BenchmarkEntry], dict]:
    """Pool SKEMPI and ATLAS, preferring ATLAS where a mutation appears in both.

    The two overlap heavily (both cover the A6 and LC13 alanine scans), so
    pooling without de-duplication would double-count and flatter any
    correlation computed afterwards.
    """
    atlas = load_atlas(atlas_path) if atlas_path and atlas_path.exists() else []
    skempi = load_skempi(skempi_path) if skempi_path and skempi_path.exists() else []

    merged: dict[tuple, BenchmarkEntry] = {}
    for entry in atlas:
        merged[entry.key] = entry
    overlap = 0
    for entry in skempi:
        if entry.key in merged:
            overlap += 1
            continue
        merged[entry.key] = entry
    entries = list(merged.values())
    provenance = {
        "atlas_entries": len(atlas),
        "skempi_entries": len(skempi),
        "overlapping_mutations": overlap,
        "pooled_unique": len(entries),
    }
    return entries, provenance


def _average_duplicates(entries: list[BenchmarkEntry]) -> list[BenchmarkEntry]:
    """Collapse repeat measurements of the same mutation to their mean ddG."""
    grouped: dict[tuple, list[BenchmarkEntry]] = {}
    for entry in entries:
        grouped.setdefault(entry.key, []).append(entry)
    merged: list[BenchmarkEntry] = []
    for group in grouped.values():
        first = group[0]
        if len(group) > 1:
            first.ddg = sum(e.ddg for e in group) / len(group)
            first.note = f"mean of {len(group)} measurements"
        merged.append(first)
    return merged


def score_entries(
    entries: list[BenchmarkEntry],
    cache_dir: Path,
    verbose: bool = True,
) -> tuple[list[BenchmarkEntry], list[str]]:
    """Attach score components to each entry, loading each structure once."""
    by_pdb: dict[str, list[BenchmarkEntry]] = {}
    for entry in entries:
        by_pdb.setdefault(entry.pdb_id, []).append(entry)

    scored: list[BenchmarkEntry] = []
    problems: list[str] = []
    for pdb_id in sorted(by_pdb):
        try:
            complex_obj = load_complex(pdb_id, cache_dir)
        except Exception as error:  # noqa: BLE001 - report and continue
            problems.append(f"{pdb_id}: could not load ({error})")
            continue

        accepted: dict[str, list[BenchmarkEntry]] = {}
        attempts: dict[str, int] = {}
        for entry in by_pdb[pdb_id]:
            attempts[entry.chain] = attempts.get(entry.chain, 0) + 1
            residue = complex_obj.residue_by_label(entry.chain, entry.number)
            if residue is None:
                problems.append(f"{pdb_id} {entry.name}: residue absent from structure")
                continue
            if residue.aa != entry.wild_type:
                problems.append(
                    f"{pdb_id} {entry.name}: structure has "
                    f"{residue.aa}{entry.number}, benchmark says {entry.wild_type}"
                )
                continue
            if entry.chain not in complex_obj.tcr_chain_ids():
                problems.append(f"{pdb_id} {entry.name}: chain is not a TCR chain here")
                continue
            entry.residue = residue
            accepted.setdefault(entry.chain, []).append(entry)

        # A chain whose residue identities mostly disagree with the benchmark is
        # numbered on a different scheme. The handful that did line up are then
        # more likely to be coincidences than real matches, so the whole chain is
        # dropped rather than trusted.
        resolved = 0
        for chain, chain_entries in accepted.items():
            tried = attempts[chain]
            rate = len(chain_entries) / tried if tried else 0.0
            if tried >= 3 and rate < NUMBERING_AGREEMENT_FLOOR:
                problems.append(
                    f"{pdb_id} chain {chain}: only {len(chain_entries)}/{tried} "
                    "mutations match the structure's residue identities; the "
                    "numbering schemes disagree, so the chain was discarded"
                )
                continue
            for entry in chain_entries:
                mutation = score_mutation(entry.residue, entry.mutant)
                entry.scored = mutation
                entry.features = tuple(getattr(mutation, name) for name in COMPONENTS)
                entry.triage_score = mutation.triage_score
                scored.append(entry)
                resolved += 1
        if verbose:
            print(f"      {pdb_id}: {resolved}/{len(by_pdb[pdb_id])} mutations scored")
    return scored, problems


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def spearman(a, b) -> float | None:
    if len(a) < 3:
        return None
    return pearson(_rankdata(a), _rankdata(b))


def pearson(a, b) -> float | None:
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if len(a) < 3 or a.std() == 0 or b.std() == 0:
        return None
    return float(np.corrcoef(a, b)[0, 1])


def _rankdata(values) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = values.argsort()
    ranks = np.empty(len(values), dtype=float)
    ranks[order] = np.arange(1, len(values) + 1, dtype=float)
    # Average ties.
    unique, inverse, counts = np.unique(values, return_inverse=True, return_counts=True)
    for index in np.where(counts > 1)[0]:
        mask = inverse == index
        ranks[mask] = ranks[mask].mean()
    return ranks


def roc_auc(scores, labels) -> float | None:
    """AUC for ranking positives above negatives, via the rank-sum identity."""
    scores = np.asarray(scores, dtype=float)
    labels = np.asarray(labels, dtype=bool)
    positives, negatives = labels.sum(), (~labels).sum()
    if positives == 0 or negatives == 0:
        return None
    ranks = _rankdata(scores)
    return float((ranks[labels].sum() - positives * (positives + 1) / 2) / (positives * negatives))


# ---------------------------------------------------------------------------
# Evaluating the existing, unfitted score
# ---------------------------------------------------------------------------

def evaluate_unfitted(entries: list[BenchmarkEntry]) -> dict:
    """Correlate the hand-weighted triage score with measured ddG.

    The score is meant to rise as binding tightens and ddG falls as binding
    tightens, so a useful score gives a negative correlation.
    """
    scores = [e.triage_score for e in entries]
    measured = [e.ddg for e in entries]

    per_structure = {}
    by_pdb: dict[str, list[BenchmarkEntry]] = {}
    for entry in entries:
        by_pdb.setdefault(entry.pdb_id, []).append(entry)
    for pdb_id, group in sorted(by_pdb.items()):
        if len(group) < 5:
            continue
        per_structure[pdb_id] = {
            "n": len(group),
            "spearman": spearman([e.triage_score for e in group], [e.ddg for e in group]),
            "pearson": pearson([e.triage_score for e in group], [e.ddg for e in group]),
        }

    stabilising = [e.ddg < -0.2 for e in entries]
    destabilising = [e.ddg > 0.5 for e in entries]
    return {
        "n": len(entries),
        "structures": len(by_pdb),
        "spearman": spearman(scores, measured),
        "pearson": pearson(scores, measured),
        "auc_find_stabilising": roc_auc(scores, stabilising),
        "auc_avoid_destabilising": roc_auc([-s for s in scores], destabilising),
        "n_stabilising": sum(stabilising),
        "n_destabilising": sum(destabilising),
        "per_structure": per_structure,
        "weights": dict(WEIGHTS),
        "convention": (
            "A working score correlates NEGATIVELY with ddG. AUC values above "
            "0.5 mean the score ranks the named class correctly."
        ),
    }


# ---------------------------------------------------------------------------
# Fitting
# ---------------------------------------------------------------------------

@dataclass
class CalibratedModel:
    """A linear model mapping score components to ddG in kcal/mol."""

    coefficients: dict[str, float]
    intercept: float
    ridge_lambda: float
    n_train: int
    held_out: dict = field(default_factory=dict)
    in_sample: dict = field(default_factory=dict)
    feature_set: str = "components"

    def predict(self, mutation: Mutation) -> float:
        """Predict ddG for a mutation scored with the component feature set."""
        if self.feature_set != "components":
            raise ValueError(
                "this model was fitted on the extended structural features; "
                "use predict_entry, which has access to the residue"
            )
        total = self.intercept
        for name, coefficient in self.coefficients.items():
            total += coefficient * getattr(mutation, name)
        return total

    def predict_entry(self, entry: "BenchmarkEntry") -> float:
        if self.feature_set == "components":
            return self.predict(entry.scored)
        vector = extended_feature_vector(entry)
        return self.intercept + sum(
            c * v for c, v in zip(self.coefficients.values(), vector)
        )

    def as_dict(self) -> dict:
        return {
            "feature_set": self.feature_set,
            "coefficients": {k: round(v, 5) for k, v in self.coefficients.items()},
            "intercept": round(self.intercept, 5),
            "ridge_lambda": self.ridge_lambda,
            "n_train": self.n_train,
            "held_out": self.held_out,
            "in_sample": self.in_sample,
            "units": "kcal/mol, positive means weaker binding",
        }


# The six hand-weighted score terms, and a richer set built from the raw
# structural measurements. The second set exists because the thing that most
# obviously drives ddG in this data - how much side chain was buried against
# the partner, and how much of it the mutation removes - is only indirectly
# visible to the six terms.
EXTENDED_FEATURES = (
    "buried_fraction",
    "delta_sasa",
    "relative_sasa_free",
    "contacts_peptide",
    "contacts_mhc",
    "void_slots",
    "growth_reaches_partner",
    "partner_charge",
    "partner_polar_atoms",
    "delta_volume",
    "delta_hydropathy",
    "delta_charge",
    "delta_hbond_capacity",
    "blosum62",
    "mutant_is_proline",
    "mutant_is_glycine",
    "wild_type_is_proline",
    "wild_type_is_glycine",
    "volume_times_burial",
    "hydropathy_times_burial",
)


def extended_feature_vector(entry: BenchmarkEntry) -> tuple[float, ...]:
    """Raw structural descriptors for one mutation.

    Everything here is measured from the structure or read from a reference
    table; none of it is a tuned quantity.
    """
    from .constants import CHARGE, HBOND_CAPACITY, HYDROPATHY, VOLUME

    residue = entry.residue
    mutation = entry.scored
    wt, mut = entry.wild_type, entry.mutant

    burial = residue.buried_fraction
    delta_volume = (VOLUME[mut] - VOLUME[wt]) / 100.0
    delta_hydropathy = (HYDROPATHY[mut] - HYDROPATHY[wt]) / 9.0
    try:
        blosum = float(BLOSUM62_MATRIX[wt, mut]) / 4.0
    except (KeyError, IndexError):
        blosum = 0.0

    return (
        burial,
        (residue.sasa_free - residue.sasa_complex) / 100.0,
        residue.relative_sasa_free,
        mutation.contacts_peptide / 10.0,
        mutation.contacts_mhc / 10.0,
        float(residue.void_slots),
        1.0 if residue.growth_reaches_partner else 0.0,
        residue.partner_charge,
        residue.partner_polar_atoms / 10.0,
        delta_volume,
        delta_hydropathy,
        CHARGE[mut] - CHARGE[wt],
        float(HBOND_CAPACITY[mut] - HBOND_CAPACITY[wt]),
        blosum,
        1.0 if mut == "P" else 0.0,
        1.0 if mut == "G" else 0.0,
        1.0 if wt == "P" else 0.0,
        1.0 if wt == "G" else 0.0,
        delta_volume * burial,
        delta_hydropathy * burial,
    )


from Bio.Align import substitution_matrices as _substitution_matrices  # noqa: E402

BLOSUM62_MATRIX = _substitution_matrices.load("BLOSUM62")


def _design_matrix(
    entries: list[BenchmarkEntry], feature_set: str = "components"
) -> tuple[np.ndarray, np.ndarray]:
    if feature_set == "extended":
        X = np.array([extended_feature_vector(e) for e in entries], dtype=float)
    else:
        X = np.array([e.features for e in entries], dtype=float)
    y = np.array([e.ddg for e in entries], dtype=float)
    return X, y


def feature_names(feature_set: str = "components") -> tuple[str, ...]:
    return EXTENDED_FEATURES if feature_set == "extended" else COMPONENTS


def _fit_ridge(X: np.ndarray, y: np.ndarray, ridge_lambda: float) -> tuple[np.ndarray, float]:
    """Ridge regression with an unpenalised intercept."""
    mean_x, mean_y = X.mean(axis=0), y.mean()
    Xc, yc = X - mean_x, y - mean_y
    n_features = Xc.shape[1]
    coefficients = np.linalg.solve(
        Xc.T @ Xc + ridge_lambda * np.eye(n_features), Xc.T @ yc
    )
    intercept = float(mean_y - mean_x @ coefficients)
    return coefficients, intercept


def fit_calibration(
    entries: list[BenchmarkEntry],
    ridge_lambda: float = 1.0,
    feature_set: str = "components",
) -> CalibratedModel:
    """Fit ddG against the score components, validated leave-one-structure-out.

    Mutations from one structure share a backbone, an epitope and usually an
    experimental method, so they are not independent. Holding out whole
    structures is the only split that answers the question anyone cares about:
    does this transfer to a receptor it has never seen?
    """
    X, y = _design_matrix(entries, feature_set)
    groups = np.array([e.pdb_id for e in entries])

    predictions = np.full(len(entries), np.nan)
    for pdb_id in np.unique(groups):
        test = groups == pdb_id
        train = ~test
        if train.sum() < 10 or test.sum() == 0:
            continue
        coefficients, intercept = _fit_ridge(X[train], y[train], ridge_lambda)
        predictions[test] = X[test] @ coefficients + intercept

    usable = ~np.isnan(predictions)
    held_out = {
        "n": int(usable.sum()),
        "structures_held_out": int(len(np.unique(groups[usable]))),
        "spearman": spearman(predictions[usable], y[usable]),
        "pearson": pearson(predictions[usable], y[usable]),
        "rmse": float(np.sqrt(np.mean((predictions[usable] - y[usable]) ** 2))),
        "mae": float(np.mean(np.abs(predictions[usable] - y[usable]))),
        "auc_find_stabilising": roc_auc(-predictions[usable], y[usable] < -0.2),
    }
    # The null model: predict the training mean for everything.
    baseline_rmse = float(np.sqrt(np.mean((y[usable] - y[usable].mean()) ** 2)))
    held_out["rmse_predicting_the_mean"] = baseline_rmse
    held_out["beats_predicting_the_mean"] = bool(held_out["rmse"] < baseline_rmse)

    coefficients, intercept = _fit_ridge(X, y, ridge_lambda)
    fitted = X @ coefficients + intercept
    in_sample = {
        "n": len(y),
        "spearman": spearman(fitted, y),
        "pearson": pearson(fitted, y),
        "rmse": float(np.sqrt(np.mean((fitted - y) ** 2))),
    }
    return CalibratedModel(
        coefficients=dict(zip(feature_names(feature_set), coefficients.tolist())),
        intercept=intercept,
        ridge_lambda=ridge_lambda,
        n_train=len(y),
        held_out=held_out,
        in_sample=in_sample,
        feature_set=feature_set,
    )


def run_calibration(
    cache_dir: Path,
    use_atlas: bool = True,
    use_skempi: bool = True,
    verbose: bool = True,
) -> tuple[dict, list[BenchmarkEntry]]:
    """Download the benchmarks, score every mutation, evaluate and fit."""
    benchmark_dir = cache_dir / "benchmark"
    skempi_path = atlas_path = None
    notes: list[str] = []

    if use_skempi:
        try:
            skempi_path = download_skempi(benchmark_dir)
        except Exception as error:  # noqa: BLE001
            notes.append(f"SKEMPI unavailable: {error}")
    if use_atlas:
        try:
            atlas_path = download_atlas(benchmark_dir)
        except Exception as error:  # noqa: BLE001
            notes.append(
                f"ATLAS unavailable: {error}. Its server omits the InCommon "
                "intermediate certificate; see build_atlas_ca_bundle."
            )

    entries, provenance = load_combined(skempi_path, atlas_path)
    if not entries:
        raise RuntimeError("no benchmark data could be loaded: " + "; ".join(notes))

    scored, problems = score_entries(entries, cache_dir / "pdb", verbose=verbose)
    if len(scored) < 30:
        raise RuntimeError(
            f"only {len(scored)} mutations could be mapped onto structures; "
            "too few to calibrate"
        )

    unfitted = evaluate_unfitted(scored)

    results = {
        "provenance": provenance,
        "scored": len(scored),
        "unmapped": len(problems),
        "problems": problems[:50],
        "notes": notes,
        "unfitted": unfitted,
    }
    for feature_set in ("components", "extended"):
        best, trace = choose_ridge_lambda(scored, feature_set=feature_set)
        model = fit_calibration(scored, ridge_lambda=best, feature_set=feature_set)
        payload = model.as_dict()
        payload["ridge_search"] = trace
        results[f"fitted_{feature_set}"] = payload
    return results, scored


def fit_excluding(
    entries: list[BenchmarkEntry],
    exclude_pdb_ids: tuple[str, ...],
    ridge_lambda: float = 1.0,
    feature_set: str = "extended",
):
    """Fit a model with specific structures removed, for a leakage-free test.

    Validating the calibrated model on a receptor that is in its training data
    would measure memorisation. This refits without the named structures so the
    test is honest.
    """
    from .calibrated import DdgModel

    excluded = {p.upper() for p in exclude_pdb_ids}
    training = [e for e in entries if e.pdb_id.upper() not in excluded]
    if len(training) < 30:
        raise RuntimeError(
            f"only {len(training)} mutations left after excluding "
            f"{sorted(excluded)}; too few to fit"
        )
    model = fit_calibration(training, ridge_lambda=ridge_lambda, feature_set=feature_set)
    return DdgModel(
        coefficients=model.coefficients,
        intercept=model.intercept,
        feature_set=feature_set,
        performance={"held_out": model.held_out, "in_sample": model.in_sample},
        provenance={
            "excluded_structures": sorted(excluded),
            "training_mutations": len(training),
            "training_structures": len({e.pdb_id for e in training}),
        },
    )


def held_out_predictions(
    entries: list[BenchmarkEntry],
    ridge_lambda: float,
    feature_set: str,
) -> np.ndarray:
    """Per-entry predictions, each made by a model that never saw its structure."""
    X, y = _design_matrix(entries, feature_set)
    groups = np.array([e.pdb_id for e in entries])
    predictions = np.full(len(entries), np.nan)
    for pdb_id in np.unique(groups):
        test = groups == pdb_id
        train = ~test
        if train.sum() < 10:
            continue
        coefficients, intercept = _fit_ridge(X[train], y[train], ridge_lambda)
        predictions[test] = X[test] @ coefficients + intercept
    return predictions


def format_calibration_report(result: dict) -> str:
    """Render the calibration run as markdown."""
    provenance = result["provenance"]
    unfitted = result["unfitted"]
    lines = [
        "# Calibrating the triage score against measured ddG",
        "",
        "## Data",
        "",
        "| source | entries |",
        "|---|---|",
        f"| ATLAS (TCR-side single mutants with ddG) | {provenance['atlas_entries']} |",
        f"| SKEMPI 2.0 (`TCR/pMHC` subset, single mutants) | {provenance['skempi_entries']} |",
        f"| overlapping between the two | {provenance['overlapping_mutations']} |",
        f"| **pooled unique mutations** | **{provenance['pooled_unique']}** |",
        f"| successfully mapped onto structures and scored | {result['scored']} |",
        f"| structures covered | {unfitted['structures']} |",
        "",
        f"Of the scored mutations, {unfitted['n_stabilising']} tighten binding "
        f"(ddG < -0.2 kcal/mol) and {unfitted['n_destabilising']} weaken it "
        "(ddG > 0.5). ddG is positive for weaker binding throughout.",
        "",
        "## 1. The hand-weighted score, as shipped",
        "",
        "| metric | value |",
        "|---|---|",
        f"| Spearman vs measured ddG | {_fmt(unfitted['spearman'])} |",
        f"| Pearson vs measured ddG | {_fmt(unfitted['pearson'])} |",
        f"| AUC, ranking stabilising mutations first | {_fmt(unfitted['auc_find_stabilising'])} |",
        f"| AUC, ranking destabilising mutations last | {_fmt(unfitted['auc_avoid_destabilising'])} |",
        "",
        "The score is built so that higher means tighter, so a working score "
        "gives a **negative** Spearman and AUCs **above 0.5**.",
        "",
        "Per structure, where at least five mutations were measured:",
        "",
        "| structure | n | Spearman |",
        "|---|---|---|",
    ]
    for pdb_id, stats in unfitted["per_structure"].items():
        lines.append(f"| {pdb_id} | {stats['n']} | {_fmt(stats['spearman'])} |")

    lines += ["", "## 2. Refitting the six weights", ""]
    lines += _model_section(result["fitted_components"])
    lines += ["", "## 3. Fitting the raw structural measurements", ""]
    lines += _model_section(result["fitted_extended"])
    return "\n".join(lines) + "\n"


def _model_section(model: dict) -> list[str]:
    held = model["held_out"]
    lines = [
        "| metric | held out | in sample |",
        "|---|---|---|",
        f"| n | {held['n']} | {model['in_sample']['n']} |",
        f"| Spearman | {_fmt(held['spearman'])} | {_fmt(model['in_sample']['spearman'])} |",
        f"| Pearson | {_fmt(held['pearson'])} | {_fmt(model['in_sample']['pearson'])} |",
        f"| RMSE (kcal/mol) | {_fmt(held['rmse'])} | {_fmt(model['in_sample']['rmse'])} |",
        f"| AUC, finding stabilising mutations | {_fmt(held['auc_find_stabilising'])} | - |",
        "",
        "The AUC row is the one that matters for design. Most of this benchmark "
        "is alanine scanning, so a model can score well overall simply by "
        "learning that deleting a large buried side chain weakens binding. "
        "Finding the rarer mutations that *improve* binding is the actual task, "
        "and it is harder.",
        "",
        f"Predicting the training mean for everything gives RMSE "
        f"{_fmt(held['rmse_predicting_the_mean'])}. This model "
        + ("**beats**" if held["beats_predicting_the_mean"] else "**does not beat**")
        + " that baseline.",
        "",
        "Fitted coefficients (kcal/mol per unit of each feature):",
        "",
        "| feature | coefficient |",
        "|---|---|",
    ]
    ordered = sorted(
        model["coefficients"].items(), key=lambda kv: -abs(kv[1])
    )
    for name, value in ordered:
        lines.append(f"| {name} | {value:+.3f} |")
    lines.append(f"| _intercept_ | {model['intercept']:+.3f} |")
    return lines


def _fmt(value) -> str:
    if value is None:
        return "n/a"
    return f"{value:+.3f}" if isinstance(value, float) else str(value)


def choose_ridge_lambda(
    entries: list[BenchmarkEntry],
    candidates: tuple[float, ...] = (0.01, 0.1, 1.0, 10.0, 100.0),
    feature_set: str = "components",
) -> tuple[float, list[dict]]:
    """Pick the ridge strength with the best held-out RMSE.

    Selecting the hyperparameter on the same held-out split that is then
    reported is mildly optimistic. The candidate grid is small and the reported
    spread across it is narrow, so the effect is small, but it is a real caveat
    rather than a clean nested cross-validation.
    """
    trace = []
    best, best_rmse = candidates[0], math.inf
    for value in candidates:
        model = fit_calibration(entries, ridge_lambda=value, feature_set=feature_set)
        rmse = model.held_out.get("rmse", math.inf)
        trace.append(
            {
                "ridge_lambda": value,
                "held_out_rmse": rmse,
                "held_out_spearman": model.held_out.get("spearman"),
            }
        )
        if rmse < best_rmse:
            best, best_rmse = value, rmse
    return best, trace
