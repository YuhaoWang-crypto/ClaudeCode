"""Build LNP_ML / LiON in-silico-screen input files from a SMILES list.

The LiON ``predict`` command expects a library folder with three aligned CSVs:

* ``<name>.csv``            : a single ``smiles`` column (the ionizable lipids).
* ``<name>_extra_x.csv``    : 36 numeric formulation/context features (X_val).
* ``<name>_metadata.csv``   : bookkeeping columns (not used by the model).

This module constructs those three files for a chosen **formulation** and
**biological context** so we can score our own enumerated libraries. It runs on
CPU with no chemprop dependency, so it is testable in this session; the actual
prediction happens on GPU/Modal (see ``modal_app/lion.py``).

The default context is the repo-recommended screening baseline: the **KK
formulation** (35:16:46.5:2.5 ionizable:DOPE:cholesterol:PEG molar) in HeLa,
mRNA cargo, in vitro — which the LNP_ML authors note gives the most stable
predictions.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors

RDLogger.DisableLog("rdApp.*")

# aza-Michael addition: an amine N-H adds across the acrylate/acrylamide C=C,
# giving a tertiary amine + an ester/amide-linked tail. Generic in the linker:
# the [*:6] captures the ester -O-R or amide -N(H)-R of the Michael acceptor.
# nucleophile is an AMINE N-H (exclude amide/carbonyl-adjacent N, else the
# acrylamide's own N-H would react and runaway-polymerise).
_MICHAEL = AllChem.ReactionFromSmarts(
    "[N;!H0;!$(N-C=O);!$(N=*):1].[CH2:2]=[CH1:3][C:4](=[O:5])[*:6]>>"
    "[N:1][CH2:2][CH2:3][C:4](=[O:5])[*:6]")


def _exhaustively_alkylate(head_smiles: str, acceptor_smiles: str) -> str | None:
    """Add the Michael acceptor to every N-H of the head (NH2 -> 2 tails, NH -> 1)."""
    cur = Chem.MolFromSmiles(head_smiles)
    acc = Chem.MolFromSmiles(acceptor_smiles)
    if cur is None or acc is None:
        return None
    for _ in range(12):  # safety cap; each pass consumes one N-H
        prods = _MICHAEL.RunReactants((cur, acc))
        if not prods:
            break
        nxt = None
        for (p,) in prods:
            try:
                Chem.SanitizeMol(p)
                nxt = p
                break
            except Exception:
                continue
        if nxt is None:
            break
        cur = nxt
    return Chem.MolToSmiles(cur)


def _acceptor(linker: str, tail_smiles: str) -> str:
    """Build the Michael acceptor reagent for a linker + tail.
    ester -> acrylate  CH2=CH-C(=O)-O-R ; amide -> acrylamide CH2=CH-C(=O)-NH-R."""
    if linker == "ester":
        return f"C=CC(=O)O{tail_smiles}"
    if linker == "amide":
        return f"C=CC(=O)N{tail_smiles}"
    raise ValueError(f"unknown linker {linker!r}")


_AMINE_NH = Chem.MolFromSmarts("[N;!H0;!$(N-C=O);!$(N=*)]")


def _n_reactive_nh(smiles: str) -> int:
    """Count reactive AMINE N-H (excludes amide/imine N-H)."""
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return 0
    return sum(m.GetAtomWithIdx(i[0]).GetTotalNumHs() for i in m.GetSubstructMatches(_AMINE_NH))


def enumerate_michael_lipids(
    heads: dict[str, str],
    tails: dict[str, str],
    linkers: tuple[str, ...] = ("ester", "amide"),
    mw_range: tuple[float, float] = (400.0, 1400.0),
) -> pd.DataFrame:
    """Combinatorially enumerate ionizable lipids by aza-Michael addition.

    ``heads`` / ``tails`` map id -> SMILES (tail = an alkyl/hetero group R; it is
    attached as the ester -O-R or amide -N(H)-R of an acrylate). Every head N-H is
    alkylated, so the tail count = the head's reactive N-H count. Products outside
    ``mw_range`` or that fail to build are dropped.

    Returns a DataFrame: head, linker, tail, n_tails, mw, smiles.
    """
    rows, seen = [], set()
    for hid, hsmi in heads.items():
        nh = _n_reactive_nh(hsmi)
        for linker in linkers:
            for tid, tsmi in tails.items():
                smi = _exhaustively_alkylate(hsmi, _acceptor(linker, tsmi))
                if not smi or smi in seen:
                    continue
                m = Chem.MolFromSmiles(smi)
                if m is None:
                    continue
                mw = Descriptors.MolWt(m)
                if not (mw_range[0] <= mw <= mw_range[1]):
                    continue
                # require a fully-substituted amine (no remaining N-H) and >=2 tails
                if _n_reactive_nh(smi) > 0 or nh < 2:
                    continue
                seen.add(smi)
                rows.append({"head": hid, "linker": linker, "tail": tid,
                             "n_tails": nh, "mw": round(mw, 1), "smiles": smi})
    return pd.DataFrame(rows)

# --- exact X_val schema from data/libraries/test_screen/*_extra_x.csv ---------

CONTINUOUS = [
    "Cationic_Lipid_to_mRNA_weight_ratio",
    "Cationic_Lipid_Mol_Ratio",
    "Phospholipid_Mol_Ratio",
    "Cholesterol_Mol_Ratio",
    "PEG_Lipid_Mol_Ratio",
]

# one-hot groups: {feature_prefix: [allowed suffixes]}
ONEHOT_GROUPS: dict[str, list[str]] = {
    "Delivery_target": ["dendritic_cell", "generic_cell", "liver", "lung",
                        "lung_epithelium", "macrophage", "muscle", "spleen"],
    "Helper_lipid_ID": ["DOPE", "DOTAP", "DSPC", "MDOA", "None"],
    "Route_of_administration": ["in_vitro", "intramuscular", "intratracheal",
                               "intravenous"],
    "Batch_or_individual_or_barcoded": ["Barcoded", "Individual"],
    "Cargo_type": ["mRNA", "pDNA", "siRNA"],
    "Model_type": ["A549", "BDMC", "BMDM", "HBEC_ALI", "HEK293T", "HeLa",
                  "IGROV1", "Mouse", "RAW264p7"],
}

# full ordered X column list, exactly matching the reference file
EXTRA_X_COLUMNS: list[str] = CONTINUOUS + [
    f"{prefix}_{suffix}"
    for prefix, suffixes in ONEHOT_GROUPS.items()
    for suffix in suffixes
]


@dataclass
class Formulation:
    """A formulation + biological context for a LiON screen.

    Molar ratios should sum to ~100. Categorical fields must be one of the
    allowed values in ``ONEHOT_GROUPS``.
    """

    # continuous
    cationic_to_mrna_wt: float = 10.0
    cationic_mol: float = 35.0
    phospholipid_mol: float = 16.0
    cholesterol_mol: float = 46.5
    peg_mol: float = 2.5
    # categorical (defaults = KK baseline in HeLa)
    delivery_target: str = "generic_cell"
    helper_lipid: str = "DOPE"
    route: str = "in_vitro"
    batch: str = "Individual"
    cargo: str = "mRNA"
    model_type: str = "HeLa"

    _CHOICE_FIELDS = {
        "Delivery_target": "delivery_target",
        "Helper_lipid_ID": "helper_lipid",
        "Route_of_administration": "route",
        "Batch_or_individual_or_barcoded": "batch",
        "Cargo_type": "cargo",
        "Model_type": "model_type",
    }

    def validate(self) -> None:
        for prefix, attr in self._CHOICE_FIELDS.items():
            val = getattr(self, attr)
            if val not in ONEHOT_GROUPS[prefix]:
                raise ValueError(
                    f"{attr}={val!r} not in allowed {ONEHOT_GROUPS[prefix]}")
        total = self.cationic_mol + self.phospholipid_mol + self.cholesterol_mol + self.peg_mol
        if abs(total - 100) > 1.0:
            raise ValueError(f"molar ratios sum to {total:.1f}, expected ~100")

    def extra_x_row(self) -> dict[str, float]:
        self.validate()
        row = {
            "Cationic_Lipid_to_mRNA_weight_ratio": self.cationic_to_mrna_wt,
            "Cationic_Lipid_Mol_Ratio": self.cationic_mol,
            "Phospholipid_Mol_Ratio": self.phospholipid_mol,
            "Cholesterol_Mol_Ratio": self.cholesterol_mol,
            "PEG_Lipid_Mol_Ratio": self.peg_mol,
        }
        for prefix, attr in self._CHOICE_FIELDS.items():
            chosen = getattr(self, attr)
            for suffix in ONEHOT_GROUPS[prefix]:
                row[f"{prefix}_{suffix}"] = int(suffix == chosen)
        return row


def build_library(
    smiles: list[str],
    out_dir: str | Path,
    name: str = "my_screen",
    formulation: Formulation | None = None,
) -> Path:
    """Write the three LiON screen CSVs for ``smiles`` into ``out_dir/name/``.

    Returns the created directory. Point ``main_script.py predict`` at it (after
    copying under ``data/libraries/``).
    """
    formulation = formulation or Formulation()
    row = formulation.extra_x_row()

    out_dir = Path(out_dir) / name
    out_dir.mkdir(parents=True, exist_ok=True)

    pd.DataFrame({"smiles": smiles}).to_csv(out_dir / f"{name}.csv", index=False)

    extra = pd.DataFrame([row] * len(smiles))[EXTRA_X_COLUMNS]
    extra.to_csv(out_dir / f"{name}_extra_x.csv", index=False)

    # minimal metadata (model ignores it, but predict expects the file)
    meta = pd.DataFrame({
        "Lipid_name": [f"cand_{i}" for i in range(len(smiles))],
        "Helper_lipid_ID": formulation.helper_lipid,
        "Cargo_type": formulation.cargo,
        "Model_type": formulation.model_type,
        "Delivery_target": formulation.delivery_target,
        "Route_of_administration": formulation.route,
    })
    meta.to_csv(out_dir / f"{name}_metadata.csv", index=False)
    return out_dir


# named baseline contexts, handy for organ-targeted screens
BASELINES: dict[str, Formulation] = {
    "KK_HeLa": Formulation(),  # repo default
    "lung_epithelium_IT": Formulation(delivery_target="lung_epithelium",
                                      route="intratracheal", model_type="Mouse"),
    "liver_IV": Formulation(delivery_target="liver", route="intravenous",
                            model_type="Mouse"),
    "spleen_IV": Formulation(delivery_target="spleen", route="intravenous",
                             model_type="Mouse"),
}
