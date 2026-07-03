"""Modal app: reproduce the LiON (Witten et al. 2024) D-MPNN LNP-potency model.

Runs the upstream LNP_ML pipeline (Chemprop 1.7.0, Python 3.8) on a Modal GPU and
persists trained models + screen results to a Modal Volume.

This session is CPU-only, so this file is authored to be run from a machine with
the Modal CLI + a token:

    pip install modal && modal token new
    modal run modal_app/lion.py::train    --split all_random_split_for_paper --epochs 30
    modal run modal_app/lion.py::analyze  --split all_random_split_for_paper
    modal run modal_app/lion.py::demo_screen

    # score your own library (JSON list of SMILES):
    modal run modal_app/lion.py::screen --smiles-json candidates.json \
        --split all_random_split_for_paper --context KK_HeLa

Notes
-----
* Pre-made splits already ship in the repo (``all_random_split_for_paper``,
  ``all_amine_split_for_paper``), so ``train`` works with no ``split`` step.
* Full-data training is ~3 h on one GPU. Use ``--epochs 3`` on a toy split first.
* The LNP_ML pipeline is orchestrated via ``scripts/main_script.py`` run with
  ``cwd=scripts`` (it references ``../data`` relatively).
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import modal

REPO = "https://github.com/jswitten/LNP_ML.git"
REPO_DIR = "/LNP_ML"
SCRIPTS = f"{REPO_DIR}/scripts"
MODELS_MNT = "/models"

# Chemprop 1.7.0 hard-pins Python >=3.7,<3.9, but Modal's image builder only
# supports >=3.10. So the Modal *function* runs in 3.10 (its client needs that),
# while Chemprop lives in a separate 3.8 venv built with `uv` inside the image.
# main_script.py is invoked with that 3.8 interpreter via subprocess.
CHEMPROP_PY = "/opt/lnp38/bin/python"

image = (
    modal.Image.debian_slim(python_version="3.10")
    .apt_install("git", "build-essential",
                 # RDKit's Draw module (imported by main_script.py) needs these
                 "libxrender1", "libxext6", "libsm6", "libglib2.0-0")
    .pip_install("uv", "pandas")  # pandas for the function's own post-processing (3.10)
    .run_commands(
        # standalone Python 3.8 + Chemprop 1.7.0 (pulls CUDA torch, rdkit, pandas…)
        "uv venv --python 3.8 /opt/lnp38",
        # uv venvs omit setuptools; hyperopt (a chemprop dep) needs pkg_resources
        "uv pip install --python /opt/lnp38/bin/python setuptools wheel chemprop==1.7.0",
        f"git clone --depth 1 {REPO} {REPO_DIR}",
        # make LiON's hardcoded 5-fold ensemble configurable via LION_CV_NUM, so a
        # lightweight run can train/predict with fewer folds (both train & predict
        # loops read cv_num). Surgical one-liner; default stays 5.
        f"""sed -i "s/cv_num = 5/cv_num = int(__import__('os').environ.get('LION_CV_NUM','5'))/g" {SCRIPTS}/main_script.py""",
    )
)

app = modal.App("lion-lnp")
vol = modal.Volume.from_name("lion-models", create_if_missing=True)


def _run(cmd: list[str], env_extra: dict | None = None) -> None:
    # cmd[0] == "python" -> run with the 3.8 Chemprop interpreter
    if cmd and cmd[0] == "python":
        cmd = [CHEMPROP_PY] + cmd[1:]
    import os
    env = os.environ.copy()
    if env_extra:
        env.update({k: str(v) for k, v in env_extra.items()})
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=SCRIPTS, check=True, env=env)


def _restore_split_from_volume(split: str) -> None:
    """Copy a trained split folder from the Volume back into the repo tree."""
    import shutil
    src = Path(MODELS_MNT) / split
    dst = Path(REPO_DIR) / "data" / "crossval_splits" / split
    if src.exists():
        if dst.exists():
            shutil.rmtree(dst)
        shutil.copytree(src, dst)


def _save_split_to_volume(split: str) -> None:
    import shutil
    src = Path(REPO_DIR) / "data" / "crossval_splits" / split
    dst = Path(MODELS_MNT) / split
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    vol.commit()


@app.function(image=image, gpu="A10G", volumes={MODELS_MNT: vol}, timeout=60 * 60 * 6)
def train(split: str = "all_random_split_for_paper", epochs: int = 30,
          folds: int = 5) -> str:
    """Train the D-MPNN ensemble on a pre-made split (full: 5 folds, 30 epochs)."""
    _run(["python", "main_script.py", "train", split, "--epochs", str(epochs)],
         env_extra={"LION_CV_NUM": folds})
    _save_split_to_volume(split)
    return f"trained '{split}' ({folds} folds, {epochs} epochs); models on volume"


@app.function(image=image, gpu="T4", volumes={MODELS_MNT: vol}, timeout=60 * 60 * 2)
def train_lite(split: str = "all_random_split_for_paper", epochs: int = 15,
               folds: int = 1) -> str:
    """Lightweight training on the FULL dataset: a single fold, fewer epochs, on a
    cheaper T4 GPU. ~5-10x cheaper/faster than the 5-fold ensemble; yields one
    usable model (no ensemble averaging / cross-val error bars). Screen it with
    ``screen(..., folds=1)``."""
    _run(["python", "main_script.py", "train", split, "--epochs", str(epochs)],
         env_extra={"LION_CV_NUM": folds})
    _save_split_to_volume(split)
    return f"lite-trained '{split}' ({folds} fold(s), {epochs} epochs, T4)"


@app.function(image=image, gpu="A10G", volumes={MODELS_MNT: vol}, timeout=60 * 60 * 2)
def analyze(split: str = "all_random_split_for_paper") -> str:
    """Compute held-out Pearson/Spearman/Kendall correlations for a trained split."""
    _restore_split_from_volume(split)
    _run(["python", "main_script.py", "analyze", split])
    # persist the results folder too
    import shutil
    res = Path(REPO_DIR) / "results" / "crossval_splits" / split
    if res.exists():
        dst = Path(MODELS_MNT) / f"{split}__results"
        if dst.exists():
            shutil.rmtree(dst)
        shutil.copytree(res, dst)
        vol.commit()
    return f"analyzed '{split}'; correlations under {split}__results on the volume"


# Hermetic copy of the extra-X schema + a few named contexts. The canonical,
# richer version (with validation and more baselines) lives in
# ``lipidlib/lion_library.py``; keep this in sync if that schema changes.
_CONTINUOUS = ["Cationic_Lipid_to_mRNA_weight_ratio", "Cationic_Lipid_Mol_Ratio",
               "Phospholipid_Mol_Ratio", "Cholesterol_Mol_Ratio", "PEG_Lipid_Mol_Ratio"]
_ONEHOT = {
    "Delivery_target": ["dendritic_cell", "generic_cell", "liver", "lung",
                        "lung_epithelium", "macrophage", "muscle", "spleen"],
    "Helper_lipid_ID": ["DOPE", "DOTAP", "DSPC", "MDOA", "None"],
    "Route_of_administration": ["in_vitro", "intramuscular", "intratracheal", "intravenous"],
    "Batch_or_individual_or_barcoded": ["Barcoded", "Individual"],
    "Cargo_type": ["mRNA", "pDNA", "siRNA"],
    "Model_type": ["A549", "BDMC", "BMDM", "HBEC_ALI", "HEK293T", "HeLa",
                   "IGROV1", "Mouse", "RAW264p7"],
}
# context -> (delivery_target, helper, route, cargo, model, batch, mol ratios, wt)
_CONTEXTS = {
    "KK_HeLa":              ("generic_cell",   "DOPE", "in_vitro",     "mRNA", "HeLa"),
    "lung_epithelium_IT":   ("lung_epithelium","DOPE", "intratracheal","mRNA", "Mouse"),
    "liver_IV":             ("liver",          "DOPE", "intravenous",  "mRNA", "Mouse"),
    "spleen_IV":            ("spleen",         "DOPE", "intravenous",  "mRNA", "Mouse"),
}


def _write_library(smiles: list[str], name: str, context: str) -> None:
    """Build the 3 LiON screen CSVs under data/libraries/<name> (KK molar ratios)."""
    import pandas as pd
    if context not in _CONTEXTS:
        raise ValueError(f"unknown context {context!r}; choose {list(_CONTEXTS)}")
    dtgt, helper, route, cargo, model = _CONTEXTS[context]
    choices = {"Delivery_target": dtgt, "Helper_lipid_ID": helper,
               "Route_of_administration": route, "Batch_or_individual_or_barcoded": "Individual",
               "Cargo_type": cargo, "Model_type": model}
    row = {"Cationic_Lipid_to_mRNA_weight_ratio": 10.0, "Cationic_Lipid_Mol_Ratio": 35.0,
           "Phospholipid_Mol_Ratio": 16.0, "Cholesterol_Mol_Ratio": 46.5, "PEG_Lipid_Mol_Ratio": 2.5}
    cols = list(_CONTINUOUS)
    for prefix, suffixes in _ONEHOT.items():
        for suf in suffixes:
            col = f"{prefix}_{suf}"
            row[col] = int(choices[prefix] == suf)
            cols.append(col)

    lib = Path(REPO_DIR) / "data" / "libraries" / name
    lib.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"smiles": smiles}).to_csv(lib / f"{name}.csv", index=False)
    pd.DataFrame([row] * len(smiles))[cols].to_csv(lib / f"{name}_extra_x.csv", index=False)
    pd.DataFrame({"Lipid_name": [f"cand_{i}" for i in range(len(smiles))],
                  "Helper_lipid_ID": helper, "Cargo_type": cargo, "Model_type": model,
                  "Delivery_target": dtgt, "Route_of_administration": route}
                 ).to_csv(lib / f"{name}_metadata.csv", index=False)


@app.function(image=image, gpu="A10G", volumes={MODELS_MNT: vol}, timeout=60 * 60 * 2)
def screen(smiles_json: str = "", split: str = "all_random_split_for_paper",
           context: str = "KK_HeLa", name: str = "user_screen", folds: int = 5) -> str:
    """Score a user library. ``smiles_json`` is a JSON list of SMILES strings.
    Set ``folds=1`` to match a ``train_lite`` model (must equal the trained folds)."""
    _restore_split_from_volume(split)
    smiles = json.loads(smiles_json) if smiles_json else []
    if not smiles:
        return "no SMILES provided"
    _write_library(smiles, name, context)
    _run(["python", "main_script.py", "predict", split, name],
         env_extra={"LION_CV_NUM": folds})

    # LiON writes screen output to "<model_folder>_preds/<library>/"
    pred = (Path(REPO_DIR) / "results" / "screen_results" / f"{split}_preds" / name / "pred_file.csv")
    out = pred.read_text() if pred.exists() else "(no pred_file.csv produced)"
    # persist predictions
    dst = Path(MODELS_MNT) / f"screen__{name}.csv"
    if pred.exists():
        dst.write_text(out)
        vol.commit()
    return out[:4000]


@app.function(image=image, gpu="A10G", volumes={MODELS_MNT: vol}, timeout=60 * 60)
def demo_screen(split: str = "small_test_split_with_ultra_held_out_for_in_silico_screen") -> str:
    """Smoke test: train 3 epochs on the toy split and screen the shipped library."""
    _run(["python", "main_script.py", "train", split, "--epochs", "3"])
    _run(["python", "main_script.py", "predict", split, "test_screen"])
    pred = (Path(REPO_DIR) / "results" / "screen_results" / f"{split}_preds" / "test_screen" / "pred_file.csv")
    if not pred.exists():
        return "(no pred_file.csv produced)"
    import pandas as pd
    df = pd.read_csv(pred)
    # persist to volume and return a compact, ANSI-safe summary
    (Path(MODELS_MNT) / "demo_pred_file.csv").write_text(pred.read_text())
    vol.commit()
    pred_col = [c for c in df.columns if "pred" in c.lower()]
    head = df[["smiles"] + pred_col].head(5).to_string(index=False) if pred_col else df.head().to_string()
    summary = f"OK: pred_file.csv has {len(df)} rows, cols={list(df.columns)}\n\ntop rows:\n{head}"
    print(summary, flush=True)  # visible in container logs on any invocation
    return summary


@app.local_entrypoint()
def main():
    """Default: run the toy smoke test end-to-end on Modal."""
    print(demo_screen.remote())
