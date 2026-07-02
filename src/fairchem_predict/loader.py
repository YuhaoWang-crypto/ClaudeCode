"""Shared helper for loading a UMA predictor and building an ASE calculator.

The UMA model is a single universal network with several *task heads*, one per
DFT dataset / level of theory. Pick the task that matches your system:

    omat  -> inorganic bulk materials (OMat24 level of theory)
    omol  -> isolated molecules / polymers (OMol25)
    odac  -> MOFs for direct air capture: CO2 / H2O adsorption (ODAC)
    omc   -> organic molecular crystals / polymorphs (OMC25) -- drug crystals
    oc20  -> catalysis: adsorbates on metal slabs

All heads share one checkpoint, so you download the weights once.
"""

from __future__ import annotations

VALID_TASKS = {"omat", "omol", "odac", "omc", "oc20", "oc22", "oc25"}


def get_calculator(task_name: str, model: str = "uma-s-1p1", device: str = "cpu"):
    """Return an ASE calculator backed by the UMA model for `task_name`.

    Downloading the checkpoint requires a HuggingFace account that has accepted
    the model license at https://huggingface.co/facebook/UMA and a token
    exported as HF_TOKEN (run `huggingface-cli login` once).
    """
    if task_name not in VALID_TASKS:
        raise ValueError(f"task_name={task_name!r} not in {sorted(VALID_TASKS)}")

    # Imported lazily so `--help` / install checks work without the heavy deps.
    from fairchem.core import FAIRChemCalculator, pretrained_mlip

    predictor = pretrained_mlip.get_predict_unit(model, device=device)
    return FAIRChemCalculator(predictor, task_name=task_name)
