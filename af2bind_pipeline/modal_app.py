"""Modal GPU app: run the AlphaFold2 forward pass and return AF2BIND features.

This is the only step that needs a GPU, jax and AlphaFold weights. It returns
the (L, 5120) pair-feature matrix rather than predictions, so the head can be
re-run locally with any seed/ensemble without paying for another AF2 pass.

    modal run af2bind_pipeline/modal_app.py --target 6w70 --chain A

Environment: needs MODAL_TOKEN_ID / MODAL_TOKEN_SECRET (or `modal setup`).
"""

from __future__ import annotations

import os

import modal

APP_NAME = "af2bind"

# ColabDesign v1.1.1 is the commit the AF2BIND paper used. jax is pinned to a
# version contemporary with it; newer jax removes APIs ColabDesign still calls.
image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("git", "wget", "curl", "tar", "unzip", "aria2")
    .pip_install(
        "numpy<2",
        "scipy",
        "pandas",
        "biopython",
        "ml-collections",
        "immutabledict",
        "joblib",
        "tqdm",
        "matplotlib",
        "dm-haiku==0.0.12",
        "optax==0.2.2",
        "chex==0.1.86",
        "dm-tree",
    )
    .pip_install(
        "jax[cuda12]==0.4.28",
        extra_options="-f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html",
    )
    # dm-haiku unconditionally imports flax from haiku.experimental, and a
    # current flax calls jax.sharding.AbstractMesh, which jax 0.4.28 does not
    # have. Downgrade flax to a release contemporary with the pinned jax.
    # --no-deps so it cannot drag jax back up.
    .run_commands("pip install --no-deps 'flax==0.8.5'")
    .run_commands(
        "pip install --no-deps git+https://github.com/sokrypton/ColabDesign.git@v1.1.1"
    )
    .env({"XLA_PYTHON_CLIENT_PREALLOCATE": "false"})
    # Modal imports the defining module inside the container, so the package
    # has to travel with it. It is added as a local source mount, which does not
    # invalidate the built image layers when the code changes.
    .add_local_python_source("af2bind_pipeline")
)

app = modal.App(APP_NAME, image=image)
volume = modal.Volume.from_name("af2bind-params", create_if_missing=True)
VOL = "/params"

#: AlphaFold parameter releases. The upstream notebook (the artifact shipped
#: alongside the trained heads) downloads 2021-07-14; the paper's README states
#: the reported experiments used 2022-03-02. The head was fit to features from
#: one specific checkpoint, so this choice shifts absolute p(bind). Default to
#: the notebook's, which is what the released heads were validated against.
AF2_PARAM_DATES = ("2021-07-14", "2022-03-02", "2022-12-06")
DEFAULT_AF2_PARAMS = "2021-07-14"


def _ensure_af2_params(date: str) -> str:
    """Download AlphaFold weights into the shared volume if absent."""
    import subprocess

    if date not in AF2_PARAM_DATES:
        raise ValueError(f"unknown AlphaFold params release {date!r}")
    root = f"{VOL}/af2-{date}"
    marker = f"{root}/params/params_model_1_ptm.npz"
    if os.path.isfile(marker):
        return root
    os.makedirs(f"{root}/params", exist_ok=True)
    url = f"https://storage.googleapis.com/alphafold/alphafold_params_{date}.tar"
    print(f"[af2bind] downloading AlphaFold params {date} ...", flush=True)
    subprocess.run(
        f"curl -fsSL {url} | tar x -C {root}/params", shell=True, check=True
    )
    if not os.path.isfile(marker):
        raise RuntimeError(f"AlphaFold params {date} did not unpack as expected")
    volume.commit()
    return root


@app.function(
    gpu=os.environ.get("AF2BIND_GPU", "A10G"),
    volumes={VOL: volume},
    timeout=3600,
    retries=modal.Retries(max_retries=1),
)
def pair_features(
    pdb_text: str,
    chain: str = "A",
    mask_sidechains: bool = True,
    mask_sequence: bool = False,
    af2_params: str = DEFAULT_AF2_PARAMS,
) -> dict:
    """Run AF2 on target + 20 bait residues; return AF2BIND pair features.

    Returns a dict with:
      features (L, 5120) float32, chain/resi/resn per target residue,
      plddt of the target from the AF2 pass, and the run's provenance.
    """
    import numpy as np
    from colabdesign import clear_mem, mk_afdesign_model
    from colabdesign.af.alphafold.common import residue_constants

    data_dir = _ensure_af2_params(af2_params)

    with open("/tmp/target.pdb", "w") as fh:
        fh.write(pdb_text)

    n_bait = 20
    bait_seq = "ACDEFGHIKLMNPQRSTVWY"

    clear_mem()
    af_model = mk_afdesign_model(protocol="binder", debug=True, data_dir=data_dir)
    af_model.prep_inputs(
        pdb_filename="/tmp/target.pdb",
        chain=chain,
        binder_len=n_bait,
        rm_target_sc=mask_sidechains,
        rm_target_seq=mask_sequence,
    )

    # Give every bait residue its own residue-index island (+50 apart) so AF2
    # treats them as 20 separate single-residue chains rather than one peptide.
    # Skipping this silently changes the features and invalidates the head.
    r_idx = af_model._inputs["residue_index"][-n_bait] + (1 + np.arange(n_bait)) * 50
    af_model._inputs["residue_index"][-n_bait:] = r_idx.flatten()

    af_model.set_seq(bait_seq)
    af_model.predict(verbose=False)

    outputs = af_model.aux["debug"]["outputs"]
    pair = np.asarray(outputs["representations"]["pair"])
    n_target = af_model._target_len

    pair_a = pair[:-n_bait, -n_bait:].reshape(n_target, -1)
    pair_b = pair[-n_bait:, :-n_bait].swapaxes(0, 1).reshape(n_target, -1)
    features = np.concatenate([pair_a, pair_b], -1).astype(np.float32)

    aa_order = {v: k for k, v in residue_constants.restype_order.items()}
    idx = af_model._pdb["idx"]
    aatype = af_model._pdb["batch"]["aatype"]
    # aux["plddt"] is (L_total,) for a single model but (n_models, L_total) when
    # ColabDesign keeps the per-model stack; collapse before slicing the target.
    plddt = np.asarray(af_model.aux["plddt"], dtype=np.float32)
    while plddt.ndim > 1:
        plddt = plddt.mean(0)
    plddt = plddt[:n_target]

    return {
        "features": features,
        "chain": [str(idx["chain"][i]) for i in range(n_target)],
        "resi": [int(idx["residue"][i]) for i in range(n_target)],
        "resn": [aa_order.get(int(aatype[i]), "X") for i in range(n_target)],
        "plddt": plddt.astype(np.float32),
        "meta": {
            "af2_params": af2_params,
            "mask_sidechains": bool(mask_sidechains),
            "mask_sequence": bool(mask_sequence),
            "chain_requested": chain,
            "n_target": int(n_target),
            "n_bait": n_bait,
            "bait_seq": bait_seq,
            "colabdesign": "v1.1.1",
        },
    }


@app.local_entrypoint()
def main(
    target: str = "6w70",
    chain: str = "A",
    out: str = "af2bind_out",
    mask_sidechains: bool = True,
    af2_params: str = DEFAULT_AF2_PARAMS,
    seeds: str = "0",
):
    """Convenience entrypoint: fetch a target, run AF2, score, write results."""
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from af2bind_pipeline.run import run_from_features
    from af2bind_pipeline.structure import fetch_structure

    out_dir = Path(out)
    out_dir.mkdir(parents=True, exist_ok=True)
    pdb_path = fetch_structure(target, out_dir)
    result = pair_features.remote(
        pdb_text=Path(pdb_path).read_text(),
        chain=chain,
        mask_sidechains=mask_sidechains,
        af2_params=af2_params,
    )
    seed_list = tuple(int(s) for s in seeds.split(",") if s.strip())
    run_from_features(result, pdb_path, out_dir, seeds=seed_list)
    print(f"[af2bind] wrote results to {out_dir}/")
