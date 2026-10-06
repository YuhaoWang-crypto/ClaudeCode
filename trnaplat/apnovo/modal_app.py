"""Run the corrected AP Novo package on Modal.

The app is split in two deliberately, because the two halves have different
licence status:

| function | weights | licence |
|---|---|---|
| `featurize` | none | Apache 2.0 code only — runnable by anyone |
| `gpu_probe` | none | Apache 2.0 code only |
| `generate` | AP Novo Generator | ⚠️ **non-commercial organizations only** |

## ⚠️ Why `generate` is gated

`WEIGHTS_TERMS_OF_USE.md` and `WEIGHTS_PROHIBITED_USE_POLICY.md` restrict the
generator parameters to non-commercial organizations, and the restriction is
tighter than it first looks:

* "only non-commercial organizations (*i.e.*, universities, non-profit
  organizations and research institutes, educational, journalism and government
  bodies) may access the AP Novo Generator … The AP Novo Generator Assets are
  **not available for any other types of organization, even if conducting
  non-commercial work**."
* such an organization "**must not process Input provided by or on behalf of any
  commercial organization**".
* `OUTPUT_TERMS_OF_USE.md` carries the same restriction forward: the designs are
  "for non-commercial use only", must not be shared with a commercial
  organization except via publication, and must not be used to train another
  design model.

So a motif prepared for a commercial platform is Input the licence excludes, and
the designs would be Output that platform could not use. `generate` therefore
refuses unless `--attest-non-commercial` is passed, and it does not ask Claude to
judge anyone's organizational status — that is the caller's statement to make.

## What `featurize` proves without any weights

Everything up to the first diffusion step: `Manifest.from_file`, then
`MotifSpec.parsed()` (CIF load, author→internal renumbering, designable-length
sampling, residue mapping) and `create_designable_protein_features` plus
`add_ligand_features` / `add_unindexed_motif_features`. That is the whole input
path, including the CCD lookup of `YLY` — the part a grammar check cannot reach,
and where a malformed motif actually bites.

```bash
modal run trnaplat/apnovo/modal_app.py::featurize
modal run trnaplat/apnovo/modal_app.py::gpu_probe
# only with the attestation above:
modal run trnaplat/apnovo/modal_app.py::generate --attest-non-commercial \\
    --job ncaa_pocket_unindexed --num-designs 2 --num-sampling-steps 50
```
"""

from __future__ import annotations

import pathlib

import modal

REPO_URL = "https://github.com/google-deepmind/alphaprotein-novo.git"
REPO_DIR = "/opt/alphaprotein-novo"
PACKAGE_DIR = "/opt/fixed"
WEIGHTS_DIR = "/weights"

GENERATOR_URL = "https://storage.googleapis.com/alphaprotein_novo/generator.bin.zst"
AF3_LA_URL = "https://storage.googleapis.com/alphafold3/af3_leaving_atom.bin.zst"

LOCAL_FIXED = pathlib.Path(__file__).parent / "fixed"

#: AF3 builds native extensions, so the image needs a toolchain. `build_data`
#: compiles the CCD into a pickle and is the step that makes a `YLY` lookup work.
image = (
    modal.Image.debian_slim(python_version="3.12")
    # ✅ zlib1g-dev is load-bearing: without it the alphafold3 wheel fails at
    # "Could NOT find ZLIB (missing: ZLIB_LIBRARY ZLIB_INCLUDE_DIR)".
    .apt_install("git", "build-essential", "cmake", "pkg-config",
                 "zlib1g-dev", "libbz2-dev", "zstd", "wget", "curl")
    .pip_install("jax[cuda12]>=0.4.30", "numpy", "pydantic", "etils[epath]",
                 "absl-py", "dm-tree", "zstandard")
    .run_commands(
        f"git clone --depth 1 {REPO_URL} {REPO_DIR}",
        "pip install git+https://github.com/google-deepmind/alphafold3.git",
        f"pip install -e {REPO_DIR}",
        "build_data",
    )
    .add_local_dir(LOCAL_FIXED, PACKAGE_DIR, copy=True)
)

app = modal.App("apnovo-pylrs")
weights = modal.Volume.from_name("apnovo-weights", create_if_missing=True)


@app.function(image=image, timeout=900)
def featurize(manifest: str = "manifest_fixed.json") -> dict:
    """Drive every stage of the input path that needs no model parameters."""
    import json
    import sys
    import traceback

    sys.path.insert(0, f"{REPO_DIR}/src")
    from alphaprotein_novo.data import design_manifest, motif_spec

    package = pathlib.Path(PACKAGE_DIR)
    results: dict = {"manifest": manifest, "jobs": {}}

    try:
        parsed = design_manifest.Manifest.from_file(package / manifest)
        results["schema"] = "PASSED"
    except Exception as err:                                  # noqa: BLE001
        results["schema"] = f"FAILED: {type(err).__name__}: {err}"
        return results

    for job in parsed.designs:
        entry: dict = {}
        try:
            spec = motif_spec.MotifSpec(
                name=job.name,
                input_file=str(package / job.input_file),
                motif_str=job.motif_str,
                is_author_naming=job.is_author_naming,
                motif_atoms=job.motif_atoms,
                reseq_residues=job.reseq_residues,
                seq_length=job.seq_length,
                unindexed_motif_residues=job.unindexed_motif_residues,
            ).parsed(seed=0)
            entry["sampled_motif_str"] = spec.sampled_motif_str
            entry["is_unindexed"] = spec.is_unindexed
            # ResidueMap's arrays are (num_tokens,) over protein residues plus
            # ligand atoms, so this is the designed protein's token count.
            residue_map = spec.residue_map
            entry["n_mapped_tokens"] = (
                int(len(residue_map.source_residue_ids))
                if residue_map is not None else None
            )
            # ⚠️ source_residue_ids is an object array whose designable slots
            # hold None, so it cannot be compared numerically. Count the slots
            # that actually come from the input structure instead.
            entry["n_fixed_from_source"] = (
                sum(1 for v in residue_map.source_residue_ids if v is not None)
                if residue_map is not None else None
            )
            entry["parsed"] = "PASSED"
        except Exception as err:                              # noqa: BLE001
            entry["parsed"] = f"FAILED: {type(err).__name__}: {err}"
            entry["traceback"] = traceback.format_exc()[-1500:]
            results["jobs"][job.name] = entry
            continue

        # The featurisation itself, via the same call `run_generator.py` makes:
        # this is where a bad motif or an unknown CCD code stops being a string
        # problem and becomes a tensor problem. Nothing here loads weights.
        try:
            dinput = spec.to_diffusion_input(max_num_res=512,
                                             max_num_res_unindexed=128)
            shapes = {}

            def walk(obj, prefix=""):
                shape = getattr(obj, "shape", None)
                if shape is not None:
                    shapes[prefix.rstrip(".")] = tuple(int(d) for d in shape)
                    return
                items = (obj.items() if isinstance(obj, dict)
                         else getattr(obj, "__dict__", {}).items())
                for key, value in items:
                    walk(value, f"{prefix}{key}.")

            walk(dinput)
            entry["featurised"] = "PASSED"
            entry["n_feature_arrays"] = len(shapes)
            entry["feature_shapes"] = dict(sorted(shapes.items())[:14])
        except Exception as err:                              # noqa: BLE001
            entry["featurised"] = f"FAILED: {type(err).__name__}: {err}"
            entry["traceback"] = traceback.format_exc()[-1500:]
        results["jobs"][job.name] = entry

    print(json.dumps(results, indent=2, default=str))
    return results


@app.function(image=image, gpu="T4", timeout=600)
def gpu_probe() -> dict:
    """Confirm the image's JAX actually sees a GPU, before anything expensive."""
    import jax

    devices = [{"platform": d.platform, "kind": d.device_kind, "id": d.id}
               for d in jax.devices()]
    out = {"jax": jax.__version__, "backend": jax.default_backend(),
           "devices": devices}
    print(out)
    return out


@app.function(image=image, volumes={WEIGHTS_DIR: weights},
              gpu="A100-80GB", timeout=60 * 60 * 4)
def generate(attest_non_commercial: bool = False,
             job: str = "ncaa_pocket_unindexed",
             manifest: str = "manifest_fixed.json",
             num_designs: int = 2,
             num_sampling_steps: int = 50) -> dict:
    """Run the weighted generator. ⚠️ Non-commercial organizations only.

    Refuses without an explicit attestation, because the licence turns on a fact
    about the caller's organization that nothing here can check: the AP Novo
    Generator Assets are unavailable to any organization that is not a
    university, non-profit, research institute, educational, journalism or
    government body, and such a body must not process Input supplied by or on
    behalf of a commercial organization.
    """
    import json
    import subprocess
    import sys

    if not attest_non_commercial:
        raise RuntimeError(
            "Refusing to download or run the AP Novo Generator parameters.\n"
            "WEIGHTS_TERMS_OF_USE.md restricts them to non-commercial "
            "organizations (universities, non-profits, research institutes, "
            "educational, journalism, government bodies) and forbids processing "
            "Input provided by or on behalf of a commercial organization. "
            "OUTPUT_TERMS_OF_USE.md carries the same restriction to the designs.\n"
            "Pass --attest-non-commercial only if that describes your "
            "organization and this input."
        )

    weights_path = pathlib.Path(WEIGHTS_DIR)
    generator = weights_path / "apnovo_generator" / "generator.bin.zst"
    af3_la = weights_path / "af3_la" / "af3_leaving_atom.bin.zst"
    for target, url in ((generator, GENERATOR_URL), (af3_la, AF3_LA_URL)):
        target.parent.mkdir(parents=True, exist_ok=True)
        if not target.exists():
            print(f"fetching {url} -> {target}")
            subprocess.run(["wget", "-q", "-O", str(target), url], check=True)
        print(f"{target.name}: {target.stat().st_size} bytes")
    weights.commit()

    # Cut the manifest down to one job at smoke-test scale, so a first run is
    # minutes rather than hours.
    package = pathlib.Path(PACKAGE_DIR)
    spec = json.loads((package / manifest).read_text())
    spec["designs"] = [d for d in spec["designs"] if d["name"] == job]
    if not spec["designs"]:
        raise RuntimeError(f"no job named {job!r} in {manifest}")
    spec["defaults"]["num_designs"] = num_designs
    spec["defaults"]["num_sampling_steps"] = num_sampling_steps
    # LigandMPNN lives in its own interpreter, so the smoke run skips it -- and
    # `folding.inputs` has to follow, because the pipeline refuses the pair:
    # "The manifest specifies folding.inputs=["resequenced"], but
    #  resequence.enabled is false."
    spec["resequence"] = {"enabled": False}
    spec.setdefault("folding", {})["inputs"] = ["generated"]
    run_manifest = pathlib.Path("/tmp/manifest_run.json")
    run_manifest.write_text(json.dumps(spec, indent=2))

    cmd = [
        sys.executable, f"{REPO_DIR}/run_pipeline.py",
        f"--manifest={run_manifest}",
        "--output_dir=/tmp/out",
        "--only_stage=generation",
        f"--apn_model_dir={generator.parent}",
    ]
    print(" ".join(cmd))
    proc = subprocess.run(cmd, cwd=PACKAGE_DIR, capture_output=True, text=True)
    out = {"returncode": proc.returncode,
           "stdout": proc.stdout[-4000:], "stderr": proc.stderr[-4000:]}
    produced = sorted(str(p.relative_to("/tmp/out"))
                      for p in pathlib.Path("/tmp/out").rglob("*")
                      if p.is_file()) if pathlib.Path("/tmp/out").exists() else []
    out["files"] = produced[:40]
    print(json.dumps({k: v for k, v in out.items() if k != "stdout"}, indent=2))
    return out


@app.local_entrypoint()
def main(manifest: str = "manifest_fixed.json"):
    print("=== featurize (no weights, licence-clean) ===")
    result = featurize.remote(manifest)
    for name, entry in result.get("jobs", {}).items():
        print(f"\n  {name}")
        for key in ("parsed", "featurised", "sampled_motif_str",
                    "n_mapped_tokens", "n_features"):
            if key in entry:
                print(f"    {key}: {entry[key]}")
        if "traceback" in entry:
            print(f"    traceback tail:\n{entry['traceback']}")
