#!/usr/bin/env python3
"""
DTU Health Tech multi-model orchestrator (BioLib cloud backend).

Runs one or more protein sequences through several DTU immunology / epitope
models hosted on BioLib, in parallel, and merges results into a single
per-residue table plus summary scores. No BioLib token required for the public
DTU apps (requires_user_identity=False).

Sequence-input models wired in:
  - DTU/ImmunoGeNN : population-level MHC-II immunogenicity risk (pIRS), per 15mer peptide
  - DTU/BepiPred-3 : linear B-cell epitope probability, per residue
  - DTU/NetSurfP-3 : surface accessibility (RSA) + secondary structure, per residue

Structure-input models (DTU/DiscoTope-3) take a PDB, handled separately.
"""
import os, json, argparse, traceback
from concurrent.futures import ThreadPoolExecutor, as_completed

os.environ.setdefault("BIOLIB_CACHE_DIR", "/tmp/biolib_cache")
os.makedirs(os.environ["BIOLIB_CACHE_DIR"], exist_ok=True)
import biolib

MODELS = {
    "immunogenn": {"uri": "DTU/ImmunoGeNN", "arg": "--fasta_file",
                   "outputs": ["pIRS.csv", "scores.csv"], "kind": "peptide"},
    "bepipred":   {"uri": "DTU/BepiPred-3", "arg": None,
                   "outputs": ["raw_output.csv"], "kind": "residue"},
    "netsurfp":   {"uri": "DTU/NetSurfP-3", "arg": "--input_data",
                   "outputs": [], "kind": "residue"},
}
BEPIPRED_ARG = os.environ.get("BEPIPRED_ARG", "-i")


def run_one(model_key, fasta_path, outroot, timeout=1200):
    m = MODELS[model_key]
    arg = m["arg"] or BEPIPRED_ARG
    outdir = os.path.join(outroot, model_key)
    os.makedirs(outdir, exist_ok=True)
    fasta_name = os.path.basename(fasta_path)
    rec = {"model": model_key, "uri": m["uri"], "status": None,
           "outdir": outdir, "files": [], "error": None}
    try:
        app = biolib.load(m["uri"], suppress_version_warning=True)
        job = app.cli(args=f"{arg} {fasta_name}", files=[fasta_path],
                      blocking=True, timeout=timeout)
        rec["status"] = job.get_status()
        if rec["status"] == "completed":
            job.save_files(outdir)
            rec["files"] = sorted(os.listdir(outdir))
        else:
            rec["error"] = "job did not complete"
    except Exception as e:
        rec["status"] = "error"
        rec["error"] = f"{type(e).__name__}: {e}"
        rec["trace"] = traceback.format_exc()
    return rec


def run_all(fasta_path, outroot, models, max_workers=4, timeout=1200):
    os.makedirs(outroot, exist_ok=True)
    results = {}
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futs = {ex.submit(run_one, k, fasta_path, outroot, timeout): k for k in models}
        for fut in as_completed(futs):
            k = futs[fut]
            results[k] = fut.result()
            print(f"[{results[k]['status']}] {k} -> {results[k]['files']}", flush=True)
    with open(os.path.join(outroot, "run_manifest.json"), "w") as f:
        json.dump(results, f, indent=2)
    return results


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--models", nargs="+",
                    default=["immunogenn", "bepipred", "netsurfp"])
    ap.add_argument("--timeout", type=int, default=1200)
    a = ap.parse_args()
    res = run_all(a.fasta, a.outdir, a.models, timeout=a.timeout)
    print(json.dumps({k: v["status"] for k, v in res.items()}, indent=2))
