"""
Shared helpers for the cell-state-space pipeline.

Network access notes (this environment):
  * EBI BioModels REST search works, but ONLY when redirects are followed
    (a bare request returns HTTP 301). requests follows by default.
  * biomodels.org / EBI *download* endpoints are proxy-blocked; SBML must come
    from the one-repo-per-model GitHub mirror `github.com/biomodels/<ID>`,
    the same route M20b established.
"""
import json
import os
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
FIGDIR = os.path.join(ROOT, "figures")
CACHE = os.path.join(FIGDIR, "_cellstate_cache")
SBML_CACHE = os.path.join(FIGDIR, "_sbml_cache")

SEARCH = "https://www.ebi.ac.uk/biomodels/search"
MIRROR = "https://raw.githubusercontent.com/biomodels/{id}/master/{id}/{id}.xml"


def _ensure(d):
    os.makedirs(d, exist_ok=True)
    return d


def biomodels_search(query, num=25, timeout=60):
    """Query the live BioModels index. Cached on disk so re-runs are offline.

    Returns the parsed JSON dict, or None if the network is unavailable AND
    nothing is cached (callers must degrade honestly rather than invent hits).
    """
    _ensure(CACHE)
    key = "".join(c if c.isalnum() else "_" for c in query)[:60]
    path = os.path.join(CACHE, f"search_{key}_{num}.json")
    if os.path.exists(path) and os.path.getsize(path) > 20:
        with open(path) as fh:
            return json.load(fh)
    url = f"{SEARCH}?query={query.replace(' ', '%20')}&format=json&numResults={num}"
    try:
        out = subprocess.run(
            ["curl", "-sL", "--max-time", str(timeout), url],
            capture_output=True, text=True, check=True).stdout
        data = json.loads(out)
    except Exception:
        return None
    if "models" not in data:
        return None
    with open(path, "w") as fh:
        json.dump(data, fh)
    return data


def fetch_sbml(biomd_id, timeout=90):
    """Download a curated SBML from the GitHub mirror. Returns path or None."""
    _ensure(SBML_CACHE)
    path = os.path.join(SBML_CACHE, f"{biomd_id}.xml")
    if os.path.exists(path) and os.path.getsize(path) > 500:
        return path
    try:
        subprocess.run(["curl", "-sL", "--max-time", str(timeout),
                        MIRROR.format(id=biomd_id), "-o", path], check=True)
    except Exception:
        return None
    if not os.path.exists(path) or os.path.getsize(path) < 500:
        return None
    return path


def load_rr(biomd_id):
    """Load a curated model into libRoadRunner. Returns (rr, path) or (None, None)."""
    path = fetch_sbml(biomd_id)
    if path is None:
        return None, None
    try:
        import roadrunner
        roadrunner.Logger.setLevel(roadrunner.Logger.LOG_ERROR)
        return roadrunner.RoadRunner(path), path
    except Exception:
        return None, path


def curated(model_id):
    """BioModels convention: BIOMD* = manually curated, MODEL* = not curated."""
    return str(model_id).startswith("BIOMD")


def banner(title):
    line = "=" * 76
    print(f"\n{line}\n{title}\n{line}")


def save_json(name, obj):
    _ensure(CACHE)
    path = os.path.join(CACHE, name)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, default=float)
    return path


def figpath(name):
    _ensure(FIGDIR)
    return os.path.join(FIGDIR, name)
