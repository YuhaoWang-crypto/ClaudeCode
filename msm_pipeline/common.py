"""
Shared utilities for the conformational-flux (MSM/TPT) pipeline.

Data source
-----------
48 x 500 ns all-atom GROMACS trajectories of TEM beta-lactamase (C-alpha
coordinates, 263 residues, 8001 frames each; 24 us aggregate), published as
Zenodo record 6565237 by the authors of

  Avery C., Baker L., Jacobs D.J. "Functional Dynamics of Substrate
  Recognition in TEM Beta-Lactamase." Entropy 24(5):729, 2022.
  doi:10.3390/e24050729

  * 32 apo systems  = {TEM-1, TEM-2, TEM-10, TEM-52} x 8 crystal starts
                      (1ERM 1ERO 1ERQ 1HTZ 1JWP 1LHY 1XPB 3JYI)
  * 16 holo systems = the same 4 variants x 4 antibiotics
                      (AIC=ampicillin, AXL=amoxicillin,
                       CAZ=ceftazidime,  CEF=cefotaxime)

Variants (Ambler numbering):
  TEM-1  wild type
  TEM-2  Q39K
  TEM-10 R164S / E240K          (extended spectrum; R164S is in the omega loop)
  TEM-52 E104K / M182T / G238S  (extended spectrum)

Residue numbering
-----------------
The trajectories carry 263 C-alpha atoms spanning Ambler residues 26-290.
Ambler numbering for class A beta-lactamases skips 239 and 253, so the
mapping to a 0-based column index is handled by `amb_to_idx` / `idx_to_amb`.
Never index these arrays with a raw Ambler number.
"""
from __future__ import annotations

import json
import os
import subprocess
import zipfile

import numpy as np

# ----------------------------------------------------------------- paths ----
CACHE = os.environ.get("MSM_CACHE", "/tmp/claude-0/tem_cache")
DATA = os.environ.get("MSM_DATA", "/tmp/claude-0/tem_data")
PDB = os.environ.get("MSM_PDB", "/tmp/claude-0/pdb")
FIG = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "figures")

ZENODO = "https://zenodo.org/records/6565237/files"
RCSB = "https://files.rcsb.org/download"

N_RES = 263
N_FRAMES = 8001
TRAJ_NS = 500.0
DT_NS = TRAJ_NS / (N_FRAMES - 1)          # 0.0625 ns between stored frames

APO_STARTS = ["1ERM", "1ERO", "1ERQ", "1HTZ", "1JWP", "1LHY", "1XPB", "3JYI"]
VARIANTS = ["TEM-1", "TEM-2", "TEM-10", "TEM-52"]
LIGANDS = ["AIC", "AXL", "CAZ", "CEF"]
LIGAND_NAMES = {"AIC": "ampicillin", "AXL": "amoxicillin",
                "CAZ": "ceftazidime", "CEF": "cefotaxime"}
# holo runs were started from these crystal structures (paper, section 2.1)
HOLO_START = {"AIC": "1XPB", "AXL": "1XPB", "CAZ": "1HTZ", "CEF": "1HTZ"}

# Reference crystal structures.
#   1BTL / 1ZG4 / 1ERM : apo TEM-1, cryptic pocket CLOSED
#   1PZO / 1PZP        : TEM-1 with a non-competitive allosteric inhibitor
#                        bound in the H11/H12 cryptic pocket (core-disrupting)
XTALS = ["1BTL", "1ZG4", "1ERM", "1PZO", "1PZP"]
XTAL_CLOSED, XTAL_OPEN = "1BTL", "1PZO"


# ------------------------------------------------------------- numbering ----
def amb_to_idx(amb: int) -> int:
    """Ambler residue number -> 0-based column index in the (263, 3) arrays."""
    if not 26 <= amb <= 290:
        raise ValueError(f"Ambler residue {amb} outside the modelled range 26-290")
    if amb in (239, 253):
        raise ValueError(f"Ambler numbering skips residue {amb}")
    seq = amb if amb <= 238 else (amb - 1 if amb <= 252 else amb - 2)
    return seq - 26


def idx_to_amb(idx: int) -> int:
    """Inverse of `amb_to_idx`."""
    seq = 26 + idx
    return seq if seq <= 238 else (seq + 1 if seq <= 251 else seq + 2)


# ------------------------------------------------------------------ data ----
def _sh(cmd, **kw):
    return subprocess.run(cmd, check=True, **kw)


def fetch_trajectories(verbose=True):
    """Download + unzip the Zenodo archives (idempotent)."""
    os.makedirs(DATA, exist_ok=True)
    for name in ("apo", "holo"):
        out = os.path.join(DATA, name)
        if os.path.isdir(out) and os.listdir(out):
            continue
        zp = os.path.join(DATA, f"{name}.zip")
        if not os.path.exists(zp):
            if verbose:
                print(f"  downloading {name}.zip from Zenodo 6565237 ...")
            _sh(["curl", "-sS", "-L", "--retry", "4", "-o", zp,
                 f"{ZENODO}/{name}.zip?download=1"])
        if verbose:
            print(f"  unpacking {name}.zip ...")
        with zipfile.ZipFile(zp) as z:
            z.extractall(out)


def fetch_crystals(codes=XTALS, verbose=True):
    """Download reference PDB structures (idempotent)."""
    os.makedirs(PDB, exist_ok=True)
    for c in codes:
        p = os.path.join(PDB, f"{c}.pdb")
        if not os.path.exists(p):
            if verbose:
                print(f"  downloading {c}.pdb ...")
            _sh(["curl", "-sS", "-L", "--retry", "4", "-o", p, f"{RCSB}/{c}.pdb"])


def _traj_files():
    import glob
    f = sorted(glob.glob(f"{DATA}/apo/**/*.txt", recursive=True))
    f += sorted(glob.glob(f"{DATA}/holo/**/*.txt", recursive=True))
    if len(f) != 48:
        raise RuntimeError(f"expected 48 trajectory files, found {len(f)}")
    return f


def _parse_one(path):
    """
    One JEDi-exported file is (789, 8001): rows are degrees of freedom ordered
    x(1..263), y(1..263), z(1..263); columns are frames.  Returns (8001,263,3)
    in angstroms.
    """
    X = np.loadtxt(path, dtype=np.float32)
    if X.shape[0] != 3 * N_RES:
        raise RuntimeError(f"{path}: expected {3*N_RES} rows, got {X.shape[0]}")
    return X.T.reshape(-1, 3, N_RES).transpose(0, 2, 1).copy()


def build_cache(verbose=True):
    """Parse all 48 trajectories into one memory-mappable array + name list."""
    os.makedirs(CACHE, exist_ok=True)
    xyz_p, nm_p = f"{CACHE}/xyz.npy", f"{CACHE}/names.json"
    if os.path.exists(xyz_p) and os.path.exists(nm_p):
        return xyz_p, nm_p
    fetch_trajectories(verbose)
    files = _traj_files()
    from concurrent.futures import ProcessPoolExecutor
    if verbose:
        print(f"  parsing {len(files)} trajectories ...")
    with ProcessPoolExecutor(max(1, os.cpu_count() or 1)) as ex:
        arrs = list(ex.map(_parse_one, files))
    names = [os.path.basename(f).replace("_AlignedCoordinates.txt", "")
             for f in files]
    np.save(xyz_p, np.stack(arrs))
    json.dump(names, open(nm_p, "w"))
    return xyz_p, nm_p


def load(mmap=True):
    """-> (X, names). X is (48, 8001, 263, 3) float32 angstrom, memory-mapped."""
    build_cache()
    X = np.load(f"{CACHE}/xyz.npy", mmap_mode="r" if mmap else None)
    names = json.load(open(f"{CACHE}/names.json"))
    return X, names


def apo_indices(names, variant=None):
    out = [k for k, n in enumerate(names) if n.split("_")[0] in APO_STARTS]
    if variant:
        out = [k for k in out if names[k].endswith(variant)]
    return out


def holo_indices(names, variant=None):
    out = [k for k, n in enumerate(names) if n.split("_")[0] in LIGANDS]
    if variant:
        out = [k for k in out if names[k].endswith(variant)]
    return out


# -------------------------------------------------------------- geometry ----
def kabsch_rmsd(P, Q):
    """Optimally superposed RMSD between two (n,3) coordinate sets."""
    P = P - P.mean(0)
    Q = Q - Q.mean(0)
    V, _, W = np.linalg.svd(P.T @ Q)
    d = np.sign(np.linalg.det(V @ W))
    R = V @ np.diag([1.0, 1.0, d]) @ W
    return float(np.sqrt(((P @ R - Q) ** 2).sum(1).mean()))


def crystal_ca(code):
    """
    -> callable(ambler_number) giving the C-alpha xyz (angstrom) of chain A.

    Some TEM PDB entries are deposited in Ambler numbering (239/253 absent),
    others in sequential numbering; this resolves both.
    """
    import mdtraj as md
    fetch_crystals([code], verbose=False)
    t = md.load(os.path.join(PDB, f"{code}.pdb"))
    ca = {}
    for a in t.topology.atoms:
        if a.name == "CA" and a.residue.is_protein and a.residue.chain.index == 0:
            ca[a.residue.resSeq] = t.xyz[0, a.index] * 10.0
    is_ambler = 239 not in ca

    def get(amb):
        if is_ambler:
            return ca[amb]
        return ca[amb if amb <= 238 else (amb - 1 if amb <= 252 else amb - 2)]
    return get


def pair_distances(x, pairs):
    """x: (T,263,3); pairs: [(amb_i, amb_j), ...] -> (T, n_pairs) distances."""
    i = [amb_to_idx(a) for a, _ in pairs]
    j = [amb_to_idx(b) for _, b in pairs]
    return np.linalg.norm(x[:, i] - x[:, j], axis=-1)


def cached(name, fn):
    """Memoise a numpy result under CACHE/<name>.npy."""
    os.makedirs(CACHE, exist_ok=True)
    p = f"{CACHE}/{name}.npy"
    if os.path.exists(p):
        return np.load(p)
    v = fn()
    np.save(p, v)
    return v


def banner(title):
    print("=" * 78)
    print(title)
    print("=" * 78)
