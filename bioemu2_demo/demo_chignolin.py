"""BioEmu2 CPU demo: chignolin (CLN025, GYDPETGTWG).

1. Predict the folded structure (mode="structure").
2. Sample an equilibrium ensemble (mode="equilibrium").
3. Compare against the experimental NMR structure PDB 5AWL (C-alpha RMSD).
"""

import sys
import time
from pathlib import Path

import mdtraj
import numpy as np
import requests
import torch

from bioemu2 import EmuSession, System, load_model, save_frames
from bioemu2.paths import EXAMPLES_DIR

ckpt = Path(sys.argv[1])
out = Path(sys.argv[2])
n_eq = int(sys.argv[3]) if len(sys.argv) > 3 else 8
torch.set_num_threads(4)

t0 = time.time()
model = load_model(checkpoint_path=ckpt, device="cpu")
print(f"[load] model loaded in {time.time() - t0:.1f}s", flush=True)

system = System.from_yaml(EXAMPLES_DIR / "systems/chignolin.yaml")
print(system, flush=True)

t0 = time.time()
pred_session = EmuSession(system, model=model, mode="structure")
pred = pred_session.sample(num_samples=1, seed=0)
save_frames(pred, out / "structure")
print(f"[structure] 1 sample in {time.time() - t0:.1f}s", flush=True)

t0 = time.time()
eq_session = EmuSession(system, model=model, mode="equilibrium", msa=pred_session.msa)
eq = eq_session.sample(num_samples=n_eq, seed=0)
save_frames(eq, out / "equilibrium")
print(f"[equilibrium] {n_eq} samples in {time.time() - t0:.1f}s", flush=True)

# Reference: experimental CLN025 structure (PDB 5AWL, first model).
ref_path = out / "5AWL.pdb"
ref_path.write_text(requests.get("https://files.rcsb.org/download/5AWL.pdb", timeout=60).text)
ref = mdtraj.load(ref_path)[0]
ref_ca = ref.atom_slice(ref.topology.select("name CA"))


def ca_rmsd(traj: mdtraj.Trajectory) -> np.ndarray:
    ca = traj.atom_slice(traj.topology.select("name CA"))
    assert ca.n_atoms == ref_ca.n_atoms, (ca.n_atoms, ref_ca.n_atoms)
    ca = mdtraj.Trajectory(ca.xyz, ref_ca.topology)
    return mdtraj.rmsd(ca, ref_ca) * 10.0  # nm -> Angstrom


pred_rmsd = ca_rmsd(pred)
eq_rmsd = ca_rmsd(eq)
print(f"[result] predicted structure CA-RMSD vs 5AWL: {pred_rmsd[0]:.2f} A")
print("[result] equilibrium samples CA-RMSD vs 5AWL (A):", np.round(eq_rmsd, 2).tolist())
print(f"[result] fraction folded (RMSD < 2 A): {(eq_rmsd < 2.0).mean():.2f}")
