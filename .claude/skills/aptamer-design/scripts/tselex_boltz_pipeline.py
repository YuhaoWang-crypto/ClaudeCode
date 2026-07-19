#!/usr/bin/env python3
"""
Orchestrator SCAFFOLD for the T-SELEX × Boltz consensus in-silico SELEX loop.
See references/t-selex-integration.md for the design.

This is a TEMPLATE, not a turnkey script: the T-SELEX half needs a provisioned Linux
env (ViennaRNA, IntaRNA, a manual HDOCK download, and Selenium+browser for RNAComposer),
and the Boltz half runs through the Boltz_API MCP server (call it from your agent, not
from plain Python). Functions that need those are stubbed with clear TODOs so the control
flow is runnable/inspectable offline.

Round loop:
  pool -> ViennaRNA 2D -> IntaRNA liability -> RNAComposer 3D -> {HDOCK, Boltz} -> consensus
       -> specificity counter-screen -> select top-k -> LM-mutate winners -> repeat
"""
import csv, json, os, random, subprocess, sys

# ---- knobs ----
ROUNDS      = 4
POOL_SIZE   = 200
KEEP_TOPK   = 20
MUT_RATE    = 0.20          # doped-pool per-position mutation on winners
RNA_ALPHA   = "ACGU"

# ---------------------------------------------------------------- pool generation
def gen_pool(round_idx, winners=None):
    """Round 0: T-SELEX gen_aptamers (+ optional LM-biased seed). Later rounds: mutate winners."""
    if round_idx == 0 or not winners:
        try:
            from T_SELEX_program import gen_aptamers
            return list(gen_aptamers(seed=1, length="randomize", aptamers_num=POOL_SIZE))
        except Exception:
            # offline fallback: uniform random RNA (replace with generate_lm.py output)
            rng = random.Random(round_idx)
            return ["".join(rng.choice(RNA_ALPHA) for _ in range(rng.randint(26, 40)))
                    for _ in range(POOL_SIZE)]
    # doped mutation of winners (LM-guided in production; uniform here)
    rng = random.Random(1000 + round_idx)
    pool = []
    while len(pool) < POOL_SIZE:
        w = rng.choice(winners)
        pool.append("".join(rng.choice(RNA_ALPHA) if rng.random() < MUT_RATE else b for b in w))
    return pool

# ---------------------------------------------------------------- T-SELEX stages
def vienna_2d(pool):
    from T_SELEX_program import fold_and_composition       # ViennaRNA
    return fold_and_composition(pool)                       # df: Aptamer, MFE, MFE structure...

def intarna_liability(df, target_seq):
    """Drop self-dimerizing / pool-cross-hybridizing candidates (IntaRNA)."""
    from T_SELEX_program import intarna
    df.to_csv("pool.csv", index=False)
    return intarna("pool.csv", "Aptamer", target_seq, "target", True)

def rnacomposer_3d(df):
    from T_SELEX_program import tertiary_structure          # RNAComposer via Selenium
    return tertiary_structure(df["Aptamer"], df["MFE structure"])

def hdock_score(df, receptor_pdb, receptor_name, ligands_dir, hdock_path):
    from T_SELEX_program import Mol_docking_calc            # HDOCK
    return Mol_docking_calc(data_frame=df, MFE_column="Minimum free Energy",
                            receptor_name=receptor_name, receptor=receptor_pdb,
                            ligands_directory=ligands_dir, directory_path=hdock_path,
                            Ap_folded=True)                  # -> HDOCK score + confidence

# ---------------------------------------------------------------- Boltz stage (via MCP)
def boltz_iptm(seq, target_aa, is_rna=True):
    """TODO(agent): call Boltz_API MCP boltz_start_structure_and_binding with
    entities=[{protein: target_aa}, {rna|dna: seq}], poll, return best_sample.metrics.iptm.
    Kept as a stub so this file stays importable without the MCP context."""
    raise NotImplementedError("Run Boltz via the Boltz_API MCP from the agent loop.")

# ---------------------------------------------------------------- consensus + select
def consensus(cands):
    """cands: [{id, chem, hdock_score, boltz_iptm, ...}] -> ordered ids (best first)."""
    with open("_cons.json", "w") as fh:
        json.dump(cands, fh)
    subprocess.run([sys.executable,
                    os.path.join(os.path.dirname(__file__), "consensus_rank.py"),
                    "_cons.json"], check=False)
    # parse/order in production; here we defer to consensus_rank.py's printed ranking

def run(target_aa, target_rna_seq, receptor_pdb, receptor_name, ligands_dir, hdock_path,
        paralog_aa=None):
    winners = None
    for r in range(ROUNDS):
        print(f"\n===== SELEX round {r} =====")
        pool = gen_pool(r, winners)
        # --- T-SELEX half (needs provisioned env) ---
        # df = vienna_2d(pool)
        # df = drop liabilities via intarna_liability(df, target_rna_seq)
        # structs = rnacomposer_3d(df)           # RNA only
        # dock = hdock_score(df, receptor_pdb, receptor_name, ligands_dir, hdock_path)
        # --- Boltz half (agent/MCP): iptm = boltz_iptm(seq, target_aa) per survivor ---
        # --- specificity: repeat HDOCK+Boltz vs paralog_aa; require on >> off ---
        # cands = merge(dock, boltz, offtarget); consensus(cands); winners = top-k
        print("  [scaffold] wire T-SELEX + Boltz per references/t-selex-integration.md")
        break   # remove once stages are wired
    print("Done. Feed consensus survivors to rank_candidates.py and the report template.")

if __name__ == "__main__":
    print(__doc__)
    print("This is a scaffold — see references/t-selex-integration.md to provision & wire it.")
