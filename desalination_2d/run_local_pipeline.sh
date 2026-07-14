#!/bin/bash
# One command to exercise everything that runs WITHOUT the heavy DFT/MD engines.
# Verifies the whole scaffold end-to-end: geometry -> DFT/NEB inputs -> data
# plumbing -> piston calibration -> NEMD analysis -> toy validation.
#
# The steps that need real compute (QE SCF, CI-NEB, DeePMD training, LAMMPS
# NEMD) are printed as the exact commands to run on Modal / HPC.
set -e
cd "$(dirname "$0")"
echo "############ 1. BUILD GEOMETRIES ############"
python3 01_build/build_graphene_pore.py --pore-radius 3.5 --out figures/graphene_pore
python3 01_build/build_mos2_pore.py     --pore-radius 3.5 --out figures/mos2_pore
python3 01_build/build_saltwater_box.py --membrane figures/graphene_pore.xyz \
        --n-water-feed 300 --n-water-perm 300 --n-nacl 10 --out figures/nemd_cell

echo; echo "############ 2. GENERATE DFT INPUTS (Milestone 1 & 2) ############"
python3 02_dft/make_qe_input.py --xyz figures/graphene_pore.xyz \
        --out 02_dft/qe/scf_graphene_ready.in --ecutwfc 60 --kpts 3 3 1
python3 02_dft/build_neb_endpoints.py --membrane figures/graphene_pore.xyz \
        --ion Na --out 02_dft/qe/neb_na_ready.in

echo; echo "############ 3. MLP DATA PLUMBING (Milestone 3) ############"
python3 03_mlp/prepare_deepmd_data.py --demo --out data/train/demo_system

echo; echo "############ 4. PISTON FORCE CALIBRATION (Milestone 4) ############"
python3 04_nemd/piston_force.py --data figures/nemd_cell.lammps-data \
        --pressure-mpa 100 --n-piston 60

echo; echo "############ 5. NEMD ANALYSIS + TOY VALIDATION ############"
python3 04_nemd/analyze_flux_rejection.py --demo
python3 05_toy_validation/toy_nemd.py --plot figures/toy_nemd.png

cat <<'EOF'

############ NEXT: the compute-heavy steps (run on Modal / HPC) ############
  M1  QE SCF        : (cd 02_dft && modal run modal_run_dft.py --infile qe/scf_graphene_ready.in)
  M2  CI-NEB        : neb.x -inp 02_dft/qe/neb_na_ready.in > neb.out
                      python 02_dft/parse_neb_barrier.py --dat neb_na.dat
  M3  DFT->data     : python 03_mlp/prepare_deepmd_data.py --from-dft <qe_runs> --out data/train/pore_crossing
      DP-GEN loop   : dpgen run 03_mlp/dpgen_param.json machine.json
      train MLP     : (cd 03_mlp && modal run modal_train_deepmd.py --input deepmd_input.json)
  M4  NEMD 1 ns     : lmp -in 04_nemd/in.desalination.lammps   (or sbatch slurm/lammps_nemd.slurm)
  M5  scans + MoS2  : loop M1-M4 over pore radius / pressure, then rerun for figures/mos2_pore.xyz
EOF
echo; echo "LOCAL PIPELINE OK."
