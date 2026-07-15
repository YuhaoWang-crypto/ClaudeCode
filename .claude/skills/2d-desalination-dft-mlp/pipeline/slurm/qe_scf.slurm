#!/bin/bash
#SBATCH --job-name=qe_scf
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32       # QE (plane-wave DFT) is CPU/MPI-bound
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --partition=cpu
#SBATCH --output=qe_scf_%j.log

# Quantum ESPRESSO SCF / relax on an HPC cluster (the on-prem alternative to
# modal_run_dft.py). Edit modules/paths for your site.

module load quantum-espresso/7.2 openmpi 2>/dev/null || true

export OMP_NUM_THREADS=1
INFILE=${1:-../02_dft/qe/scf_graphene.in}

srun --mpi=pmix pw.x -in "$INFILE" > "${INFILE%.in}.out"

echo "Done. Grep the final energy with:  grep '!    total energy' ${INFILE%.in}.out"
