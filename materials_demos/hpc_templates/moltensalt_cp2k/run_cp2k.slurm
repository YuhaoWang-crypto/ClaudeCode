#!/bin/bash
#SBATCH --job-name=ms_aimd
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32        # <<adjust to your node>>
#SBATCH --time=48:00:00
#SBATCH --partition=<<PARTITION>>
#SBATCH --output=%x-%j.out

module load cp2k/2024.1              # <<your cluster's module>>

# one job per temperature point; loop or submit as an array
export OMP_NUM_THREADS=1
srun cp2k.psmp -i aimd.inp -o aimd.out

# Deliverable: aimd-pos-1.xyz + aimd-frc-1.xyz + aimd-1.ener
# -> parse with dpdata (fmt="cp2k/aimd_output") into deepmd/npy for training.
