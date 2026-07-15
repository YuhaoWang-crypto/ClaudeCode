#!/bin/bash
#SBATCH --job-name=nemd
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1               # LAMMPS + DeePMD inference is GPU-accelerated
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH --output=nemd_%j.log

# Production NEMD desalination run (on-prem alternative to a Modal GPU launch).

module load cuda/12.1 2>/dev/null || true
source activate deepmd 2>/dev/null || true   # env with lammps + deepmd pair style

cd ../04_nemd
lmp -in in.desalination.lammps

echo "Analysing flux + rejection..."
python analyze_flux_rejection.py --demo   # replace with a dump-driven call:
# python -c "import analyze_flux_rejection as a; \
#   zw,zi,t=a.load_lammps_dump('traj.lammpstrj',water_types=[3],ion_types=[4,5]); \
#   print(a.summarize(zw,zi,z_membrane=0.0,times_ns=t*5e-4,area_A2=150.,\
#         n_ions_feed_initial=20,pressure_MPa=100.))"
