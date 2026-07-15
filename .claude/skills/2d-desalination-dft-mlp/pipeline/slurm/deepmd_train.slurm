#!/bin/bash
#SBATCH --job-name=dp_train
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1               # MLP training IS GPU-bound
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --partition=gpu
#SBATCH --output=dp_train_%j.log

# DeePMD-kit training on an HPC GPU node (on-prem alternative to
# modal_train_deepmd.py).

module load cuda/12.1 2>/dev/null || true
source activate deepmd 2>/dev/null || true   # conda env with deepmd-kit[gpu]

cd ../03_mlp
dp train deepmd_input.json
dp freeze -o graph.pb
dp compress -i graph.pb -o graph-compress.pb

echo "Trained model: 03_mlp/graph-compress.pb"
echo "Check accuracy:  dp test -m graph.pb -s ../data/valid/held_out -n 1000"
