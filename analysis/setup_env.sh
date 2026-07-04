set -e
echo "[$(date +%T)] installing torch (CPU) + fair-esm + deps"
pip install --quiet --no-input torch --index-url https://download.pytorch.org/whl/cpu 2>&1 | tail -2
pip install --quiet --no-input fair-esm numpy biopython 2>&1 | tail -2
echo "[$(date +%T)] torch/esm install done"
python3 -c "import torch,esm; print('torch',torch.__version__,'esm ok, threads',torch.get_num_threads())"
echo "[$(date +%T)] cloning ProteinMPNN"
cd /home/user/ClaudeCode/analysis
[ -d ProteinMPNN ] || git clone --depth 1 https://github.com/dauparas/ProteinMPNN.git 2>&1 | tail -2
ls ProteinMPNN/vanilla_model_weights/ 2>/dev/null | head
echo "[$(date +%T)] SETUP COMPLETE"
