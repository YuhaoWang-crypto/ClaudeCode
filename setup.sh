#!/usr/bin/env bash
# 一键安装 FairChem（Meta FAIR Chemistry）用于材料 / MOF 吸附 / 分子晶体预测
#
# 前提：
#   1) 能访问 huggingface.co（模型权重从 HF 下载，且 UMA/OMat24/ODAC 均为 gated 仓库）
#   2) 已在 https://huggingface.co/facebook/UMA 接受许可证，并生成 HF token
#   3) Python 3.11–3.13
#
# 用法：
#   export HF_TOKEN=hf_xxx        # 你的 HuggingFace token
#   bash setup.sh                 # 默认装 CPU 版 torch；GPU 见下方注释
set -euo pipefail

PYTHON="${PYTHON:-python3}"
echo ">> 使用 Python: $($PYTHON --version)"

# 建议用虚拟环境
if [[ "${USE_VENV:-1}" == "1" ]]; then
  "$PYTHON" -m venv .venv
  # shellcheck disable=SC1091
  source .venv/bin/activate
  PYTHON=python
fi

pip install --upgrade pip

# --- 1) 先装 torch（按硬件二选一）---
# CPU 版（默认，能跑但慢）：
pip install "torch~=2.8.0" --index-url https://download.pytorch.org/whl/cpu
# GPU 版（示例 CUDA 12.4，取消注释并注掉上面一行）：
# pip install "torch~=2.8.0" --index-url https://download.pytorch.org/whl/cu124

# --- 2) 装 fairchem 及其余依赖 ---
pip install -r requirements.txt

# --- 3) 登录 HuggingFace 以拉取 gated 模型 ---
if [[ -n "${HF_TOKEN:-}" ]]; then
  "$PYTHON" - <<'PY'
from huggingface_hub import login
import os
login(token=os.environ["HF_TOKEN"])
print(">> HuggingFace 登录成功")
PY
else
  echo ">> 警告：未设置 HF_TOKEN。请先 export HF_TOKEN=hf_xxx，否则无法下载 gated 模型。"
fi

echo ">> 安装完成。用 'python examples/01_material_relax_omat.py' 做冒烟测试。"
