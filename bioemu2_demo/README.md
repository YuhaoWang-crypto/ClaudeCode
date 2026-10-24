# BioEmu2 CPU 演示（chignolin）

在无 GPU 的机器（4 核 CPU、15 GB 内存）上运行 [microsoft/bioemu2](https://huggingface.co/microsoft/bioemu2)
的最小演示。官方只支持 NVIDIA GPU；本演示通过一处配置改动让它在 CPU 上跑通小体系。

## 结果

体系：chignolin CLN025（`GYDPETGTWG`，10 个残基），与实验 NMR 结构 PDB 5AWL 对比 Cα-RMSD。

| 步骤 | 耗时（CPU） | 结果 |
|---|---|---|
| 加载模型（约 5 亿参数，2 GB） | 29 s | — |
| 结构预测，1 个样本 | 55 s | RMSD 1.12 Å |
| 平衡采样，4 个样本 | 104 s | RMSD 0.79 / 1.21 / 0.92 / 0.85 Å，4/4 折叠（< 2 Å） |

输出在 `outputs/`：每个目录包含 `topology.cif` + `samples.xtc`（可用 PyMOL、VMD、ChimeraX 或 mdtraj 打开），
运行日志为 `outputs/demo.log`。

## 复现步骤

```bash
# 1. 获取代码并安装（CPU 版 PyTorch，不装 cuda extra）
git clone https://github.com/microsoft/bioemu2.git
python3 -m venv venv && . venv/bin/activate
pip install "torch>=2.11,<2.14" --index-url https://download.pytorch.org/whl/cpu
pip install ./bioemu2

# 2. 下载权重，并把 GPU 专用内核换成纯 PyTorch 实现
mkdir ckpt
curl -L -o ckpt/checkpoint.ckpt https://huggingface.co/microsoft/bioemu2/resolve/main/checkpoints/bioemu2-preview/checkpoint.ckpt
curl -L https://huggingface.co/microsoft/bioemu2/resolve/main/checkpoints/bioemu2-preview/config.yaml \
  | sed 's/trimul_variant: cueq/trimul_variant: full/' > ckpt/config.yaml

# 3. 运行（参数：权重路径、输出目录、平衡采样数）
python demo_chignolin.py ckpt/checkpoint.ckpt out 4
```

## 注意事项

- 唯一的改动是 `trimul_variant: cueq → full`。`cueq` 依赖 NVIDIA 的 cuEquivariance 内核，无法在 CPU 上运行；
  `full` 是同一运算的纯 PyTorch 实现，两者权重布局相同。
- MSA 会自动从公共 ColabFold 服务器（api.colabfold.com）获取，因此需要联网。
- CPU 只适合很小的体系。官方测试过的最小硬件是 80 GB 的 A100；在 CPU 上跑较大的蛋白（例如 214 个残基的
  1AKE）会非常慢。
- 样本数只有 4 个，只能说明模型能跑通并给出合理结构，不能作为对构象分布的定量评估。
