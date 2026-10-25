# RFD4-Proteina 在无 GPU 机器上跑 demo（实测可行）

实测环境：Linux x86_64，4 核 CPU，15 GB 内存，**无 GPU**，Python 3.12。
结论：官方要求 NVIDIA GPU，但 CPU-only 路径能跑通，需要改一行代码（见第 3 步）。

## 1. 克隆 + 建环境（CPU 版 torch）

```bash
git clone --depth 1 https://github.com/RosettaCommons/RFD4-Proteina.git
cd RFD4-Proteina
# 系统库（openbabel 等需要）
apt-get install -y --no-install-recommends libxrender1 libx11-6 libxext6
uv venv --python 3.12 .venv && source .venv/bin/activate
uv pip install --no-sources --only-binary tmol --torch-backend cpu -r pyproject.toml
uv pip install --no-sources --no-deps .
uv pip install --only-binary tmol "tmol>=0.2.0,<0.3"
```

不要用 `make setup`（默认走 CUDA）。装完约 5 GB，耗时几分钟。
`torch` 应为 `2.14.1+cpu`。

## 2. 权重

第一次运行会自动从 Hugging Face 拉 `nvidia/RFD4-Proteina`（约 1.7 GB，含 autoencoder）。
设 `HF_HOME` 控制缓存位置。公开仓库，无需 token。

## 3. 必须的一行补丁

`src/rfproteina/generate.py` 里 Lightning Trainer 把 accelerator 写死成 GPU，
导致配置项 `hardware.accelerator` 被忽略，CPU 上直接报
`MisconfigurationException: No supported gpu backend found!`：

```python
# src/rfproteina/generate.py, 约第 709 行
-        accelerator="gpu",
+        accelerator=cfg.hardware.accelerator,
```

训练代码（`src/rfproteina/train.py:535`）读的就是 `cfg_exp.hardware.accelerator`，
所以这个改法与仓库自身约定一致。改完重装：`uv pip install --no-sources --no-deps .`

模型本身对 CPU 是友好的：`rfproteina/__init__.py` 只在 `torch.cuda.is_available()`
时才加载 cuEquivariance 融合算子，否则自动回退到纯 torch 实现。

## 4. 跑 demo

无条件单体设计（输入样例仓库里自带）：

```bash
rfd4-generate inference_experiment=quickstart \
  inputs=data/benchmarks/monomer/monomer-short.cif.zst \
  hardware.accelerator=cpu precision=32-true num_workers=0 \
  out_dir=./runs/monomer
```

条件化设计（PD-L1 binder，README 的例子，把设计长度固定成 60 以省时间）；
存为 `binder.yaml`：

```yaml
binder:
  input: docs/examples/inputs/5o45_cropped.cif.zst
  contig: "60-60, /0, A17-131"
  conditions: {C_CRDGRP: {1: [{select: "sel('A/*/*/[N,CA,C,O]')"}]},
               C_SEQ: {true: [{select: "sel('A')"}]}}
```

```bash
rfd4-generate inference_experiment=quickstart inputs=binder.yaml \
  hardware.accelerator=cpu precision=32-true num_workers=0 \
  out_dir=./runs/binder
```

要点：
- `inference_experiment=quickstart` 跳过打分和 `torch.compile` 预热。
- `precision=32-true`：CPU 上比默认的 `bf16-mixed` 稳。
- `num_workers=0`：核少时省内存。
- 不要装 metrics/refolding 环境（AlphaFold2/3），CPU 上不现实；
  `quickstart` 预设本身就不跑打分。

## 5. 实测耗时与结果

| demo | 规模 | 耗时（4 核 CPU） |
|---|---|---|
| 单体 | 100 残基，60 步 | 3 分 38 秒 |
| PD-L1 binder | 60 设计 + 115 靶标 = 175 残基 | 12 分 16 秒 |

内存峰值远低于 15 GB，没有触发交换。

单体输出 `runs/monomer/inference_samples/monomer-short/monomer-short_generated.cif`
几何合理：100 个残基，相邻 CA–CA 平均 3.82 Å（理论 3.8），回转半径 14.9 Å，
序列与骨架共设计：

```
MTDEEILEEALRMIAEGKSEEEVNDFLEENNASSEVKLKVNLVFKVLKDEEVDIDDLDTLKSLLSKAESLELTDLASMLSEKISELETKKERLEEEEEEE
```

Binder 输出 `runs/binder/inference_samples/binder_0/binder_0_generated.cif`
两条链，条件化按预期生效：
- 链 B：新设计的 60 残基 binder，`SMATDYFNKGNTYYEEGKYYEALSEYMKALALAHENEDWDLASKIAENRNLVLSKIDELR`
- 链 A：115 残基靶标，序列保持原生 PD-L1（`NAFTVTVPKDLYVVEYGSNMTIECKFP...`），即 `C_SEQ` + `C_CRDGRP` 固定住了靶标
- 60 个 binder 残基里有 19 个落在靶标 5 Å 以内 —— 是个真实的结合界面，不是飘在旁边

## 6. CPU 的现实限制

- 能力只到「小规模试跑、验证输入/配置、学 API」。真实设计活动要 GPU：
  要出成百上千个设计、要 refolding 打分（AF2/AF3）、要 reward-guided beam search。
- 长度放大后时间按注意力的平方增长，几百残基的复合体在 CPU 上不实用。
- 想无痛体验：官方 Hugging Face Space（ZeroGPU 或专用 GPU 池）不用装环境。
