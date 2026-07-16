# MatterGen 安装与 Demo 测试

对 [microsoft/mattergen](https://github.com/microsoft/mattergen)（Nature 2025，跨元素周期表的无机材料生成扩散模型）做了**安装 + 无条件生成 demo**，分别在 **CPU** 和 **Modal 云端 GPU** 上跑通并保存了结果。

## 结论

| 环境 | 状态 | batch | 1000 步采样耗时 | 速度 | 产出 |
|------|------|-------|----------------|------|------|
| 本地 CPU（无 GPU） | ✅ 跑通 | 2 | ~9 分 22 秒 | ~1.8 it/s | 2 个结构 |
| Modal GPU（NVIDIA T4） | ✅ 跑通 | 16 | ~6 分 16 秒 | ~2.6 it/s | 16 个结构 |

> T4 单次处理 batch=16，单位结构吞吐约为 CPU（batch=2）的 **8 倍**。

生成的样例结构（见 `results/`）：
- **CPU**：`Ce₃ZrSc`、`Sb₂OF₃`
- **GPU**：`BaHoGa₄`、`Li₂In₂AgI₁₁`、`Ba₂LiH₆Rh`、`LiZrRh₂`、`CsBi₂F₆`、`Ba₂Sn₂Pt`、`RbSi(HgS₂)₂`、`AgGeAu₂`、`CrNiRh₂`、`Dy₂GaNi`、`YbLuCoGeIr`、`Ce₃Pr₂Nd₃Se₈`、`Er₄H`、`HfSc₂Cd`、`ErGeRh₃` 等

每个结构都是带周期性边界（`pbc: T T T`）的合法晶体：晶格矢量 + 分数坐标 + 元素占位。

## 一、本地安装（CPU / GPU 通用）

MatterGen 要求 **Python 3.10**（pyg 扩展只有 cp310 wheel），推荐用 `uv`。

```bash
git clone https://github.com/microsoft/mattergen.git
cd mattergen
uv venv .venv --python 3.10
source .venv/bin/activate
uv pip install -e .
```

### 踩过的坑

1. **无 GPU 也能装能跑**：`pyproject.toml` 固定了 `torch==2.2.1+cu118` 及 cu118 的 pyg wheel，但这些 wheel 在无 GPU 机器上照样能 `import`，MatterGen 的 `get_device()` 会自动回退到 CPU。所以 CPU 环境无需改任何依赖。
2. **权重是 Git LFS 指针**：仓库里 `checkpoints/<model>/checkpoints/last.ckpt` 只是 134 字节的 LFS 指针（真实 ~461MB）。本机没装 git-lfs 也没关系——用 `--pretrained-name=<model>` 时会**自动从 HuggingFace 下载**权重。

## 二、CPU 运行 demo（无条件生成）

```bash
source .venv/bin/activate
mattergen-generate results/mattergen_base_cpu \
    --pretrained-name=mattergen_base \
    --batch_size=2 --num_batches=1
```

产出（`results/cpu/`）：
- `generated_crystals_cif.zip`：每个结构一个 `.cif`
- `generated_crystals.extxyz`：所有结构合并为多帧
- `generated_trajectories.zip`：完整去噪轨迹（体积较大，本目录未收录）

CPU 上 1000 步 ancestral sampling 约 9 分钟。想更快可减小采样步数（`--sampling_config_overrides="sampler_partial.N=250"`），但会牺牲质量。

## 三、Modal 云端 GPU 运行 demo

用 `mattergen_modal.py` 把安装 + 生成整个流程搬到 Modal 的 GPU 容器上，结果自动拉回本地。

```bash
pip install "modal[api-proxy-support]"   # 见下方“代理”说明
export MODAL_TOKEN_ID=...   MODAL_TOKEN_SECRET=...
modal run mattergen_modal.py --batch-size 16 --num-batches 1
# 条件生成示例（磁密度）：
modal run mattergen_modal.py --model-name dft_mag_density --batch-size 16 \
    --extra "--properties_to_condition_on={'dft_mag_density':0.15} --diffusion_guidance_factor=2.0"
```

脚本要点：
- 镜像 = `python3.10` + `git clone mattergen` + `uv pip install -e .`（全部在 Modal 云端完成，cu118 torch 正好匹配 GPU）。
- `@app.function(gpu="T4", timeout=3600)` 里跑 `mattergen-generate`，把 `generated_crystals_cif.zip` / `.extxyz` 的字节读回，`local_entrypoint` 保存到本地 `modal_results/<model>/`。
- 结果同时持久化到 Modal Volume `mattergen-results`。

产出见 `results/gpu/`（16 个结构）。

### 在受限网络（CONNECT 代理）下连 Modal

若本机出站流量走 HTTP CONNECT 代理（`HTTPS_PROXY` 已设），Modal 的 gRPC 客户端默认直连会被拒。Modal 支持读取 `HTTPS_PROXY` 走代理，但需要额外的包：

```bash
pip install "python-socks[asyncio]"       # 或 pip install "modal[api-proxy-support]"
export SSL_CERT_FILE=/path/to/ca-bundle.crt   # 让 gRPC TLS 信任代理 CA
```

装上后 `modal app list` / `modal run` 即可正常连通。

## 可用的预训练模型

`mattergen_base`（无条件基座）、`chemical_system`、`space_group`、`dft_mag_density`、`dft_band_gap`、`ml_bulk_modulus`、`dft_mag_density_hhi_score`、`chemical_system_energy_above_hull`、`mp_20_base`。把上面命令里的 `--pretrained-name` / `--model-name` 换掉即可做属性条件生成。
