# 材料与药物结晶预测 (fairchem / UMA / OMat24)

用 Meta FAIR Chemistry 的机器学习原子间势 (MLIP) 做**材料**、**MOF 吸附**、
**药物分子晶体 / 多晶型**的能量与结构预测。

## 结论：可行，而且非常匹配

你提到的两个链接背后是同一套东西：

- **[fairchem](https://github.com/facebookresearch/fairchem)** —— FAIR Chemistry 的库，
  现在主推 **UMA (Universal Model for Atoms)** 通用模型。
- **[facebook/OMAT24](https://huggingface.co/facebook/OMAT24)** —— 无机材料数据集
  (1.1 亿条 DFT 计算) 及其 EquiformerV2 模型；现已被 UMA 覆盖并统一。

UMA 是**一个网络 + 多个任务头 (task)**，每个 task 对应一套 DFT 训练数据。
你的三类需求正好一一对应到官方 task，不需要自己训练：

| 你的需求 | 对应 task | 训练数据 |
|---|---|---|
| MOF 吸水 / CO₂ 吸附能检测 | **`odac`** | Open Direct Air Capture（MOF + H₂O/CO₂）|
| 药物结晶 / 多晶型稳定性排序 | **`omc`** | Open Molecular Crystals 2025（有机分子晶体）|
| 无机材料筛选 / 结构弛豫 | **`omat`** | OMat24 |
| 单分子 / 聚合物 | `omol` | OMol25 |
| 催化 (吸附在金属表面) | `oc20` | Open Catalyst |

- **MOF 吸附**：`odac` 就是为「MOF 加/不加 H₂O、CO₂ 的能量变化」而训练的，
  可直接算吸附能 → 筛选吸水性 / 直接空气捕集效果。
- **药物结晶**：`omc` + FastCSP 工作流 (arXiv:2508.02641) 已能做晶体结构预测，
  实验多晶型通常被排到全局最低或前 10，精度接近含色散校正的 DFT。

### 注意点
1. **模型是 gated 的**：需要 HuggingFace 账号，先到
   <https://huggingface.co/facebook/UMA>（及 `facebook/OMAT24`）**接受许可**，
   再 `huggingface-cli login` 或设 `HF_TOKEN`。许可为**研究用途**。
2. **GPU**：本环境是 **CPU-only**。小体系可跑（慢），几百原子以上的 MOF / 晶胞
   建议上 GPU。
3. **药物结晶最匹配 `omc`**（有机分子晶体），不要用 `omat`（无机）。真正的多晶型
   预测 (CSP) 是一整套采样流程，这里给的是「候选结构弛豫 + 能量排序」这一核心步骤。

## 安装

```bash
python -m pip install --upgrade pip setuptools wheel   # 旧 setuptools 会导致 antlr4 依赖编译失败
python -m pip install -r requirements.txt
```

已在本仓库 CPU 环境**实测安装并导入成功**：`fairchem-core 2.21.0` + `torch 2.8.0`
（CPU 模式）+ `ase 3.29`，UMA API (`FAIRChemCalculator`, `pretrained_mlip`) 可正常导入，
三个示例脚本均可正常解析运行（模型权重需你自备 HF token 才能下载，见下）。

安装注意（本仓库 CPU 环境踩坑）：
- 装依赖前**先升级 `setuptools`/`wheel`**，否则 `antlr4-python3-runtime`（fairchem 依赖）
  会以 `AttributeError: install_layout` 编译失败。
- **不要用** `--index-url https://download.pytorch.org/whl/cpu` 装 torch：该主机在本
  受限网络里被代理拦截 (403)。直接从 PyPI `pip install torch` 即可（PyPI 已放行）。
- 确认 `pip` 与你运行脚本的 `python` 是**同一个解释器**（本环境两者曾不一致，
  用 `python -m pip ...` 最保险）。

## 在 Claude Code cloud 上运行（重要）

云端容器默认是 **Trusted** 网络级别，只放行了包管理器（pip 能用），**huggingface.co
被代理拦截 (403)**。所以下载 gated 权重前，必须在环境设置里做两件事，然后**新开一个
session** 生效：

1. **加 HF token（作为环境变量，别贴进聊天）**
   编辑环境 → Environment variables（`.env` 格式，**不要加引号**）：
   ```
   HF_TOKEN=hf_你的token
   ```
   `huggingface_hub` 会自动读取 `HF_TOKEN`。注意：云端目前**无专用 secrets 保险库**，
   环境变量对可编辑该环境的人可见。

2. **把 HuggingFace 加入网络白名单**
   编辑环境 → Network access 选 **Custom** → Allowed domains 每行一个，并勾选
   "Also include default list of common package managers"：
   ```
   huggingface.co
   *.huggingface.co
   *.hf.co
   ```

3. **接受模型许可**：用同一 HF 账号在 <https://huggingface.co/facebook/UMA>
   （及 `facebook/OMAT24`）点接受（研究用途）。

三步齐了、开新 session 后，`python examples/mof_adsorption.py ...` 会自动下载权重并运行。

## 快速检查（无需下载模型）

```bash
python examples/check_install.py
```

打印依赖版本、是否有 GPU、HF token 是否就绪。

## 用法

```bash
# 1) MOF 吸附能（吸水 / CO2）
python examples/mof_adsorption.py --mof your_mof.cif --guest H2O --relax
python examples/mof_adsorption.py --mof your_mof.cif --guest CO2 --relax

# 2) 药物 / 有机分子晶体多晶型排序
python examples/drug_crystal.py --cif form_I.cif form_II.cif form_III.cif --relax-cell

# 3) 无机材料弛豫（OMat24）
python examples/material_relax.py --cif structure.cif --relax-cell
```

需要你自备结构文件（CIF）：MOF 可取自 CoRE MOF / QMOF 数据库，
药物多晶型可取自 CSD / 文献。

## 目录

```
requirements.txt                依赖
src/fairchem_predict/loader.py  加载 UMA、按 task 建 ASE calculator
examples/check_install.py       环境自检（不下载模型）
examples/mof_adsorption.py      MOF 吸附能 (odac)
examples/drug_crystal.py        分子晶体 / 多晶型排序 (omc)
examples/material_relax.py      无机材料弛豫 (omat)
```

## 参考

- fairchem: <https://github.com/facebookresearch/fairchem>
- UMA 论文: <https://arxiv.org/abs/2506.23971>
- OMat24: <https://arxiv.org/abs/2410.12771>
- FastCSP（分子晶体预测）: <https://arxiv.org/abs/2508.02641>
