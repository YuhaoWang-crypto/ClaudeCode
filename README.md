# FairChem 用于材料 / MOF 吸附 / 药物结晶预测

评估「用 [fairchem](https://huggingface.co/fairchem) 与 [facebook/OMAT24](https://huggingface.co/facebook/OMAT24)
模型做材料和药物结晶预测、以及 MOF 吸水/吸附检测」的可行性，并给出可直接使用的安装脚本与示例。

## 一句话结论

**技术路线可行**，但要用对模型；**在当前云环境里跑不通**——模型权重必须从
HuggingFace 下载，而本会话的组织出网策略把 `huggingface.co` 拦截了（403/407）。
请在**能访问 HuggingFace 的机器**上按本仓库的 `setup.sh` 安装运行。

## 两个关键前提（缺一不可）

1. **网络**：能访问 `huggingface.co`。当前 Claude Code 云会话被组织策略拦截
   （已验证 `huggingface.co`、`hf-mirror.com` 等均返回 403/407），代理规则明确
   要求不得绕过。PyPI 可用，所以**库能装、模型下不了**。
2. **许可证 + Token**：OMAT24 / UMA 等是 **gated 仓库**。需先在模型页
   （如 https://huggingface.co/facebook/UMA ）接受许可证，并用 HF token 登录。

## 选对模型（重要）

OMAT24 只是其中一个数据集/模型。FairChem v2 的 **UMA** 通用模型用一个权重、
通过 `task_name` 覆盖多个领域，正好对上你的三个目标：

| 你的目标 | 该用的任务 | 说明 |
|---|---|---|
| 无机材料预测（稳定性/形成能/力学） | `omat`（OMat24） | OMAT24 的主场，无机晶体 |
| **MOF 吸水 / 吸附** | `odac`（Open DAC） | 专为 MOF + CO2/H2O 吸附训练，比 OMAT24 合适 |
| **药物结晶 / 有机分子晶体** | `omc`（分子晶体）+ `omol`（孤立分子） | OMAT24 是无机训练，**不适合**有机药物晶体 |

> 直接拿 OMAT24 去做药物结晶或 MOF 吸附是选型错误：药物晶体是靠范德华力+氢键
> 的有机分子晶体，MOF 是有机-无机杂化并涉及客体分子吸附，两者都不在 OMAT24 的
> 训练分布内。

## 能力边界（别期望过高）

- **MOF 吸水量/等温线**：MLIP 给的是单点**吸附能**。要得到 water uptake、
  工作容量、等温线，还需 **GCMC**（如 RASPA）或热力学积分，MLIP 只是提供能量的一环。
- **药物结晶预测**：完整的晶体结构预测（CSP）= 生成候选结构 + 能量排序。
  本仓库示例只做第二步（用势给多晶型排序）；生成候选需另配 CSP 流程。弱相互作用
  对色散敏感，必要时叠加 D3 校正。
- **算力**：UMA 是等变 Transformer，**CPU 能跑但很慢**，大规模筛选需要 GPU。
  本云会话无 GPU。

## 快速开始

```bash
# 1) 在模型页接受许可证，拿到 token
export HF_TOKEN=hf_xxx
# 2) 安装（默认装 CPU 版 torch，GPU 见脚本注释）
bash setup.sh
# 3) 冒烟测试
python examples/01_material_relax_omat.py
```

## 示例

- `examples/01_material_relax_omat.py` — 无机材料能量/弛豫（`omat`）
- `examples/02_mof_adsorption_odac.py` — MOF 吸附能（`odac`）
- `examples/03_drug_molecular_crystal_omc.py` — 药物多晶型晶格能排序（`omc`/`omol`）

> 示例是**模板**，具体 API 以你安装的 `fairchem-core` 版本为准（基于 v2.x UMA 接口编写）。
> 由于本环境无法访问 HF/无 GPU，示例未在此实机验证过。

## 如果一定要在 Claude Code 云端跑

需要你的组织把 `huggingface.co` 及其 CDN 加入出网白名单（网络策略在环境创建时设定，
见 https://code.claude.com/docs/en/claude-code-on-the-web ）。否则请在本地/自有 GPU 机器运行。
