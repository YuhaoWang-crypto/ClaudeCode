# Sanofi-Public 三个开源库 —— 安装 + Demo 结论

**结论：三个都可以学习、安装、并成功跑出 demo。** 全部在本机（CPU，4 核，无 GPU）实测通过。

环境根目录：`/tmp/claude-0/-home-user-ClaudeCode/7c540e2c-2fab-54d9-9da4-4e45967b1e8b/scratchpad`

| 库 | 做什么 | 环境 | Demo 结果 |
|---|---|---|---|
| **SNAC-DB** | 从 PDB/CIF 结构文件抽取、整理抗体/纳米抗体-抗原复合物 | conda(micromamba) `snacdb` + ANARCI/HMMER/Foldseek | ✅ 单元测试 3/3 过；7 个真实 PDB → 6 个复合物；Foldseek 聚类正确 |
| **RiboNN** | 从 mRNA 序列预测翻译效率(TE) | venv `ribonn` (torch1.13 CPU, PL1.8.5) | ✅ 30 条转录本 × 78 个人类细胞/组织 TE 预测 |
| **CodonBERT** | mRNA 语言模型，密码子优化/表达预测 | venv `codonbert` (torch1.13, transformers4.29) | ✅ 87M 参数模型加载；mRFP 表达基准 Pearson 0.66 |

---

## 1. SNAC-DB （抗体/纳米抗体复合物数据库流水线）

- 安装：`micromamba create -n snacdb python=3.10 anarci hmmer foldseek`，再 `pip install -r requirements.txt .`，并运行官方补丁 `snacdb-patch-anarci`（让 ANARCI 容忍未知残基 `X`）。
- 单元测试：`pytest unit_tests/` → **3 passed**。
- Demo：新下载 7 个 PDB（1MLC 溶菌酶、1N8Z HER2/trastuzumab、5JDS PD-L1 纳米抗体、6WAQ SARS-CoV-1 RBD 纳米抗体、7BWJ SARS-CoV-2 抗体、7KGK 合成纳米抗体、6M0J ACE2-RBD 无抗体），跑 `data_curation.sh`：
  - 正确识别出 **6 个 Ab/Nb-抗原复合物**（6M0J 无抗体被正确排除），自动标注 VH/VL/VHH 链、抗原链、相互作用强度。
  - Foldseek `easy-multimercluster` 按抗原结构聚类，**把三个冠状病毒 RBD 结合体（6WAQ / 7BWJ / 7KGK）正确归为同一簇**，溶菌酶/HER2/PD-L1 各自单独成簇 —— 生物学上合理。
- 产物：`snac_demo/pdbs_curated/`（6 个整理好的复合物 PDB+npy）、`pdbs_curation_summary.csv`、`pdbs_curated_antigen_all_cluster.tsv`。

## 2. RiboNN （翻译效率预测，Nature Biotech 2025）

- 安装：pytorch 1.13.1(cpu) + pytorch-lightning 1.8.5 等固定版本（README 的 mamba 环境换成 uv venv 装等价包）。权重从 Zenodo 17258709 自动下载（`weights.zip`，含 90 个 human + 90 个 mouse 嵌套交叉验证模型）。
- Demo：`python -m src.main --predict human`，对 `data/prediction_input1.txt` 的示例转录本预测。
  - 输出 `results/human/prediction_output.txt`：**30 条转录本 × 78 个细胞系/组织**的 TE 预测 + 跨组织平均 TE。耗时 ~75s（CPU）。
  - 预测值分布合理（mean TE 从 -0.15 到 1.32，可排序），模型区分开了高/低翻译效率转录本。

## 3. CodonBERT （mRNA 大语言模型）

- 安装：`poetry` 的依赖较老，改用 uv venv 装 torch1.13 + transformers4.29 + peft/accelerate 等。模型（BertForPreTraining, 69 词表 = 64 密码子 + 5 特殊符, 12 层, 87M 参数）从 Sanofi CDN 下载解压。
- Demo（冻结特征探针，自写 `codonbert_demo.py`）：用 CodonBERT 对论文基准数据 `mRFP_Expression.csv`（1459 条，官方 train/test split）提取嵌入，训练 ridge 回归预测蛋白表达：

  | 特征 | Pearson | Spearman | R² |
  |---|---|---|---|
  | 密码子频率基线 (64d) | 0.384 | 0.513 | 0.139 |
  | CodonBERT CLS (768d) | 0.658 | 0.662 | 0.430 |
  | CodonBERT mean-pool (768d) | **0.661** | **0.713** | **0.430** |

  冻结的预训练嵌入把相关性从 0.38 提到 0.66（接近翻倍），说明模型确实学到了与表达相关的序列信号。论文里带微调的 CodonBERT 会更高，这里只做轻量 CPU 探针。

---

## 复现命令速查

本目录（`sanofi-demos/`）里每个库都有一键脚本，会自动克隆上游仓库、建独立环境、
下载权重并运行。产物默认落在各自的 `work/`（已 gitignore）。

```bash
# SNAC-DB  —— 需要 micromamba/conda（ANARCI/HMMER/Foldseek 只有 conda 装）
bash snac-db/run_snacdb_demo.sh
# RiboNN   —— 需要 uv（或 python3.10），从 Zenodo 下 ~200MB 权重
bash ribonn/run_ribonn_demo.sh
# CodonBERT —— 需要 uv（或 python3.10），从 Sanofi CDN 下 ~950MB 模型
bash codonbert/run_codonbert_demo.sh
```

`outputs/` 里收录了几个小体积示例输出作为运行证据（SNAC-DB 整理汇总 + 抗原聚类、
RiboNN 预测结果）。上游仓库、权重、环境、大体积产物不入库，由脚本重新生成。
