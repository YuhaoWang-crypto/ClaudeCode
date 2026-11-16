# mRNA-LM DEMO

学习并复现 Sanofi 的 **[mRNA-LM](https://github.com/Sanofi-Public/mRNA-LM)**:每个 mRNA 功能区
(`5'UTR — CDS — 3'UTR`)各用一个 BERT 编码,再融合预测分子性质(半衰期、表达量等)。

本目录提供三样东西:

| 文件 | 说明 |
|------|------|
| `mrna_lm_explainer.html` | **交互式网页图解**:可动手的分词器 + 三塔结构图 + 实测结果。直接浏览器打开即可 |
| `mrna_lm_demo.py` | **忠实复现、CPU 可跑**:三塔结构在真实半衰期数据上从零训练并评估 |
| `codonbert_real_demo.py` | 加载 Sanofi **真·CodonBERT** 预训练权重,跑掩码预测 + 嵌入提取 |

> 网页图解无需任何依赖,双击 `mrna_lm_explainer.html` 即可查看。下面两个 Python demo 才需要环境和数据。

---

## 架构一图流

```
  5'UTR ──[5UTRBERT]──┐  (单核苷酸分词)
   CDS  ──[CodonBERT]──┼─ mean-pool → concat(768×3) ─[融合头]→ 预测值
  3'UTR ──[3UTRBERT]──┘  (CDS 按密码子分词)
```

- **分词**:5'/3'UTR 按单核苷酸,CDS 按密码子(3-mer);先 T→U(同原仓库 `dataload.py`)
- **融合**:三区域 mean-pool 向量拼接 → Dense+GELU+LayerNorm → 回归/分类头(同 `FullModel.py`)
- **微调**:LoRA 适配器冻结主干;可选 CLIP 对比损失对齐 UTR↔CDS 表示

---

## 环境

```bash
python -m venv venv && source venv/bin/activate
pip install torch transformers scikit-learn scipy pandas numpy
```

## 准备数据(两个 Python demo 需要)

数据和权重**体积大,不入库**,按需拉取到本目录:

```bash
# 1) 真实数据:从原仓库取 CSV 放到 ./data/
git clone --depth 1 https://github.com/Sanofi-Public/mRNA-LM.git /tmp/mRNA-LM
mkdir -p data && cp /tmp/mRNA-LM/data/mrna_half-life.csv data/

# 2) 真·CodonBERT 权重(约 911MB)放到 ./weights/codonbert/
mkdir -p weights && cd weights
curl -O "https://cdn.prod.accelerator.sanofi/llm/CodonBERT.zip" && unzip CodonBERT.zip
cd ..
```

---

## 运行

### ① 三塔模型预测半衰期(约 2–3 分钟,CPU)

```bash
python mrna_lm_demo.py                 # 默认读 data/mrna_half-life.csv
python mrna_lm_demo.py --csv path.csv --epochs 4 --train_n 4000
```

实测(CPU):测试集 Spearman 从随机的 **+0.009 → +0.348**。证明"三区域分词 + 三塔融合"能从序列读出半衰期信号。
(小型从零训练,非论文 SOTA;论文用大预训练权重 + LoRA 相关性更高。)

为 CPU 可跑做了简化:塔小、随机初始化、序列截断(5'UTR≤100nt、CDS≤200密码子、3'UTR≤200nt)、训练集下采样。

### ② 真·CodonBERT 推理

```bash
python codonbert_real_demo.py                       # 默认 weights/codonbert + data/mrna_half-life.csv
python codonbert_real_demo.py --ckpt DIR --data CSV
```

- **掩码密码子预测**:60 条基因汇总 top-1 **15.5%** / top-5 **38.3%**(随机仅 1.6%),~10 倍于随机;top-5 常是真实或**同义**密码子(例:真实 `UUA`,top-5 里 `CUU/UUU/CUG/UUG/UUA` 四个都是亮氨酸)→ 确实学到了遗传密码。
- **嵌入提取**:演示喂给融合头的 768 维 CDS 向量。

> 注:发布的 CodonBERT 词表为 69(5 特殊符 + 64 密码子,基于 AUGC),依据原版
> `benchmarks/utils/tokenizer.py`,而非 mRNA-LM 仓库 `FullModel.build_tokenizer` 里构造的 125 密码子表。
> 本脚本已按正确词表加载。

---

## 四个下游任务(原仓库)

| task | 含义 | 类型 | 数据 |
|------|------|------|------|
| `halflife` | mRNA 半衰期 | 回归 | `mrna_half-life.csv` |
| `tr` | 翻译速率 | 回归 | `translation_rate.csv` |
| `liver` | 肝脏转录本表达 | 回归 | `transcript_expression_liver.csv` |
| `5class` | 蛋白表达分级 | 5 分类 | `protein_expression_5class.csv` |

原仓库训练入口:`python finetune_all.py --task halflife --output out/`(需 GPU + 预训练权重)。

---

*复现自开源项目 Sanofi-Public / mRNA-LM,用于学习与演示。*
