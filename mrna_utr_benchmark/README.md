# 小核酸与 mRNA / UTR 构建设计:公开数据实测

检验矩阵这一行的承诺与边界。
这一行里 UTR-LM/Optimus 负责"5′UTR→平均核糖体负载(MRL)/翻译效率",RiboNN 负责全长 mRNA TE,
Saluki 负责半衰期,LinearDesign 负责 CDS 优化,ViennaRNA 负责结构/可及性。
第 05 栏承诺"零样本/基线直接输出 TE 与稳定性排序";第 06 栏边界:**"翻译效率(TE)不等于累计蛋白产量或持续表达曲线;单构建数据无法微调出客户专用序列规律"**。

本机 CPU 可跑 Optimus 5-Prime 预训练 CNN(Sample 2019)和 ViennaRNA。这是 5′UTR→MRL 领域的标准公开基准,结论可直接检验。

结论先行:**这一行是少数"基本兑现"的——预训练 TE 模型在随机文库 r²=0.93,泛化到从没见过的真实人源 5′UTR 仍有 r²=0.75–0.80,并正确学到 uAUG 降低翻译的机制**。
但有明确边界:它预测的是核糖体负载(TE),不是蛋白产量/持续表达;它是 50nt、HEK293-eGFP 语境专用,换细胞/换长度/换读出需要重训。

## 数据(全部公开)

| 项目 | 数据 | 规模 |
|---|---|---|
| 模型 | Optimus 5-Prime CNN 预训练权重(Sample 2019,pjsample/human_5utr_modeling) | 零样本,未重训 |
| 随机文库真值 | GEO GSE114002 GSM3130435,随机 50nt 5′UTR 多聚体分析 MRL | 测试集 2 万(模型训练用其余 26 万) |
| 内源真值 | 同系列 GSM3130443 设计文库中的 **真实人源 5′UTR** | 41,862 条,模型从未见过 |
| 基线 | ViennaRNA 结构(MFE、3′端可及性)+ uAUG + GC + 碱基组成 → 梯度提升 | 随机文库 6 万训练 |

## 结果

### 1. TE 预测:兑现(承诺成立)

| 数据 | Optimus CNN r² | ViennaRNA 基线 r² |
|---|---|---|
| 随机文库·留出测试(2 万) | **0.934** | 0.542 |
| 真实人源 5′UTR(41,862) | **0.746** | 0.561 |
| 真实人源·高测序深度子集(2.1 万) | **0.804** | — |

- 预训练模型在它设计的任务上几乎完美(r²=0.93,Spearman 0.95),远超简单生物物理基线(0.54)。
- 泛化到从未见过的真实人源 5′UTR,r² 从 0.93 降到 0.75(高质量子集 0.80):有下降,但仍可用——比前面几行(祖源/细胞外推)稳健得多。
- 对设计/变异序列也成立:sergii 设计 UTR r²=0.87,SNV 变异 r²=0.65。可用于设计排序。

### 2. 机制自洽(已知答案)

上游 AUG(uAUG)会抑制翻译。实测:无 uAUG 的 UTR 平均 MRL 7.7,有 uAUG 降到 4.86;人源 UTR 6.07→4.27。
Optimus 的预测同向(无 uAUG 0.72 → 有 uAUG −1.0)。模型学到了核心调控语法,不是在拟合噪声。

### 3. 边界(第 06 栏)

- **TE ≠ 蛋白产量/持续表达**:本基准的真值是核糖体负载(翻译效率代理)。累计蛋白量、表达持续时间受 mRNA 半衰期、降解、免疫原性等影响,不在该模型范围内——这正是矩阵标注的边界,本次未直接检验(需 Saluki/半衰期数据与蛋白读出)。
- **语境专用**:模型是 50nt 固定长度、HEK293 + eGFP 报告系统训练的。换细胞系、换 UTR 长度、换修饰核苷酸,序列-TE 规律会变,需要重训;"单构建客户数据"不足以微调出客户专用规律(与第 06 栏一致)。

## 对矩阵这一行的修改建议

1. 第 05 栏"TE 排序"可按兑现度保留,但注明:交付的是翻译效率排序(r² 内源约 0.75),不是蛋白产量预测;稳定性/半衰期是另一套模型(Saluki),需单独验证。
2. 注明模型的语境(长度、细胞系、报告系统);跨语境使用须重训或至少校准。
3. 第 06 栏的边界完全正确,保留。

## 限制

- **只测了 5′UTR→MRL 一个任务一个模型**:RiboNN(全长 mRNA TE)、Saluki(半衰期)、LinearDesign(CDS/MFE-CAI 联合优化)、小核酸(siRNA/ASO 效力与离靶)均未测;"稳定性排序、小核酸效力初筛"这部分承诺未检验。
- Optimus 的输出是训练集标准化后的 MRL,这里用 Pearson r²(=相关系数平方)和 Spearman 作尺度无关比较,未做绝对标定。
- 真实人源 UTR 被截断/补齐到 50nt(与模型输入一致),长 UTR 的远端结构信息丢失,可能低估也可能高估真实难度。
- 没有独立的蛋白产量读出,无法量化"TE vs 蛋白"的偏差,只能作为限制陈述。

## 复现

```bash
# svenv: numpy<2 + tensorflow-cpu==2.15.1 (加载 Keras 预训练模型);主 venv: ViennaRNA scikit-learn
git clone https://github.com/pjsample/human_5utr_modeling repo   # 含 saved_models/main_MRL_model.hdf5
# GEO: samples/GSM3130nnn/GSM3130435/suppl/GSM3130435_egfp_unmod_1.csv.gz(随机文库)
#      samples/GSM3130nnn/GSM3130443/suppl/GSM3130443_designed_library.csv.gz(含 human_utrs)
export NB_DATA=...
svenv/bin/python 01_cnn_predict.py && python 02_vienna_baseline.py && python 03_eval.py
```
数字在 `results/summary.json`,图在 `results/fig1_mrl.png`。
