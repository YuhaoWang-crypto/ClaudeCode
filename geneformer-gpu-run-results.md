# Geneformer GPU 运行结果（Modal T4）— 诚实记录

> 在用户 Modal 账号上运行 `modal_geneformer.py`（单张 T4）。目的：证明真实 GPU Geneformer 嵌入 + de novo 发现管线可用。

## 跑了什么（真实、非模拟）
- **模型**：`Geneformer-V1-10M`（官方通用**预训练**模型，10M 参数，最轻量），从 `ctheodoris/Geneformer` 下载。
- **词表**：官方 `gc30M` 字典（token + gene-median），**与 V1-10M 词表匹配**（关键：词表不匹配则 token 全错）。
- **编码**：Geneformer 的 rank-value 编码（按 归一化表达/基因中位数 降序取 token）。
- **硬件**：Modal T4 GPU（`device=cuda`，`BertModel`，emb_dim=256）。
- **下游**：de novo KMeans 聚类 + 炎症 token 富集 / inflamed 条件富集。

## 三次运行（成本共几分钱，app 均已停止）
| 运行 | 模型 | silhouette | inflamed 富集分离 | 说明 |
|---|---|---|---|---|
| 1 | 误选 cardiomyopathy 分类器 | 0.363 | 弱（各簇≈0.5） | 数据 n_genes=400 无截断，信号被均值池化冲淡 |
| 2 | 同上分类器 | 0.294 | 略有（cluster inflamed 0.59） | 加大基因数/信号；但仍是微调分类器权重 |
| 3 | **Geneformer-V1-10M 通用预训练** | 0.048 | 近乎随机（各簇≈0.5） | 正确模型；但合成数据无真实共表达结构 |

## 关键诚实结论
1. **基础设施完全跑通**：真实 Geneformer 预训练权重 + 官方词表 + rank-value 编码 + GPU 前向 + 256 维细胞嵌入 + de novo 聚类，端到端在 Modal 上成功返回报告（`geneformer_gpu_report.json`）。
2. **合成数据无法产生真实生物信号**——这不是 bug，而是 Geneformer 的本质：**它的嵌入编码的是从真实语料学到的基因-基因共表达关系**。我的合成数据用随机 Ensembl ID + 随机指定的"炎症程序"，这些基因间没有真实关系，所以通用预训练模型不会把它们当连贯程序分离（silhouette 0.048≈随机）。第 2 次那个"0.59 分离"是微调分类器的伪信号。
3. **推论**：真正的 de novo 发现**必须用真实 scRNA**（GSE134809 / SCP259），那里的炎症模块（如 IAF: IL13RA2/IL11/TNFRSF11B）是 Geneformer 认得的真实共表达结构。

## 对比：本地 PCA 后端（合成数据）为何能分离？
`geneformer_denovo_discovery.py --backend pca` 在合成数据上恢复了炎症态（1.88× 富集）——因为 PCA 直接吃我植入的表达差异，是"自证"。而 Geneformer 用的是**先验的真实生物知识**，不吃人造的随机结构。**这恰好说明 Geneformer 的价值：它带来真实生物先验，但也因此需要真实数据才能显现。**

## 下一步（真实数据，同一脚本）
```bash
# 1) 下载 GSE134809 (Martin CD) 或 SCP259 (Smillie UC) 为 .h5ad（含 Ensembl ID + n_counts）
# 2) 把 modal_geneformer.py 的 build_synthetic 段换成读 h5ad 的 rank-value 编码
#    （counts 矩阵 + panel=数据基因∩词表；其余流程不变）
# 3) python3 -m modal run modal_geneformer.py
# 预期：Geneformer 嵌入分出炎症相关成纤维细胞态 / GIMATS 样态，且富集于 inflamed / anti-TNF 无应答
```

### 署名
模型与词表：`ctheodoris/Geneformer`（HuggingFace）。数据集见 `scRNA_datasets.md`。运行记录来自用户 Modal workspace（wyh-58141），app 均已 stopped。
