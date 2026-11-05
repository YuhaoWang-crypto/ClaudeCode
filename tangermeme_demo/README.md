# tangermeme 安装与 Demo

[tangermeme](https://github.com/jmschrei/tangermeme)（[文档](https://tangermeme.readthedocs.io/en/latest/)）是一个对**已训练好的基因组序列深度学习模型**（输入 DNA one-hot，输出活性/信号）做事后分析的工具包：批量预测、DeepLIFT/SHAP 归因、饱和突变 (ISM)、motif 插入/间距实验、变异效应、seqlet 调用、motif 注释、序列设计等。它不训练模型，任何 PyTorch 模型都能接进来。

## 安装（已在本环境验证）

```bash
pip install torch --index-url https://download.pytorch.org/whl/cpu   # 有 GPU 就装 CUDA 版
```

```bash
pip install tangermeme logomaker
```

验证环境：Python 3.13，torch 2.15.0+cpu，tangermeme 1.6.0（会顺带装上 memelite、numba、pandas 等）。纯 CPU 即可运行。

## Demo 思路

真实的预训练模型（如 BPNet/ChromBPNet）需要下载大文件和参考基因组，所以这里用一个**答案已知**的合成任务：

1. 自己写一套"调控语法"，在 200 bp 随机 DNA 中埋入 motif，按规则算出活性：
   | motif | 序列 | 作用 |
   |---|---|---|
   | GATA | `GATAAG` | 激活 +1.0 |
   | EBOX | `CACGTG` | 激活 +1.5 |
   | AP1 | `TGACTCA` | 抑制 −1.0 |
   | GATA→EBOX 协同 | 间隔 10 bp 时最强 | 额外 +2.0（高斯形，σ=2 bp） |
2. 训练一个小 CNN（3 层卷积，约 7.7 万参数）去拟合它。
3. 只把模型交给 tangermeme，看它能否把上面这张表**反推出来**。

## 运行

```bash
python3 tangermeme_demo/demo.py                               # 从头训练 + 全部分析，约 5 分钟（4 核 CPU）
```

```bash
CKPT=tangermeme_demo/model.pt python3 tangermeme_demo/demo.py # 用已存权重跳过训练，约 30 秒
```

输出：`figures/*.png`、`results.json`、完整日志 `run_log.txt`。

## 结果

| 步骤 | tangermeme 函数 | 真值 | 恢复结果 |
|---|---|---|---|
| 1 预测 | `predict` | — | 测试集 Pearson r = **0.995** |
| 2 单 motif 效应 | `marginalize` | GATA +1.0 / EBOX +1.5 / AP1 −1.0 / 对照 0 | **+1.007 / +1.509 / −0.989 / +0.006** |
| 3 间距协同 | `space` | 最佳间隔 10 bp，峰值 +2.0 | 最佳间隔 **10 bp**，峰值 +1.71 |
| 4 归因 + seqlet | `deep_lift_shap` → `recursive_seqlets` → `annotate_seqlets` | 1416 个真实 motif 位点 | 召回 **98.4%**；正 seqlet 为 GATA/EBOX，负 seqlet 全是 AP1 |
| 5 饱和突变 | `saturation_mutagenesis` | — | 与 DeepLIFT 的相关 r = 0.93 |
| 6 变异效应 | `substitution_effect` | 破坏 GATA ≈ −1.0 | GATA 内 18 个 SNV 平均 **−1.03** |
| 7 序列设计 | `greedy_substitution` | — | 预测值 0 → 12；自发放了两组**间隔正好 10 bp 的 GATA→EBOX** |

### 图

- `figures/1_predict.png`：预测 vs 真值
- `figures/2_marginalize.png`：各 motif 单独效应 vs 真值
- `figures/3_space.png`：GATA–EBOX 间距扫描，红线（模型）与黑虚线（真值）几乎重合
- `figures/4_attribution_ism.png`：一条测试序列的 DeepLIFT/SHAP logo（带自动标注的 seqlet）和 ISM logo
- `figures/5_design.png`：设计出的序列及其归因

## 注意事项

- **tangermeme 反映的是模型学到的东西，不是生物学真相。** 第一版模型（感受野较稀疏）把协同峰抹平了，`space` 只测到约 1.0 的协同。查测试集残差后确认是模型本身在 10 bp 附近系统性低估约 0.9，tangermeme 的测量没错。换成更密的卷积并增加成对样本后，恢复到 1.71。
- **设计出的序列超出了训练分布。** 模型预测 12.0，但按真实规则算只有 10.0：真实规则里协同奖励只算一次，模型却把两组协同都算上了。设计结果需要用实验或独立模型复核。
- 只有 3 个 motif 时，TOMTOM 会警告 p 值依赖于 motif 集合大小（少于 25 个）。这里只用它来打标签，不解读 p 值。真实数据应使用 JASPAR/HOCOMOCO 等完整数据库。
- 接真实模型时：多任务模型要显式传 `target=`，`deep_lift_shap` 要固定 `random_state`，`greedy_substitution` 必须给正的 `max_iter`（默认 −1 会一轮都不跑）。tangermeme 包里自带一份说明（`tangermeme/_skills/data/SKILL.md`），列出了这些坑。
