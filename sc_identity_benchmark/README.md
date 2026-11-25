# 单细胞身份 / 干性 / 衰老分型:开源模型 + 未见数据的实测

问题:方案里"参考映射 + 预训练权重"这条路线(scVI/scANVI、CellTypist、CytoTRACE 2、SenCID),
用公开数据预训练后,在**没参与训练**的数据上,能否给出它承诺的四类输出?

结论先行:**思路可行,但四项输出的可靠程度差别很大**,其中"未知群预警"不能用分类器置信度来做。

## 实验设计

| 角色 | 数据 | 说明 |
|---|---|---|
| 参考集(预训练) | CELLxGENE / HCA 骨髓图谱(Cell 2023),去掉 Tabula Sapiens 后 52,891 细胞、13 个供体、34 个细亚群 | scVI 30 轮 → scANVI 15 轮,4 核 CPU 约 10 分钟 |
| 测试集 A | Triana 2021,BD Rhapsody 全转录组,13,165 细胞 | 不同实验室、**不同平台** |
| 测试集 B | Tabula Sapiens 骨髓 10x,23,981 细胞 | 整个研究从参考集中剔除 |
| 衰老测试 | GSE115301(Teo 2019),IMR90 成纤维细胞,增殖 207 / RAS 诱导衰老 240 / 旁观者继发衰老 33 | 有实验真值 |

防泄漏措施:参考集与测试集无供体重叠;高变基因只在参考集上选;测试集标签从不输入模型;
浆细胞和 pDC 两个谱系在训练前从参考集**整体删除**,用来检验"未知群预警"。

## 结果

### 1. 细亚群注释 — 大谱系可用,细亚群不稳

| 已知谱系的粗分类 | 测试集 A 准确率 / macro-F1 | 测试集 B 准确率 / macro-F1 |
|---|---|---|
| scANVI + scArches | **0.872 / 0.794** | **0.652 / 0.631** |
| CellTypist(用同一参考集训练) | 0.790 / 0.638 | 0.501 / 0.539 |
| CellTypist 官方预训练权重(全部谱系) | 0.719 / 0.709 | 0.454 / 0.460 |

- T、B、NK、单核、红系/巨核谱系召回率 0.89–1.00。
- 失败集中在两处。其一是中性粒谱系(A 0.25,B 0.43):剔除 Tabula Sapiens 后参考集只剩 114 个成熟中性粒细胞,
  参考集覆盖不到的群体会被强行归到邻近类(单核)。测试集 B 去掉成熟中性粒后准确率为 0.787。
  其二是测试集 B 的 HSPC(0.20):Tabula Sapiens 的"造血前体细胞"与参考集的 HSC/MPP 定义不一致,这部分误差无法区分是模型错还是标签口径不同。
- 细亚群(测试集 A 中 16 个与参考集同名的标签):准确率 0.717,macro-F1 0.510。
  初始 B(0.92)、pro-B(0.94)、单核亚型(0.96–1.00)很好;CD8 记忆亚型(0.02–0.05)、γδ T(0.00)、MEP(0.14)基本失败。

### 2. 未知群预警 — 分类器置信度无效,潜空间距离部分有效

| 打分方式(AUROC:被删谱系 vs 已知谱系) | 测试集 A | 测试集 B |
|---|---|---|
| scANVI 1 − 最大概率 | 0.671 | 0.364 |
| CellTypist 1 − 最大概率 | 0.623 | 0.859 |
| **到参考集的 kNN 潜空间距离** | **0.947** | **0.815** |

被删除的浆细胞有 80%–94% 被**高置信度**地判成成熟 B 细胞,pDC 被判成 cDC。
用 kNN 距离、在 5% 假阳性率下,能召回 77%(A)和 34%(B)的未知细胞。即:预警必须用距离类指标,且仍会漏掉与已知类相近的新群。

### 3. 无速率分化潜能排序(CytoTRACE 2 预训练权重)— 效果最好

| | 测试集 A | 测试集 B |
|---|---|---|
| 与三级潜能真值的 Spearman | 0.816 | 0.690 |
| AUROC 多能干/祖细胞 vs 成熟细胞 | 0.996 | 0.979 |
| AUROC 多能 vs 定向祖细胞 | 0.809 | 0.655 |
| 基线(检出基因数)Spearman | 0.644 | 0.413 |

干/祖细胞与成熟细胞几乎完全分开;HSC/MPP 与定向祖细胞之间的细排序较弱。每个细胞类型抽样 ≤300 细胞。

### 4. 衰老候选打分(SenCID 预训练 SVM)— 排序有效,阈值未校准

RAS 诱导衰老 vs 增殖的 AUROC:SID3 0.889、SID2 0.888、SID4 0.838、SID1 0.761;单基因基线 CDKN1A 0.731、CDKN2A 0.715。
但按官方 0.5 阈值,SID4 把 64% 的增殖细胞判为衰老,SID1 为 42%,SID3 为 18%。只能当候选排序,不能直接当阳性判定。

## 未做到 / 需要知道的限制

- **scGPT / Geneformer 没有测**:本环境只有 4 核 CPU、无 GPU。
- **SenCID 不是官方完整流程**:未做 DCA 去噪(依赖 TensorFlow < 2.5);仓库未附 SID 推荐模型和 SID5/6 权重,所以四个 SID 全部报告而非自动选型。
- **预训练权重的训练集无法完全核实**:CellTypist Immune_All_Low、CytoTRACE 2、SenCID 的训练数据可能包含本测试集或同一研究的数据;只有 scANVI 和自训 CellTypist 是严格无泄漏的。
- 只测了健康骨髓和一个体外成纤维细胞衰老模型;疾病样本、实体组织、体内衰老未验证。
- "真值"是原作者注释,自身带噪声;单次运行、单一随机种子,没有置信区间。
- 为绕开 cytotrace2-py 与新版 pandas 的一处不兼容,改了已安装包 `common/gen_utils.py` 中一行(把 `i is not np.nan` 换成 `isinstance(i, str)`)。

## 复现

```bash
pip install scanpy scvi-tools celltypist cytotrace2-py "setuptools<81" matplotlib scikit-learn
git clone --depth 1 https://github.com/JackieHanLab/SenCID ext/sencid
# 按 common.py 中 DATASETS 的链接下载三个 h5ad 到 $SCB_DATA,GSE115301 的两个 Growing_Sen_10x 文件到 $SEN_DATA
export SCB_DATA=... SCB_WORK=... SENCID_DIR=ext/sencid/SenCID SEN_DATA=...
for s in 01_prepare 02_train_map 03_celltypist 04_cytotrace2 05_sencid 06_evaluate; do python $s.py; done
```

全部数字在 `results/summary.json`,图在 `results/fig1–fig4`。CytoTRACE 2 仅限非商业使用,商用须另行取得许可。
