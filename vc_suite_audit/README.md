# VC 套件模型可用性审计

对 BioVenic / Creative Diagnostics 干实验套件（`vc_antiviral_suite`）中已上线 QSAR 端点的
独立复算与外部审计。目标是回答两个问题：这些模型现在能不能用于药物筛选服务的指导与快速筛选，
以及应该用什么测试集给它们定级。

发布版报告：<https://claude.ai/code/artifact/90b78eab-7814-4b7a-8398-fb92c912bf11>

## 做了什么

审计对象的原始报告只给出**骨架划分交叉验证**指标。骨架划分对类似物系列几乎不设防，
因此本审计从 ChEMBL 独立重建六个端点的数据集，并在四种严苛度递增的划分下重训模型：

| 划分 | 含义 |
|---|---|
| `random` | 随机 5 折，乐观上界 |
| `scaffold` | Bemis-Murcko 骨架分组 5 折（原报告使用的协议） |
| `cluster0.4` | 球面排除聚类分组 5 折，簇心间 Tanimoto < 0.4 |
| `temporal` | 按文献年份，早期 70% 训练、后期 30% 测试；最接近前瞻筛选 |

同时复核了原报告的 q-RASAR 相似性特征增益、适用域边界曲线，以及跨靶点/跨菌种迁移能力。

## 复算的端点

NA 流感神经氨酸酶 pIC50、PA 核酸内切酶 pIC50、泛细胞毒性 pCC50、
流感嗜血杆菌 pMIC、阴沟肠杆菌 pMIC、hERG 抑制 pIC50。

数据溯源核对通过：六个端点的化合物数、骨架数与基线骨架划分 Spearman
均落在原报告声明值的几个百分点内。

## 主要结论

- 抗菌 MIC 是唯一接近可投产的一支：聚类划分下前 5% 富集仍有 4.8–10.9 倍。
- PA 靶点模型在新化学型上前 5% 富集为 0，R² 转负，不具备筛选能力。
- 六个端点一致显示：最近邻 Tanimoto 0.4 以下基本没有排序能力。
- q-RASAR 增益只在最大的数据集（泛细胞毒性）上复现，且只有声称量级的四成。
- MIC 绝对值不可用：落在一个倍比稀释度内的比例最高只有 54%。

## 运行

```bash
pip install numpy scipy pandas scikit-learn rdkit
cd scripts
python fetch.py            # 从 ChEMBL REST API 取活性数据
python curate.py           # 策展为建模表
python run_eval.py         # 四级划分 × 基线/q-RASAR
python run_cc50.py         # 泛细胞毒性（大数据集，单独跑，支持断点续算）
python transfer.py         # 跨端点迁移
python utility.py          # 指标 → 筛选决策量换算
python make_report_data.py # 汇总并注入报告页
```

`fetch.py` 与 `run_cc50.py` 会跳过已完成的部分，可安全重跑。
泛细胞毒性端点数据量大，`evalkit._topk_stable` 采用分块实现以控制内存。

## 目录

- `scripts/` — 取数、策展、评估、迁移、汇总与出图脚本
- `results/` — 全部实测指标（CSV）与报告页数据（JSON）
- `report/` — 报告页模板与注入数据后的成品页

## 限制

复算的是「同数据、同配置重训的等价模型」，不是套件页面内嵌的那份权重
（含模型的 HTML 未提供）。结论指向数据与方法的边界，不针对某一份权重的实现质量。
VC-AMR 与 VC-PTM 不是 QSAR，本轮未测。
