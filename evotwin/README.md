# evotwin — EvoTwin 的位置/结构/遗传图谱扩展

```bash
pip install numpy scipy
python3 -m evotwin.run_all          # E1-E4（纯理论，无需数据）
python3 -m evotwin.e5_real_kmer     # E5（需要 data/ 下的真实序列）
```

| 模块 | 内容 |
|---|---|
| `e1_target_supply` | 功能突变靶大小 L_eff = L·3k·4^-k；U_gain/U_loss；等待时间；受选择边比例反演 |
| `e2_liquidity` | 3'UTR 位点的标记点过程；Λ_liq = L·4^-k；精确 CTMC 的边保留 vs 坐标保留 |
| `e3_recombination` | 两位点选择-重组递推；r < s_ε 阈值；共依赖模块的最大遗传跨度 |
| `e4_avoidance` | 回避的系统发育似然比与所需物种数功效分析 |
| `e5_real_kmer` | **用真实 3'UTR + HKY+CpG 突变谱重标定 Λ_liq**；案例基因逐条报告 |

## data/ 的来源

- `data/utr3/*.fa` — Ensembl BioMart，人类 canonical protein-coding 转录本的 3'UTR
  （FASTA header = `基因名|转录本ID`）。`fetch.sh` 可重新拉取；BioMart 会限流，
  脚本带重试与分段。
- `data/miR_Family_Info.txt` — TargetScan v8.0 `miR_Family_Info.txt.zip` 解压。

两者都不入 git 之外的处理；重跑 `e5_real_kmer` 会自动读取 `data/utr3/` 下所有 `.fa`。
