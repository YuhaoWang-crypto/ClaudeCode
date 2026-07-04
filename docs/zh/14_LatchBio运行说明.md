# 在 LatchBio 上补跑打分与脱靶(运行说明)

这一步补齐 `13_引导设计实跑结果.md` 里标记为待办的两项:**Rule Set 3 on-target 打分**
与 **全基因组脱靶(GuideScan2)**。输入已备好、脚本已就绪、已本地自测通过——**唯一
的拦路点是 Latch 账号还没有默认 workspace**。

---

## 0. 当前状态:workspace 已就绪 ✅

`list_workspaces` 现返回 **`default_workspace_id: 42942`(名称 "Claude")**,Latch 连接
完全可用。

### 但:Latch 目录里**没有引导设计/脱靶工作流**

我按 `guidescan / off-target / sgrna / crispr / screen` 逐个搜索,结论:

| 搜索 | 结果 |
|------|------|
| guidescan / off-target / sgrna | **无**——Latch 没有 GuideScan2 或引导脱靶设计工作流 |
| crispr | CRISPResso2、nf-core/crisprseq(**编辑结局分析**,非引导设计) |
| mageck | **MAGeCK Count / Test / MLE / Pathway / Plot**(**筛选数据分析**) |

**结论:** 你想在 Latch 上跑的**脱靶打分(GuideScan2)在 Latch 目录里没有**,MCP 也
不能启动任意 Pod/脚本——所以这一步仍走**路线 B**(Pod/本地 hg38,见下),且注意:
**2,024 条 Brunello 引导原库已做脱靶过滤,真正需要新算脱靶的只有 30 条 de-novo 引导**。

Latch 真正能帮上的是**下游筛选分析**(MAGeCK),那是 `06/09` 的分析阶段——需要**测序
FASTQ**(本设计阶段还没有)。为让它即插即用,我已生成 MAGeCK 所需的库文件(见下)。

---

## 1. 已经为你备好的东西(无需等 Latch)

| 文件 | 内容 | 状态 |
|------|------|------|
| `data/all_spacers.txt` | 2,054 条 spacer(20 nt),脱靶工具输入 | ✅ 已生成 |
| `data/rs3_context.tsv` | 每条引导的 30-mer 上下文,Rule Set 3 输入 | ✅ 已生成(0 缺失) |
| `scripts/score_guides.py` | on-target(rs3)+ 脱靶合并 + CFD 过滤 | ✅ 已自测通过 |
| `scripts/run_scoring.sh` | 一键编排(rs3 → GuideScan2 → 合并) | ✅ 就绪 |

`selected_guides.tsv`(2,024)+ `gap_genes_guides.tsv`(30)= **2,054 条**全部纳入。

---

## 1b. Latch 分析阶段:MAGeCK 库文件已备好

Latch 的 MAGeCK Count 需要两样输入:**库文件**(sgRNA/序列/基因)+ **样本 FASTQ**。
库文件我已从本设计生成:

- **`data/mageck_library.csv`** —— **2,304 条**(511 激酶 + 250 NTC),列
  `sgRNA,sequence,gene`,即 MAGeCK Count 的 `list_seq_file`。
- **`data/control_sgrnas.txt`** —— 250 条 NTC 的 id 列表(供 MAGeCK/drugZ 归一)。

**MAGeCK Count 工作流(id 85605)关键参数**(已核对 schema):
`list_seq_file`=上传的库文件,`sample_fastqs`=各样本 FASTQ,`sample_labels`=样本名,
`day0_label`=`T0`(打开负选 QC),`output_location`=`latch://42942.account/...`。

> **前提:** MAGeCK 需要**真实测序 FASTQ**——本设计阶段尚未产生,待实际筛选测序后才能跑。
> 那时把 `mageck_library.csv` 和 FASTQ 传到 Latch Data,我就能 `launch_workflow`
> 启动 Count → Test/MLE(A 线必需性),drugZ(B 线合成致死)在本地或 Pod 跑。
>
> **上传说明:** Latch 的 MCP 工具只有列出/下载,**没有上传**;文件需你在 console.latch.bio
> 的 Data 页上传(或用 `latch` CLI `latch cp`)。传好后告诉我路径,我来启动工作流。

## 2. 补跑脱靶/打分:在 Latch Pod 或任何有 hg38 的机器上跑脚本(路线 B)

```bash
pip install rs3
conda install -c bioconda guidescan
bash scripts/run_scoring.sh /path/to/hg38.fa
# 产出 data/guides_scored.tsv:
#   id, gene, spacer, pam, source, rs3_score, cfd_specificity, n_offtargets, pass_filter
```

脚本三步:①rs3 打 on-target 分 → ②GuideScan2 枚举脱靶+特异性 → ③合并并按
**CFD 特异性 ≥ 0.2** 过滤。脱靶工具也可换 CRISPOR(同样输出 CFD)。

---

## 3. 产出与收尾

`data/guides_scored.tsv` 出来后的最终整理:

1. **每基因在 PASS 引导中按 rs3_score 取前 6**;
2. 标记仍 **<6 条**的基因(可能因脱靶过滤淘汰过多)做人工复核或放宽阈值;
3. 更新最终文库清单,连同对照(~250 NTC + ~150 其他)定稿 → 送芯片合成。

到这一步,`13` 里两个 `pending` 字段(`rs2_score`/`offtarget`)就全部落地,
整个 gRNA 设计闭环完成。
