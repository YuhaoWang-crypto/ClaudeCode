# 在 LatchBio 上补跑打分与脱靶(运行说明)

这一步补齐 `13_引导设计实跑结果.md` 里标记为待办的两项:**Rule Set 3 on-target 打分**
与 **全基因组脱靶(GuideScan2)**。输入已备好、脚本已就绪、已本地自测通过——**唯一
的拦路点是 Latch 账号还没有默认 workspace**。

---

## 0. 当前状态:被账号设置卡住(需你操作)

我本会话多次调用 `list_workspaces`,返回始终是:

```json
{"default_workspace_id": null, "workspaces": []}
```

所有 Latch 工具(`list_workflows`/`launch_workflow`/…)都因此报错:

```
No default workspace is configured. Please make sure you have completed
setting up your Latch account.
```

**这一步我无法代做**——Latch 的 MCP 工具里**没有"创建 workspace"能力**,只能列出/
启动。需要**你在 Latch 控制台完成账号设置**:

1. 登录 <https://console.latch.bio>
2. 完成账号/团队设置,使账号拥有一个**默认 workspace**(新建或加入一个团队工作区)
3. 完成后告诉我——我会重新 `list_workspaces` 确认拿到 `default_workspace_id`,即可继续

> 会话是非交互的,我无法在这里替你走 OAuth/控制台流程;这一步必须在浏览器里由你完成。

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

## 2. workspace 就绪后,两条路任选

### 路线 A —— 我用 MCP 直接在 Latch 上跑(推荐)

你把默认 workspace 建好后,我会:

1. `list_workflows`(search `guidescan` / `crispr` / `rule set`)找现成工作流;
2. `get_workflow_schema` 读参数;
3. 把 `data/all_spacers.txt`(或 `rs3_context.tsv`)传到 Latch Data;
4. `launch_workflow` 启动,`get_execution` 轮询状态,失败时 `get_task_logs` 排查;
5. 结果拉回,跑 `score_guides.py merge` 生成 `data/guides_scored.tsv`。

> 若 workspace 里没有现成的 GuideScan2/引导设计工作流,我会改走路线 B,或帮你把脚本
> 封装成一个可复用的 Latch 工作流(需要你的 workspace 有注册工作流的权限)。

### 路线 B —— 在 Latch Pod 或任何有 hg38 的机器上跑脚本

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
