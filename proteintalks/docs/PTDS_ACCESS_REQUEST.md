# PTDS 完整矩阵申请材料

申请入口：<https://db.prottalks.com> → Download

想要的文件（属于 **PTV1 dataset** 分组，`requireApproval: true`）：

| 文件 | 规模 | 大小 |
|---|---|---|
| `02_protein_intensity_matrix.xlsx` | 16,311 样本（含 1,238 技术重复）× 5,583 protein groups | 325 MB |
| `01_sample_info.xlsx` | 15,073 样本的元数据 | 675 KB |

---

## 为什么这份申请必须由你本人提交

表单的每一个剩余空缺都指向你本人的真实身份，我无法代填，也不应该代填：

1. **机构邮箱**。选择 academic 路径时，前端会调用 `/pupy/email/validate` 做服务端校验。
   实测结果：`gmail.com → isEducational: false`，`qq.com` 被单独硬拒绝，`163.com` 同样不通过；
   `stanford.edu`、`westlake.edu.cn` 通过。你在本会话中已知的邮箱是 Gmail，会被直接驳回。
2. **PI 姓名、所属机构、职称、实验室主页**。这些是关于真实个人与真实实验室的声明。
   为拿到受限数据而编造任何一项属于身份不实陈述，无论出发点如何。
3. **许可协议勾选**（`agreedToTerms`）。这是一份需要你本人接受的法律协议
   （Guomics Lab Software and Data Usage License Agreement，站内 `/agreement.html`）。
4. **4 位图形验证码**（`/pupy/captcha` 返回图片与 token）。它的存在本身就是要求人工提交。

表单里还有一个 `website` 蜜罐字段（`tabindex="-1"`、关闭自动填充）和一个 `_ts` 时间戳，
说明站方明确在检测自动化提交。

**所以：请你自己打开页面填写。下面是可直接粘贴的内容。**

---

## 表单字段与约束（从前端校验规则逐条提取）

| 字段 | 提交键名 | 必填 | 约束 |
|---|---|---|---|
| Name | `name` | 是 | 最长 100 字符 |
| Usage type | `usageType` | 是 | `academic` 或 `non-academic` |
| Email | `email` | 是 | 合法邮箱；academic 路径须通过教育域校验，QQ 邮箱硬拒 |
| Institution | `institution` | 是 | — |
| PI name | `piName` | 是 | — |
| Job title | `jobTitle` | 是 | — |
| Lab web URL | `labWebUrl` | 是 | 必须是合法 URL |
| Country | `country` | 是 | 下拉选择 |
| Intended use | `usage` | 是 | **至少 20 个字符** |
| Verification code | `captchaCode` | 是 | 恰好 4 个字符 |
| License agreement | `agreedToTerms` | 是 | 必须勾选 |

提交后端点为 `POST /pupy/download/request`，随请求一并发送 `datasetTitles`（所选数据集分组标题）。
站方提示：部分数据集需审稿人批准，通过后邮件发送链接。

**如果没有机构邮箱**：选 `non-academic / commercial use`，请求会转发至
`aivc@westlakeomics.com`，不走教育域校验，但按商业用途处理。

---

## Intended use 草稿（英文，可直接粘贴）

选其一。两份都超过 20 字符下限，且如实描述用途。

### A. 独立复现与基准（贴合我们已经做的工作）

```
We are conducting an independent computational reproduction of the ProteinTalks
virtual cell model (Nature 2026, doi:10.1038/s41586-026-11001-9). We have already
verified our port of the released ppODE architecture against best_checkpoint.pth
to a residual of 1.8e-7, and reproduced the reported leave-one-cell-line-out and
leave-one-drug-out comparisons from Supplementary Tables S5 and S6.

We request the full PTDS protein intensity matrix in order to (1) retrain the
published model from scratch and confirm the reported AUROC on the drug-response
task, (2) run an ablation isolating the contribution of the neural-ODE dynamics
module, which is not reported in the paper, and (3) evaluate the model against
drug-identity and cell-line-identity control predictors on the proteomic features
rather than on labels alone.

No redistribution of the data is intended. Any resulting analysis code will cite
the original publication and the PTDS resource. We are happy to share our
findings with the authors prior to any public write-up.
```

### B. 方法学基准（若你的用途是更广的模型比较）

```
We are benchmarking perturbation-response models across modalities and require a
large-scale perturbation proteomics corpus as the protein-readout arm of the
comparison. The PTDS matrix is, to our knowledge, the only perturbation
proteomics dataset of sufficient scale for pretraining.

We intend to use it to compare proteome-based and transcriptome-based perturbation
prediction under matched evaluation protocols, including trivial and linear
control baselines. Results will be reported with full attribution to the original
publication, and the data will not be redistributed.
```

---

## 另外两点值得先想清楚

**一、缺的可能不只是矩阵。** 已发布 checkpoint 的第一层形状要求 **5,585** 个输入通道，
而门户标注的是 **5,583** 个 protein groups，差 2。同时，checkpoint 期望的蛋白排序索引
（README 里的 `node_Index.csv`）在任何渠道都没有发布。拿到矩阵之后，若列顺序与索引对不上，
**已训练权重仍然无法直接使用**，只能从头重训。建议在申请备注或后续邮件里一并索要
`node_Index.csv` 及 5583/5585 的口径说明。

**二、很多事情不需要这份矩阵。** 已经完全公开、无需登录的部分包括：模型代码与权重（MIT）、
1,116 条疗效标签、63 个药物 SMILES、2,487 × 3,631 的多时间点蛋白矩阵（Table S2）、
501 例患者矩阵与生存数据（Table S13）、以及全部逐细胞系与逐药物基准结果（Table S5、S6）。
本仓库的三项验证全部建立在这些之上。矩阵真正解锁的是**重训模型**与**在真实数据上做动力学模块消融**。
