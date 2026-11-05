# FineDNA demo

[`HuggingFaceBio/FineDNA`](https://huggingface.co/datasets/HuggingFaceBio/FineDNA) 是 Carbon-X（把 Qwen3-14B
训练成能"读 DNA"的模型）三个训练阶段的全部数据：

| 阶段 | 配置 | 规模 |
|---|---|---|
| Mid-training | `mid_*`，23 个来源（基因组、mRNA、含 DNA 的生物文献、英文网页、指令数据） | 514B tokens，约 950 GB |
| SFT | `sft_mix` + 8 个领域子集 | 203,923 条对话 |
| RL | `rl_mix` | 98,433 条可用代码判分的题 |

数据集约 1 TB，不能全量下载。这个 demo 用 `datasets` 的**流式读取**，每个配置只拉几百到几千行
（总共约 150 MB 网络读取，2–3 分钟），然后做三件事。

## 运行

```bash
pip install -r finedna_demo/requirements.txt
```

```bash
python finedna_demo/finedna_demo.py
```

```bash
python finedna_demo/finedna_demo.py --n-dna 100 --n-sft 500 --n-rl 2000   # 更快的小样本
```

抽到的行缓存在 `finedna_demo/cache/`（不入库），重跑只重做分析。

## 1. 预训练 DNA：组成是否符合生物学常识

5 个 DNA 配置各抽 300 条（流式 + 1000 行洗牌缓冲）：

![DNA composition](figures/dna_composition.png)

| 组 | 序列长度中位数 | GC% 中位数 | CpG 观测/期望 | 长度为 6 的倍数 | 带物种/基因类型标签 |
|---|---:|---:|---:|---:|---:|
| 哺乳动物 | 7,455 bp | 45.5 | **0.31** | 100% | 48% |
| 无脊椎动物 | 4,614 bp | 37.4 | 0.91 | 100% | 50% |
| 植物 | 2,484 bp | 38.3 | 0.56 | 100% | 47% |
| 原核生物 | 7,887 bp | 50.4 | **1.00** | 100% | 0% |
| mRNA | 1,698 bp | 48.0 | 0.71 | 100% | 0% |

- **CpG 缺失梯度**符合预期：哺乳动物 DNA 甲基化导致 CpG 强烈缺失（≈0.3），植物居中，
  很多无脊椎动物（如昆虫）几乎不甲基化（≈0.9），细菌没有缺失（≈1.0）。这说明数据是真实基因组序列，没有被打乱。
- 原核生物 GC% 呈双峰（约 35% 和 65%），反映不同菌种的基因组 GC 差异。
- README 的格式声明都对得上：所有序列长度都是 6 的倍数（模型按不重叠的 6-mer 分词），
  真核序列约一半带标签，如 `<mammalian_species><protein_coding_region><dna>ATG…</dna>`。

分词示例（`summary.json` → `tokenization_example`）：
`ATTTTTTTTTAGACAG…` → `ATTTTT | TTTTAG | ACAGAA | …`

## 2. SFT / RL 的任务构成

![SFT and RL mix](figures/sft_rl_mix.png)

SFT 以序列续写、变异致病性判断、DNA 计算三类为主（共约 78%）。每个领域的一条完整对话样例在
[`results/sft_examples.md`](results/sft_examples.md)。

## 3. 独立复核 RL 标准答案

`rl_mix` 中 `basic_dna` 组的题都是可算出确定答案的小题（GC 含量、反向互补、酶切片段数、翻译、
Cas9 切点、Tm 值……共 60 多种）。[`verifiers.py`](verifiers.py) 用纯 Python 为其中 24 种题型重写了
求解器，从题面解析序列、自己算答案，再和数据集给的 `answer_text` 比对：

| 抽样 5,000 行 RL 数据 | 行数 |
|---|---:|
| `basic_dna` 题 | 1,581 |
| 有我们求解器的题型 | 720 |
| 与标准答案一致 | **717（100%）** |
| 不一致 | 0 |
| 题面措辞没解析出来 | 3 |

逐条结果在 [`results/rl_check_details.csv`](results/rl_check_details.csv)。这些求解器也可以直接当作 RL
的奖励函数使用。

## 4. 在 Hugging Face 上考一考线上的 Carbon-X

Carbon-X（14.8B 参数）没有官方托管的推理接口，但社区有 Space 在 ZeroGPU（A10G）上跑它。
[`eval_carbonx_space.py`](eval_carbonx_space.py) 通过 `gradio_client` 远程调用，本地不需要 GPU。
三个运行 Carbon-X 的 Space：

| Space | 能否用脚本调用 | token 上限 | 模型版本 |
|---|---|---:|---|
| [`mrfakename/carbon-x`](https://huggingface.co/spaces/mrfakename/carbon-x)（默认） | ✓ `/chat` | 1536 | 固定在 `1fac52c` |
| [`mbcredits/carbon-x-demo`](https://huggingface.co/spaces/mbcredits/carbon-x-demo) | ✓ `/ask_carbon` | 1024 | 最新 `main`（`17de891`） |
| [`fffiloni/Carbon-X-Playground`](https://huggingface.co/spaces/fffiloni/Carbon-X-Playground) | ✗ 关闭了 API，只能在网页上用 | — | `1fac52c` |

两个版本的模型权重和聊天模板完全相同，只差一次 tokenizer 修复和结束符配置调整。两边都关闭 thinking，
temperature 0.6（作者报告结果时用的设置）。每条记录都写明是哪个 Space 答的。

```bash
pip install gradio_client
```

```bash
export HF_TOKEN=hf_...   # 用你自己账号的 ZeroGPU 配额
python finedna_demo/eval_carbonx_space.py --per-type 3 --skip-types gc_content
```

```bash
python finedna_demo/eval_carbonx_space.py --per-type 3 --space mbcredits/carbon-x-demo --skip-types gc_content
```

**ZeroGPU 配额按 Hugging Face 账号计算，所有 Space 共用**，免费账号每天约 5 分钟 GPU 时间，换 Space 不会多出额度。

每种题型出两类题：

- **seen**：`rl_mix` 原题。Carbon-X 的 RL 阶段就是用这些题训练的，所以它见过。
- **fresh**：题面不变，把 DNA 换成同长度的随机序列，标准答案由 `verifiers.py` 重新计算。模型没见过，考的是算，不是记。

回答缓存在 `results/carbonx_space_eval.jsonl`，配额用完后第二天重跑会接着问。被截断的回答只有在换用更高的
token 上限时才会重问。

**目前结果（用了两天的免费配额，共 7 道题）：**

| 题型 | seen | fresh |
|---|---|---|
| 反向互补 | 0/1（步骤对，结尾多写了 `TG`） | 1/1 |
| 互补链 | 1/1 | 1/1 |
| 转录（DNA→mRNA） | 1/1 | — |
| GC 含量（118 bp） | 截断 ×2 | 截断 ×2 |

**GC 含量题在这两个 Space 上答不完。** 模型的做法是逐个碱基累加计数（`97. C (G=22, C=24)` 这样一行一个碱基），
每个碱基约 15 个 token，118 bp 的序列需要约 1,900 个 token。在 `mbcredits`（512 token）停在第 28 个碱基，
在 `mrfakename`（1536 token）停在第 97–99 个碱基，而且每题耗掉约 70 秒 GPU 时间。截断的回答不计对错；
继续跑时建议加 `--skip-types gc_content`。

7 题太少，还不能下结论。剩下 61 题（不含 GC）大约需要 8–10 分钟 GPU 时间：免费账号要跑 2–3 天，
PRO 账号（每天 40 分钟）一次就能跑完。

## 局限

- 这是抽样，不是全量统计。真核配置的 300 条只来自 25–41 个不同的基因组 contig（流式读取集中在一个分片），
  表中数字是示意性的，不代表整个配置的分布。
- 只看了数据，没有加载 Carbon-X 模型做推理（14B 参数需要 GPU）。
- 数据集 Apache-2.0；mid-training 部分来自第三方来源，各有自己的许可。
