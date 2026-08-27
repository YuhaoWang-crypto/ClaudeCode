"""
E5 — 用真实的人类 3'UTR k-mer 组成 + 真实突变谱重标定 M2 的 Λ_liq。

E2 里的 Λ_liq = L·4^-k 有两个理想化假设，两个都错：
  (1) 3'UTR 是 25% 等频组成  → 实际强烈 AU 偏倚
  (2) 所有替换等速率          → 实际有 transition 偏倚与 CpG 超突变

本模块换成：
  · 真实序列：Ensembl BioMart 全部人类 canonical protein-coding 3'UTR
  · 真实靶位点：TargetScan miR_Family_Info 的 human seed+m8（7mer-m8 位点 = 其反向互补）
  · 突变模型：HKY(κ) + CpG 超突变，按全 UTR 平均归一化

得到逐 (miRNA family × 基因) 的
      λ_gain = Σ_{单突变邻居 n} count(n)·rate(n→t)      （每单位中性分歧度 d）
      λ_loss = Σ_{位点内 7 个位置} 总替换率                （每位点每单位 d）
      Λ_liq  = λ_gain / λ_loss
并与 L·4^-7 的朴素零模型比较。
"""
import glob
import os
import numpy as np

from .e2_liquidity import branch_stats_rates

K = 7
KAPPA = 2.0          # transition / transversion
CPG_MULT = 6.0       # CpG 位点的 C->T / G->A 超突变倍数
B2I = {"A": 0, "C": 1, "G": 2, "T": 3}
I2B = "ACGT"
TRANSITION = {(0, 2), (2, 0), (1, 3), (3, 1)}

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
UTR_DIR = os.path.join(DATA, "utr3")
FAM_FILE = os.path.join(DATA, "miR_Family_Info.txt")


# ---------------------------------------------------------------- 数据载入
def load_utrs(min_len=50):
    """返回 {gene: 编码后的 int8 数组}（同基因取最长的一条）。"""
    best = {}
    for path in sorted(glob.glob(os.path.join(UTR_DIR, "*.fa"))):
        name, buf = None, []
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if line.startswith(">"):
                    if name:
                        _keep(best, name, "".join(buf))
                    name, buf = line[1:].split("|")[0], []
                else:
                    buf.append(line)
        if name:
            _keep(best, name, "".join(buf))
    return {g: encode(s) for g, s in best.items()
            if len(s) >= min_len and g and "Sequence unavailable" not in s}


def _keep(best, name, seq):
    if name not in best or len(seq) > len(best[name]):
        best[name] = seq.upper()


def encode(seq):
    a = np.frombuffer(seq.encode(), dtype=np.uint8)
    out = np.full(a.shape, 255, dtype=np.uint8)
    for b, i in B2I.items():
        out[a == ord(b)] = i
    return out


def load_families(min_cons=2):
    """human (9606) 的 (family, seed+m8) 去重表。min_cons: 2=广泛保守, 1=保守。"""
    fams = {}
    with open(FAM_FILE, errors="ignore") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6 or f[2] != "9606":
                continue
            try:
                cons = int(f[5])
            except ValueError:
                continue
            if cons < min_cons:
                continue
            seed = f[1].strip().replace("U", "T").upper()
            if len(seed) == K and set(seed) <= set("ACGT"):
                fams.setdefault(f[0].strip(), seed)
    return fams


# ---------------------------------------------------------------- k-mer 工具
def kmer_index(arr, k=K):
    """int8 序列 → 每个位置起始的 k-mer 整数编码（含 N 的窗口标 -1）。"""
    n = len(arr) - k + 1
    if n <= 0:
        return np.empty(0, dtype=np.int64)
    idx = np.zeros(n, dtype=np.int64)
    bad = np.zeros(n, dtype=bool)
    for j in range(k):
        col = arr[j:j + n].astype(np.int64)
        bad |= col > 3
        idx = idx * 4 + np.where(col > 3, 0, col)
    idx[bad] = -1
    return idx


def to_int(kmer):
    v = 0
    for c in kmer:
        v = v * 4 + B2I[c]
    return v


def to_str(v, k=K):
    return "".join(I2B[(v >> (2 * (k - 1 - i))) & 3] for i in range(k))


def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


# ---------------------------------------------------------------- 突变模型
def sub_rate(x, y, pi, left=None, right=None, p_left_C=0.0, p_right_G=0.0):
    """HKY(κ) + CpG 的 x→y 相对替换率。侧翼未知时用其期望概率。"""
    r = pi[y] * (KAPPA if (x, y) in TRANSITION else 1.0)
    if x == 1 and y == 3:                                   # C->T，需右侧 G
        f = (CPG_MULT if right == 2 else 1.0) if right is not None \
            else 1 + (CPG_MULT - 1) * p_right_G
        r *= f
    elif x == 2 and y == 0:                                 # G->A，需左侧 C
        f = (CPG_MULT if left == 1 else 1.0) if left is not None \
            else 1 + (CPG_MULT - 1) * p_left_C
        r *= f
    return r


def _bases(v, k=K):
    return [(v >> (2 * (k - 1 - i))) & 3 for i in range(k)]


def kmer_loss_rate(t_int, pi, k=K):
    """一个已有位点每单位时间的总丢失率（未归一化）。两端侧翼用期望组成。"""
    b, tot = _bases(t_int, k), 0.0
    for i in range(k):
        left = b[i - 1] if i > 0 else None
        right = b[i + 1] if i < k - 1 else None
        for y in range(4):
            if y != b[i]:
                tot += sub_rate(b[i], y, pi, left, right,
                                p_left_C=pi[1], p_right_G=pi[2])
    return tot


def neighbor_gain_rates(t_int, pi, k=K):
    """返回 (邻居 k-mer 编码, 邻居→目标 的特定替换率)。"""
    b, idxs, rates = _bases(t_int, k), [], []
    for i in range(k):
        for x in range(4):
            if x == b[i]:
                continue
            nb = list(b)
            nb[i] = x
            v = 0
            for c in nb:
                v = v * 4 + c
            left = nb[i - 1] if i > 0 else None
            right = nb[i + 1] if i < k - 1 else None
            idxs.append(v)
            rates.append(sub_rate(x, b[i], pi, left, right,
                                  p_left_C=pi[1], p_right_G=pi[2]))
    return np.array(idxs, dtype=np.int64), np.array(rates, dtype=np.float64)


def mean_position_rate(utrs, pi):
    """全 UTR 实测的平均每位点总替换率 → 把速率归一化到"每单位中性分歧度 d"。"""
    base_tot = np.array([sum(sub_rate(x, y, pi) for y in range(4) if y != x)
                         for x in range(4)])
    cpg_c = KAPPA * pi[3] * (CPG_MULT - 1)      # C->T 在 CpG 中的增量
    cpg_g = KAPPA * pi[0] * (CPG_MULT - 1)      # G->A 在 CpG 中的增量
    tot, n = 0.0, 0
    for arr in utrs.values():
        a = arr.astype(np.int16)
        ok = a <= 3
        left = np.concatenate(([-1], a[:-1]))
        right = np.concatenate((a[1:], [-1]))
        r = np.where(ok, base_tot[np.clip(a, 0, 3)], 0.0)
        r += np.where(ok & (a == 1) & (right == 2), cpg_c, 0.0)
        r += np.where(ok & (a == 2) & (left == 1), cpg_g, 0.0)
        tot += float(r.sum())
        n += int(ok.sum())
    return tot / n, n


# ---------------------------------------------------------------- 主分析
def analyse(min_cons=2):
    utrs = load_utrs()
    fams = load_families(min_cons)
    comp = np.zeros(4)
    for arr in utrs.values():
        comp += np.bincount(arr[arr <= 3], minlength=4)
    comp /= comp.sum()
    rbar, n_pos = mean_position_rate(utrs, comp)

    names = sorted(fams)
    targets = np.array([to_int(revcomp(fams[f])) for f in names])
    loss = np.array([kmer_loss_rate(t, comp) for t in targets]) / rbar
    nb_idx = np.zeros((len(names), 3 * K), dtype=np.int64)
    nb_rate = np.zeros((len(names), 3 * K))
    for i, t in enumerate(targets):
        nb_idx[i], r = neighbor_gain_rates(t, comp)
        nb_rate[i] = r / rbar

    genes = sorted(utrs)
    lens = np.array([len(utrs[g]) for g in genes], dtype=float)
    n_obs = np.zeros((len(genes), len(names)))
    n_a1 = np.zeros((len(genes), len(names)))
    lam_g = np.zeros((len(genes), len(names)))
    a1_idx = np.array([to_int(revcomp(fams[f][:K - 1]) + "A") for f in names])
    kmer_tot = np.zeros(4 ** K)
    for gi, g in enumerate(genes):
        idx = kmer_index(utrs[g])
        idx = idx[idx >= 0]
        if idx.size == 0:
            continue
        cnt = np.bincount(idx, minlength=4 ** K).astype(np.float64)
        kmer_tot += cnt
        n_obs[gi] = cnt[targets]
        n_a1[gi] = cnt[a1_idx]
        lam_g[gi] = (cnt[nb_idx] * nb_rate).sum(axis=1)
    lam_liq = lam_g / loss[None, :]
    return dict(genes=genes, names=names, fams=fams, lens=lens, comp=comp,
                targets=targets, loss=loss, n_obs=n_obs, n_a1=n_a1, lam_g=lam_g,
                lam_liq=lam_liq, kmer_tot=kmer_tot, rbar=rbar, n_pos=n_pos,
                utrs=utrs)


def canonical_sites(seq, seed_dna):
    """8mer / 7mer-m8 / 7mer-A1 计数与并集（TargetScan 的三类规范位点）。"""
    m8 = revcomp(seed_dna)
    a1 = revcomp(seed_dna[:K - 1]) + "A"     # revcomp(miRNA nt2-7) + A
    m8a1 = m8 + "A"
    n8, nm8, na1 = seq.count(m8a1), seq.count(m8), seq.count(a1)
    return n8, nm8, na1, nm8 + na1 - n8      # 并集（8mer 同时被前两者计入）


def run(min_cons=2):
    print("=" * 78)
    print("E5 — 用真实 3'UTR 组成 + HKY+CpG 突变谱重标定 Λ_liq")
    print("=" * 78)
    R = analyse(min_cons)
    genes, names, lens = R["genes"], R["names"], R["lens"]
    lam_liq, n_obs, loss = R["lam_liq"], R["n_obs"], R["loss"]
    null = np.repeat(lens[:, None], len(names), axis=1) * 0.25 ** K

    print(f"\n[0] 数据：{len(genes)} 个基因的 3'UTR，{len(names)} 个广泛保守 miRNA family")
    print(f"    UTR 总碱基 {R['n_pos']/1e6:.1f} Mb；组成 A/C/G/T = "
          + "/".join(f"{x:.3f}" for x in R["comp"]))
    print(f"    UTR 长度中位数 {np.median(lens):.0f} bp，均值 {lens.mean():.0f} bp")
    print(f"    平均每位点总替换率 R̄ = {R['rbar']:.3f}（相对单位；用于归一化到 d）")

    print("\n[1] 真实 vs 朴素零模型（总量口径，不受短 UTR 主导）")
    agg = lam_liq.sum() / null.sum()
    print(f"    Σ Λ_liq(真实) / Σ Λ_liq(L·4^-7) = {agg:.2f}")
    frac_r, frac_n = (lam_liq > 1).mean(), (null > 1).mean()
    print(f"    处于 liquid 区 (Λ_liq>1) 的 (基因×family) 对："
          f"真实 {frac_r:.3%}  vs  朴素 {frac_n:.3%}  →  放大 {frac_r/max(frac_n,1e-12):.1f}×")
    print(f"    真实 liquid 对的绝对数量 = {int((lam_liq>1).sum()):,} "
          f"（覆盖 {int((lam_liq>1).any(axis=1).sum()):,} 个基因）")

    print("\n[2] 位点丰度：观测 7mer-m8 计数 vs 均匀组成预期")
    tot = R["kmer_tot"]
    exp_uniform = tot.sum() * 0.25 ** K
    obs_fam = tot[R["targets"]]
    print(f"    全 UTR 7mer 总数 {tot.sum()/1e6:.1f} M；均匀预期每种 {exp_uniform:,.0f}")
    print(f"    保守 family 位点中位数 {np.median(obs_fam):,.0f} → 富集 "
          f"{np.median(obs_fam)/exp_uniform:.2f}×，"
          f"跨 family 变异 {obs_fam.min():,.0f}–{obs_fam.max():,.0f}（{obs_fam.max()/obs_fam.min():.0f}×）")

    # family 的内禀量：每 kb 的 Λ_liq，以及达到 liquid 所需的临界 UTR 长度
    per_kb = 1000.0 * (lam_liq.sum(axis=0) / lens.sum())
    Lstar = 1000.0 / per_kb
    order = np.argsort(-per_kb)
    print("\n[3] family 内禀量：临界 UTR 长度 L*（Λ_liq=1 所需长度），朴素模型给 16.4 kb")
    print(f"{'family':>26} {'位点':>9} {'AU%':>5} {'CpG':>4} {'λ_loss':>7} "
          f"{'Λ_liq/kb':>9} {'L*(kb)':>8} {'位点总数':>9}")
    for tag, sel in (("最液态 TOP6", order[:6]), ("最晶态 TOP6", order[-6:])):
        print(f"  --- {tag} ---")
        for i in sel:
            site = to_str(R["targets"][i])
            au = 100 * sum(c in "AT" for c in site) / K
            print(f"{names[i][:26]:>26} {site:>9} {au:>5.0f} "
                  f"{'Y' if 'CG' in site else '-':>4} {loss[i]:>7.2f} "
                  f"{per_kb[i]:>9.4f} {Lstar[i]/1000:>8.1f} {tot[R['targets'][i]]:>9,.0f}")

    cg = np.array(["CG" in to_str(t) for t in R["targets"]])
    print(f"\n    含 CpG 的 family (n={cg.sum():2d}): λ_loss {np.median(loss[cg]):5.2f}，"
          f"L* 中位 {np.median(Lstar[cg])/1000:5.1f} kb")
    print(f"    不含 CpG      (n={(~cg).sum():2d}): λ_loss {np.median(loss[~cg]):5.2f}，"
          f"L* 中位 {np.median(Lstar[~cg])/1000:5.1f} kb")
    print(f"    → CpG 家族的位点丢失快 {np.median(loss[cg])/np.median(loss[~cg]):.2f}×，"
          f"进入 liquid 所需 UTR 长 {np.median(Lstar[cg])/np.median(Lstar[~cg]):.1f}×")

    print("\n[4] 重标定后的 conservation-scan 漏检率（人鼠 d=0.5, α=1e-4）")
    print(f"{'(基因×family) 对分层':>22} {'Λ_liq':>8} {'边保留':>8} {'坐标保留':>9} {'漏检率':>8}")
    flat = lam_liq.ravel()
    lo = np.median(loss)
    for tag, val in (("全部中位数", np.median(flat)),
                     ("P90", np.percentile(flat, 90)),
                     ("P99", np.percentile(flat, 99)),
                     ("P99.9", np.percentile(flat, 99.9)),
                     ("最高 1 对", flat.max())):
        e, p = branch_stats_rates(val * lo, lo, 1e-4, 0.5)
        print(f"{tag:>22} {val:>8.3f} {e:>8.3f} {p:>9.3f} {1 - p/e:>8.1%}")
    print("\n    → 漏检率是 Λ_liq 的单调函数；真正被 conservation scan 系统性漏掉的，")
    print("      是 Λ_liq 高的那个尾巴（长 UTR × AU 富集 seed），不是全基因组平均。")

    print("\n[5] 案例基因逐条报告（用于 biomarker 案例）")
    print("    位点数 = 8mer / 7mer-m8 / 7mer-A1 / 并集；Λ_liq 与保留率按 7mer-m8 类算")
    print(f"{'基因':>10} {'family':>22} {'L(bp)':>7} {'8m':>3} {'7m8':>4} {'7A1':>4} "
          f"{'并集':>4} {'Λ_liq':>7} {'边保留':>7} {'坐标保留':>8} {'漏检率':>7} {'判定':>12}")
    for gene, key in CASES:
        gene_report(R, gene, key)

    print("\n[6] 单位点变异的预测效应量（案例 6 的定律）")
    print("    ΔR/R = 1/n_union（冗余缓冲） × 1/(1+Λ_liq)（可补偿性）")
    print(f"{'基因×family':>34} {'n_union':>8} {'Λ_liq':>7} {'预测 ΔR/R':>11} {'判定':>18}")
    rows = []
    for gene, key in CASES:
        if gene not in R["genes"]:
            continue
        hits = [i for i, n in enumerate(R["names"]) if key in n]
        if not hits:
            continue
        fi, gi = hits[0], R["genes"].index(gene)
        seq = "".join(I2B[b] if b < 4 else "N" for b in R["utrs"][gene])
        _, _, _, nu = canonical_sites(seq, R["fams"][R["names"][fi]])
        ll = R["lam_liq"][gi, fi]
        if nu == 0:
            eff, verd = 0.0, "无规范位点→预测无效应"
        else:
            eff = (1.0 / nu) / (1.0 + ll)
            verd = "可检测" if eff > 0.4 else ("边缘" if eff > 0.2 else "预测检测不到")
        rows.append((f"{gene}×{R['names'][fi][:14]}", nu, ll, eff, verd))
    for name, nu, ll, eff, verd in sorted(rows, key=lambda r: -r[3]):
        print(f"{name:>34} {nu:>8} {ll:>7.3f} {eff:>11.3f} {verd:>18}")

    print("\n[7] 全基因组：有多少边是'单位点、无冗余'的（变异敏感边）")
    nm8 = R["n_obs"]
    na1 = R["n_a1"]
    nu_all = nm8 + na1                                  # 上界（未扣 8mer 重叠）
    has = nu_all > 0
    print(f"    有 >=1 个规范位点的 (基因×family) 对：{int(has.sum()):,} "
          f"（占 {has.mean():.1%}）")
    for n in (1, 2, 3):
        sel = nu_all == n
        print(f"    恰好 {n} 个位点：{int(sel.sum()):>7,} 对 "
              f"（占有位点者 {sel.sum()/has.sum():>5.1%}），预测单变异效应 ~{1/n:.2f}")
    print(f"    >=4 个位点：{int((nu_all>=4).sum()):>7,} 对 "
          f"（占 {(nu_all>=4).sum()/has.sum():.1%}），预测单变异效应 <0.25")
    print("\n    → 约一半的边是单位点边：它们才是 UTR 变异 biomarker 应该找的地方；")
    print("      多位点边（HMGA2/IGF2BP1 型）对单变异免疫，只对整段 UTR 丢失敏感。")
    return R


CASES = [
    ("HMGA2", "let-7"), ("LIN28B", "let-7"), ("KRAS", "let-7"),
    ("IGF2BP1", "let-7"), ("MYC", "let-7"),
    ("CCND1", "miR-15-5p"), ("BCL2", "miR-15-5p"), ("CCND1", "miR-17-5p"),
    ("E2F1", "miR-17-5p"), ("CDKN1A", "miR-17-5p"),
    ("SIRT1", "miR-34-5p"), ("MYC", "miR-34-5p"),
    ("SOCS1", "miR-155-5p"), ("INPP5D", "miR-155-5p"),
    ("PAK1", "miR-17-5p"), ("NUDT21", "let-7"),
    ("GAPDH", "miR-125-5p"), ("FOXA1", "miR-1-3p"), ("NFKBIZ", "miR-155-5p"),
]


def gene_report(R, gene, key, alpha=1e-4, d=0.5):
    if gene not in R["genes"]:
        print(f"{gene:>10} {'(UTR 未取到)':>22}")
        return
    gi = R["genes"].index(gene)
    hits = [i for i, n in enumerate(R["names"]) if key in n]
    if not hits:
        print(f"{gene:>10} {key + ' (family 未匹配)':>22}")
        return
    fi = hits[0]
    L, ll = R["lens"][gi], R["lam_liq"][gi, fi]
    lo = R["loss"][fi]
    seq = "".join(I2B[b] if b < 4 else "N" for b in R["utrs"][gene])
    n8, nm8, na1, nu = canonical_sites(seq, R["fams"][R["names"][fi]])
    e, p = branch_stats_rates(ll * lo, lo, alpha, d)
    verdict = "liquid" if ll > 1 else ("marginal" if ll > 0.3 else "crystalline")
    print(f"{gene:>10} {R['names'][fi][:22]:>22} {L:>7.0f} {n8:>3} {nm8:>4} {na1:>4} "
          f"{nu:>4} {ll:>7.3f} {e:>7.3f} {p:>8.3f} {1 - p/e:>7.1%} {verdict:>12}")


if __name__ == "__main__":
    run()
