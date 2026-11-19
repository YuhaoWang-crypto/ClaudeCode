# -*- coding: utf-8 -*-
r"""abproto 矩阵教学演示：从 PDB 坐标到 Twisted-SMC，一步步用小矩阵看懂全流程。

纯 numpy、确定性、维度极小（6 个残基，特征宽度 D=4），所有数字都能按计算器核对。
迷你复合物（与配套 docx 一致）：
    1~3 号 = 固定抗原表位  E, W, K
    4~6 号 = 待生成抗体 HCDR3 残基  N, Y, S
词表：{MASK:0, E:1, W:2, K:3, N:4, Y:5, S:6, G:7}

运行：  python -m abproto.demo.matrix_walkthrough
"""

from __future__ import annotations

import numpy as np

np.set_printoptions(precision=3, suppress=True, floatmode="fixed")
rng = np.random.default_rng(0)

AA = {0: "MASK", 1: "E", 2: "W", 3: "K", 4: "N", 5: "Y", 6: "S", 7: "G"}
D = 4            # 特征宽度（真实代码里是 256）
N_ANTIGEN = 3   # 前 3 行是抗原
N_CDR = 3       # 后 3 行是待设计 CDR


def banner(t):
    print("\n" + "=" * 72 + f"\n{t}\n" + "=" * 72)


def show(name, M):
    print(f"\n{name}  shape={tuple(np.shape(M))}")
    print(np.array(M))


# ==========================================================================
def stage1_pdb_to_matrices():
    banner("STAGE 1 · 从 PDB 坐标到可计算矩阵（interpolants / ctmc）")
    # 真实终点坐标 X1（抗原中心对齐后），与 docx 完全一致。
    X1 = np.array([
        [0.0, 0.0, 0.0],   # E  抗原
        [3.8, 0.0, 0.0],   # W  抗原
        [1.9, 3.3, 0.0],   # K  抗原
        [1.0, 1.0, 4.0],   # N  CDR
        [3.5, 1.2, 3.8],   # Y  CDR
        [2.1, 3.0, 4.2],   # S  CDR
    ])
    seq1 = np.array([1, 2, 3, 4, 5, 6])          # E W K N Y S
    # 起点噪声 X0：抗原固定不动（=X1），CDR 注入噪声。
    X0 = X1.copy()
    X0[3:] = np.array([[5.0, 5.0, 10.0], [7.5, 5.2, 9.8], [6.1, 7.0, 10.2]])

    t = 0.5
    Xt = t * X1 + (1 - t) * X0                   # 半成品坐标（流匹配插值）
    U = X1 - X0                                   # 真实速度靶子（线性插值下为常数）

    # 结构先行调度 κ(t)=t²=0.25：只保留一部分 CDR 真值，其余掩码为 MASK。
    gen_mask = np.array([False, False, False, True, True, True])   # 哪些位要生成
    kappa = t * t
    keep = rng.random(6) < kappa                 # 被保留（未掩码）的概率 = κ
    keep[~gen_mask] = True                        # 抗原永远可见
    seq_in = np.where(keep, seq1, 0)             # 0 = MASK

    show("X1  真实终点坐标 (t=1.0)", X1)
    show("X0  起点噪声 (t=0.0, 抗原不动)", X0)
    show("Xt  半成品输入坐标 (t=0.5)", Xt)
    show("U = X1 - X0  速度靶子（只对 CDR 行有意义）", U)
    print("\n序列掩码 κ(t)=t²=0.25 后 seq_in:",
          [AA[i] for i in seq_in], " (gen_mask=", gen_mask.tolist(), ")")
    return dict(Xt=Xt, U=U, seq1=seq1, seq_in=seq_in, gen_mask=gen_mask)


# ==========================================================================
def stage2_embedding(st):
    banner("STAGE 2 · 4 股信息查表相加 → 输入特征矩阵 H0 [6×4]")
    Xt, seq_in = st["Xt"], st["seq_in"]
    # ① 氨基酸身份查表 nn.Embedding(8, D)
    E_tok = rng.normal(scale=0.3, size=(8, D))
    E_tok[0] = 0.0                                # MASK 行置零
    H_id = E_tok[seq_in]
    # ② 链角色查表 nn.Embedding(2, D)： 0=抗原 1=CDR
    E_chain = rng.normal(scale=0.2, size=(2, D))
    role = np.array([0, 0, 0, 1, 1, 1])
    H_chain = E_chain[role]
    # ③ 3D 坐标投影 nn.Linear(3, D)
    W_coord = rng.normal(scale=0.1, size=(3, D))
    H_coord = Xt @ W_coord
    # ④ 时间步编码（这里简化为把 t 放进最后一维）
    H_time = np.zeros((6, D)); H_time[:, -1] = 0.5

    H0 = H_id + H_chain + H_coord + H_time        # 按位相加
    print("四股分别查表/投影得到 [6×4]，再逐元素相加：")
    show("① H_id   氨基酸身份", H_id)
    show("② H_chain 链角色", H_chain)
    show("③ H_coord 坐标投影 = Xt @ W_coord", H_coord)
    show("④ H_time  时间编码", H_time)
    show("H0 = ①+②+③+④  输入骨干特征", H0)
    return dict(H0=H0, Xt=Xt, **{k: st[k] for k in ("U", "seq1", "gen_mask")})


# ==========================================================================
def stage3_attention(st):
    banner("STAGE 3 · 表位距离偏置自注意力 → 高阶特征 H1 [6×4]")
    H0, Xt = st["H0"], st["Xt"]
    Wq = rng.normal(scale=0.3, size=(D, D))
    Wk = rng.normal(scale=0.3, size=(D, D))
    Wv = rng.normal(scale=0.3, size=(D, D))
    Q, K, V = H0 @ Wq, H0 @ Wk, H0 @ Wv

    # 成对 3D 距离矩阵（埃），用于距离衰减偏置（IPA 的几何偏置精神）
    Dist = np.linalg.norm(Xt[:, None, :] - Xt[None, :, :], axis=-1)
    gamma = 0.15
    scores = (Q @ K.T) / np.sqrt(D) - gamma * Dist
    # softmax 按行
    A = np.exp(scores - scores.max(axis=1, keepdims=True))
    A = A / A.sum(axis=1, keepdims=True)
    H1 = A @ V + H0                               # 残差连接

    show("Q = H0 @ Wq", Q)
    show("成对距离矩阵 Dist (Å)", Dist)
    show("注意力矩阵 A = softmax(QKᵀ/√d − γ·Dist)  [6×6]", A)
    print("\n看第 5 行(CDR Y)：它把多少注意力投给每个位置——离得近的抗原表位权重更高：")
    for j in range(6):
        print(f"   Y → {AA[st['seq1'][j]]:>4}  权重 {A[4, j]:.3f}   距离 {Dist[4, j]:.2f} Å")
    show("H1 = A·V + H0  更新后的高阶特征", H1)
    return dict(H1=H1, **{k: st[k] for k in ("U", "seq1", "gen_mask")})


# ==========================================================================
def stage4_heads(st):
    banner("STAGE 4 · 两个输出头：速度头 [4×3] 与 序列头 [4×8]")
    H1, gen = st["H1"], st["gen_mask"]
    H_cdr = H1[gen]                               # 只对 CDR 行做预测  [3×4]
    W_vel = rng.normal(scale=0.3, size=(D, 3))   # 速度头：4 特征 → 3 个空间方向
    V_pred = H_cdr @ W_vel                        # [3×4]×[4×3] = [3×3]
    W_seq = rng.normal(scale=0.3, size=(D, 8))   # 序列头：4 特征 → 8 类氨基酸
    logits = H_cdr @ W_seq
    P = np.exp(logits - logits.max(1, keepdims=True))
    P = P / P.sum(1, keepdims=True)

    print("速度头权重 W_vel 的物理含义：3 列 = Vx, Vy, Vz 三个移动方向的“翻译表”。")
    show("H_cdr  待设计 3 个 CDR 残基的特征 [3×4]", H_cdr)
    show("V_pred = H_cdr @ W_vel  预测速度 [3×3] (Å/单位时间)", V_pred)
    show("序列头 softmax 概率 P [3×8]（列 = MASK,E,W,K,N,Y,S,G）", P)
    return dict(H_cdr=H_cdr, W_vel=W_vel, V_pred=V_pred,
                U_cdr=st["U"][gen])


# ==========================================================================
def stage5_training(st):
    banner("STAGE 5 · 反馈优化（训练）：对速度头做梯度下降，看 loss 下降")
    H, U = st["H_cdr"], st["U_cdr"]
    W = st["W_vel"].copy()
    lr, n = 0.02, H.shape[0]
    print("流匹配损失 L = ‖H_cdr·W_vel − U‖² / n；闭式梯度 dL/dW = (2/n)·Hᵀ(HW − U)")
    print("（这正是反向传播在一个线性头上的样子——一次矩阵乘法就是一次“反馈”）\n")
    print(" step |   loss   | W_vel 第一列 (Vx 方向权重)")
    for step in range(8):
        Vp = H @ W
        resid = Vp - U
        loss = (resid ** 2).mean()
        grad = (2 / n) * H.T @ resid             # [4×3] 梯度
        if step in (0, 1, 2, 4, 7):
            print(f"  {step:3d}  | {loss:8.3f} | {W[:,0].round(3)}")
        W -= lr * grad                            # 参数更新（反馈）
    print("\nloss 稳定下降 = 速度头学会了把特征翻译成正确的移动指令。")
    print("真实训练里同样的链式法则会一路把梯度传回注意力层、Embedding，整网一起更新。")


# ==========================================================================
def stage6_twisted_smc():
    banner("STAGE 6 · 推理：Twisted SMC (FK-Steering) + ESS 动态淘汰克隆")
    # 新抗原口袋中心（教学用）：CDR 应当扎进这个口袋。
    pocket = np.array([2.0, 2.0, 4.0])
    K = 6                                          # 粒子数（候选设计）
    T = 0.5                                         # 温度（越小重采样越尖锐）
    # 初始化 K 个粒子：每个是一个 CDR 残基的 3D 位置（简化为单点代表一个设计）
    particles = pocket + rng.normal(scale=3.0, size=(K, 3))
    names = [f"P{i}" for i in range(K)]

    def score(p):                                  # 打分 = 越贴口袋越高（界面置信度代理）
        return -np.linalg.norm(p - pocket, axis=1)

    print(f"口袋中心 = {pocket}，粒子数 K={K}，温度 T={T}")
    print("每轮：打分 → 权重 w∝exp(s/T) → ESS=1/Σw² → ESS 太低就按权重重采样（克隆好的/淘汰差的）→ 扰动前进\n")
    for rnd in range(4):
        s = score(particles)
        w = np.exp((s - s.max()) / T); w /= w.sum()
        ess = 1.0 / np.sum(w ** 2)                  # 有效样本量
        print(f"--- 第 {rnd} 轮 ---")
        print("  打分 s   :", s.round(2))
        print("  权重 w   :", w.round(3), " (高权重=更贴口袋)")
        print(f"  ESS = {ess:.2f}  (阈值 K/2 = {K/2:.1f})", end="")
        if rnd == 3:
            print("  → 结束")
            break
        if ess < K / 2:
            # 系统重采样：高权重粒子被克隆，低权重被淘汰
            pos = (rng.random() + np.arange(K)) / K
            cdf = np.cumsum(w)
            idx = np.searchsorted(cdf, pos)
            print("  → ESS 偏低，重采样！存活索引:", idx.tolist(),
                  " 淘汰:", sorted(set(range(K)) - set(idx.tolist())))
            particles = particles[idx]
        else:
            print("  → ESS 足够，不重采样")
        # 传播（propagate）：对每个粒子做一步小扰动（相当于 generator.refine）
        particles = particles + rng.normal(scale=0.6, size=particles.shape)
    best = particles[np.argmax(score(particles))]
    print(f"\n最终胜出粒子位置 {best.round(2)}，离口袋 {np.linalg.norm(best-pocket):.2f} Å")
    print("算力被集中到“贴口袋”的高分分支上——这就是测试时算力扩展(test-time scaling)的增益来源。")


# ==========================================================================
def main():
    print("abproto 矩阵教学演示：6 残基迷你复合物（抗原 E,W,K + 待设计 CDR N,Y,S）")
    s1 = stage1_pdb_to_matrices()
    s2 = stage2_embedding(s1)
    s3 = stage3_attention(s2)
    s4 = stage4_heads(s3)
    stage5_training(s4)
    stage6_twisted_smc()
    banner("全流程回顾")
    print("PDB坐标 → [6×3]矩阵 →(4股查表相加) H0[6×4] →(距离偏置注意力) H1[6×4]")
    print("        →(速度头[4×3]/序列头[4×8]) 预测 →(梯度反馈) 训练模型")
    print("        →(新抗原) Twisted SMC + ESS 动态克隆/淘汰 → 高分结合构象")


if __name__ == "__main__":
    main()
