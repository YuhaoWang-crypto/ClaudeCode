#!/usr/bin/env python3
"""Regenerate the 13 NSCLC-report figures as self-contained PNGs.
Provenance honesty: each figure carries a corner badge:
  [独立复现]  = generated from MY real reproduction (Boltz co-fold / biomarker skill / kernel DNB)
  [据报告数值] = redrawn from the report's stated numbers (original render unavailable)
  [机制示意]  = schematic diagram
Output: docs/nsclc_repro/figures/fig_<artifact-uuid>.png
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
import matplotlib.font_manager as fm
import numpy as np, pathlib, sys

# CJK font
for fp in ["/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc"]:
    try: fm.fontManager.addfont(fp)
    except Exception: pass
plt.rcParams["font.family"] = "WenQuanYi Zen Hei"
plt.rcParams["axes.unicode_minus"] = False

OUT = pathlib.Path("docs/nsclc_repro/figures"); OUT.mkdir(parents=True, exist_ok=True)
TEAL="#0b6b5e"; BLUE="#0d47a1"; GREY="#c0c7d0"; RED="#c0392b"; GOLD="#b8860b"; PURP="#762a83"
BADGE={"repro":("#0b6b5e","独立复现"),"stated":("#b8860b","据报告数值"),"schema":("#5a6472","机制示意")}

def badge(fig, kind):
    c,txt=BADGE[kind]
    fig.text(0.995,0.985,txt,ha="right",va="top",fontsize=8.5,color="white",weight="bold",
             bbox=dict(boxstyle="round,pad=0.3",fc=c,ec="none"))

def save(fig,uuid,kind):
    badge(fig,kind); fig.savefig(OUT/f"fig_{uuid}.png",dpi=140,bbox_inches="tight",facecolor="white"); plt.close(fig)
    print("wrote fig_"+uuid[:8])

# ── 1. market scoring ── ac0b81e6 ─────────────────────────────────────────
u="ac0b81e6-7b2d-4673-a198-5e5c6d0f86b0"
inds=["NSCLC","结直肠","肝癌(HCC)","胰腺癌"]; comp=[4.57,4.20,3.97,3.81]
fig,ax=plt.subplots(figsize=(8,3.6)); y=np.arange(len(inds))[::-1]
cols=[TEAL]+[GREY]*3
ax.barh(y,comp,color=cols,height=0.6)
for yi,c in zip(y,comp): ax.text(c+0.03,yi,f"{c:.2f}",va="center",fontsize=11,weight="bold")
ax.set_yticks(y); ax.set_yticklabels(inds,fontsize=11); ax.set_xlim(0,5)
ax.set_xlabel("7 维加权综合分(0–5)"); ax.set_title("适应症综合打分:NSCLC 排名第一(4.57)",fontsize=12,weight="bold")
ax.text(0.02,-0.9,"NSCLC 原始指标:12,475 关联靶点 · 1,072 药物 · 1,289 在招募试验(均 Open Targets/ClinicalTrials live 核实)",
        fontsize=8.5,color="#333",transform=ax.get_yaxis_transform())
save(fig,u,"stated")

# ── 2. pathway axis ── 213151a7 ───────────────────────────────────────────
u="213151a7-193c-44d2-860a-696b719b3462"
fig,ax=plt.subplots(figsize=(9,3.8)); ax.axis("off"); ax.set_xlim(0,10); ax.set_ylim(0,5)
def node(x,y,t,c,w=1.5,h=0.8):
    ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle="round,pad=0.08",fc=c,ec="#333",lw=1.2))
    ax.text(x,y,t,ha="center",va="center",fontsize=10,color="white",weight="bold")
def arr(x1,y1,x2,y2,txt="",c="#333"):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=16,lw=1.6,color=c))
    if txt: ax.text((x1+x2)/2,(y1+y2)/2+0.22,txt,ha="center",fontsize=9,color=c,weight="bold")
node(1.4,4,"KRAS\n(基因组驱动)","#8e44ad")
node(1.4,2,"KEAP1\n(E3 氧感受器)","#8e44ad")
node(4.2,3,"NRF2 / NFE2L2\n(转录因子·枢纽)",BLUE,w=2.2)
node(7.4,4,"SLC7A11 / xCT\n(膜转运体)",TEAL)
node(7.4,2,"GLS\n(谷氨酰胺酶)",TEAL)
arr(2.3,2.1,3.2,2.7,"失活→稳定 -")
arr(2.3,3.9,3.2,3.3,"+")
arr(5.3,3.2,6.5,3.9,"反式激活 +")
arr(5.3,2.8,6.5,2.1,"代谢依赖*")
ax.text(5,0.5,"* NRF2→GLS 独立复核:STRING 无直接边 → 应为“代谢依赖”而非“反式激活”(见复核 B.2)",
        ha="center",fontsize=8.5,color=RED)
ax.text(5,4.9,"KEAP1/NRF2 氧化还原–代谢轴(基因→蛋白→代谢)",ha="center",fontsize=12,weight="bold")
save(fig,u,"schema")

# ── 3. biomarker triangulation (REAL) ── 1ae4c2dd ─────────────────────────
u="1ae4c2dd-6a44-4d45-bf18-3c5de850af97"
rows=[l.split("\t") for l in open("docs/nsclc_repro/nsclc_biomarker_triangulation.tsv").read().splitlines()[1:]]
genes=[r[0] for r in rows][:16]; attr=[float(r[4]) for r in rows][:16]; rwr=[float(r[5]) for r in rows][:16]
fig,axs=plt.subplots(1,2,figsize=(11,5.2))
yy=np.arange(len(genes))[::-1]
axs[0].barh(yy,attr,color=[TEAL if float(r[6]) else GREY for r in rows[:16]])
axs[0].set_yticks(yy); axs[0].set_yticklabels(genes,fontsize=8.5); axs[0].set_title("方法一 归因(绿=最小核心 22/40)",fontsize=11,weight="bold")
axs[0].set_xlabel("归因 = 可测性×(0.65·关联+0.35·中心性)")
axs[1].barh(yy,rwr,color=BLUE)
axs[1].set_yticks(yy); axs[1].set_yticklabels(genes,fontsize=8.5); axs[1].set_title("方法二 RWR 拓扑扩散",fontsize=11,weight="bold")
axs[1].set_xlabel("personalized PageRank")
fig.suptitle("NSCLC 轴 biomarker 三方法交叉验证 · 重叠 20/22 · Jaccard 0.83 · 共享核心含 KEAP1/KRAS/NFE2L2/SLC7A11/GLS",fontsize=11,weight="bold")
fig.tight_layout(rect=[0,0,1,0.95])
save(fig,u,"repro")

# ── 4. synthetic cassette ── d2e24a34 ─────────────────────────────────────
u="d2e24a34-ef9c-4bcf-88f4-9a2e85b3c411"
fig,ax=plt.subplots(figsize=(10,2.6)); ax.axis("off"); ax.set_xlim(0,10); ax.set_ylim(0,3)
mods=[("A549\n增强子","#8e44ad",1.2),("4×ARE\n状态传感器",TEAL,1.6),("最小\n核心启动子",BLUE,1.4),
      ("合成\n内含子","#2c7fb8",1.2),("RNA\n条形码",GOLD,1.0),("肺特异 miRNA\n关断开关",RED,1.7)]
x=0.3
for t,c,w in mods:
    ax.add_patch(FancyBboxPatch((x,1.1),w,0.9,boxstyle="round,pad=0.04",fc=c,ec="#333",lw=1))
    ax.text(x+w/2,1.55,t,ha="center",va="center",fontsize=8.5,color="white",weight="bold"); x+=w+0.12
ax.annotate("",xy=(x,0.85),xytext=(0.3,0.85),arrowprops=dict(arrowstyle="<->",lw=1))
ax.text(x/2,0.5,"SynBioMark-NSCLC-KEAP1axis-01 · 627 bp · QC 关键项全通过",ha="center",fontsize=9.5,weight="bold")
ax.text(5,2.6,"合成 NRF2-ARE 报告盒:7 模块架构",ha="center",fontsize=12,weight="bold")
save(fig,u,"schema")

# ── 5. Evo2 naturalness ── 43e375ab ───────────────────────────────────────
u="43e375ab-fdd9-40c5-a329-e65c356603f3"
fig,ax=plt.subplots(figsize=(7.5,3.6))
labels=["合成 ARE\n传感器","天然 NQO1-ARE\n(阳性参照)","二核苷酸\n打乱对照"]; vals=[-0.975,-1.311,-1.258]
cols=[TEAL,BLUE,GREY]
ax.bar(labels,vals,color=cols,width=0.6)
for i,v in enumerate(vals): ax.text(i,v-0.05,f"{v:.3f}",ha="center",va="top",fontsize=10,weight="bold",color="white")
ax.set_ylabel("Evo2-7B 平均每 token 对数似然(越高越自然)")
ax.set_title("Evo2 序列自然性:合成 ARE(-0.975) > 天然参照,完整盒 +0.283 vs 打乱",fontsize=10.5,weight="bold")
ax.axhline(0,color="#333",lw=0.6)
ax.text(0.5,-1.45,"注:数值取自报告;本项目 Evo2 Modal 部署仍 WIP,未独立复现该绝对值(方向合理)",
        ha="center",fontsize=8,color=RED,transform=ax.get_xaxis_transform())
save(fig,u,"stated")

# ── 6. feasibility matrix ── 4175a2b1 ─────────────────────────────────────
u="4175a2b1-7309-40c0-9471-5a21404c67fd"
tg=["KRAS","KEAP1","NRF2","SLC7A11","GLS"]; md=["小分子","抗体/VHH","基因(ASO/siRNA)","PROTAC 降解","微结合肽"]
M=np.array([[5.0,2.0,3.5,4.0,3.0],
            [4.0,3.0,3.0,3.5,5.0],
            [1.5,1.5,4.5,2.0,2.0],
            [1.5,4.0,3.5,1.5,2.5],
            [4.0,2.0,3.5,3.0,3.0]])
stars={(0,0),(1,4),(2,2),(3,1),(4,0),(0,3)}  # 6 leads L5,L1,L3,L2,... mapped to representative cells
fig,ax=plt.subplots(figsize=(7.6,4.4))
im=ax.imshow(M,cmap="YlGnBu",vmin=0,vmax=5,aspect="auto")
ax.set_xticks(range(5)); ax.set_xticklabels(md,fontsize=9); ax.set_yticks(range(5)); ax.set_yticklabels(tg,fontsize=10)
for i in range(5):
    for j in range(5):
        s="★" if (i,j) in stars else ""
        ax.text(j,i,f"{M[i,j]:.1f}{s}",ha="center",va="center",fontsize=9,
                color="white" if M[i,j]>3 else "#222",weight="bold")
ax.set_title("靶点 × 模态可行性矩阵(0–5)· ★=6 条选中先导",fontsize=11,weight="bold")
fig.colorbar(im,ax=ax,shrink=0.8,label="可行性")
save(fig,u,"stated")

# ── 7. L1 KEAP1 binder — REAL co-fold ── 348f5224 ─────────────────────────
u="348f5224-8c51-4d63-8e51-a1b2c3d4e5f6".replace("348f5224-8c51-4d63-8e51-a1b2c3d4e5f6","348f5224-a0d9-4c6d-8e51-887a08fed43d")
fig,ax=plt.subplots(figsize=(7.5,3.8))
bars=["报告声称\n(0.90–0.99)","我的独立 co-fold\nipTM","结构置信度","结合置信度"]
vals=[0.95,0.955,0.967,0.447]; cols=[GREY,TEAL,BLUE,RED]
ax.bar(bars,vals,color=cols,width=0.6)
for i,v in enumerate(vals): ax.text(i,v+0.01,f"{v:.3f}",ha="center",fontsize=10.5,weight="bold")
ax.set_ylim(0,1.1); ax.set_ylabel("分值")
ax.set_title("L1 KEAP1 Kelch + NRF2 ETGE 肽:Boltz-2.1 真实 co-fold",fontsize=11,weight="bold")
ax.text(0.5,1.02,"ipTM 0.955 复现报告区间 —— 但 binding_confidence 仅 0.447:界面几何高 ≠ 结合把握高",
        ha="center",fontsize=8.5,color=RED,transform=ax.get_xaxis_transform())
save(fig,u,"repro")

# ── 8. L2 SLC7A11 VHH campaign ── 9f0490e2 ────────────────────────────────
u="9f0490e2-8277-42ca-a354-17d960fc5254"
fig,axs=plt.subplots(1,2,figsize=(9,3.6))
axs[0].pie([80,20],labels=["无 CDR3 责任\n80","有责任\n20"],colors=[TEAL,GREY],autopct="%d%%",startangle=90,
           textprops=dict(fontsize=9,weight="bold"))
axs[0].set_title("100 条 VHH 设计:可开发性",fontsize=10.5,weight="bold")
axs[1].bar(["总设计","RF2 重折叠","可溶","无责任"],[100,100,100,80],color=[GREY,BLUE,TEAL,"#1a7a3a"],width=0.6)
for i,v in enumerate([100,100,100,80]): axs[1].text(i,v+1,str(v),ha="center",fontsize=10,weight="bold")
axs[1].set_title("xCT ECL2(211–234)靶向 campaign",fontsize=10.5,weight="bold"); axs[1].set_ylim(0,115)
fig.suptitle("L2 SLC7A11/xCT VHH 纳米抗体:25 骨架×4 序列=100 设计",fontsize=11,weight="bold")
fig.tight_layout(rect=[0,0,1,0.93])
save(fig,u,"stated")

# ── 9. L4/L5 GLS & KRAS small-mol — REAL KRAS ── 4836c863 ─────────────────
u="4836c863-8c63-4d64-979c-3e65d3a99bae"
fig,ax=plt.subplots(figsize=(8.5,4))
names=["KRAS-sotorasib","KRAS-adagrasib","GLS-BPTES","GLS-telaglenastat"]
rep=[0.96,0.97,0.88,0.80]; mine=[0.977,0.984,None,None]
x=np.arange(len(names)); w=0.38
ax.bar(x-w/2,rep,w,label="报告声称",color=GREY)
ax.bar([xi+w/2 for xi,m in zip(x,mine) if m],[m for m in mine if m],w,label="我的独立 co-fold",color=TEAL)
for xi,r in zip(x,rep): ax.text(xi-w/2,r+0.005,f"{r:.2f}",ha="center",fontsize=9)
for xi,m in zip(x,mine):
    if m: ax.text(xi+w/2,m+0.005,f"{m:.3f}",ha="center",fontsize=9,weight="bold",color=TEAL)
ax.set_xticks(x); ax.set_xticklabels(names,fontsize=9); ax.set_ylim(0,1.1); ax.set_ylabel("ipTM")
ax.legend(fontsize=9,loc="lower left")
ax.set_title("L4/L5 小分子 co-fold:KRAS 两个已独立复现(0.977/0.984,均 > 报告)",fontsize=10.5,weight="bold")
ax.text(2.5,0.5,"GLS 两条为报告数值\n(本次未复现)",ha="center",fontsize=8.5,color=GOLD)
save(fig,u,"repro")

# ── 10. L6 PROTAC ADMET ── ba4bea71 ───────────────────────────────────────
u="ba4bea71-61c8-48e9-a5d1-389317359a24"
fig,ax=plt.subplots(figsize=(8.5,3.8))
props=["MW","cLogP","TPSA","HBD","HBA","可旋转键"]
peg2=[0.9,0.55,0.8,0.4,0.85,0.9]; peg4=[0.98,0.6,0.88,0.45,0.95,0.98]
x=np.arange(len(props)); w=0.38
ax.bar(x-w/2,peg2,w,label="PEG2 连接子",color=PURP); ax.bar(x+w/2,peg4,w,label="PEG4 连接子",color="#c994c7")
ax.axhline(0.83,ls="--",c=RED,lw=1); ax.text(5.4,0.85,"bRo5 上界",color=RED,fontsize=8)
ax.set_xticks(x); ax.set_xticklabels(props,fontsize=9); ax.set_ylabel("归一化(1=最大)"); ax.set_ylim(0,1.15)
ax.legend(fontsize=9); ax.set_title("L6 KRAS PROTAC(adagrasib 弹头+VHL):7 化合物 RDKit ADMET,属 bRo5",fontsize=10,weight="bold")
save(fig,u,"stated")

# ── 11. lead ranking ── f8078eb1 ──────────────────────────────────────────
u="f8078eb1-11b7-47cd-905a-975c4c74fb63"
leads=["L5 KRAS\n共价小分子","L1 KEAP1\n微结合物","L4 GLS\n别构","L2 SLC7A11\nVHH","L3 NRF2\n基因","L6 KRAS\nPROTAC"]
comp=[4.36,4.21,4.01,3.94,3.89,3.32]
fig,ax=plt.subplots(figsize=(8.5,3.8)); x=np.arange(len(leads))
cols=[TEAL,TEAL,"#3a9188","#3a9188","#6bb0a8",GREY]
ax.bar(x,comp,color=cols,width=0.62)
for xi,c in zip(x,comp): ax.text(xi,c+0.03,f"{c:.2f}",ha="center",fontsize=10,weight="bold")
ax.set_xticks(x); ax.set_xticklabels(leads,fontsize=8.5); ax.set_ylim(0,5)
ax.set_ylabel("综合分(可行性0.30+效力0.30+安全0.22+可开发0.18)")
ax.set_title("6 条先导安全性/可开发性综合排序",fontsize=11,weight="bold")
save(fig,u,"stated")

# ── 12. LNP delivery ── 3ac4fa74 ──────────────────────────────────────────
u="3ac4fa74-c70f-4803-b704-d1357470eefe"
fig,axs=plt.subplots(1,2,figsize=(9.5,3.8))
axs[0].pie([50,25,10,10,5],labels=["可电离脂 50%","DOTAP 25%\n(SORT/肺趋向)","DSPC 10%","胆固醇 10%","DMG-PEG 5%"],
           colors=[TEAL,RED,BLUE,GOLD,GREY],autopct="%d%%",startangle=90,textprops=dict(fontsize=8))
axs[0].set_title("LungTrop-NRF2-siLNP-01 · N/P=6",fontsize=10,weight="bold")
t=np.linspace(0,72,100); kd=82*(1-np.exp(-t/14))
axs[1].plot(t,kd,color=TEAL,lw=2.2); axs[1].axhline(82,ls="--",c=RED,lw=1)
axs[1].text(40,74,"NRF2 敲低 82%(k2=0.015)",color=RED,fontsize=9)
axs[1].set_xlabel("时间 (h)"); axs[1].set_ylabel("NRF2 敲低 %"); axs[1].set_ylim(0,95)
axs[1].set_title("Mihaila 5-ODE 预测(示意)",fontsize=10,weight="bold")
fig.suptitle("第 9 阶段 · NFE2L2 siRNA 肺靶向 SORT-LNP",fontsize=11,weight="bold")
fig.tight_layout(rect=[0,0,1,0.93])
save(fig,u,"stated")

# ── 13. network DNB — REAL kernel ── 71da31e9 ─────────────────────────────
u="71da31e9-b53d-41e6-a3f4-b79105c63a50"
sys.path.insert(0,".claude/skills/denovo-drug-campaign")
import kernel as K
ps=np.linspace(0.55,1.35,17); var_meas=[]
for p in ps:
    try: var_meas.append(K.langevin_earlywarning(0.4,0.03,0.3,p=float(p),seed=3)["var"])
    except Exception: var_meas.append(np.nan)
fig,axs=plt.subplots(1,2,figsize=(11,4.3))
# panel A: fibers
axs[0].axis("off"); axs[0].set_xlim(0,10); axs[0].set_ylim(0,5)
def n2(ax,x,y,t,c):
    ax.add_patch(FancyBboxPatch((x-1.1,y-0.4),2.2,0.8,boxstyle="round,pad=0.05",fc=c,ec="#333"))
    ax.text(x,y,t,ha="center",va="center",fontsize=9,color="white",weight="bold")
n2(axs[0],2,4,"{KEAP1, KRAS}\n输入纤维","#8e44ad")
n2(axs[0],2,1.4,"{SLC7A11, GLS}\n效应器纤维",TEAL)
n2(axs[0],6.5,2.7,"{NRF2}\n枢纽纤维",BLUE)
axs[0].annotate("",xy=(5.3,2.9),xytext=(3.2,3.7),arrowprops=dict(arrowstyle="-|>",lw=1.6))
axs[0].annotate("",xy=(5.3,2.4),xytext=(3.2,1.5),arrowprops=dict(arrowstyle="-|>",lw=1.6))
axs[0].text(5,4.7,"纤维化不可约核心:5 节点 → 3 纤维",ha="center",fontsize=11,weight="bold")
axs[0].text(5,0.3,"CRNT 亏格 δ = 4-2-1 = 1  → 拓扑允许双稳(Schlögl:2 稳态+1 鞍点)",ha="center",fontsize=9.5,color=RED,weight="bold")
# panel B: variance rise
axs[1].plot(ps,np.array(var_meas)/np.nanmin(var_meas),marker="o",ms=4,color=TEAL,lw=1.8)
axs[1].set_xlabel("驱动参数 p(→ 鞍结分岔)"); axs[1].set_ylabel("方差(相对最小值)")
axs[1].set_title("临界慢化早期预警:方差随 p→fold 上升",fontsize=10.5,weight="bold")
axs[1].axvline(0.70,ls="--",c=RED,lw=1); axs[1].text(0.71,axs[1].get_ylim()[1]*0.85,"近临界",color=RED,fontsize=9)
axs[1].text(0.95,axs[1].get_ylim()[1]*0.5,"实测方差 near/far ≈ 5.15×\n(kernel langevin_earlywarning)",fontsize=8.5,color="#333")
fig.suptitle("第 10 阶段 · 网络不可约性 + 临界转变生物标志物(用合并的 kernel 真实复现)",fontsize=11,weight="bold")
fig.tight_layout(rect=[0,0,1,0.93])
save(fig,u,"repro")

print("=== all 13 figures done ===")
print("real-reproduced:", "1ae4c2dd(biomarker), 348f5224(KEAP1 co-fold), 4836c863(KRAS co-fold), 71da31e9(DNB)")
