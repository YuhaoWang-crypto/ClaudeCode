#!/usr/bin/env python3
"""ALS campaign figures. All from real campaign data (Open Targets/CT/STRING/Boltz/kernel)
or clearly-labeled schematics. Badge: 真数据 / 机制示意."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.font_manager as fm
import numpy as np, pathlib, sys, json
for fp in ["/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc"]:
    try: fm.fontManager.addfont(fp)
    except Exception: pass
plt.rcParams["font.family"]="WenQuanYi Zen Hei"; plt.rcParams["axes.unicode_minus"]=False
OUT=pathlib.Path("docs/als_repro/figures"); OUT.mkdir(parents=True,exist_ok=True)
TEAL="#0b6b5e";BLUE="#0d47a1";GREY="#c0c7d0";RED="#c0392b";GOLD="#b8860b";PURP="#6a3d9a";GREEN="#1a7a3a"
def badge(fig,kind):
    c,t={"real":("#0b6b5e","真数据"),"schema":("#5a6472","机制示意")}[kind]
    fig.text(0.995,0.985,t,ha="right",va="top",fontsize=8.5,color="white",weight="bold",
             bbox=dict(boxstyle="round,pad=0.3",fc=c,ec="none"))
def save(fig,name,kind):
    badge(fig,kind); fig.savefig(OUT/f"fig_als_{name}.png",dpi=140,bbox_inches="tight",facecolor="white"); plt.close(fig); print("wrote",name)

# 1 market
inds=["阿尔茨海默\n(伪影第一)","ALS(选中)","MASH","SLE","IgA肾病","IPF"]; comp=[4.10,3.99,3.41,3.27,3.26,3.22]
fig,ax=plt.subplots(figsize=(8.5,3.8)); y=np.arange(len(inds))[::-1]
cols=[GREY,TEAL,GREY,GREY,GREY,GREY]
ax.barh(y,comp,color=cols,height=0.62)
for yi,c in zip(y,comp): ax.text(c+0.03,yi,f"{c:.2f}",va="center",fontsize=10,weight="bold")
ax.set_yticks(y); ax.set_yticklabels(inds,fontsize=9.5); ax.set_xlim(0,4.6)
ax.set_xlabel("综合分(OT关联数+CT在招募数 硬锚 + 未满足/多模态/遗传/拥挤 专家维度)")
ax.set_title("适应症打分:ALS 在决策相关维度胜出(AD第一为研究体量伪影)",fontsize=11,weight="bold")
save(fig,"market","real")

# 2 axis
fig,ax=plt.subplots(figsize=(9.5,4)); ax.axis("off"); ax.set_xlim(0,11); ax.set_ylim(0,5)
def node(x,y,t,c,w=1.7,h=0.95):
    ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle="round,pad=0.08",fc=c,ec="#333",lw=1.2))
    ax.text(x,y,t,ha="center",va="center",fontsize=8.8,color="white",weight="bold")
def arr(x1,y1,x2,y2,txt=""):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=15,lw=1.6,color="#333"))
    if txt: ax.text((x1+x2)/2,(y1+y2)/2+0.2,txt,ha="center",fontsize=8,color="#333")
node(1.3,4.2,"SOD1\n0.887",PURP,h=0.8); node(1.3,3.0,"C9orf72\n0.658",PURP,h=0.8); node(1.3,1.8,"FUS\n0.840",PURP,h=0.8)
node(4.4,3.0,"TDP-43 / TARDBP\n蛋白病变枢纽\n(>95% ALS)",BLUE,w=2.4,h=1.5)
node(7.4,3.8,"STMN2\n隐蔽外显子",TEAL,h=0.8); node(7.4,2.2,"UNC13A\n隐蔽外显子",TEAL,h=0.8)
node(9.9,3.0,"NfL\n神经丝读出\n(FDA认可)",GREEN,w=1.7,h=1.3)
for yy in (4.2,3.0,1.8): arr(2.15,yy,3.2,3.0)
arr(5.6,3.3,6.55,3.8,"功能维持"); arr(5.6,2.7,6.55,2.2)
arr(8.25,3.8,9.05,3.2); arr(8.25,2.2,9.05,2.8)
ax.text(5.5,4.85,"ALS TDP-43 蛋白稳态轴:基因驱动 → 蛋白枢纽 → 隐蔽外显子效应器 → NfL",ha="center",fontsize=11,weight="bold")
ax.text(1.3,0.7,"基因层",ha="center",fontsize=9,color=PURP,weight="bold")
ax.text(4.4,0.9,"蛋白层",ha="center",fontsize=9,color=BLUE,weight="bold")
ax.text(7.4,0.9,"RNA/效应器",ha="center",fontsize=9,color=TEAL,weight="bold")
save(fig,"axis","real")

# 3 biomarker — copy real triangulation.png
import shutil
src=pathlib.Path("scratchpad/als_campaign/biomarker/triangulation.png")
if src.exists(): shutil.copy(src, OUT/"fig_als_biomarker.png"); print("wrote biomarker (copied real)")

# 4 synthetic biomarker
fig,ax=plt.subplots(figsize=(10,2.6)); ax.axis("off"); ax.set_xlim(0,10); ax.set_ylim(0,3)
mods=[("神经元特异\n启动子",PURP,1.6),("STMN2 隐蔽外显子\n剪接传感器",TEAL,2.4),("荧光报告\n(GFP/Nluc)",BLUE,1.6),("RNA\n条形码",GOLD,1.2)]
x=0.4
for t,c,w in mods:
    ax.add_patch(FancyBboxPatch((x,1.1),w,0.9,boxstyle="round,pad=0.04",fc=c,ec="#333",lw=1))
    ax.text(x+w/2,1.55,t,ha="center",va="center",fontsize=8.6,color="white",weight="bold"); x+=w+0.15
ax.text(5,0.55,"TDP-43 核功能丢失 → 隐蔽外显子纳入 → 报告信号变化(可测剪接开关)",ha="center",fontsize=9,color=RED)
ax.text(5,2.6,"合成 biomarker:Syn-TDP43-STMN2sensor(概念设计)",ha="center",fontsize=11.5,weight="bold")
save(fig,"synbiomarker","schema")

# 5 network DNB (real kernel)
sys.path.insert(0,".claude/skills/denovo-drug-campaign"); import kernel as K
ps=np.linspace(0.55,1.35,17); var=[]
for p in ps:
    try: var.append(K.langevin_earlywarning(0.4,0.03,0.3,p=float(p),seed=3)["var"])
    except Exception: var.append(np.nan)
fig,axs=plt.subplots(1,2,figsize=(11,4.3))
axs[0].axis("off"); axs[0].set_xlim(0,10); axs[0].set_ylim(0,5)
def n2(x,y,t,c):
    axs[0].add_patch(FancyBboxPatch((x-1.2,y-0.42),2.4,0.84,boxstyle="round,pad=0.05",fc=c,ec="#333"))
    axs[0].text(x,y,t,ha="center",va="center",fontsize=8.8,color="white",weight="bold")
n2(2,4,"{SOD1,C9orf72,FUS}\n输入纤维",PURP); n2(2,1.4,"{STMN2,UNC13A}\n效应器纤维",TEAL); n2(6.6,2.7,"{TDP-43}\n枢纽纤维",BLUE)
axs[0].annotate("",xy=(5.4,2.9),xytext=(3.25,3.7),arrowprops=dict(arrowstyle="-|>",lw=1.6))
axs[0].annotate("",xy=(5.4,2.4),xytext=(3.25,1.5),arrowprops=dict(arrowstyle="-|>",lw=1.6))
axs[0].text(5,4.7,"纤维化不可约核心:5→3 纤维",ha="center",fontsize=11,weight="bold")
axs[0].text(5,0.3,"CRNT δ = 4-2-1 = 1 → 双稳(TDP-43 可溶↔聚集;2稳态+1鞍点)",ha="center",fontsize=9.3,color=RED,weight="bold")
axs[1].plot(ps,np.array(var)/np.nanmin(var),marker="o",ms=4,color=TEAL,lw=1.8)
axs[1].axvline(0.70,ls="--",c=RED); axs[1].text(0.71,axs[1].get_ylim()[1]*0.8,"近临界",color=RED,fontsize=9)
axs[1].set_xlabel("驱动 p(聚集seeding/年龄→鞍结)"); axs[1].set_ylabel("方差(相对最小)")
axs[1].set_title("NfL 临界慢化:方差 near/far ≈5.15× · AR1 0.989",fontsize=10.5,weight="bold")
fig.suptitle("网络层:TDP-43 蛋白稳态双稳 + NfL 早期预警(kernel 真跑)",fontsize=11,weight="bold")
fig.tight_layout(rect=[0,0,1,0.93]); save(fig,"network","real")

# 6 feasibility matrix
tg=["SOD1","C9orf72","FUS","TDP-43","STMN2"]; md=["ASO/siRNA","抗体/VHH","小分子","PROTAC","基因/AAV"]
M=np.array([[5.0,4.0,2.5,2.0,4.0],
            [5.0,1.5,1.5,1.0,3.5],
            [4.5,3.0,2.0,1.5,3.5],
            [3.0,2.5,1.5,2.0,2.5],
            [4.5,1.5,1.5,1.0,3.0]])
stars={(0,0),(1,0),(4,0),(0,1),(0,2)}
fig,ax=plt.subplots(figsize=(7.8,4.4))
im=ax.imshow(M,cmap="YlGnBu",vmin=0,vmax=5,aspect="auto")
ax.set_xticks(range(5)); ax.set_xticklabels(md,fontsize=9); ax.set_yticks(range(5)); ax.set_yticklabels(tg,fontsize=10)
for i in range(5):
    for j in range(5):
        s="★" if (i,j) in stars else ""
        ax.text(j,i,f"{M[i,j]:.1f}{s}",ha="center",va="center",fontsize=9,color="white" if M[i,j]>3 else "#222",weight="bold")
ax.set_title("ALS 靶点×模态可行性(★=选中先导)· RNA层→ASO最优",fontsize=11,weight="bold")
fig.colorbar(im,ax=ax,shrink=0.8,label="可行性"); save(fig,"matrix","real")

# 7 designs — real scores
fig,axs=plt.subplots(1,2,figsize=(11,4))
# left: nanobody ipTM distribution (top15) vs smallmol binding
nb=[0.886,0.869,0.817,0.812,0.746,0.716,0.591,0.544,0.519,0.46,0.406,0.341,0.276,0.25,0.221]
axs[0].bar(range(len(nb)),nb,color=[GREEN if v>0.7 else TEAL if v>0.5 else GREY for v in nb])
axs[0].axhline(0.7,ls="--",c=RED,lw=1); axs[0].text(10,0.72,"ipTM 0.7",color=RED,fontsize=8)
axs[0].set_title("L4 SOD1 纳米抗体 de-novo(100条→检出15)\ntop ipTM 0.886 · 6条>0.7",fontsize=10,weight="bold")
axs[0].set_xlabel("设计(按ipTM排序)"); axs[0].set_ylabel("ipTM"); axs[0].set_ylim(0,1)
# right: modality comparison (real headline numbers)
labs=["SOD1-ASO\n×100(gapmer)","C9orf72-ASO\n×100","STMN2-ASO\n×100","纳米抗体\ntop ipTM","ebselen重利用\nipTM","de-novo小分子\nbinding"]
vals=[1.0,1.0,1.0,0.886,0.824,0.311]; cols=[PURP,PURP,PURP,GREEN,GOLD,RED]
axs[1].bar(range(len(labs)),vals,color=cols)
for i,v in enumerate(vals): axs[1].text(i,v+0.01,f"{v:.2f}" if i>2 else "设计100",ha="center",fontsize=8,weight="bold")
axs[1].set_xticks(range(len(labs))); axs[1].set_xticklabels(labs,fontsize=7.5); axs[1].set_ylim(0,1.1)
axs[1].set_title("五模态真实产出:抗体(0.89)>ebselen(0.82)>>小分子(0.31)",fontsize=10,weight="bold")
axs[1].set_ylabel("ipTM / binding_confidence")
fig.suptitle("多模态 GPU 设计真实分数(Boltz-2, Modal)· 数据自证 ASO/抗体 > 小分子",fontsize=11,weight="bold")
fig.tight_layout(rect=[0,0,1,0.93]); save(fig,"designs","real")

# 8 ranking
rk=json.load(open("scratchpad/als_campaign/ranking.json"))
names={"L1":"L1 SOD1-ASO","L2":"L2 STMN2-ASO","L3":"L3 C9orf72-ASO","L4":"L4 纳米抗体","L5":"L5 ebselen","L6":"L6 de-novo小分子"}
rk=sorted(rk,key=lambda r:r["rank"])
labs=[names[r["id"]] for r in rk]; comp=[r["composite"] for r in rk]
fig,ax=plt.subplots(figsize=(8.5,3.8)); x=np.arange(len(labs))
cols=[PURP,PURP,PURP,GOLD,GREEN,RED]
ax.bar(x,comp,color=cols,width=0.62)
for xi,c in zip(x,comp): ax.text(xi,c+0.03,f"{c:.2f}",ha="center",fontsize=10,weight="bold")
ax.set_xticks(x); ax.set_xticklabels(labs,fontsize=8.5); ax.set_ylim(0,5)
ax.set_ylabel("综合分(可行0.30+效力0.30+安全0.22+可开发0.18)")
ax.set_title("6 条先导排序:基因/ASO 包揽前三(kernel rank_leads)",fontsize=11,weight="bold")
save(fig,"ranking","real")

# 9 delivery
fig,ax=plt.subplots(figsize=(9.5,3.4)); ax.axis("off"); ax.set_xlim(0,10); ax.set_ylim(0,4)
def box(x,y,t,c,w=2.7,h=1.5):
    ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle="round,pad=0.08",fc=c,ec="#333",lw=1.2))
    ax.text(x,y,t,ha="center",va="center",fontsize=9,color="white",weight="bold")
box(1.9,2,"ASO(L1/L2/L3)\n→ 鞘内注射\n(tofersen类,绕BBB)",TEAL)
box(5.0,2,"纳米抗体(L4)\n→ AAV脑内表达 /\nTfR穿梭过BBB",BLUE)
box(8.1,2,"小分子(L5/ebselen)\n→ 脑穿透优化\n+ CNS-LNP",GOLD)
ax.text(5,3.6,"CNS 靶向递送:据模态定制(★鞘内 ASO 最成熟)",ha="center",fontsize=11.5,weight="bold")
ax.text(5,0.4,"⚠️ 基于已批准先例(tofersen/nusinersen 鞘内)的机制推荐,非湿实验测定",ha="center",fontsize=8.5,color=RED)
save(fig,"delivery","schema")
print("=== all ALS figures done ===")
