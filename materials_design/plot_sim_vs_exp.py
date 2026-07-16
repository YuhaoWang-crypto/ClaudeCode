"""
模拟(FairChem UMA / Brus) vs 真实实验数据 —— 同图对比。
四个面板: (A)量子点 (B)上转换基质声子 (C)超导声子 (D)晶格常数总 parity。
配色: Okabe-Ito 色盲安全调色板。
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ---- Okabe-Ito 色盲安全 ----
OI = {"black":"#000000","orange":"#E69F00","skyblue":"#56B4E9","green":"#009E73",
      "yellow":"#F0E442","blue":"#0072B2","vermillion":"#D55E00","purple":"#CC79A7"}
plt.rcParams.update({"font.size":11,"axes.grid":True,"grid.alpha":0.25,
    "axes.axisbelow":True,"figure.dpi":140,"savefig.dpi":160,
    "font.sans-serif":["WenQuanYi Zen Hei"],"axes.unicode_minus":False})

# ================== 数据 ==================
# --- (A) 量子点: Brus 模拟 + 实验尺寸曲线 ---
h=6.626e-34; me=9.109e-31; e=1.602e-19; eps0=8.854e-12
def brus(Eg,mes,mhs,epsr,d_nm):
    R=d_nm/2*1e-9
    return Eg+(h**2/(8*R**2))*(1/(mes*me)+1/(mhs*me))/e-1.8*e/(4*np.pi*eps0*epsr*R)
qd={
 "CdSe":  dict(par=(1.74,0.13,0.45,10.2),c=OI["vermillion"],
              exp_d=[2.0,2.5,3.0,3.5,4.0,4.5,5.0,6.0],
              exp_l=[458,495,528,552,573,590,605,628]),   # Yu&Peng 2003 / Jasieniak 2009
 "CsPbBr3":dict(par=(2.36,0.15,0.14,7.3),c=OI["green"],
              exp_d=[3.5,4.5,6.0,8.0,11.0],
              exp_l=[462,480,498,508,515]),                # Dong 2015 / Brennan 2017
 "PbS":   dict(par=(0.41,0.085,0.085,17.0),c=OI["blue"],
              exp_d=[2.5,3.0,4.0,5.0,6.0,7.0,8.0],
              exp_l=[830,920,1090,1230,1350,1470,1580]),   # Moreels 2009
}
# --- (B) 上转换基质: UMA Γ 最大声子 vs 实验最大(LO)声子 (cm^-1) ---
uc=[  # name, type, UMA计算, 实验
 ("BaF2","F",306,335),("SrF2","F",341,380),("CaF2","F",378,463),
 ("ZnO","O",511,591),("MgO","O",529,720),("TiO2","O",769,830),
]
# --- (C) 超导: 特征声子 (cm^-1) UMA vs 实验 ---
sc=[("Nb3Sn",178,200),("Nb",213,260),("MgB2 (E2g)",686,660)]  # MgB2 Raman E2g~620-660
# --- (D) 晶格常数 parity: UMA vs 实验 a(Å) ---
latt=[  # name, project, UMA a, 实验 a
 ("CsPbBr3","QD",6.004,5.874),("CdSe","QD",4.397,4.300),("PbS","QD",6.011,5.936),("Si","QD",5.453,5.431),
 ("CaF2","UC",5.515,5.463),("SrF2","UC",5.863,5.800),("BaF2","UC",6.284,6.200),
 ("MgO","UC",4.250,4.212),("ZnO","UC",3.291,3.250),("TiO2","UC",4.646,4.593),
 ("MgB2","SC",3.068,3.086),("Nb3Sn","SC",5.338,5.289),("Nb","SC",3.322,3.301),
]
proj_c={"QD":OI["vermillion"],"UC":OI["blue"],"SC":OI["orange"]}

# ================== 作图 ==================
fig,ax=plt.subplots(2,2,figsize=(13,10.5))
fig.suptitle("模拟预测 vs 真实实验数据 —— 三个材料项目的验证 (线/柱=模拟, 点=实验)",
             fontsize=14,fontweight="bold")

# (A) 量子点
a=ax[0,0]
dgrid=np.linspace(2.5,9,120)
for name,d in qd.items():
    a.plot(dgrid,[1240/brus(*d["par"],x) for x in dgrid],"-",color=d["c"],lw=2,label=f"{name} Brus模拟")
    a.plot(d["exp_d"],d["exp_l"],"o",color=d["c"],ms=7,mec="white",mew=0.8)
a.set_xlabel("量子点直径 (nm)"); a.set_ylabel("发射/首激子波长 (nm)")
a.set_title("(A) 量子点: 尺寸-发光 (Brus 模拟 vs 实验尺寸曲线)")
a.set_ylim(350,1750)
handles=[Line2D([],[],color=qd[k]["c"],lw=2,label=k) for k in qd]
handles+=[Line2D([],[],color="gray",lw=2,ls="-",label="模拟(线)"),
          Line2D([],[],color="gray",marker="o",ls="",label="实验(点)")]
a.legend(handles=handles,fontsize=8,loc="upper right",ncol=2)

# (B) 上转换基质声子
b=ax[0,1]
x=np.arange(len(uc)); w=0.38
comp=[u[2] for u in uc]; exp=[u[3] for u in uc]
b.bar(x-w/2,comp,w,color=OI["skyblue"],label="UMA 模拟 (Γ 最大声子)")
b.bar(x+w/2,exp,w,color=OI["black"],alpha=0.75,label="实验 (最高 LO 声子)")
b.set_xticks(x); b.set_xticklabels([f"{u[0]}\n({'氟化物' if u[1]=='F' else '氧化物'})" for u in uc],fontsize=8)
b.set_ylabel("最大声子能量 (cm$^{-1}$)")
b.set_title("(B) 上转换基质: 氟化物 < 氧化物 (声子越软猝灭越弱)")
b.axvline(2.5,color=OI["vermillion"],ls="--",lw=1.2)
b.text(1.0,760,"氟化物",color=OI["blue"],fontsize=9,ha="center")
b.text(4.0,760,"氧化物",color=OI["vermillion"],fontsize=9,ha="center")
b.legend(fontsize=8,loc="upper left")

# (C) 超导声子
c=ax[1,0]
x=np.arange(len(sc));
c.bar(x-w/2,[s[1] for s in sc],w,color=OI["skyblue"],label="UMA 模拟 (Γ 最大声子)")
c.bar(x+w/2,[s[2] for s in sc],w,color=OI["black"],alpha=0.75,label="实验 (拉曼/声子)")
c.set_xticks(x); c.set_xticklabels([s[0] for s in sc],fontsize=9)
c.set_ylabel("特征声子频率 (cm$^{-1}$)")
c.set_title("(C) 超导: MgB$_2$ 的 686 cm$^{-1}$ ≈ 驱动 39K 的 E$_{2g}$ 硼声子")
c.legend(fontsize=8,loc="upper left")

# (D) 晶格常数 parity
d=ax[1,1]
lo,hi=3.0,6.5
d.plot([lo,hi],[lo,hi],"--",color="gray",lw=1.3,label="y = x (完美一致)")
for proj in ["QD","UC","SC"]:
    xs=[L[3] for L in latt if L[1]==proj]; ys=[L[2] for L in latt if L[1]==proj]
    d.scatter(xs,ys,s=55,color=proj_c[proj],edgecolor="white",lw=0.7,label=proj,zorder=3)
# 统计
exp_all=np.array([L[3] for L in latt]); umd_all=np.array([L[2] for L in latt])
mae=np.mean(np.abs(umd_all-exp_all)); mape=np.mean(np.abs(umd_all-exp_all)/exp_all)*100
d.set_xlabel("实验晶格常数 a (Å)"); d.set_ylabel("UMA 模拟晶格常数 a (Å)")
d.set_title(f"(D) 晶格常数 parity (13 材料)  MAE={mae:.03f} Å, {mape:.1f}%")
d.set_xlim(lo,hi); d.set_ylim(lo,hi); d.legend(fontsize=8,loc="upper left")

plt.tight_layout(rect=[0,0,1,0.97])
plt.savefig("/home/user/ClaudeCode/materials_design/sim_vs_experiment.png",bbox_inches="tight")
print("已保存图: materials_design/sim_vs_experiment.png")
print(f"晶格 parity: MAE={mae:.3f} Å ({mape:.1f}%)")
print("上转换验证: 氟化物", [u[2] for u in uc[:3]], "< 氧化物", [u[2] for u in uc[3:]], "(UMA cm^-1)")
EOF_MARKER = None
