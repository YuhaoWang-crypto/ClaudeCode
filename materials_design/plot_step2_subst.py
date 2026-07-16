"""① 同族替换 QD(Pb→Sn) 与 SC(Nb₃Sn→V₃Si) 的性质变化, UMA vs 实验。"""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
OI={"skyblue":"#56B4E9","black":"#000000","vermillion":"#D55E00","blue":"#0072B2","green":"#009E73"}
plt.rcParams.update({"font.size":11,"axes.grid":True,"grid.alpha":0.25,"axes.axisbelow":True,
    "font.sans-serif":["WenQuanYi Zen Hei"],"axes.unicode_minus":False,"savefig.dpi":160})

fig,ax=plt.subplots(2,3,figsize=(15,8.5))
fig.suptitle("① 同族金属替换 → 性质变化 (UMA 模拟, 晶格另叠实验值)",fontsize=14,fontweight="bold")

def panel(a, labels, uma, ylab, ttl, exp=None, notes=None):
    x=np.arange(len(labels)); w=0.5
    bars=a.bar(x,uma,w,color=OI["skyblue"],edgecolor="white",label="UMA 模拟")
    if exp is not None:
        a.plot(x,exp,"D",color=OI["black"],ms=9,label="实验",zorder=5)
    a.set_xticks(x); a.set_xticklabels(labels,fontsize=10); a.set_ylabel(ylab); a.set_title(ttl)
    for xi,v in zip(x,uma): a.text(xi,v,f"{v:g}",ha="center",va="bottom",fontsize=9)
    if notes:
        for xi,t in zip(x,notes): a.text(xi,a.get_ylim()[1]*0.06,t,ha="center",fontsize=8.5,color=OI["vermillion"])
    a.legend(fontsize=8,loc="upper right")

# ---- QD: CsPbBr3 -> CsSnBr3 ----
qd=["CsPbBr₃","CsSnBr₃"]
panel(ax[0,0],qd,[6.004,5.900],"晶格常数 a (Å)","(QD-A) 晶格: Pb→Sn 收缩",exp=[5.874,5.804])
panel(ax[0,1],qd,[77.4,50.3],"最大声子 (cm⁻¹)","(QD-B) 声子软化")
panel(ax[0,2],qd,[18.6,19.5],"体模量 (GPa)","(QD-C) 体模量(卤化物钙钛矿都很软~19GPa)",
      notes=["带隙~2.3eV*","带隙~1.75eV*"])
ax[0,2].text(0.5,-4.2,"*带隙为实验值(UMA不直接给); Sn 更窄带隙",fontsize=8,ha="center",color="gray",transform=ax[0,2].transData)

# ---- SC: Nb3Sn -> V3Si ----
sc=["Nb₃Sn","V₃Si"]
panel(ax[1,0],sc,[5.336,4.694],"晶格常数 a (Å)","(SC-A) 晶格: Nb→V, Sn→Si 收缩",exp=[5.289,4.722])
panel(ax[1,1],sc,[178,322],"最大声子 (cm⁻¹)","(SC-B) 声子变硬(更轻原子)")
panel(ax[1,2],sc,[152.9,186.3],"体模量 (GPa)","(SC-C) 体模量升高(更刚)",
      notes=["Tc~18.3K*","Tc~17.1K*"])
ax[1,2].text(0.5,-40,"*Tc为实验值(UMA不直接给); 两者均为A15超导",fontsize=8,ha="center",color="gray",transform=ax[1,2].transData)

fig.tight_layout(rect=[0,0,1,0.95])
fig.savefig("/home/user/ClaudeCode/materials_design/step2_substitution.png",bbox_inches="tight")
print("saved step2_substitution.png")
