"""#3 QE DFPT: MgB2 Γ点声子(E2g驱动模) + Allen-Dynes Tc(λ) 演示。"""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
OI={"skyblue":"#56B4E9","black":"#000000","vermillion":"#D55E00","green":"#009E73","orange":"#E69F00","blue":"#0072B2"}
plt.rcParams.update({"font.size":11,"axes.grid":True,"grid.alpha":0.25,"axes.axisbelow":True,
    "font.sans-serif":["WenQuanYi Zen Hei"],"axes.unicode_minus":False,"savefig.dpi":160})

# DFPT (QE/LDA) MgB2 Γ点声子 (cm^-1) 与对称性
modes=[(25,"E1u"),(25,"E1u"),(327,"E1u"),(327,"E1u"),(389,"A2u"),(485,"E2g"),(485,"E2g"),(702,"B1g")]
# ω_log 代理(光学模对数平均) -> K
opt=np.array([327,327,389,485,485,702.0])
wlog_cm=np.exp(np.log(opt).mean()); wlog_K=wlog_cm*1.4388
def allen_dynes(lam, wlogK, mustar=0.13):
    return wlogK/1.2*np.exp(-1.04*(1+lam)/(lam-mustar*(1+0.62*lam)))

fig,ax=plt.subplots(1,2,figsize=(14,5.6))
fig.suptitle("#3 真实 DFPT (Quantum ESPRESSO) 电声耦合: MgB₂ 声子 + Allen-Dynes Tc",fontsize=13,fontweight="bold")

# 左: DFPT 声子模 (E2g 高亮), 对比 UMA 有限差分估计
a=ax[0]
xs=np.arange(len(modes))
cols=[OI["green"] if m[1]=="E2g" else OI["skyblue"] for m in modes]
a.bar(xs,[m[0] for m in modes],color=cols,edgecolor="white")
a.set_xticks(xs); a.set_xticklabels([m[1] for m in modes],fontsize=8,rotation=30)
a.axhline(686,color=OI["vermillion"],ls="--",lw=1.5); a.text(0.1,700,"UMA有限差分估计 686",color=OI["vermillion"],fontsize=8.5)
a.annotate("E₂g = 485 cm⁻¹\n驱动超导的模(Raman)",xy=(6,485),xytext=(2.2,560),fontsize=9,color=OI["green"],
           arrowprops=dict(arrowstyle="->",color=OI["green"]))
a.set_ylabel("声子频率 (cm⁻¹)"); a.set_title("(左) DFPT MgB₂ Γ点声子谱\nE₂g 严格值 485(谐波LDA); 实验~585(强非谐)")

# 右: Allen-Dynes Tc(λ) 曲线
b=ax[1]
lam=np.linspace(0.3,1.3,200)
for mu,ls in [(0.10,"-"),(0.13,"--"),(0.16,":")]:
    b.plot(lam,allen_dynes(lam,wlog_K,mu),ls,color=OI["blue"],lw=2,label=f"μ*={mu}")
b.axhline(39,color=OI["green"],ls="-",lw=1.5); b.text(0.32,41,"实验 Tc(MgB₂)=39 K",color=OI["green"],fontsize=9)
# 标出 λ≈0.87 (文献)
b.axvline(0.87,color=OI["orange"],ls=":",lw=1.3); b.text(0.88,5,"文献 λ≈0.87",color=OI["orange"],fontsize=8.5,rotation=90)
b.set_xlabel("电声耦合 λ"); b.set_ylabel("Allen-Dynes Tc (K)")
b.set_ylim(0,60)
b.set_title(f"(右) Allen-Dynes Tc(λ),  ω_log(DFPT代理)={wlog_K:.0f} K\n39K 对应 λ≈0.9-1.0, 与MgB₂强σ-E₂g耦合一致")
b.legend(fontsize=9,loc="upper left")
fig.tight_layout(rect=[0,0,1,0.93])
fig.savefig("/home/user/ClaudeCode/materials_design/step3b_qe_epc.png",bbox_inches="tight")
print(f"saved. ω_log代理={wlog_cm:.0f} cm⁻¹={wlog_K:.0f} K; Tc(λ=1.0,μ*=0.13)={allen_dynes(1.0,wlog_K):.1f}K")
