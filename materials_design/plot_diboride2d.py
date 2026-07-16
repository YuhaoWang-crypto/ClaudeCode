"""#2 多描述符筛选: 声子 × E_F处B-p(σ)占比。诚实版: 缩小到σ家族但仍需λ。"""
import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
OI={"skyblue":"#56B4E9","black":"#000000","vermillion":"#D55E00","green":"#009E73","orange":"#E69F00"}
plt.rcParams.update({"font.size":11,"axes.grid":True,"grid.alpha":0.25,"axes.axisbelow":True,
    "font.sans-serif":["WenQuanYi Zen Hei"],"axes.unicode_minus":False,"savefig.dpi":160})
d={r['name'].split('|')[1].split('_')[0]:r for r in json.load(open("/home/user/diboride_dft.json"))}
ph={r['name'].split('|')[1].split('_')[0]:r['phonon_max_cm-1'] for r in json.load(open("/home/user/screen_sc_results.json")) if 'phonon_max_cm-1' in r}
Tc={"MgB2":39,"NbB2":9,"TaB2":2,"ZrB2":0,"TiB2":0,"VB2":0,"CrB2":0,"AlB2":0,"BeB2":0,"ScB2":0,"YB2":0,"CaB2":0}
names=list(d)
x=np.array([ph[n] for n in names]); y=np.array([d[n]['Bp_fraction'] for n in names]); tc=np.array([Tc[n] for n in names])

fig,ax=plt.subplots(1,2,figsize=(14,6))
fig.suptitle("#2 多描述符筛选: 电子描述符(E_F处B-p σ占比)把范围缩小到 σ 家族, 但仍需 λ 定最终",fontsize=12.5,fontweight="bold")
a=ax[0]
for n,xi,yi in zip(names,x,y):
    sc=tc[names.index(n)]
    col=OI["green"] if n=="MgB2" else (OI["orange"] if yi>0.15 else OI["black"])
    a.scatter(xi,yi,s=140 if n=="MgB2" else 70,color=col,edgecolor="white",zorder=3)
    a.annotate(f"{n}"+(f"\nTc={sc}K" if (sc>0 or yi>0.15) else ""),(xi,yi),
               textcoords="offset points",xytext=(6,4),fontsize=8.2)
a.axhspan(0.15,0.8,color=OI["orange"],alpha=0.08)
a.text(830,0.6,"σ 家族\n(sp硼化物: B-p在E_F)",color=OI["orange"],fontsize=9)
a.text(830,0.03,"d 主导\n(过渡金属硼化物)",color=OI["black"],fontsize=9)
a.set_xlabel("最大声子 (cm⁻¹) — 描述符1(振动)"); a.set_ylabel("E_F 处 B-p(σ) 占比 — 描述符2(电子)")
a.set_title("(左) 二维描述符: 缩小到 MgB₂/YB₂/CaB₂ σ家族\n但 YB₂(Tc=0) 也在其中 → 单靠这两个仍不够")

b=ax[1]
labels="MgB2 YB2 CaB2 NbB2 TaB2".split()
for n in names:
    b.scatter(y[names.index(n)],tc[names.index(n)],s=130 if n=="MgB2" else 55,
              color=OI["green"] if n=="MgB2" else OI["black"],edgecolor="white",zorder=3)
    if n in labels: b.annotate(n,(y[names.index(n)],tc[names.index(n)]),textcoords="offset points",xytext=(6,3),fontsize=8.5)
b.set_xlabel("E_F 处 B-p(σ) 占比"); b.set_ylabel("实验 Tc (K)")
b.set_title("(右) B-p占比 vs Tc: MgB₂与YB₂同为高B-p但Tc天差地别\n→ 最终区分靠电声耦合 λ (见 #3 DFPT)")
b.annotate("MgB₂ 与 YB₂ 都高B-p,\n但只有 MgB₂ 的 σ空穴\n强耦合到 E₂g 声子(λ大)",
           xy=(0.67,39),xytext=(0.30,28),fontsize=8.5,
           bbox=dict(boxstyle="round",fc="#E9F7EF",ec=OI["green"]),
           arrowprops=dict(arrowstyle="->",color=OI["green"]))
fig.tight_layout(rect=[0,0,1,0.93])
fig.savefig("/home/user/ClaudeCode/materials_design/screening_2descriptor.png",bbox_inches="tight")
print("saved screening_2descriptor.png")
