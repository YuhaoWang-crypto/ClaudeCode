"""筛选 demo: 元素替换 → 算描述符 → 排序。
(A) 上转换基质: 描述符(最大声子)有效, 正确指向低声子氟化物。
(B) 超导硼化物: 声子描述符 vs 实验Tc → 失败, 说明需电子结构描述符。"""
import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
OI={"skyblue":"#56B4E9","black":"#000000","vermillion":"#D55E00","blue":"#0072B2",
    "green":"#009E73","orange":"#E69F00","purple":"#CC79A7"}
plt.rcParams.update({"font.size":11,"axes.grid":True,"grid.alpha":0.25,"axes.axisbelow":True,
    "font.sans-serif":["WenQuanYi Zen Hei"],"axes.unicode_minus":False,"savefig.dpi":160})

uc=json.load(open("/home/user/screen_uc_results.json"))
sc=json.load(open("/home/user/screen_sc_results.json"))

fig=plt.figure(figsize=(16,6.8))
# ---- (A) 上转换基质筛选 ----
a=fig.add_subplot(1,3,1)
rows=sorted([(r['name'].split('|')[1].split('_')[0], r['phonon_max_cm-1'],
              'F' if 'F' in r['name'].split('|')[1] and 'O' not in r['name'].split('|')[1].split('_')[0] else 'O')
             for r in uc if 'phonon_max_cm-1' in r], key=lambda z:z[1])
names=[x[0] for x in rows]; vals=[x[1] for x in rows]
cols=[OI["blue"] if x[2]=='F' else OI["vermillion"] for x in rows]
a.barh(range(len(names)),vals,color=cols,edgecolor="white")
a.set_yticks(range(len(names))); a.set_yticklabels(names,fontsize=9); a.invert_yaxis()
a.set_xlabel("最大声子能量 (cm⁻¹) — 越低越好"); a.set_title("(A) 上转换基质筛选\n描述符有效 ✓")
a.axvline(350,color=OI["green"],ls="--",lw=1.5); a.text(355,len(names)-1.5,"NaYF₄\n~350",color=OI["green"],fontsize=8)
for i,v in enumerate(vals): a.text(v+5,i,f"{v:.0f}",va="center",fontsize=8)
a.legend(handles=[plt.Rectangle((0,0),1,1,color=OI["blue"]),plt.Rectangle((0,0),1,1,color=OI["vermillion"])],
         labels=["氟化物","氧化物"],fontsize=8,loc="lower right")

# ---- (B) 硼化物筛选: 声子描述符排序 ----
Tc={"MgB2":39,"NbB2":9,"TaB2":2,"ZrB2":0,"TiB2":0,"VB2":0,"CrB2":0,"AlB2":0,
    "BeB2":0,"ScB2":0,"YB2":0,"CaB2":0}  # 实验Tc(K), 多数不超导
srows=sorted([(r['name'].split('|')[1].split('_')[0], r['phonon_max_cm-1'])
              for r in sc if 'phonon_max_cm-1' in r], key=lambda z:-z[1])
sn=[x[0] for x in srows]; sv=[x[1] for x in srows]
b=fig.add_subplot(1,3,2)
cols2=[OI["green"] if n=="MgB2" else OI["skyblue"] for n in sn]
b.barh(range(len(sn)),sv,color=cols2,edgecolor="white")
b.set_yticks(range(len(sn))); b.set_yticklabels(sn,fontsize=9); b.invert_yaxis()
b.set_xlabel("最大声子 (cm⁻¹)"); b.set_title("(B) 硼化物 MB₂ 按声子排序\nMgB₂(唯一39K超导)排倒数第3 ✗")
for i,n in enumerate(sn):
    if n=="MgB2": b.text(sv[i]+8,i,"← 真正的超导体!",va="center",fontsize=8.5,color=OI["green"])

# ---- (C) 声子 vs 实验Tc: 描述符失败 ----
c=fig.add_subplot(1,3,3)
for n,v in zip(sn,sv):
    t=Tc.get(n,0); col=OI["green"] if n=="MgB2" else OI["black"]
    c.scatter(v,t,s=80 if n=="MgB2" else 45,color=col,zorder=3,edgecolor="white")
    if t>1 or n=="MgB2": c.annotate(n,(v,t),textcoords="offset points",xytext=(6,4),fontsize=8.5)
c.set_xlabel("最大声子 (cm⁻¹, UMA 描述符)"); c.set_ylabel("实验 Tc (K)")
c.set_title("(C) 描述符 vs 实测Tc: 无相关\n→ 需电子结构描述符(σ带空穴)")
c.text(700,25,"声子高的 BeB₂/AlB₂/TiB₂\n都不超导;\nMgB₂声子中等却Tc=39K",
       fontsize=8.5,bbox=dict(boxstyle="round",fc="#FFF6E5",ec=OI["orange"]))

fig.suptitle("元素替换筛选 demo: 把性质推向更好 —— 描述符能行(左) 也会失效(右)",fontsize=14,fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.94])
fig.savefig("/home/user/ClaudeCode/materials_design/screening.png",bbox_inches="tight")
print("saved screening.png")
print("UC 最佳(最低声子):",names[0],vals[0],"cm⁻¹")
print("SC 声子最高:",sn[0],sv[0],"; 但真正超导体 MgB2 声子=",dict(srows)["MgB2"])
