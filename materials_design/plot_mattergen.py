"""#1 生成式逆向设计: MatterGen 按'体模量=300GPa'生成 -> UMA 独立验证体模量。"""
import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
OI={"skyblue":"#56B4E9","black":"#000000","vermillion":"#D55E00","green":"#009E73","orange":"#E69F00"}
plt.rcParams.update({"font.size":11,"axes.grid":True,"grid.alpha":0.25,"axes.axisbelow":True,
    "font.sans-serif":["WenQuanYi Zen Hei"],"axes.unicode_minus":False,"savefig.dpi":160})

res=[r for r in json.load(open("/home/user/mattergen_verify_results.json")) if 'bulk_modulus_GPa' in r]
res=sorted(res,key=lambda z:z['bulk_modulus_GPa'])
names=[r['name'].split('|')[1].split('_gen')[0] for r in res]
B=[r['bulk_modulus_GPa'] for r in res]
# 硼化物 vs 锇族合金 上色
cols=[OI["orange"] if any(e in n for e in ["B4","B8","Re3B"]) else OI["skyblue"] for n in names]

fig,ax=plt.subplots(figsize=(11,6.2))
bars=ax.barh(range(len(names)),B,color=cols,edgecolor="white")
ax.set_yticks(range(len(names))); ax.set_yticklabels(names,fontsize=9.5)
for i,v in enumerate(B): ax.text(v+4,i,f"{v:.0f}",va="center",fontsize=9)
# 目标线 + 已知材料参考
ax.axvline(300,color=OI["green"],ls="--",lw=2); ax.text(303,0.2,"MatterGen 目标 300",color=OI["green"],fontsize=9)
for x,lab in [(160,"钢~160"),(395,"锇~395"),(443,"金刚石~443")]:
    ax.axvline(x,color="gray",ls=":",lw=1); ax.text(x,len(names)-0.5,lab,color="gray",fontsize=8,rotation=90,va="top")
ax.set_xlabel("UMA 独立验证的体模量 (GPa)")
ax.set_title("#1 生成式逆向设计: MatterGen 按目标体模量生成 → UMA 验证\n"
             "只给'300 GPa'一个条件, 生成的全是超硬硼化物(橙)+锇族致密合金(蓝), 全部 ≥ 目标",
             fontsize=12)
ax.legend(handles=[plt.Rectangle((0,0),1,1,color=OI["orange"]),plt.Rectangle((0,0),1,1,color=OI["skyblue"])],
          labels=["重金属硼/碳化物(超硬类)","Os/Ir/Ru/Re 致密合金"],fontsize=9,loc="lower right")
fig.tight_layout()
fig.savefig("/home/user/ClaudeCode/materials_design/mattergen_inverse_design.png",bbox_inches="tight")
print("saved. UMA体模量范围:",min(B),"-",max(B),"GPa, 均值",round(np.mean(B)))
