"""③ DFT(GPAW/PBE) 态密度: 纯 SrTiO₃ (绝缘,带隙) vs Nb 掺杂 (n型, E_F 进导带)。"""
import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
OI={"skyblue":"#56B4E9","black":"#000000","vermillion":"#D55E00","blue":"#0072B2",
    "green":"#009E73","orange":"#E69F00","purple":"#CC79A7"}
plt.rcParams.update({"font.size":11,"axes.grid":True,"grid.alpha":0.25,"axes.axisbelow":True,
    "font.sans-serif":["WenQuanYi Zen Hei"],"axes.unicode_minus":False,"savefig.dpi":160})

res={r["name"]:r for r in json.load(open("/home/user/dft_results.json"))}
fig,ax=plt.subplots(1,2,figsize=(14,5.4),sharey=False)
fig.suptitle("③ 真实 DFT (GPAW/PBE) 态密度: 掺杂前后 n 型能级 (能量相对 E_F, 虚线=E_F)",
             fontsize=13,fontweight="bold")

def draw(a,r,title):
    e=np.array(r["energy_rel_Ef"]); tot=np.array(r["dos_total"])
    a.fill_between(e,tot,color=OI["skyblue"],alpha=0.45,label="总 DOS")
    a.plot(e,tot,color=OI["blue"],lw=1)
    for key,lab,c in [("pdos_Ti_d","Ti 3d",OI["vermillion"]),
                      ("pdos_Nb_d","Nb 4d",OI["green"]),
                      ("pdos_O_p","O 2p",OI["orange"])]:
        if r.get(key):
            a.plot(e,np.array(r[key]),lw=1.8,color=c,label=lab)
    a.axvline(0,color=OI["black"],ls="--",lw=1.5)
    a.text(0.1,a.get_ylim()[1]*0.9 if a.get_ylim()[1]>0 else 1,"E_F",fontsize=10)
    a.set_xlim(-8,6); a.set_xlabel("能量 E − E_F (eV)"); a.set_ylabel("态密度 (states/eV)")
    a.set_title(title); a.legend(fontsize=9,loc="upper left")

pure=res.get("pure_SrTiO3"); dope=res.get("Nbdoped_SrTiO3_2x2x1")
draw(ax[0],pure,"(左) 纯 SrTiO₃: 绝缘体, E_F 落在带隙中 (导带底=Ti 3d, 价带顶=O 2p)")
draw(ax[1],dope,"(右) Nb:SrTiO₃ (x=0.25): E_F 进入导带 → n 型 (每个Nb供1电子)")
fig.tight_layout(rect=[0,0,1,0.94])
fig.savefig("/home/user/ClaudeCode/materials_design/step3_dft_dos.png",bbox_inches="tight")
print("saved step3_dft_dos.png")
print("纯 STO Ef=",pure["E_fermi"],"; 掺杂 Ef=",dope["E_fermi"])
