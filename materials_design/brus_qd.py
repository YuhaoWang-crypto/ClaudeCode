"""量子点 Brus 方程带隙/发射波长计算。
E_QD(R) = Eg_bulk + (h^2/8R^2)(1/me*+1/mh*) - 1.8 e^2/(4 pi eps0 epsr R)
"""
import numpy as np
h=6.626e-34; me=9.109e-31; e=1.602e-19; eps0=8.854e-12
def brus(Eg, mes, mhs, epsr, R_nm):
    R=R_nm*1e-9
    conf=(h**2/(8*R**2))*(1/(mes*me)+1/(mhs*me))/e
    coul=1.8*e/(4*np.pi*eps0*epsr*R)
    return Eg+conf-coul
mats={"CdSe":(1.74,0.13,0.45,10.2),"PbS":(0.41,0.085,0.085,17.0),
      "CsPbBr3":(2.36,0.15,0.14,7.3),"Si":(1.12,0.26,0.39,11.7)}
if __name__=="__main__":
    for name,(Eg,mes,mhs,epsr) in mats.items():
        print(name, [f"d={d}nm:{brus(Eg,mes,mhs,epsr,d/2):.2f}eV" for d in (3,4,6,8)])
