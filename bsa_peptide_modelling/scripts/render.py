#!/usr/bin/env python3
"""Render a BSA-peptide complex figure (matplotlib 3D backbone trace).
Usage: render.py <pdb> <out_png> <title> [interface_resnums_commaseparated]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def parse(pdb):
    A=[]; B=[]; Ball=[]
    Ares={}
    for line in open(pdb):
        if not line.startswith(('ATOM','HETATM')): continue
        atom=line[12:16].strip(); ch=line[21]; resn=line[17:20].strip()
        try: resi=int(line[22:26])
        except: continue
        x=float(line[30:38]); y=float(line[38:46]); z=float(line[46:54])
        if ch=='A':
            if atom=='CA': A.append((resi,x,y,z)); Ares[resi]=(x,y,z)
        elif ch=='B':
            Ball.append((x,y,z))
            if atom=='CA': B.append((resi,x,y,z))
    return A,B,Ball,Ares

def main():
    pdb,out,title=sys.argv[1],sys.argv[2],sys.argv[3]
    iface=set()
    if len(sys.argv)>4 and sys.argv[4]:
        iface={int(x) for x in sys.argv[4].split(',')}
    A,B,Ball,Ares=parse(pdb)
    A.sort(); B.sort()
    Axyz=np.array([(x,y,z) for _,x,y,z in A])
    Bxyz=np.array([(x,y,z) for _,x,y,z in B]) if B else np.empty((0,3))
    Ballxyz=np.array(Ball)
    fig=plt.figure(figsize=(7.5,6.0))
    ax=fig.add_subplot(111,projection='3d')
    # BSA backbone trace (light grey)
    ax.plot(Axyz[:,0],Axyz[:,1],Axyz[:,2],color='#c8ccd4',lw=1.0,alpha=0.75)
    # interface BSA residues highlighted
    if iface:
        ipts=np.array([Ares[r] for r in iface if r in Ares])
        if len(ipts):
            ax.scatter(ipts[:,0],ipts[:,1],ipts[:,2],color='#e8833a',s=42,
                       depthshade=True,edgecolors='#a5561f',linewidths=0.4,label='BSA interface residues')
    # peptide backbone (blue tube) + all atoms
    if len(Ballxyz):
        ax.scatter(Ballxyz[:,0],Ballxyz[:,1],Ballxyz[:,2],color='#2e6fb0',s=10,alpha=0.5)
    if len(Bxyz):
        ax.plot(Bxyz[:,0],Bxyz[:,1],Bxyz[:,2],color='#16436e',lw=4.0,label='Peptide (chain B)')
        ax.scatter(Bxyz[0,0],Bxyz[0,1],Bxyz[0,2],color='#16436e',s=60,marker='o')
    ax.set_title(title,fontsize=11)
    ax.set_axis_off()
    # equalize aspect
    allp=np.vstack([Axyz,Ballxyz]) if len(Ballxyz) else Axyz
    c=allp.mean(0); r=(allp.max(0)-allp.min(0)).max()/2
    ax.set_xlim(c[0]-r,c[0]+r); ax.set_ylim(c[1]-r,c[1]+r); ax.set_zlim(c[2]-r,c[2]+r)
    ax.legend(loc='upper right',fontsize=8,framealpha=0.8)
    ax.view_init(elev=18,azim=35)
    plt.tight_layout()
    plt.savefig(out,dpi=150,bbox_inches='tight')
    print('wrote',out)

if __name__=='__main__': main()
