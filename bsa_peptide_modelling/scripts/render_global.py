#!/usr/bin/env python3
"""Global view: three peptides' binding sites mapped onto one common BSA frame.
Aligns every complex on BSA (chain A), draws a single BSA cartoon, and overlays
each peptide (chain B) in a distinct color. Two outputs: best-model-only and all-models.
"""
import sys, pymol
from pymol import cmd
pymol.finish_launching(['pymol', '-qc'])

BASE='results'
PEPS=[
 ('orig','original_boltz_CFAGTPSILMLAGGGS', 5, 'forest',  'Original CFAGTPSILMLAGGGS'),
 ('pepA','peptideA_CFAGTPSILKKNGGGS',        4, 'orange',  'Peptide A CFAGTPSILKKNGGGS'),
 ('pepB','peptideB_KKAGTPSILMLAGGGS',        4, 'purple',  'Peptide B KKAGTPSILMLAGGGS'),
]

def setup_scene():
    cmd.reinitialize()
    cmd.bg_color('white')
    cmd.set('ray_opaque_background', 1)
    cmd.set('cartoon_fancy_helices', 1)
    cmd.set('ambient', 0.4); cmd.set('specular', 0.2)
    cmd.set('ray_shadows', 0); cmd.set('antialias', 2)

def load_ref():
    # reference BSA = original best model complex
    ref_pdb=f"{BASE}/{PEPS[0][1]}/structures/model_{PEPS[0][2]}.pdb"
    cmd.load(ref_pdb,'ref')
    cmd.hide('everything','ref')
    cmd.show('cartoon','ref and chain A')
    cmd.color('gray80','ref and chain A')
    cmd.set('cartoon_transparency', 0.35, 'ref and chain A')
    # subtle domain tint to orient (BSA domains I/II/III)
    cmd.color('lightblue', 'ref and chain A and resi 1-197')
    cmd.color('palegreen', 'ref and chain A and resi 198-388')
    cmd.color('lightpink', 'ref and chain A and resi 389-583')

def add_peptide(tag, folder, model, color, obj):
    pdb=f"{BASE}/{folder}/structures/model_{model}.pdb"
    cmd.load(pdb, obj)
    # align this complex's BSA onto ref BSA, moving the whole object (incl. peptide)
    cmd.align(f'{obj} and chain A', 'ref and chain A')
    cmd.hide('everything', obj)
    pep=f'{obj}_pep'
    cmd.create(pep, f'{obj} and chain B')
    cmd.delete(obj)
    cmd.show('cartoon', pep)
    cmd.set('cartoon_transparency', 0.0, pep)
    cmd.color(color, pep)
    return pep

def finalize(out, view_sel='ref and chain A'):
    cmd.set('cartoon_fancy_helices', 1)
    cmd.orient(view_sel)
    cmd.zoom(view_sel, 3)
    cmd.set('ray_trace_mode', 0)
    cmd.ray(1800, 1350)
    cmd.png(out, dpi=200)
    print('wrote', out)

mode = sys.argv[1] if len(sys.argv)>1 else 'best'
out  = sys.argv[2]

setup_scene()
load_ref()
if mode=='best':
    for tag,folder,model,color,label in PEPS:
        pep=add_peptide(tag,folder,model,color,tag)
        cmd.show('sticks', pep)
        cmd.set('stick_radius',0.2,pep)
        cmd.set('cartoon_transparency',0.0,pep)
else:  # all 5 models per peptide, best opaque, rest transparent
    for tag,folder,model,color,label in PEPS:
        for m in range(1,6):
            obj=f'{tag}_{m}'
            pep=add_peptide(tag,folder,m,color,obj)
            if m==model:
                cmd.set('cartoon_transparency',0.0,pep)
            else:
                cmd.set('cartoon_transparency',0.55,pep)
finalize(out)
cmd.quit()
