#!/usr/bin/env python3
"""Ray-traced PyMOL cartoon of a BSA-peptide complex, matching CCB_BSA_v2 style.
BSA (chain A): gray70 cartoon, 15% transparency.
Peptide (chain B): orange cartoon + sticks.
BSA interface residues (within 5.0 A of chain B): marine sticks.
Usage: render_pymol.py <pdb> <out_png> [caption_off]
"""
import sys, pymol
from pymol import cmd

pymol.finish_launching(['pymol', '-qc'])  # quiet, no GUI

pdb, out = sys.argv[1], sys.argv[2]

cmd.load(pdb, 'cplx')
cmd.remove('solvent')
cmd.hide('everything')

# BSA cartoon (cartoon stays gray regardless of atom/stick colors below)
cmd.show('cartoon', 'chain A')
cmd.set('cartoon_transparency', 0.15, 'chain A')
cmd.color('gray70', 'chain A')
cmd.set('cartoon_color', 'gray70', 'chain A')

# Peptide cartoon + sticks
cmd.show('cartoon', 'chain B')
cmd.show('sticks', 'chain B')
cmd.set('cartoon_transparency', 0.0, 'chain B')
cmd.set('cartoon_color', 'orange', 'chain B')
cmd.color('orange', 'chain B')
cmd.util.cnc('chain B and not elem C')  # color non-carbon by element on peptide sticks

# BSA interface residues within 5 A of peptide -> marine sticks (cartoon left gray)
cmd.select('iface', 'byres (chain A within 5.0 of chain B)')
cmd.show('sticks', 'iface')
cmd.color('marine', 'iface and elem C')
cmd.util.cnc('iface and not elem C')

# nice presentation
cmd.bg_color('white')
cmd.set('ray_opaque_background', 1)
cmd.set('cartoon_fancy_helices', 1)
cmd.set('cartoon_highlight_color', 'grey50')
cmd.set('ambient', 0.35)
cmd.set('specular', 0.25)
cmd.set('ray_shadows', 1)
cmd.set('ray_shadow', 0)
cmd.set('antialias', 2)
cmd.set('stick_radius', 0.18)
cmd.set('cartoon_side_chain_helper', 1)

# orient on the interface/peptide for a focused, framed view
cmd.orient('chain B or iface')
cmd.zoom('chain B or iface', 4)
cmd.turn('y', 8)

cmd.set('ray_trace_mode', 0)  # smooth ray trace, no black outline
cmd.ray(1600, 1200)
cmd.png(out, dpi=200)
print('wrote', out)
cmd.quit()
