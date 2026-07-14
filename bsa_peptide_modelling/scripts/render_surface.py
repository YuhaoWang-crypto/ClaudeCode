#!/usr/bin/env python3
"""BSA molecular surface with the peptide nestled in its pocket.
BSA (chain A): semi-transparent gray surface; the pocket (surface within 6 A of the
peptide) highlighted marine. Peptide (chain B): orange sticks + cartoon.
Usage: render_surface.py <pdb> <out_png> [peptide_color]
"""
import sys, pymol
from pymol import cmd
pymol.finish_launching(['pymol', '-qc'])

pdb, out = sys.argv[1], sys.argv[2]
pepcol = sys.argv[3] if len(sys.argv) > 3 else 'orange'

cmd.load(pdb, 'cplx')
cmd.remove('solvent')
cmd.hide('everything')

# BSA surface
cmd.show('surface', 'chain A')
cmd.set('surface_quality', 1)
cmd.set('transparency', 0.25)                 # surface transparency
cmd.color('gray85', 'chain A')
# pocket = BSA surface contacting the peptide -> marine (atom-level for a tight footprint)
cmd.select('pocket', 'chain A within 4.5 of chain B')
cmd.color('marine', 'pocket')

# peptide: cartoon + sticks (sits inside/against the surface)
cmd.show('cartoon', 'chain B')
cmd.show('sticks', 'chain B')
cmd.set('cartoon_transparency', 0.0, 'chain B')
cmd.color(pepcol, 'chain B')
cmd.util.cnc('chain B and not elem C')
cmd.set('stick_radius', 0.2, 'chain B')

# presentation
cmd.bg_color('white')
cmd.set('ray_opaque_background', 1)
cmd.set('ambient', 0.4); cmd.set('specular', 0.2)
cmd.set('ray_shadows', 0); cmd.set('antialias', 2)

cmd.set('ray_interior_color', 'gray60')       # softer interior where surface is clipped
cmd.orient('chain B')
cmd.zoom('chain B', 13)                        # pull back to keep the front surface intact
cmd.clip('slab', 200)                          # widen depth slab to avoid clipping the pocket
cmd.turn('y', 10)
cmd.set('ray_trace_mode', 0)
cmd.ray(1600, 1200)
cmd.png(out, dpi=200)
print('wrote', out)
cmd.quit()
