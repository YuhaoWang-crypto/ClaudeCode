# QM cluster refinement of a Cd(2+) coordination sphere

Turns a docked pose into a defensible binding geometry + energy.

## Cut the cluster
- Cd2+ + all first-shell donor atoms (from demo #5 site) with their side chains.
- Cap severed Cα–Cβ bonds with H (or use link atoms in QM/MM).
- Include any coordinating waters.

## ORCA input sketch (B3LYP-D3, def2 basis, Cd ECP)
```
! B3LYP D3BJ def2-TZVP def2/J RIJCOSX TightSCF Opt
%basis  NewGTO Cd "def2-TZVP" end
        NewECP Cd "def2-ECP" end
end
* xyz 2 1                 # charge = 2 (Cd2+) + ligand charges; spin: Cd2+ is d10 -> singlet
  Cd  ...
  ... (donor atoms + capped side chains)
*
```

## What to report
- optimised Cd–L bond lengths (expect 2.3-2.6 Å) and coordination number,
- interaction energy (BSSE-corrected if you quote absolute affinity),
- comparison across the shortlisted sites -> which site binds Cd2+ best,
- effect of the designed mutations (re-cut cluster, re-optimise).

For a full protein context use QM/MM (ORCA + optional xtb MM, or CP2K QM/MM).
