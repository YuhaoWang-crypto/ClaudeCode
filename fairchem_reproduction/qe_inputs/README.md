# Fe(110) Δρ + Bader with Quantum ESPRESSO

1. Get PAW pseudos (SSSP efficiency) into qe_inputs/pseudos/.
2. Run SCF for slab, mol, complex:
   `mpirun pw.x -in espresso.pwi > scf.out` (each folder)
3. Export each charge density to cube with pp.x (plot_num=0).
4. Δρ = ρ(complex) − ρ(slab) − ρ(mol): subtract the three cubes on the
   same FFT grid (use `pp.x` `plot_num=0` with identical `nr1,nr2,nr3`,
   or ASE/pymatgen to combine).
5. Bader: `bader complex_CHGCAR` (Henkelman code) → ACF.dat → sum the
   molecule atoms vs slab atoms for the net charge transfer.
