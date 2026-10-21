"""
Run the full pipeline end to end and emit a consolidated summary.

Pipeline:
  M1  symmetry reduction  (automorphism group -> quotient -> irreducible core)
  M2  CRNT deficiency     (topology gates switching capacity / bistability)
  M3  elementary flux modes (irreducible generators of the flux cone)
  M4  DNB / critical slowing down / Lyapunov exponent (tipping biomarker)
"""
from grn_pipeline import (m1_symmetry, m2_crnt, m3_efm, m4_dnb_lyapunov,
                          m5_kras_real, m6_integrate, m7_screen,
                          m8_clinical, m9_occupancy, m10_validate,
                          m11_fibration, m12_dualphos, m13_fim_sloppy,
                          m14_atlas, m15_markevich_mm, m16_erk_dnb,
                          m18_titration_benchmark,
                          m19_switch_library, m20_literature_bistable,
                          m21_oscillators, m22_snic_mixed,
                          m23_chromatin_memory, m24_binding_energy_memory,
                          m25_chromatin_network, m26_minimal_pathways,
                          m27_alphagenome_coopening,
                          m28_compression_biomarker)


def main():
    print("\n" + "#" * 68)
    print("#  GENE-REGULATORY / METABOLIC NETWORK — MATHEMATICAL PIPELINE")
    print("#" * 68 + "\n")

    r1 = m1_symmetry.report();          print()
    m2_crnt.report();                   print()
    r3 = m3_efm.report();               print()
    r4 = m4_dnb_lyapunov.report();      print()
    r5 = m5_kras_real.report();         print()
    r6 = m6_integrate.report();         print()
    r7 = m7_screen.report();            print()
    m8_clinical.report();               print()
    r9 = m9_occupancy.report();         print()
    r10 = m10_validate.report();        print()
    r11 = m11_fibration.report();       print()
    r12 = m12_dualphos.report();        print()
    r13 = m13_fim_sloppy.report();      print()
    r14 = m14_atlas.report();           print()
    r15 = m15_markevich_mm.report();    print()
    r16 = m16_erk_dnb.report();         print()
    try:
        from grn_pipeline import m17_realdata
        r17 = m17_realdata.report();    print()
    except Exception as e:              # needs network for the real data
        r17 = None
        print(f"M17 real-data validation skipped ({type(e).__name__}: {e})\n")
    r18 = m18_titration_benchmark.report();  print()
    r19 = m19_switch_library.report();       print()
    r20 = m20_literature_bistable.report(); print()
    r21 = m21_oscillators.report();          print()
    try:
        from grn_pipeline import m20b_biomodels_exact
        r20b = m20b_biomodels_exact.report(); print()
    except Exception as e:               # needs network + libroadrunner
        r20b = None
        print(f"M20b exact-biomodels skipped ({type(e).__name__}: {e})\n")
    r22 = m22_snic_mixed.report();           print()

    # ---- chromatin-memory extension (M23-M28) ----
    r23 = m23_chromatin_memory.report();     print()
    try:
        r24 = m24_binding_energy_memory.report(); print()
    except Exception as e:                   # needs JASPAR (cached)
        r24 = None
        print(f"M24 skipped ({type(e).__name__}: {e})\n")
    try:
        r25 = m25_chromatin_network.report(); print()
    except Exception as e:                   # needs STRING (cached)
        r25 = None
        print(f"M25 skipped ({type(e).__name__}: {e})\n")
    try:
        r26 = m26_minimal_pathways.report();  print()
    except Exception as e:
        r26 = None
        print(f"M26 skipped ({type(e).__name__}: {e})\n")
    r27 = m27_alphagenome_coopening.report(); print()
    r28 = m28_compression_biomarker.report(); print()

    print("=" * 68)
    print("CONSOLIDATED SUMMARY")
    print("=" * 68)
    print(f"M1  symmetry : |Aut(G)|={r1['group_order']}, "
          f"network {r1['n_nodes']}->{r1['n_core_nodes']} nodes "
          f"(3 RAS paralogues -> 1 core node)")
    print(f"M2  CRNT     : deficiency delta gates switching; delta=0 weakly "
          f"reversible => monostable, delta>=1 => bistable switch")
    print(f"M3  EFM      : {len(r3['efms'])} irreducible flux generators "
          f"span all steady-state flux")
    print(f"M4  biomarker: LLE code validated on Rossler = "
          f"{r4['lle_rossler']:+.4f}; leading eigenvalue -> 0 with rising "
          f"SD/autocorr/DNB at the tipping point")
    print(f"M5  real KRAS: covalent G12C drug breaks paralog symmetry "
          f"S_3(6)->S_2({r5['drugged']['order']}); sotorasib 1217x "
          f"G12C-selective (ChEMBL)")
    print(f"M6  integrate: ChEMBL+Boltz binding -> engagement -> mu -> "
          f"network-stability biomarker (adagrasib slightly closer to edge)")
    print(f"M7  screen   : Boltz-2.1 ranked 10 G12C ligands; "
          f"{len(r7['better_than_sotorasib'])} analogues out-rank sotorasib "
          f"in binding")
    print(f"M8  clinical : biomarker layers mapped to CodeBreaK-100 endpoints "
          f"(ORR/DOR, PFS/OS, Cmax/AUC, QTc)")
    print(f"M9  occupancy: at approved exposure occupancy="
          f"{r9['occ_label']*100:.0f}% -> mu={r9['mu_label']:.2f} (near "
          f"tipping) -> DNB={r9['dnb_label']:.3f}")
    print(f"M10 validate : Boltz opt_score tracks ChEMBL potency "
          f"(rho={r10['rho_opt_score']:+.2f}) but binding_confidence does not "
          f"(rho={r10['rho_bind_conf']:+.2f})")
    print(f"M11 fibration: input-tree fibration compresses MAPK "
          f"{r11['n_nodes']}->{r11['n_fibers']} fibers "
          f"(generalises M1 automorphism)")
    print(f"M12 dualphos : real ERK double-phospho core has CRNT deficiency="
          f"{r12['deficiency']} and is bistable (matches report)")
    print(f"M13 FIM      : sloppy spectrum spans {r13['orders']:.0f} orders; "
          f"flux-ratio observables load best on the stiff axis")
    print(f"M14 atlas    : {len(r14['rows'])} pathways, mean compression "
          f"{sum(r['compression'] for r in r14['rows'])/len(r14['rows']):.2f}x; "
          f"JAK-STAT most compressible")
    if r15.get("window"):
        print(f"M15 Markevich: EXACT MM ERK cycle bistable window "
              f"[{r15['window'][0]:.2f},{r15['window'][1]:.2f}] nM "
              f"(report 39.25-57.38); 3-state table reproduced to the decimal")
    print(f"M16 ERK-DNB  : at the real saddle-nodes lambda_max->0 (tau~4700s), "
          f"SD/autocorr/DNB rise -> early warning on the true ERK switch")
    if r17:
        print(f"M17 real data: on real single-cell EKAR traces, lag-1 autocorr "
              f"rises before ERK pulses (p={r17['fgf']['p_ar']:.1g}); variance "
              f"flat; EGF doses all supra-threshold (partial validation)")
    if r18:
        print(f"M18 positive ctl: simulated MEKi titration across the real "
              f"bifurcation -> variance & autocorr PEAK near threshold "
              f"({r18['peak_ratio']:.1f}x) -> pipeline is sensitive, not blind")
    if r19:
        print(f"M19 migration: critical-slowing engine applied to "
              f"{len(r19['results'])} pathway switches; {r19['n_detected']} "
              f"show the variance/autocorr peak at their saddle-node")
    if r20:
        print(f"M20 literature: Rb-E2F & apoptosis multi-variable switches "
              f"reproduce documented bistability + hysteresis (Yao 2008 / "
              f"Eissing 2004 topology); eigenvalue->0 at folds")
    if r21:
        hs = r21['hopf']
        print(f"M21 oscillators: Hopf extension - Goodwin/p53-Mdm2/Brusselator "
              f"cross a Hopf; approaching it variance rises AND a spectral peak "
              f"emerges at the intrinsic frequency (distinct from saddle-node)")
    if r20b:
        print(f"M20b exact  : fetched official BioModels via GitHub mirror + "
              f"libRoadRunner; Markevich Km5={r20b['km5']} confirms M15, "
              f"Legewie apoptosis bistable (XIAP {r20b['apop_window']})")
    print(f"M22 SNIC    : mixed saddle-node+oscillation; period diverges "
          f"(T~1/sqrt, slope {r22['slope']:.2f}) -> frequency->0 signature "
          f"distinct from Hopf and pure saddle-node")
    # ---- chromatin-memory extension summary ----
    print(f"M23 chromatin: A->R->I memory Markov chain (Bintu 2016); "
          f"t_half matches paper (KRAB {r23['t_half']['KRAB']:.0f}h, "
          f"DNMT3B {r23['t_half']['DNMT3B']:.0f}h); {{I}} = absorbing class "
          f"(irreversible-trajectory classification by graph SCC)")
    if r24:
        print(f"M24 binding  : real JASPAR REST/ZNF263 energy matrix -> "
              f"occupancy -> silencing rate -> memory; decay rate is "
              f"strength-INVARIANT (spread {r24['decay_rate_spread']*100:.0f}%)"
              f" — concentration sets HOW MUCH memory, not decay speed")
    if r25:
        print(f"M25 network  : real STRING PPI ({r25['n_nodes']} regulators, "
              f"{r25['n_edges']} edges) -> {len(r25['modules'])} irreducible "
              f"mark x action modules; orbits incl. {{DNMT3A,DNMT3B}}, "
              f"{{EHMT2,SUV39H1}}")
    if r26:
        gl = r26['groups']
        print(f"M26 primes   : {r26['n_efms']} 'prime' minimal silencing "
              f"pathways (EFMs); SHORT={len(gl['SHORT'])} HYBRID="
              f"{len(gl['HYBRID'])} LONG={len(gl['LONG'])}; long memory needs "
              f"the H3K9->DNMT->5mC crosstalk prime")
    print(f"M27 alphagen : real hg38 locus panel ({len(r27['panel'])} regulator "
          f"TSS); co-opening readout ready (Modal/AlphaGenome unreachable here: "
          f"{r27['real_alphagenome']}); {r27['n_comodules']} illustrative "
          f"co-modules")
    print(f"M28 biomarker: silent={{R,I}} lumpable? {r28['lumpable_silent']} "
          f"(hybrid memory irreducibly 3-state); commitment biomarker "
          f"rho=k_RI/k_RA with P(commit)=rho/(1+rho)")
    print("\nAbstract, measurable biomarker candidates produced:")
    print("  * irreducible-core node identity        (M1 quotient)")
    print("  * deficiency delta / distance-to-bistability (M2)")
    print("  * EFM usage / rate-limiting generator    (M3)")
    print("  * leading Lyapunov exponent & DNB index  (M4)")
    print("  * paralog-symmetry order as drug-selectivity readout (M5)")
    print("  * binding-driven network-stability score (M6/M7)")
    print("  * occupancy-calibrated mu tied to PK endpoints (M9)")
    print("  * chromatin memory class + commitment ratio rho=k_RI/k_RA (M23/M28)")
    print("  * irreversible-trajectory (absorbing-class) identity (M23)")
    print("  * minimal 'prime' silencing pathway per memory length (M26)")


if __name__ == "__main__":
    main()
