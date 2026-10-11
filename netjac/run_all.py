"""Run the three network-Jacobian demos end to end (figures -> figures/netjac/).

    python3 -m netjac.run_all            # all three (D2-C needs network)
    python3 -m netjac.run_all --offline  # skip the real CGMacros download
    python3 -m netjac.run_all --tcga     # also D4 patient-specific Jacobians
                                         # (TCGA-KIRC, ~15 min first time)
"""
import sys

from netjac.demos import d1_edge_biomarker, d2_perturbation_relaxation, \
    d3_fragility


def main(offline=False, tcga=False):
    r1 = d1_edge_biomarker.report(); print()
    r2 = d2_perturbation_relaxation.report(with_real=not offline); print()
    r3 = d3_fragility.report(); print()
    r4 = None
    if tcga and not offline:
        from netjac.demos import d4_patient_jacobian
        r4 = d4_patient_jacobian.report(); print()
    print("=" * 68)
    print("NETJAC SUMMARY")
    print("=" * 68)
    C1 = r1["C"]
    print(f"D1 edge   : identical metabolome (max |dln x| "
          f"{r1['B']['max_dlnx']:.0e}); 18-analyte panel AUC "
          f"{C1['panel_auc']:.2f} vs pre-specified edge readout AUC "
          f"{C1['edge']['auc']:.2f}")
    A2, C2 = r2["A"], r2["C"]
    print(f"D2 relax  : slowest λ from a 2-analyte pulse {A2['lam_est']:.2f} "
          f"(true {A2['lam_true']:.2f}); Hopf ringing read from NADH alone")
    if C2:
        s = C2["corr"]["k"][0]
        print(f"            REAL CGM: meal recovery rate vs HbA1c Spearman "
              f"{s.statistic:+.2f} (p={s.pvalue:.0g}); adds nothing beyond "
              f"fasting glucose in n=45 (partial "
              f"{C2['corr']['k'][1].statistic:+.2f})")
    print(f"D3 fragile: total fragility predicts Hopf crossing "
          f"{r3['C']['agree']*100:.0f} % (local-only "
          f"{r3['C']['agree_local']*100:.0f} %); best probe follows mode "
          f"participation (ρ={r3['D']['rho_part']:+.2f}), not abundance "
          f"(ρ={r3['D']['rho_abund']:+.2f})")
    if r4:
        sv = r4["surv"]["uni"]
        print(f"D4 TCGA   : {len(r4['tab'])} KIRC samples instantiated; "
              f"glycolytic flux [pre-specified] adj. HR "
              f"{sv['glycolytic flux [pre]']['hr_adj']:.2f} (Holm p="
              f"{sv['glycolytic flux [pre]']['p_holm']:.2g})")


if __name__ == "__main__":
    main(offline="--offline" in sys.argv, tcga="--tcga" in sys.argv)
