# Global view — where the three peptides bind BSA

All Boltz-2.1 complexes were superimposed on the shared BSA chain (chain A) and the
peptides overlaid, so the three binding sites can be compared in one frame. BSA is
tinted by its three structural domains: I (1–197), II (198–388), III (389–583).

- `global_best_pose.png` — the single best-scoring pose of each peptide.
- `global_all_models_footprint.png` — all 5 models per peptide (best opaque, rest faded),
  showing each peptide's full site distribution.

## Site summary (Boltz-2.1 + PRODIGY)

| Peptide | Best pose site | All-model distribution |
|---|---|---|
| Original `CFAGTPSILMLAGGGS` | domain **IIIA** (ΔG −14.4) | 4/5 in domain II (IIA), best in III (IIIA) |
| Peptide A `CFAGTPSILKKNGGGS` | domain **IIIA/IIIB** (ΔG −11.1) | 3/5 in domain II, 2/5 (incl. best) in domain III |
| Peptide B `KKAGTPSILMLAGGGS` | domain **IIA** (ΔG −9.8) | 3/5 in domain **I** (IB) — a site the others never use — 2/5 in domain II |

**Takeaways**
- Original and Peptide A share the classic domain II/III hydrophobic pockets (Sudlow site I
  in IIA and site II in IIIA); Peptide A's strongest pose sits in the same IIIA patch as the
  original peptide.
- Peptide B is the outlier: its two N-terminal lysines let it uniquely engage **domain I**
  (subdomain IB, residues ~114–189) in the majority of models, in addition to domain II.
- All three are somewhat promiscuous — expected, since BSA is a multi-site carrier protein.
