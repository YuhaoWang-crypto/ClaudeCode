# Anti-Tirzepatide antibody VH/VL (scFv) — Boltz-2 de novo design

Designed with Boltz-2 protein_design (curated **human antibody** frameworks), CDRs designed against the consensus epitope (F22/V23/L26/I27 + E3) of the Aib-containing Tirzepatide backbone. Chains: **ABH1 = VH, ABL1 = VL**.

> Design ran against the Aib→Ala backbone (non-standard Aib broke the design folder); top scFvs are re-validated against the FULLY-modified drug (Aib + K20 lipid) — see scFv validation section.

## Top 5 designs (ranked by binding_confidence)

| ID | ipTM | bind_conf | fold | light | CDR-H3 | CDR-L3 |
|---|---|---|---|---|---|---|
| ab1 | 0.96 | 0.143 | 0.81 | lambda | `ARLEVSRDGNTINVYLDT` | `MTFSPHHGLV` |
| ab2 | 0.92 | 0.141 | 0.64 | lambda | `ARTAPASLPAAHHVMWL` | `ATNDPSNAAGVV` |
| ab3 | 0.82 | 0.126 | 0.46 | lambda | `ARSAHRLTTSGIAALTL` | `MTASPSSNGLPV` |
| ab4 | 0.82 | 0.121 | 0.48 | lambda | `ARGASPANANGLGVLYL` | `QVNPLSNPNLLV` |
| ab5 | 0.93 | 0.099 | 0.77 | kappa | `GLWEGDYFSH` | `QTSYEQPT` |

## Sequences (top 2)

### ab1  (light: lambda)
- **VH** (125 aa): `EVQLVESGGGLVKPGGSLRLSCAASGFTFSEYSFAWVRQAPGKGLEWVSTISADGSLTGYAPRVAGRFTISRDNAKNSLYLQMNSLRAEDTAVYFCARLEVSRDGNTINVYLDTWGQGTMVTVSS`
- **VL** (108 aa): `ESVLTQPPSVSGAPGQRVTISCTGDADGIGSFGVSWYQQLPGTAPKLLISNNKRPAGVPDRFSGSKSGTSASLAITGLQAEDEADYYCMTFSPHHGLVFGGGTKLTVL`
- **scFv** (VH–(G4S)₃–VL): `EVQLVESGGGLVKPGGSLRLSCAASGFTFSEYSFAWVRQAPGKGLEWVSTISADGSLTGYAPRVAGRFTISRDNAKNSLYLQMNSLRAEDTAVYFCARLEVSRDGNTINVYLDTWGQGTMVTVSSGGGGSGGGGSGGGGSESVLTQPPSVSGAPGQRVTISCTGDADGIGSFGVSWYQQLPGTAPKLLISNNKRPAGVPDRFSGSKSGTSASLAITGLQAEDEADYYCMTFSPHHGLVFGGGTKLTVL`

### ab2  (light: lambda)
- **VH** (123 aa): `VQLVQSGAEVKKPGASVKVSCKASGYTFTTDSGIAWVRQAPGQGLEWMGRIDGSGNTIINPKYRGRVTMTTDTSTSTAYMELRSLRSDDTAVYYCARTAPASLPAAHHVMWLWGRGTLVTVSS`
- **VL** (106 aa): `YVLTQPPSVSVAPGKTARITCTGVAEGQPVSWYQQKPGQAPVLVIYMGERPPGIPERFSGSNSGNTATLTISRVEAGDEADYYCATNDPSNAAGVVFGGGTKLTVL`
- **scFv** (VH–(G4S)₃–VL): `VQLVQSGAEVKKPGASVKVSCKASGYTFTTDSGIAWVRQAPGQGLEWMGRIDGSGNTIINPKYRGRVTMTTDTSTSTAYMELRSLRSDDTAVYYCARTAPASLPAAHHVMWLWGRGTLVTVSSGGGGSGGGGSGGGGSYVLTQPPSVSVAPGKTARITCTGVAEGQPVSWYQQKPGQAPVLVIYMGERPPGIPERFSGSNSGNTATLTISRVEAGDEADYYCATNDPSNAAGVVFGGGTKLTVL`

## Formats & humanized constant regions
- **scFv-Fc (peptibody / minibody):** scFv + human **IgG1 Fc** (hinge-CH2-CH3). Sequences per design in `antibody_constructs.fasta` (field `scFv_Fc_IgG1`).
- **Full IgG:** Heavy = VH + human **IgG1 CH1-hinge-CH2-CH3**; Light = VL + human **Cκ** (ab5) or **Cλ** (ab1–ab4). Fields `IgG_heavy_chain` / `IgG_light_chain`.
- **Effector options:** IgG1 (WT, ADCC/CDC-competent), **IgG1-LALA** (L234A/L235A, effector-silenced — for pure neutralization/PK), **IgG4-S228P** (stabilized, low effector). All in `human_constant_regions.fasta`.
- All designed regions were filtered against N-glycosylation sequons (NxS/NxT).

## Honest read
High ipTM (0.82–0.97) = confident docking geometry; binding_confidence 0.10–0.14 (top) is modest — expected for a small, flexible, modified-peptide antigen. These are **triage-ranked hypotheses for wet-lab test**, not affinity guarantees. Recommend expressing top 3–4 as scFv-Fc, test by SPR/BLI vs synthetic (modified) tirzepatide; affinity-mature CDR-H3 on the best.
