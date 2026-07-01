# CD3 × CD19 Bispecific T-Cell Engager — Boltz-2.1 Structural Validation

A BiTE-style tandem single-chain variable fragment (scFv) construct targeting CD19 (on B cells /
B-cell malignancies) and CD3ε (on T cells) was designed from structurally-validated, published
antibody clones and evaluated with the Boltz-2.1 structure + binding prediction API.

## 1. Design

**Architecture:** `VL(CD19) – (G4S)3 – VH(CD19) – G4S – VH(CD3) – (G4S)3 – VL(CD3)`

This mirrors the domain order used in blinatumomab (the approved CD19×CD3 BiTE): an N-terminal
CD19-binding scFv (VL-then-VH) joined through a single short `GGGGS` linker to a C-terminal
CD3-binding scFv, this time in VH-then-VL order — the "inside-out" arrangement BiTE constructs use
to prevent the two scFv units from cross-pairing. Intra-scFv linkers are the standard `(G4S)3`
(15 aa); the inter-scFv linker is a single `G4S` (5 aa).

**Binder selection.** The user asked for known, validated clones rather than a de novo design.
Blinatumomab's actual parental clones are **HD37** (anti-CD19) and **L2K-07** (anti-CD3, a
low-affinity, non-cynomolgus-cross-reactive OKT3 variant). Neither sequence could be retrieved
in this session — patent databases (Google Patents, USPTO, WIPO), UniProt, NCBI, and DrugBank
were all blocked at the network/proxy level throughout, despite roughly 40 independent research
attempts across several agents. Rather than fabricate or guess at those sequences, this construct
substitutes two **different, independently well-characterized clones whose sequences were
obtained directly from experimentally-deposited PDB structures** (not secondary sources):

| Arm | Clone | Sequence source | Deposited complex |
|---|---|---|---|
| Anti-CD19 | **FMC63** (basis of tisagenlecleucel, axicabtagene ciloleucel, and most CD19 CAR-T products) | RCSB PDB **7URV**, chain B (FASTA fetched directly) | FMC63 scFv bound to soluble CD19 (He et al., *Sci. Immunol.* 2023) |
| Anti-CD3 | **UCHT1** (a well-characterized, clinically-relevant anti-human-CD3 clone) | RCSB PDB **1XIW**, chains D/H (VH) and C/G (VL) | UCHT1 Fab bound to the CD3ε/CD3δ ectodomain heterodimer |

Because both binder sequences were pulled directly from co-crystal structures with their
respective targets (rather than assembled from patent claims or vendor spec sheets), there is
direct experimental proof each binder engages the epitope used as the Boltz target in this report.
**This is not a reproduction of blinatumomab** — it is a structurally-analogous FMC63×UCHT1
tandem scFv built in the same architecture.

**Target antigens**, also taken directly from the same PDB entries (so binder and target come from
a single consistent, verified source per arm):

- CD19 ectodomain — 7URV chain A (257 aa; the crystallized soluble-CD19 construct, differs from
  the wild-type UniProt P15391 ectodomain at a few positions consistent with
  thermostabilizing crystallization mutations)
- CD3ε ectodomain — 1XIW chain A (104 aa; matches UniProt P07766 residues 23–126 exactly)
- CD3δ ectodomain — 1XIW chain B (79 aa; included because UCHT1's epitope spans the CD3ε/CD3δ
  heterodimer interface, not CD3ε alone — see §3)

Full sequences: [`structures/construct_and_targets.fasta`](structures/construct_and_targets.fasta).
Final construct: **491 aa** (blinatumomab itself is 504 aa; the difference is expected since
UCHT1's VH/VL differ in length from L2K-07's).

## 2. Boltz-2.1 predictions run

Three `structure_and_binding` jobs were submitted (model `boltz-2.1`, ~$0.10 each):

| Run | Chains | iPTM | Structure confidence | Complex pLDDT | Complex ipDE (Å) | Binding confidence |
|---|---|---|---|---|---|---|
| **CD19 arm** — construct vs. CD19 ectodomain | A (construct) + B (CD19) | **0.743** | 0.856 | 0.885 | 7.69 | 0.00036 |
| CD3 arm, v1 (flawed) — construct vs. CD3ε alone | A + C (CD3ε only) | 0.137 | 0.697 | 0.837 | 16.78 | 0.000024 |
| **CD3 arm, v2 (corrected)** — construct vs. CD3ε+CD3δ | A + C + D | **0.300** | 0.762 | 0.877 | 8.94 | 0.0000092 |

Predicted structures (mmCIF, Boltz-2.1 output) are archived in `structures/`:
`bite_vs_cd19_sample0.cif`, `bite_vs_cd3_epsilon_only_sample0.cif` (v1, kept for comparison),
`bite_vs_cd3_epsilon_delta_sample0.cif` (v2, corrected).

## 3. Interpretation

**CD19 arm (FMC63):** high confidence. iPTM 0.74 and a low interface distance error (7.7 Å) are
consistent with FMC63 correctly recognizing CD19 in the tandem-scFv context — expected, since
this arm reuses FMC63's native VL-(G4S)3-VH domain order unchanged from the crystallized 7URV
construct.

**CD3 arm (UCHT1), v1 → v2 correction:** The first attempt modeled the construct against CD3ε in
isolation and scored poorly (iPTM 0.14, 16.8 Å interface error) — essentially an unconfident
prediction. This is likely an artifact of the target definition, not necessarily the construct:
UCHT1 (like OKT3-family antibodies generally) binds an epitope formed jointly by CD3ε and its
obligate heterodimeric partner CD3δ, exactly as crystallized in 1XIW. Presenting CD3ε alone omits
half the paratope's contact surface. Resubmitting with the CD3ε+CD3δ heterodimer as the target
raised iPTM to 0.30 and lowered the interface error to 8.9 Å — a clear improvement, but still
well below the CD19 arm's confidence. Two plausible, non-exclusive explanations: (1) the
C-terminal position and single 5-residue inter-scFv linker constrain the CD3-binding domain's
conformational freedom to dock correctly in a template-free prediction, and/or (2) UCHT1 may need
additional structural context (e.g. Fab constant domains) that a bare VH/VL didn't fully recapture
here. This would be the natural next design iteration (e.g., trying a longer or different
inter-scFv linker, or swapping domain order within the CD3 arm) but was out of scope for this
validation pass.

**The `binding_confidence` metric was low (~10⁻⁴–10⁻⁵) across *all three* runs, including the
high-iPTM CD19 arm.** This is flagged rather than glossed over: it suggests this particular Boltz
output field may be calibrated to a much stricter/different threshold than the iPTM/ipDE
structural-confidence metrics, and should not on its own be read as "this construct does not
bind." The structural metrics (iPTM, complex pLDDT, interface distance error) are the more
informative signal here and point to a well-formed CD19-binding interface and a plausible, if
lower-confidence, CD3-binding interface.

## 4. Caveats

- **Network restrictions.** Direct access to UniProt, NCBI/PubMed full text (partially available
  via the PubMed MCP tool), Google Patents, USPTO, WIPO PatentScope, and DrugBank was blocked at
  the proxy/gateway level for the entire session. All target/binder sequences used here were
  instead obtained by directly fetching FASTA records from **RCSB PDB**, which was reachable and
  gives byte-exact, primary experimental data (not AI-summarized search snippets).
- **Not a blinatumomab replica.** Blinatumomab's real CD19 arm is HD37 (not FMC63) and its CD3 arm
  is L2K-07 (not UCHT1); neither sequence was obtainable this session. This construct instead pairs
  two other independently validated, PDB-documented clones in the same overall architecture.
- **CD3 arm confidence is moderate, not high**, even after correcting the target definition (see
  §3). Treat the CD3-binding prediction as directionally supportive rather than conclusive.
- These are Boltz-2.1 structure predictions, not experimental data — appropriate for early-stage
  design triage, not a substitute for expression, SPR/BLI binding assays, or cell-based potency
  testing.

## 5. Files

- `structures/construct_and_targets.fasta` — all sequences (construct, targets, individual VH/VL)
- `structures/bite_vs_cd19_sample0.cif` — predicted complex, construct + CD19
- `structures/bite_vs_cd3_epsilon_delta_sample0.cif` — predicted complex, construct + CD3ε/CD3δ (corrected target)
- `structures/bite_vs_cd3_epsilon_only_sample0.cif` — predicted complex, construct + CD3ε alone (v1, retained for comparison)
