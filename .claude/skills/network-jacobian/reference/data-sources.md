# Public data sources for network-Jacobian analysis

No database stores Jacobians. A Jacobian is a property of a state, so you
always **build** it from three kinds of data:

1. a network with kinetics (what the equations are),
2. state-specific values (concentrations, enzyme levels, fluxes),
3. optionally a perturbation/time series (to read eigenvalues directly).

Reachability below was tested from this sandbox (agent proxy) on 2026-10-11.

## 1. Kinetic models: direct Jacobians

| Source | Access that works | Notes |
|---|---|---|
| **BioModels** (curated SBML+ODE) | `raw.githubusercontent.com/biomodels/{id}/master/{id}/{id}.xml` ✅ (`core.fetch_biomodel`) | biomodels.org / EBI download are proxy-blocked. Used: BIOMD64 Teusink glycolysis, BIOMD206 Wolf oscillations; also BIOMD27 Markevich ERK, BIOMD102 Legewie apoptosis, BIOMD42 Nielsen, BIOMD172 Pritchard |
| **SABIO-RK** (kinetic constants) | REST `sabiork.h-its.org/sabioRestWebServices` (301 → https; untested beyond) | Km/kcat with conditions; an alternative to BRENDA for parameterising a GEM |
| **BRENDA** | website 200; the SOAP/REST API needs free registration | Use for Km, Ki, kcat and **allosteric effector** lists, i.e. the candidate EDGES |

## 2. Genome-scale networks: topology / sparsity of ε

| Source | Access | Use |
|---|---|---|
| **BiGG** (Recon3D, iML1515…) | `https://bigg.ucsd.edu/api/v2/models` ✅ | stoichiometry N + gene–reaction rules; gives the sparsity mask for the inverse Jacobian (COVRECON does exactly this) |
| **Human-GEM** (SysBioChalmers) | GitHub raw ✅ | current human consensus GEM (YAML/SBML) |
| **KEGG** REST | `rest.kegg.jp` ✅ | pathway membership, reaction–compound links (license terms apply) |

## 3. Patient / cohort omics: state-specific values

| Source | Access | What you get |
|---|---|---|
| **Metabolomics Workbench** | REST `metabolomicsworkbench.org/rest/...` ✅ | thousands of studies with group labels; covariance per group → Lyapunov inverse (D1-D) |
| **MetaboLights** | `ebi.ac.uk/metabolights/ws` ✅ | same, EBI |
| **TCGA / GDC** | `api.gdc.cancer.gov` ✅ | tumour/normal transcriptome → enzyme-level scaling of Vmax per patient (patient-specific Jacobian) |
| **CPTAC / PDC** proteomics | `pdc.cancer.gov/graphql` (400 = endpoint alive, needs a query) | enzyme protein abundance (better than mRNA for Vmax) |
| **GTEx** | `gtexportal.org/api/v2` ✅ | normal-tissue enzyme expression reference |
| HMDB | 403 through proxy ✗ | use Metabolomics Workbench `refmet` instead |

## 4. Perturbation-relaxation time series (eigenvalues directly)

| Source | Access | Notes |
|---|---|---|
| **CGMacros** (PhysioNet, CC BY-NC-SA 4.0) | ✅ via HTTP range reads of the 627 MB zip (`netjac.data.cgmacros`, ~45 MB of CSV) | 45 adults, healthy/prediabetes/T2D, 1-min CGM + logged meals = ~10–30 natural perturbations per person. Das et al., Sci Data 12:1557 (2025) |
| PhysioNet (other) | ✅ | ICU vital-sign recoveries, etc. Check each dataset's access tier (open vs credentialed) |
| Hall 2018 "glucotypes" CGM | PLoS Biol supplement (not tested) | normoglycaemic people with CGM |
| Pertz-lab single-cell ERK (EKAR) | GitHub raw ✅ (see `network-biomarker` skill, M17) | single-cell signalling relaxations |

## What to use for which direction

- **Edge biomarkers:** BioModels (exact kinetics) or GEM + SABIO/BRENDA to
  build ε. Then group covariance from Metabolomics Workbench (passive, needs
  very large N, see D1-D), or better, design a small challenge and measure
  two time points.
- **Perturbation-relaxation:** CGMacros / OGTT / any challenge time course,
  analysed with `relax.fit_modes`. In humans, compare against the clinical
  standard (fasting glucose, HbA1c), including the partial correlation.
- **Fragility / probe placement:** needs a calibrated kinetic model
  (BioModels) or a patient-specific GEM-kinetic hybrid (TCGA/CPTAC enzyme
  levels × SABIO constants). The probe panel = species with the largest
  critical-mode participation.

## Key citations (verified via PubMed)

- Steuer, Kurths, Fiehn, Weckwerth 2003 *Bioinformatics* 19:1019 (PMID 12761066): metabolite correlations ← network via stochastic fluctuations
- Nägele … Weckwerth 2014 *PLoS One* 9:e96188 (PMID 24695071): differential Jacobian from metabolomics covariance
- Sun, Länger, Weckwerth 2015 *Front Bioeng Biotechnol* 3:188 (PMID 26636075): ill-conditioning of the inverse problem
- Li, Waldherr, Weckwerth 2023 *Bioinformatics* 39:btad397 (PMID 37402625): COVRECON
- Liu … Chen 2016 *NAR* 44:e164 (PMID 27596597): sample-specific networks; Zhang … Chen 2015 *J Mol Cell Biol* 7:231 (PMID 25917597): edge biomarkers
- Chen … Aihara 2012 *Sci Rep* 2:342 (PMID 22461973): dynamical network biomarkers
- Olde Rikkert … Scheffer 2016 *Crit Care Med* 44:601 (PMID 26765499); Scheffer … 2018 *PNAS* 115:11883 (PMID 30373844): slowing of recovery as resilience marker
- Teusink … Snoep 2000 *Eur J Biochem* 267:5313 (PMID 10951190); Bergman … Cobelli 1979 *Am J Physiol* 236:E667 (PMID 443421)
- Kacser & Burns 1973 (PMID 4148886); Heinrich & Rapoport 1974 *Eur J Biochem* 42:89 (PMID 4830198): MCA theorems
- Sugimoto … Kuroda 2025 *Commun Med* 5:103 (PMID 40263561): CGM autocorrelation index tracks clamp disposition index (closest prior art for D2-C)
