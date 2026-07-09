import pandas as pd, numpy as np, re
df = pd.read_pickle('/home/user/ClaudeCode/df.pkl')

EARLY = ['Preclinical testing','Phase 1 (or 0)','Lead discovery/optimization','Target identification/validation']
P12   = ['Phase 1 (or 0)','Phase 2']
LATE  = ['Phase 3','Pre registration','Registered','On market','Commercial (industrial biotech)']

modalities = {
 'T-cell engager (TCE/bispecific)': r'\btce\b|t.?cell engager|bispecific|bi-specific|trispecific|cd3[\s\-x]',
 'ADC / conjugate': r'\badc\b|antibody.?drug conjugate|drug conjugate|payload',
 'CAR-T / cell therapy': r'\bcar[\s\-]?t\b|car[\s\-]?nk|chimeric antigen|\btil\b|tcr[\s\-]?t\b',
 'Protein degrader (PROTAC/glue)': r'protac|degrader|molecular glue|targeted protein degrad',
 'mRNA / RNA therapeutic': r'\bmrna\b|messenger rna|self.?amplif|circular rna|circrna|\bsarna\b',
 'siRNA/ASO/oligonucleotide': r'sirna|antisense|oligonucleotide|\baso\b|\brnai\b|gapmer|\bsnorna\b',
 'Gene editing (CRISPR/base/prime)': r'crispr|base edit|prime edit|cas9|cas12|zinc finger|talen|gene edit',
 'AAV / gene therapy': r'\baav\b|adeno.?associated|gene therapy|lentivir',
 'Radiopharmaceutical/theranostic': r'radioligand|radiopharm|theranostic|actinium|lutetium|177lu|225ac|alpha.?emitter|radionuclide|radiotherapeut',
 'GLP-1/incretin (metabolic)': r'glp.?1|\bgip\b|glucagon|incretin|amylin',
 'Obesity/weight': r'obesity|weight loss|adipos',
 'KRAS/RAS': r'\bkras\b|g12c|g12d',
 'PSMA/FAP tumor stroma': r'\bpsma\b|\bfap\b|fibroblast activation',
 'Microbiome': r'microbiome|microbiota|gut bacteria|probiotic',
 'AI/ML drug discovery': r'machine learning|deep learning|generative model|in silico|computational design|ai-driven|ai-based|ai platform',
 'Nanobody/VHH/single-domain': r'nanobody|\bvhh\b|single.?domain antibody',
 'Myeloid/macrophage/TME': r'myeloid|macrophage|\btreg\b|tumor microenvironment',
 'NK-cell engager/therapy': r'nk cell|nk-cell|natural killer',
 'Senescence/longevity': r'senolytic|senescen|longevity',
 'Neurodegen (tau/synuclein/TDP)': r'synuclein|tdp.?43|tauopath|neurodegener|\btau protein',
 'GPCR/allosteric': r'gpcr|allosteric|g protein.?coupled',
 'Fibrosis': r'fibrosis|fibrotic',
 'Complement': r'complement (system|pathway|c[35])',
 'Exosome/EV delivery': r'exosome|extracellular vesicle',
 'LNP/nanoparticle delivery': r'\blnp\b|lipid nanoparticle|nanoparticle deliver',
 'Antibiotic/AMR': r'antibiotic|antimicrobial resist|\bamr\b|multidrug.?resist',
 'Cancer/neoantigen vaccine': r'cancer vaccine|therapeutic vaccine|neoantigen',
 'Cytokine/IL engineering': r'\bil-?\d|interleukin|cytokine',
 'Autoimmune/I&I': r'autoimmune|inflammatory bowel|lupus|rheumatoid|psorias|atopic derm',
}

recs=[]
N=len(df)
for name,pat in modalities.items():
    m = df['text'].str.contains(pat, regex=True, na=False)
    sub=df[m]; tot=len(sub)
    if tot<15: continue
    ph=sub['Development phase']
    early=ph.isin(EARLY).mean()
    late=ph.isin(LATE).mean()
    p12=ph.isin(P12).mean()
    acad=sub['is_academic'].mean()
    part=(sub['Asset partnered']=='Partnered').mean()
    acad_early_n=int((sub['is_academic'] & ph.isin(EARLY)).sum())
    recs.append(dict(theme=name,total=tot,early=early,late=late,p12=p12,acad=acad,part=part,acad_early_n=acad_early_n))

R=pd.DataFrame(recs)
# ---- Seller's-market-corrected NEXT-WAVE score ----
# Rationale (weights explained in report):
#  + early-phase enrichment  (leading indicator of science direction)
#  + academic origin         (upstream of commercial wave, less hype-driven)
#  + "immaturity gap" = early - late  (rising vs peaking)
#  - crowding penalty        (very high volume = current/peaking wave, seller glut)
#  - partnered rate          (already buyer-validated = today's wave, not next)
def z(s):
    return (s-s.mean())/(s.std()+1e-9)
R['crowd'] = np.log1p(R['total'])
R['gap']   = R['early']-R['late']
R['nextwave'] = (1.1*z(R['early']) + 1.0*z(R['acad']) + 1.2*z(R['gap'])
                 - 0.7*z(R['crowd']) - 0.6*z(R['part']))
R=R.sort_values('nextwave',ascending=False)
pd.set_option('display.width',200)
print('%-34s%6s %5s %5s %5s %5s %6s  %6s'%('THEME','N','early','late','acad','part%','acadE','SCORE'))
for _,r in R.iterrows():
    print('%-34s%6d %5.0f %5.0f %5.0f %5.1f %6d  %+6.2f'%(
        r.theme,r.total,r.early*100,r.late*100,r.acad*100,r.part*100,r.acad_early_n,r.nextwave))
R.to_pickle('/home/user/ClaudeCode/scored.pkl')
