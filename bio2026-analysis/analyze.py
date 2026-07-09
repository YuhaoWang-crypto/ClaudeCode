import pandas as pd, re
F='/root/.claude/uploads/56fa63aa-90a0-506f-a6f6-b3ad6138e2de/a45c0e99-AssetsExport.xlsx'
df = pd.read_excel(F, sheet_name='Asset')

acad_kw = re.compile(r'univ|institut|college|hospital|research foundation|yissum|ramot|yeda|school of|medical center|clinic\b|\bnih\b|academ|centre national|cnrs|max planck|foundation for|technology development|tech manage|commercializ', re.I)
df['is_academic'] = df['Company name'].fillna('').str.contains(acad_kw)
print('Academic assets:', int(df['is_academic'].sum()), '/', len(df),
      '(%.1f%%)' % (df['is_academic'].mean()*100))

df['text'] = (df['Asset name'].fillna('')+' '+df['Description'].fillna('')+' '+
              df['Mechanisms of action'].fillna('')+' '+df['Technologies'].fillna('')+' '+
              df['Clinical indication(s)'].fillna('')).str.lower()

EARLY = ['Preclinical testing','Phase 1 (or 0)','Lead discovery/optimization','Target identification/validation']

modalities = {
 'T-cell engager (TCE/bispecific)': r'\btce\b|t.?cell engager|bispecific|bi-specific|trispecific|cd3[\s\-x]',
 'ADC / conjugate': r'\badc\b|antibody.?drug conjugate|drug conjugate|payload|conjugate',
 'CAR-T / cell therapy': r'\bcar[\s\-]?t\b|car[\s\-]?nk|chimeric antigen|\btil\b|tcr[\s\-]?t\b',
 'Protein degrader (PROTAC/glue)': r'protac|degrader|molecular glue|targeted protein degrad',
 'mRNA / RNA therapeutic': r'\bmrna\b|messenger rna|self.?amplif|circular rna|circrna|\bsarna\b',
 'siRNA/ASO/oligonucleotide': r'sirna|antisense|oligonucleotide|\baso\b|\brnai\b|gapmer|\bsnorna\b',
 'Gene editing (CRISPR/base/prime)': r'crispr|base edit|prime edit|cas9|cas12|zinc finger|talen|gene edit',
 'AAV / gene therapy': r'\baav\b|adeno.?associated|gene therapy|lentivir|gene replace',
 'Radiopharmaceutical/theranostic': r'radioligand|radiopharm|theranostic|actinium|lutetium|177lu|225ac|alpha.?emitter|radionuclide|radiotherapeut',
 'GLP-1/incretin/metabolic': r'glp.?1|\bgip\b|glucagon|incretin|amylin',
 'Obesity/weight': r'obesity|weight loss|adipos',
 'KRAS/RAS': r'\bkras\b|g12c|g12d|\bras\b onco',
 'PSMA/FAP tumor stroma': r'\bpsma\b|\bfap\b|fibroblast activation',
 'Peptide/macrocycle': r'macrocycl|cyclic peptide|stapled peptide',
 'Microbiome': r'microbiome|microbiota|gut bacteria|probiotic',
 'AI/ML drug discovery': r'machine learning|deep learning|generative model|in silico|computational design|\bai-\b|ai-driven|ai-based|ai platform',
 'Nanobody/VHH/single-domain': r'nanobody|\bvhh\b|single.?domain antibody|heavy.?chain antibody',
 'Myeloid/macrophage/TME': r'myeloid|macrophage|\btreg\b|tumor microenvironment',
 'NK-cell engager/therapy': r'\bnk cell|nk-cell|natural killer',
 'Senolytic/aging': r'senolytic|senescen|longevity|aging',
 'Neurodegen (tau/synuclein/TDP)': r'\btau\b|synuclein|tdp.?43|amyloid|tauopath|neurodegener',
 'GPCR/allosteric': r'gpcr|allosteric|g protein.?coupled',
 'Fibrosis': r'fibrosis|fibrotic',
 'Complement': r'complement (system|pathway|c[35])|\bc3\b|\bc5\b',
 'Exosome/EV': r'exosome|extracellular vesicle',
 'LNP/nanoparticle delivery': r'\blnp\b|lipid nanoparticle|nanoparticle deliver|nanoparticle-based',
 'Autoimmune/IBD/I&I': r'autoimmune|inflammatory bowel|lupus|rheumatoid|psorias|atopic',
 'Antibiotic/AMR': r'antibiotic|antimicrobial resist|\bamr\b|multidrug.?resist',
 'Vaccine (therapeutic/cancer)': r'cancer vaccine|therapeutic vaccine|neoantigen',
 'Cytokine/IL engineering': r'\bil-?\d|interleukin|cytokine',
 'Bispecific/multispecific antibody': r'bispecific antibody|multispecific|dual.?target',
 'Antibody-oligo / novel conjugate': r'aoc\b|antibody.?oligonucleotide|immune.?stimulat.*conjugate|isac',
}
rows=[]
for name,pat in modalities.items():
    m = df['text'].str.contains(pat, regex=True, na=False)
    tot=int(m.sum()); acad=int((m&df['is_academic']).sum())
    early=int((m & df['Development phase'].isin(EARLY)).sum())
    p12=int((m & df['Development phase'].isin(['Phase 1 (or 0)','Phase 2'])).sum())
    acad_early=int((m & df['is_academic'] & df['Development phase'].isin(EARLY)).sum())
    rows.append((name,tot,acad,early,p12,acad_early))
rows.sort(key=lambda x:-x[1])
print('\n%-42s%7s%7s%7s%7s%9s'%('Modality/Theme','Total','Acad','Early','Ph1-2','AcadErly'))
for r in rows:
    print('%-42s%7d%7d%7d%7d%9d'%r)
df.to_pickle('/home/user/ClaudeCode/df.pkl')
