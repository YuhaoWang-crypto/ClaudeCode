import pandas as pd, re
from collections import Counter
df = pd.read_pickle('/home/user/ClaudeCode/df.pkl')
EARLY=['Preclinical testing','Phase 1 (or 0)','Lead discovery/optimization','Target identification/validation']

# ---- Convergence intersections (the real alpha: modality x new domain) ----
def has(pat): return df['text'].str.contains(pat, regex=True, na=False)
tce=has(r'\btce\b|t.?cell engager|bispecific|cd3[\s\-x]')
deg=has(r'protac|degrader|molecular glue')
rad=has(r'radioligand|radiopharm|theranostic|177lu|225ac|actinium|lutetium|alpha.?emitter')
adc=has(r'\badc\b|antibody.?drug conjugate|drug conjugate')
ii =has(r'autoimmune|inflammatory|lupus|rheumatoid|psorias|atopic|colitis|\bibd\b|fibrosis')
cns=has(r'\bcns\b|brain|neuro|alzheimer|parkinson|glioblastoma')
onco=has(r'tumou?r|cancer|neoplasm|carcinoma|oncolog')
oralbio=has(r'oral (delivery|bioavail|peptide|formulation)|orally (available|bioavail)')
glue=has(r'molecular glue')
alpha=has(r'alpha.?emitter|225ac|actinium|astatine|211at|thorium|227th')

inter = {
 'TCE x non-oncology (I&I/autoimmune)': tce & ii & ~onco,
 'Degrader x I&I (beyond onco)': deg & ii,
 'Molecular glue (any)': glue,
 'Radiopharma x alpha-emitter': rad & alpha,
 'ADC x non-oncology': adc & ii & ~onco,
 'Oral biologic / oral peptide': oralbio,
 'CNS-penetrant biologic/BBB shuttle': cns & has(r'blood.?brain|bbb|brain.?penetr|shuttle|transferrin receptor'),
 'TCE x solid tumor': tce & onco & has(r'solid tumou?r'),
 'In-vivo CAR / in-vivo cell engineering': has(r'in.?vivo car|in vivo cell|in.?vivo engineer|in.?vivo reprogram'),
 'Bispecific/masked/conditional activation': has(r'conditional|masked|prodrug antibody|logic.?gated|tumou?r.?activat'),
}
print('=== CONVERGENCE INTERSECTIONS ===')
print('%-42s%6s%7s%7s'%('Intersection','N','early%','acad%'))
for n,mask in inter.items():
    s=df[mask]
    if len(s)==0:
        print('%-42s%6d'%(n,0)); continue
    print('%-42s%6d%7.0f%7.0f'%(n,len(s),s['Development phase'].isin(EARLY).mean()*100,s['is_academic'].mean()*100))

# ---- Emerging specific targets in early academic subset ----
sub = df[(df['is_academic']) & (df['Development phase'].isin(EARLY))]
print('\n=== Early-stage ACADEMIC subset: %d assets ==='%len(sub))
txt=' '.join(sub['text'].tolist())
# gene-like tokens
genes=Counter(re.findall(r'\b[A-Z]{2,5}[0-9]{0,2}\b', ' '.join((sub['Asset name'].fillna('')+' '+sub['Description'].fillna('')).tolist())))
stop=set('DNA RNA FDA USA IND ADC CAR NIH THE AND FOR NEW CD PD IL TGF WHO EU US OF IN TO A AN IS AS BY II III IV IP CEO'.split())
print('\nTop recurring uppercase target/gene tokens (early-academic):')
for k,v in genes.most_common(60):
    if k in stop or len(k)<2: continue
    print('  %4d  %s'%(v,k))
