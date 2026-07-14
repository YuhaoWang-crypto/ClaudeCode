#!/usr/bin/env python3
"""Replicate CCB100725-YW01 pipeline: strip BTN/NH2 caps, run Prodigy, analyze interface.
Usage: analyze.py <peptide_label> <dir_with_model_1..5.pdb> <outdir>
Produces per-model output.csv + ic.csv, plus analysis.json summarizing the report tables.
"""
import sys, os, csv, json, subprocess, collections

THREE2ONE = {
 'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
 'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
 'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
CAPS = {'BTN','NH2'}

def strip_caps(inp, out):
    with open(inp) as f, open(out,'w') as o:
        for line in f:
            if line.startswith(('ATOM','HETATM')) and line[17:20].strip() in CAPS:
                continue
            o.write(line)

def run_prodigy(pdb, workdir):
    # prodigy writes <base>.ic in cwd
    base = os.path.splitext(os.path.basename(pdb))[0]
    r = subprocess.run(['prodigy', pdb, '--selection','A','B','--contact_list'],
                       cwd=workdir, capture_output=True, text=True)
    out = r.stdout
    metrics = {}
    def val(line): return line.split(':',1)[1].strip()
    for line in out.splitlines():
        line=line.strip()
        if 'intermolecular contacts' in line: metrics['contacts']=float(val(line))
        elif 'charged-charged' in line: metrics['cc']=float(val(line))
        elif 'charged-polar' in line: metrics['cp']=float(val(line))
        elif 'charged-apolar' in line: metrics['ca']=float(val(line))
        elif 'apolar-apolar' in line: metrics['aa']=float(val(line))
        elif 'apolar-polar' in line: metrics['ap']=float(val(line))
        elif 'polar-polar' in line: metrics['pp']=float(val(line))
        elif 'apolar NIS' in line: metrics['nis_apolar']=float(val(line))
        elif 'charged NIS' in line: metrics['nis_charged']=float(val(line))
        elif 'binding affinity' in line: metrics['dG']=float(val(line))
        elif 'dissociation constant' in line: metrics['Kd']=val(line)
    ic_file = os.path.join(workdir, base+'.ic')
    return metrics, ic_file, out

def write_output_csv(metrics, path, temp=25.0):
    rows = [
      ('No. of intermolecular contacts', metrics['contacts']),
      ('No. of charged-charged contacts', metrics['cc']),
      ('No. of charged-polar contacts', metrics['cp']),
      ('No. of charged-apolar contacts', metrics['ca']),
      ('No. of polar-polar contacts', metrics['pp']),
      ('No. of apolar-polar contacts', metrics['ap']),
      ('No. of apolar-apolar contacts', metrics['aa']),
      ('Percentage of apolar NIS residues', metrics['nis_apolar']),
      ('Percentage of charged NIS residues', metrics['nis_charged']),
      ('Predicted binding affinity (kcal.mol-1)', metrics['dG']),
      (f'Predicted dissociation constant (M) at {temp}˚C', metrics['Kd']),
    ]
    with open(path,'w',newline='') as f:
        w=csv.writer(f); w.writerow(['Metric','Value'])
        for k,v in rows: w.writerow([k,v])

def parse_ic(ic_file, csv_out):
    contacts=[]
    with open(ic_file) as f:
        for line in f:
            p=line.split()
            if len(p)!=6: continue
            r1,n1,c1,r2,n2,c2=p
            contacts.append((r1,int(n1),c1,r2,int(n2),c2))
    with open(csv_out,'w',newline='') as f:
        w=csv.writer(f)
        w.writerow(['Residue1','Residue1_Number','Chain1','Residue2','Residue2_Number','Chain2'])
        for c in contacts: w.writerow(c)
    return contacts

def main():
    label, srcdir, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
    os.makedirs(outdir, exist_ok=True)
    all_models=[]
    for m in range(1,6):
        src=os.path.join(srcdir, f'model_{m}.pdb')
        if not os.path.exists(src):
            print('MISSING', src); continue
        mdir=os.path.join(outdir, f'model_{m}'); os.makedirs(mdir, exist_ok=True)
        strip=os.path.join(mdir, f'strip_{m}.pdb'); strip_caps(src, strip)
        metrics, ic_file, raw = run_prodigy(f'strip_{m}.pdb', mdir)
        write_output_csv(metrics, os.path.join(mdir,'output.csv'))
        contacts=parse_ic(ic_file, os.path.join(mdir,'ic.csv'))
        # BSA residue ranking (chain A): count distinct peptide residues touched
        bsa=collections.defaultdict(set); pep=collections.defaultdict(set)
        pep_count=collections.Counter()
        for r1,n1,c1,r2,n2,c2 in contacts:
            # chain A = BSA, chain B = peptide
            if c1=='A' and c2=='B':
                bsa[(r1,n1)].add((r2,n2)); pep[(r2,n2)].add((r1,n1)); pep_count[(r2,n2)]+=1
        bsa_rank=sorted(bsa.items(), key=lambda kv:(-len(kv[1]), kv[0][1]))
        pep_rank=sorted(pep_count.items(), key=lambda kv:(-kv[1], kv[0][1]))
        all_models.append({
            'model':m,'metrics':{**metrics},
            'bsa_rank':[{'res':f'{r}{n}','count':len(s),
                         'peptide':sorted([f'{pr}{pn}' for pr,pn in s], key=lambda x:int(''.join(ch for ch in x if ch.isdigit())))}
                        for (r,n),s in bsa_rank],
            'pep_rank':[{'res':f'{r}{n}','count':c} for (r,n),c in pep_rank],
        })
    # consensus across models: BSA residues appearing anywhere in each model's interface
    appear=collections.Counter()
    for md in all_models:
        for e in md['bsa_rank']:
            appear[e['res']]+=1
    consensus={'ge3':sorted([r for r,c in appear.items() if c>=3]),
               'eq2':sorted([r for r,c in appear.items() if c==2])}
    summary={'label':label,'models':all_models,'consensus':consensus,
             'appear':dict(appear)}
    with open(os.path.join(outdir,'analysis.json'),'w') as f:
        json.dump(summary,f,indent=2)
    # ranking table
    print(f'\n===== {label} =====')
    print('Model  dG     Kd         Contacts  Apolar-apolar%  Hydrophobic%')
    for md in all_models:
        mt=md['metrics']; c=mt['contacts']
        hyd=100*(mt['ca']+mt['ap']+mt['aa'])/c if c else 0
        apo=100*mt['aa']/c if c else 0
        print(f"  {md['model']}   {mt['dG']:6.1f} {mt['Kd']:>10}  {int(c):5d}     {apo:5.1f}          {hyd:5.1f}")
    print('Consensus >=3:', consensus['ge3'])

if __name__=='__main__':
    main()
