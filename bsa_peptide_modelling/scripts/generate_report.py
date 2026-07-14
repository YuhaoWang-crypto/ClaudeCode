#!/usr/bin/env python3
"""Generate the BSA-peptide modelling PDF report (CCB-format) for the new peptides.
Usage: generate_report.py <out_pdf>
Reads results/<label>/analysis.json and results/<label>/model_5/... figures.
"""
import json, os, collections
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib import colors
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
                                Image, PageBreak, HRFlowable)
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_CENTER, TA_LEFT

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
NAVY = colors.HexColor('#16436e')
STEEL = colors.HexColor('#2e6fb0')
ORANGE = colors.HexColor('#e8833a')
LIGHT = colors.HexColor('#eef2f7')
GREY = colors.HexColor('#5b6472')

styles = getSampleStyleSheet()
def S(name,**kw):
    base=kw.pop('parent','Normal')
    parent = base if isinstance(base, ParagraphStyle) else styles[base]
    return ParagraphStyle(name,parent=parent,**kw)
H1=S('H1',parent='Heading1',fontSize=17,textColor=NAVY,spaceAfter=8,spaceBefore=4)
H2=S('H2',parent='Heading2',fontSize=13,textColor=NAVY,spaceAfter=5,spaceBefore=10)
H3=S('H3',parent='Heading3',fontSize=11,textColor=STEEL,spaceAfter=3,spaceBefore=6)
BODY=S('BODY',fontSize=9.5,leading=13,textColor=colors.HexColor('#222'),spaceAfter=5)
SMALL=S('SMALL',fontSize=8,leading=10,textColor=GREY)
CAP=S('CAP',fontSize=8,leading=10,textColor=GREY,alignment=TA_CENTER,spaceAfter=8)
BULLET=S('BULLET',parent=BODY,leftIndent=12,bulletIndent=2,spaceAfter=2)

PEPTIDES=[
 {'label':'peptideA_CFAGTPSILKKNGGGS','seq':'CFAGTPSILKKNGGGS','name':'Peptide A',
  'construct':'Biotin–CFAGTPSILKKNGGGS–NH2',
  'diff':'Relative to the original (CFAGTPSILMLAGGGS): M10K, L11K, A12N — three central substitutions introducing two lysine residues (net +2 charge) into the C-terminal hydrophobic segment.'},
 {'label':'peptideB_KKAGTPSILMLAGGGS','seq':'KKAGTPSILMLAGGGS','name':'Peptide B',
  'construct':'Biotin–(K)KAGTPSILMLAGGGS–NH2',
  'diff':'Relative to the original (CFAGTPSILMLAGGGS): C1K, F2K — two N-terminal substitutions replacing Cys/Phe with two lysines (net +2 charge). As in the original protocol the residue-1 position is represented by the biotin cap, so the modeled construct is Biotin–K–AGTPSILMLAGGG–NH2 (first Lys represented by biotin).'},
]

def hyd_pct(m):
    c=m['contacts'];
    return 100*(m['ca']+m['ap']+m['aa'])/c if c else 0
def apo_pct(m):
    c=m['contacts']
    return 100*m['aa']/c if c else 0

def th(t):  # header style helper
    return [('BACKGROUND',(0,0),(-1,0),NAVY),('TEXTCOLOR',(0,0),(-1,0),colors.white),
            ('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),('FONTSIZE',(0,0),(-1,-1),8),
            ('GRID',(0,0),(-1,-1),0.4,colors.HexColor('#c3ccd8')),
            ('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,LIGHT]),
            ('VALIGN',(0,0),(-1,-1),'MIDDLE'),('ALIGN',(1,0),(-1,-1),'CENTER'),
            ('LEFTPADDING',(0,0),(-1,-1),4),('RIGHTPADDING',(0,0),(-1,-1),4),
            ('TOPPADDING',(0,0),(-1,-1),3),('BOTTOMPADDING',(0,0),(-1,-1),3)]

def ranking_table(models):
    data=[['Model','ΔG (kcal/mol)','Kd (M @25°C)','Interfacial\ncontacts','Hydrophobic\ncontacts %','Apolar–apolar %']]
    best=sorted(models,key=lambda m:m['metrics']['dG'])
    for m in models:
        mt=m['metrics']
        data.append([f"Model {m['model']}", f"{mt['dG']:.1f}", f"{mt['Kd']}",
                     f"{int(mt['contacts'])}", f"{hyd_pct(mt):.1f}", f"{apo_pct(mt):.1f}"])
    t=Table(data,colWidths=[24*mm,26*mm,30*mm,22*mm,26*mm,24*mm])
    st=th(t)
    # highlight best (lowest dG) row
    bi=models.index(best[0])+1
    st.append(('BACKGROUND',(0,bi),(-1,bi),colors.HexColor('#d9ead0')))
    st.append(('FONTNAME',(0,bi),(-1,bi),'Helvetica-Bold'))
    t.setStyle(TableStyle(st)); return t,best[0]

def bsa_table(model):
    data=[['Rank','BSA residue','Contacts','Peptide residues contacted']]
    for i,e in enumerate(model['bsa_rank'][:20],1):
        data.append([str(i),e['res'],str(e['count']),', '.join(e['peptide'])])
    t=Table(data,colWidths=[12*mm,26*mm,20*mm,110*mm])
    st=th(t); st.append(('ALIGN',(3,1),(3,-1),'LEFT')); st.append(('FONTSIZE',(0,0),(-1,-1),7.5))
    t.setStyle(TableStyle(st)); return t

def build_peptide_section(story, pep):
    with open(os.path.join(RESULTS,pep['label'],'analysis.json')) as f:
        d=json.load(f)
    models=d['models']
    story.append(Paragraph(f"{pep['name']} &nbsp;&mdash;&nbsp; {pep['seq']}", H1))
    story.append(HRFlowable(width='100%',thickness=1.4,color=NAVY,spaceAfter=8))
    story.append(Paragraph('Construct modeled', H3))
    story.append(Paragraph(f"<b>{pep['construct']}</b> (16-mer, N-biotinylated, C-amidated) bound to bovine serum albumin (BSA, 66.5 kDa, 583 aa).", BODY))
    story.append(Paragraph(pep['diff'], BODY))
    story.append(Paragraph('Method', H3))
    story.append(Paragraph("Structure: Boltz-2.1 co-folding of the BSA–peptide complex (5 sampled models), biotin represented by CCD ligand BTN at residue 1 and C-terminal amidation by CCD NH2 at residue 16. Binding affinity / interface contacts: PRODIGY (contact-based predictor, 5.5 Å cutoff, 25 °C). Identical protocol to project CCB100725-YW01.", BODY))

    # 1) ranking
    story.append(Paragraph('1) Model ranking and interface metrics', H2))
    t,best=ranking_table(models)
    story.append(t)
    story.append(Paragraph("Hydrophobic contacts % = (charged–apolar + apolar–polar + apolar–apolar)/total. Best predicted binder (lowest ΔG) highlighted.", SMALL))

    # figure of best model
    fig=os.path.join(RESULTS,pep['label'],f"model_{best['model']}",'figure.png')
    if os.path.exists(fig):
        story.append(Spacer(1,4))
        story.append(Image(fig,width=135*mm,height=108*mm))
        story.append(Paragraph(f"Predicted BSA–peptide complex, best model (Model {best['model']}, ΔG {best['metrics']['dG']:.1f} kcal/mol). BSA backbone in grey; peptide (chain B) in blue; highlighted spheres = BSA interface residues.", CAP))

    # 2) core interface residues (best + one more)
    story.append(PageBreak())
    story.append(Paragraph('2) Core interface residues (top-20 per model)', H2))
    story.append(Paragraph('Columns: BSA residue (chain A), number of contacts to the peptide, and which peptide residues it touches.', SMALL))
    for m in models:
        story.append(Paragraph(f"Model {m['model']} &nbsp; (ΔG {m['metrics']['dG']:.1f} kcal/mol, {int(m['metrics']['contacts'])} contacts)", H3))
        story.append(bsa_table(m))
        story.append(Spacer(1,6))

    # 3) consensus
    story.append(Paragraph('3) BSA residues recurring across models', H2))
    appear=d['appear']
    ge4=sorted([r for r,c in appear.items() if c>=4], key=lambda r:int(''.join(ch for ch in r if ch.isdigit())))
    ge3=sorted([r for r,c in appear.items() if c==3], key=lambda r:int(''.join(ch for ch in r if ch.isdigit())))
    eq2=sorted([r for r,c in appear.items() if c==2], key=lambda r:int(''.join(ch for ch in r if ch.isdigit())))
    story.append(Paragraph(f"<b>Seen in all/≥4 of 5 models (strong consensus):</b> {', '.join(ge4) if ge4 else '—'}", BODY))
    story.append(Paragraph(f"<b>Seen in 3 of 5 models:</b> {', '.join(ge3) if ge3 else '—'}", BODY))
    story.append(Paragraph(f"<b>Seen in 2 of 5 models:</b> {', '.join(eq2) if eq2 else '—'}", BODY))

    # hotspot bins
    nums=collections.Counter()
    binhits=collections.defaultdict(list)
    for r,c in appear.items():
        if c>=2:
            n=int(''.join(ch for ch in r if ch.isdigit()))
            lo=(n//30)*30
            binhits[lo].append((n,r))
    story.append(Paragraph('Interface hotspot regions on BSA (residue-number bins, residues seen in ≥2 models):', H3))
    for lo in sorted(binhits):
        rs=sorted(binhits[lo])
        story.append(Paragraph(f"• {lo}–{lo+29}: {', '.join(r for _,r in rs)}", BULLET))

    # 4) peptide residues
    story.append(Paragraph('4) Which peptide residues contact BSA most', H2))
    story.append(Paragraph('Top peptide residues by interface contact count within each model.', SMALL))
    pdata=[['Model']+[f'#{i}' for i in range(1,7)]]
    for m in models:
        row=[f"Model {m['model']}"]
        for e in m['pep_rank'][:6]:
            row.append(f"{e['res']} ({e['count']})")
        while len(row)<7: row.append('—')
        pdata.append(row)
    pt=Table(pdata,colWidths=[24*mm]+[24*mm]*6)
    pt.setStyle(TableStyle(th(pt))); story.append(pt)
    story.append(Spacer(1,10))

def main():
    import sys
    out=sys.argv[1]
    doc=SimpleDocTemplate(out,pagesize=A4,topMargin=16*mm,bottomMargin=16*mm,leftMargin=15*mm,rightMargin=15*mm,
                          title='BSA peptide binding modelling report')
    story=[]
    # cover
    story.append(Spacer(1,30))
    story.append(Paragraph('PROJECT REPORT', S('cover',parent='Title',fontSize=26,textColor=NAVY,alignment=TA_CENTER)))
    story.append(Spacer(1,6))
    story.append(Paragraph('Modelling of binding of two 16-amino-acid biotinylated peptides to BSA',
                           S('sub',parent='Normal',fontSize=13,alignment=TA_CENTER,textColor=STEEL)))
    story.append(Spacer(1,14))
    story.append(HRFlowable(width='60%',thickness=1.2,color=ORANGE,hAlign='CENTER'))
    story.append(Spacer(1,14))
    meta=[['Peptides', 'A: CFAGTPSILKKNGGGS   |   B: KKAGTPSILMLAGGGS'],
          ['Target','Bovine serum albumin (BSA), 583 aa, 66.5 kDa'],
          ['Modifications','N-biotin (BTN, res 1), C-amide (NH2, res 16)'],
          ['Structure prediction','Boltz-2.1 co-folding, 5 models per peptide'],
          ['Binding affinity','PRODIGY (contact-based ΔG / Kd, 25 °C)'],
          ['Reference project','CCB100725-YW01 (identical protocol)'],
          ['Date','2026-07-14']]
    mt=Table(meta,colWidths=[45*mm,120*mm])
    mt.setStyle(TableStyle([('FONTSIZE',(0,0),(-1,-1),9.5),('FONTNAME',(0,0),(0,-1),'Helvetica-Bold'),
        ('TEXTCOLOR',(0,0),(0,-1),NAVY),('BOTTOMPADDING',(0,0),(-1,-1),5),('TOPPADDING',(0,0),(-1,-1),5),
        ('LINEBELOW',(0,0),(-1,-1),0.3,colors.HexColor('#d5dbe4')),('VALIGN',(0,0),(-1,-1),'MIDDLE')]))
    story.append(mt)
    story.append(Spacer(1,20))
    story.append(Paragraph('Executive summary', H2))
    # summary from data
    summ=[]
    for pep in PEPTIDES:
        with open(os.path.join(RESULTS,pep['label'],'analysis.json')) as f:
            d=json.load(f)
        best=min(d['models'],key=lambda m:m['metrics']['dG'])
        dgs=[m['metrics']['dG'] for m in d['models']]
        summ.append((pep,best,dgs))
    for pep,best,dgs in summ:
        story.append(Paragraph(f"<b>{pep['name']} ({pep['seq']}):</b> predicted ΔG ranges {max(dgs):.1f} to {min(dgs):.1f} kcal/mol across 5 models; best model = Model {best['model']} (ΔG {best['metrics']['dG']:.1f} kcal/mol, Kd {best['metrics']['Kd']} M, {int(best['metrics']['contacts'])} interfacial contacts).", BODY))
    story.append(Paragraph("Full per-model tables, interface-residue rankings, consensus hotspots and structure figures follow. A cross-comparison with the original peptide is given at the end.", BODY))
    story.append(PageBreak())

    for pep in PEPTIDES:
        build_peptide_section(story, pep)
        story.append(PageBreak())

    # comparison section built by caller-added file if present
    comp=os.path.join(RESULTS,'comparison.json')
    if os.path.exists(comp):
        add_comparison(story, comp)

    doc.build(story)
    print('wrote',out)

def add_comparison(story, comp):
    d=json.load(open(comp))
    story.append(Paragraph('5) Cross-comparison: new peptides vs. original', H1))
    story.append(HRFlowable(width='100%',thickness=1.4,color=NAVY,spaceAfter=8))
    data=[['Peptide','Sequence','Best ΔG\n(kcal/mol)','Best Kd (M)','Mean ΔG','Max contacts']]
    for r in d['rows']:
        data.append([r['name'],r['seq'],f"{r['best_dg']:.1f}",r['best_kd'],f"{r['mean_dg']:.1f}",str(r['max_contacts'])])
    t=Table(data,colWidths=[26*mm,44*mm,22*mm,26*mm,20*mm,22*mm])
    t.setStyle(TableStyle(th(t))); story.append(t)
    story.append(Spacer(1,8))
    for para in d.get('notes',[]):
        story.append(Paragraph(para, BODY))

if __name__=='__main__': main()
