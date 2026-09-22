"""Rule-based structuring of FDA 510(k) decision summaries (OIVD templates, 2003-2026).

Two heading templates are recognised:
  old  (≈2003-2018): A. 510(k) Number … H. Intended Use … K. Standard/Guidance … M. Performance Characteristics
                     1. Analytical performance a. Precision … f. Assay cut-off; 2. Comparison studies; 3. Clinical studies;
                     4. Clinical cut-off; 5. Expected values/Reference range
  new  (≈2018-2026): I Background … III Intended Use/Indications … VI Standards/Guidance … VII Performance Characteristics
                     A Analytical … 7. Assay Cut-Off … B Comparison … C Clinical Studies … D Clinical Cut-Off … E Expected Values/Reference Range

For each summary we emit a dict of fields (text of the section, trimmed) plus derived fields:
  measurand, type_of_test, intended_use, indications, special_conditions, instrument, predicate, standards (list),
  precision, linearity, traceability, detection_limit, analytical_specificity, assay_cutoff, method_comparison,
  matrix_comparison, clinical_sensitivity, clinical_specificity, clinical_studies, clinical_cutoff, reference_range,
  specimen_types (derived), clsi_codes (derived), cutoff_numbers (derived), sample_n (derived), conclusion
Confidence is 'section' when taken from a recognised heading, 'regex' when derived.
"""
import re, json, sys, os

# heading label -> canonical field. Matched case-insensitively on the heading text after the enumerator.
HEADINGS = [
    (r'measurand', 'measurand'),
    (r'type of test', 'type_of_test'),
    (r'purpose for submission', 'purpose'),
    (r'intended use\(?s?\)?:?$|intended use/indications for use', 'intended_use'),
    (r'indication\(?s?\)? for use', 'indications'),
    (r'special conditions? for use', 'special_conditions'),
    (r'special instrument requirements?', 'instrument'),
    (r'device description', 'device_description'),
    (r'predicate device name', 'predicate'),
    (r'method comparison', 'method_comparison'),
    (r'^comparison with predicate', 'predicate_comparison'),
    (r'standards?/guidance documents? referenced', 'standards'),
    (r'test principle|principle of operation', 'test_principle'),
    (r'precision ?/ ?reproducibility|precision', 'precision'),
    (r'linearity', 'linearity'),
    (r'traceability|stability|expected values \(controls', 'traceability'),
    (r'detection limit|limit of detection', 'detection_limit'),
    (r'analytical specificity', 'analytical_specificity'),
    (r'assay cut-? ?off', 'assay_cutoff'),
    (r'assay reportable range|reportable range', 'reportable_range'),
    (r'matrix comparison', 'matrix_comparison'),
    (r'clinical sensitivity', 'clinical_sensitivity'),
    (r'clinical specificity', 'clinical_specificity'),
    (r'other clinical supportive data|clinical studies', 'clinical_studies'),
    (r'clinical cut-? ?off', 'clinical_cutoff'),
    (r'expected values ?/ ?reference range|expected values|reference range', 'reference_range'),
    (r'instrument name', 'instrument_name'),
    (r'conclusion', 'conclusion'),
]
ENUM = r'^\s*(?:[A-Z]|[IVX]{1,4}|\d{1,2}|[a-z])[\.\)]?\s+'
HEAD_RE = re.compile(ENUM + r'([A-Za-z][A-Za-z/ \-\(\)\.,&]{2,80}?):?\s*$')

def split_docs(txt):
    """Return (decision_summary_text, summary_text) from a stored txt file."""
    ds, sm = '', ''
    for block in re.split(r'^===== ', txt, flags=re.M):
        if block.startswith('DECISION_SUMMARY'): ds = block.split('\n', 1)[1] if '\n' in block else ''
        elif block.startswith('510K_SUMMARY') or block.startswith('DENOVO_SUMMARY'): sm = block.split('\n', 1)[1] if '\n' in block else ''
    return ds, sm

def detect_template(ds):
    if re.search(r'^\s*VII\s+Performance Characteristics', ds, flags=re.M): return 'new'
    if re.search(r'^\s*M\.\s+Performance Characteristics', ds, flags=re.M): return 'old'
    if re.search(r'^\s*L\.\s+Performance Characteristics', ds, flags=re.M): return 'denovo'
    if re.search(r'Performance Characteristics', ds): return 'other'
    return 'none'

def sections(ds):
    """Yield (field, text) by scanning heading lines; the text runs to the next recognised heading."""
    lines = ds.split('\n')
    marks = []
    for i, l in enumerate(lines):
        s = re.sub(r'\s+', ' ', l.strip())
        m = HEAD_RE.match(s)
        if not m: continue
        label = m.group(1).strip().lower()
        for pat, field in HEADINGS:
            if re.search(pat, label):
                marks.append((i, field)); break
    out = {}
    for j, (i, field) in enumerate(marks):
        end = marks[j+1][0] if j + 1 < len(marks) else len(lines)
        body = '\n'.join(lines[i+1:end]).strip()
        body = re.sub(r'\n?K\d{6} - Page \d+ of \d+\n?', '\n', body)
        body = re.sub(r'\n{3,}', '\n\n', body)
        if field in out and out[field]: out[field] += '\n' + body
        else: out[field] = body
    return out

SPECIMEN_KW = [
    ('serum', r'\bser(?:um|a)\b'), ('EDTA plasma', r'\bEDTA\b'), ('heparin plasma', r'heparin'), ('citrate plasma', r'citrat'),
    ('whole blood', r'whole blood'), ('capillary/fingerstick', r'capillary|finger ?stick'), ('urine', r'\burine\b'),
    ('stool', r'\bstool\b|\bfec(?:al|es)\b'), ('CSF', r'cerebrospinal|\bCSF\b'), ('swab', r'\bswab'), ('saliva/oral fluid', r'saliva|oral fluid'),
    ('plasma', r'\bplasma\b'), ('dried blood spot', r'dried blood spot'), ('tissue/FFPE', r'\bFFPE\b|formalin'), ('sputum', r'sputum'), ('synovial fluid', r'synovial'),
]

def derive(fields, ds):
    d = {}
    iu = ' '.join(filter(None, [fields.get('intended_use', ''), fields.get('indications', '')]))
    found = []
    for name, pat in SPECIMEN_KW:
        if re.search(pat, iu, flags=re.I): found.append(name)
    if 'plasma' in found and any(x in found for x in ('EDTA plasma', 'heparin plasma', 'citrate plasma')): found.remove('plasma')
    d['specimen_types'] = '; '.join(found)
    std = fields.get('standards', '') or ''
    codes = sorted(set(re.findall(r'\b(?:CLSI|NCCLS)?\s*((?:EP|C|H|M|GP|I/LA|ILA|MM|POCT)\d{1,2}[A-Za-z0-9\-]*)', std)))
    iso = sorted(set(re.findall(r'\bISO\s*[0-9]{4,5}(?:-\d+)?', std)))
    d['clsi_codes'] = '; '.join(codes + iso)
    cut = ' '.join(filter(None, [fields.get('assay_cutoff', ''), fields.get('clinical_cutoff', '')]))
    nums = re.findall(r'(?:[<>≤≥=]\s*)?\d+(?:[.,]\d+)?\s*(?:ng/mL|ng/L|pg/mL|µg/L|ug/L|mg/L|mg/dL|g/dL|U/mL|IU/mL|kIU/L|kU/L|mIU/L|µIU/mL|uIU/mL|pmol/L|nmol/L|µmol/L|umol/L|mmol/L|µg/g|ug/g|mcg/g|%|index|AU/mL|RU/mL|CU|LIU|S/CO|copies/mL|IU/L|U/L)', cut, flags=re.I)
    d['cutoff_numbers'] = '; '.join(dict.fromkeys(n.strip() for n in nums))[:500]
    clin = ' '.join(filter(None, [fields.get('clinical_studies', ''), fields.get('clinical_sensitivity', ''), fields.get('clinical_specificity', '')]))
    ns = re.findall(r'\b[nN]\s*=\s*(\d{2,5})\b', clin)
    d['sample_n'] = '; '.join(dict.fromkeys(ns))[:200]
    sens = re.findall(r'(?:sensitivity|PPA|positive percent agreement)[^%\n]{0,90}?(\d{1,3}(?:\.\d+)?)\s*%', clin, flags=re.I)
    spec = re.findall(r'(?:specificity|NPA|negative percent agreement)[^%\n]{0,90}?(\d{1,3}(?:\.\d+)?)\s*%', clin, flags=re.I)
    d['sens_pct'] = '; '.join(dict.fromkeys(sens))[:200]; d['spec_pct'] = '; '.join(dict.fromkeys(spec))[:200]
    mc = ' '.join(filter(None, [fields.get('method_comparison', ''), fields.get('predicate_comparison', '')]))
    sl = re.findall(r'slope[^0-9\-\n]{0,40}?(-?[01]\.\d{2,4})\b', mc, flags=re.I) + re.findall(r'y\s*=\s*(-?[01]\.\d{2,4})\s*\(?x', mc, flags=re.I)
    rr = re.findall(r'(?:\br\b|\bR\b|correlation coefficient|\br2\b|\bR²)[^0-9\n]{0,25}?(0\.\d{2,4}|1\.0+)\b', mc)
    d['mc_slope'] = '; '.join(dict.fromkeys(sl))[:120]; d['mc_r'] = '; '.join(dict.fromkeys(rr))[:120]
    # regulatory bits from the header block
    pc = re.search(r'Product [Cc]ode\(?s?\)?:?\s*\n?\s*([A-Z]{3})\b', ds) or re.search(r'\b([A-Z]{3})\s*[–-]\s*[A-Z][a-z]', ds[:3000])
    d['product_code_in_doc'] = pc.group(1) if pc else ''
    reg = re.search(r'21 CFR\s*§?\s*(\d{3}\.\d{4})', ds)
    d['regulation_in_doc'] = reg.group(1) if reg else ''
    return d

def parse_file(path):
    txt = open(path, encoding='utf-8', errors='ignore').read()
    ds, sm = split_docs(txt)
    tpl = detect_template(ds)
    fields = sections(ds) if ds else {}
    derived = derive(fields, ds) if ds else {}
    return tpl, fields, derived

if __name__ == '__main__':
    import glob, collections
    files = sorted(glob.glob(os.path.join(sys.argv[1], '*.txt')))
    cov = collections.Counter(); tpls = collections.Counter(); n = 0
    for f in files:
        tpl, fields, derived = parse_file(f)
        tpls[tpl] += 1
        if tpl == 'none': continue
        n += 1
        for k, v in fields.items():
            if v and v.strip() and not re.match(r'^(N/?A|Not applicable|None)\.?$', v.strip(), flags=re.I): cov[k] += 1
        for k, v in derived.items():
            if v: cov['~' + k] += 1
    print('files', len(files), 'templates', dict(tpls), 'parsed', n)
    for k, v in sorted(cov.items(), key=lambda x: -x[1]): print(f'  {k:28s} {v:5d} {100*v/max(n,1):5.1f}%')
    if len(sys.argv) > 2:
        tpl, fields, derived = parse_file(sys.argv[2]); print(json.dumps({'template': tpl, **{k: v[:300] for k, v in fields.items()}, **derived}, ensure_ascii=False, indent=1))
