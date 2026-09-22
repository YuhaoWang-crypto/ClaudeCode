"""Load the hand-curated 25-analyte atlas (docs/fda_510k_sections/*.md) into curated_target / curated_section,
and link each target to catalog marker_ids by label search and to product codes stated in its §1.

Usage: python3 db/load_curated.py [db path]
"""
import sys, os, re, sqlite3, glob
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
DB = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'fda_ivd_markers.sqlite')
SEC_DIR = os.path.join(ROOT, 'docs', 'fda_510k_sections')
ATLAS = os.path.join(ROOT, 'docs', 'FDA_510k_IVD_cutoff_atlas.md')
GROUPS = {'A': '心脏标志物', 'B': '感染 / 炎症', 'C': '肿瘤标志物', 'D': '内分泌 / 代谢', 'E': '自身免疫', 'F': '感染性疾病血清学'}
# target title keyword -> marker label search terms (SQL LIKE, case-insensitive)
MARKER_TERMS = [  # (title keyword, label regexes (case-insensitive, word-bounded), key in the §2.1 archetype table)
    ('肌钙蛋白', [r'troponin'], 'hs-cTnI'), ('利钠肽', [r'natriuretic', r'\bBNP\b'], 'NT-proBNP'), ('C 反应蛋白', [r'C[- ]reactive'], 'hsCRP'),
    ('降钙素原', [r'procalcitonin'], 'PCT'), ('钙卫蛋白', [r'calprotectin'], '钙卫蛋白'), ('D-二聚体', [r'D[- ]dimer'], 'D-dimer'),
    ('癌胚抗原', [r'carcinoembryonic', r'\bCEA\b'], 'CEA'), ('糖类抗原 125', [r'CA[- ]?125'], 'CA 125'), ('糖类抗原 19-9', [r'CA[- ]?19-9'], 'CA 19-9'),
    ('前列腺特异抗原', [r'prostate[- ]specific antigen', r'\bPSA\b'], 'PSA'), ('甲胎蛋白', [r'alpha[- ]?fetoprotein', r'\bAFP\b'], 'AFP'),
    ('人附睾蛋白 4', [r'\bHE4\b', r'epididym'], 'HE4'), ('甲状腺球蛋白', [r'thyroglobulin'], 'Tg'), ('糖化血红蛋白', [r'\bA1c\b', r'glycated', r'glycosylated h'], 'HbA1c'),
    ('促甲状腺激素', [r'thyroid[- ]stimulating', r'thyrotropin', r'\bTSH\b'], 'TSH'), ('25-羟维生素 D', [r'25[- ]?(?:hydroxy|OH)', r'hydroxyvitamin D'], '25-OH VitD'),
    ('铁蛋白', [r'ferritin'], 'Ferritin'), ('胱抑素 C', [r'cystatin'], 'Cystatin C'), ('抗环瓜氨酸肽抗体', [r'citrullinated', r'\bCCP\b'], 'anti-CCP'),
    ('抗核抗体', [r'antinuclear', r'\bANA\b'], 'ANA'), ('可提取核抗原', [r'dsDNA', r'double[- ]stranded DNA', r'extractable nuclear', r'\bSm\b', r'SS-?A\b', r'SS-?B\b', r'Scl-?70', r'Jo-?1\b', r'centromere', r'\bRNP\b', r'\bENA\b'], 'ENA'),
    ('单纯疱疹病毒 2', [r'herpes simplex'], 'HSV-2'), ('莱姆病', [r'borrelia', r'\blyme'], 'Lyme'), ('梅毒', [r'treponema', r'syphilis'], '梅毒'), ('巨细胞病毒', [r'cytomegalovirus', r'\bCMV\b'], 'CMV'),
]

def archetype_table(md):
    """Parse §2.1 table: 靶点 | 原型 | 代表性 cutoff | 备注 -> list of tuples."""
    rows = []
    m = re.search(r'### 2\.1[^\n]*\n(.*?)\n\n', md, flags=re.S)
    if not m: return rows
    for line in m.group(1).split('\n'):
        if line.startswith('|') and not re.match(r'^\|\s*-', line) and '靶点' not in line:
            cells = [c.strip() for c in line.strip('|').split('|')]
            if len(cells) >= 3: rows.append(cells)
    return rows

def main():
    con = sqlite3.connect(DB); cur = con.cursor()
    cur.execute('DELETE FROM curated_target'); cur.execute('DELETE FROM curated_section')
    atlas = open(ATLAS, encoding='utf-8').read() if os.path.exists(ATLAS) else ''
    arche = archetype_table(atlas)
    labels = cur.execute('SELECT marker_id, label FROM marker').fetchall()
    n = 0
    for f in sorted(glob.glob(os.path.join(SEC_DIR, '*.md'))):
        letter = os.path.basename(f)[0]
        t = open(f, encoding='utf-8').read()
        parts = re.split(r'^## ([^\n]+)\n', t, flags=re.M)
        k = 0
        for j in range(1, len(parts), 2):
            title = parts[j].strip().replace('*', ''); body = parts[j+1]
            if title.startswith('附'): continue
            k += 1; tid = f't-{letter}-{k}'
            subs = re.split(r'^### ([^\n]+)\n', body, flags=re.M)
            sec1 = subs[2] if len(subs) > 2 else ''
            codes = sorted(set(re.findall(r'\*\*([A-Z]{3})\*\*', sec1)) | set(re.findall(r'[Pp]roduct [Cc]odes?[：:]?\s*\**([A-Z]{3})', sec1)))
            for kk in range(1, len(subs), 2):
                st = subs[kk].strip(); num = re.match(r'(\d+)\.', st)
                cur.execute('INSERT OR REPLACE INTO curated_section VALUES (?,?,?,?)', (tid, int(num.group(1)) if num else kk, st, subs[kk+1].strip()))
            # archetype / representative cutoff from the overview table (match by keyword)
            arch, rep, mids = '', '', []
            for kw, terms, akey in MARKER_TERMS:
                if kw in title:
                    for cells in arche:
                        if akey in cells[0].replace('*', ''):
                            arch, rep = cells[1], cells[2]; break
                    for mid, label in labels:
                        if any(re.search(t, label, flags=re.I) for t in terms): mids.append(mid)
                    break
            mids = sorted(set(mids))
            cur.execute('INSERT OR REPLACE INTO curated_target VALUES (?,?,?,?,?,?,?,?,?)',
                        (tid, letter, GROUPS[letter], title, '; '.join(codes), arch, rep, '; '.join(mids), 'atlas_2026-09-22'))
            n += 1
    con.commit()
    print('curated targets', n, 'sections', cur.execute('SELECT COUNT(*) FROM curated_section').fetchone()[0])
    for r in cur.execute('SELECT target_id, title, product_codes, archetype, marker_ids FROM curated_target').fetchall():
        print('  ', r[0], r[1][:34], '|', r[2], '|', r[3], '| markers', len(r[4].split('; ')) if r[4] else 0)
    con.close()

if __name__ == '__main__':
    main()
