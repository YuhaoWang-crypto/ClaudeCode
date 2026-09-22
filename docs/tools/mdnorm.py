"""Normalise hand-written markdown so python-markdown (and GitHub) parse tables and lists as blocks.
- blank line before a table block that follows a paragraph / heading line
- blank line before a list block that follows a paragraph / heading line
- two/three-space nested list indentation -> four spaces
"""
import re
LIST = re.compile(r'^(\s*)([-*]|\d+\.)\s+')
TABLE = re.compile(r'^\s*\|')

def normalise(text):
    out = []
    lines = text.split('\n')
    in_fence = False
    for i, l in enumerate(lines):
        if l.strip().startswith('```'):
            in_fence = not in_fence
            out.append(l); continue
        if in_fence:
            out.append(l); continue
        prev = out[-1] if out else ''
        m = LIST.match(l)
        if m:
            ind = m.group(1)
            if 1 <= len(ind) <= 3:
                l = '    ' + l[len(ind):]
                m = LIST.match(l)
            # blank line before a top-level list that follows prose/heading
            if not m.group(1) and prev.strip() and not LIST.match(prev) and not TABLE.match(prev) and not prev.startswith('    '):
                out.append('')
        elif TABLE.match(l):
            if prev.strip() and not TABLE.match(prev):
                out.append('')
        out.append(l)
    return '\n'.join(out)

if __name__ == '__main__':
    import sys
    sys.stdout.write(normalise(open(sys.argv[1], encoding='utf-8').read()))
