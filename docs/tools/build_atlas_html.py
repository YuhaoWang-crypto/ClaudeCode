"""Build a single-page HTML reader from docs/FDA_510k_IVD_cutoff_atlas.md."""
import re, html, os, sys, json
import markdown

SRC = '/home/user/ClaudeCode/docs/FDA_510k_IVD_cutoff_atlas.md'
OUT = '/home/user/ClaudeCode/docs/FDA_510k_IVD_cutoff_atlas.html'
GROUP_COLORS = {'A': '#B5473A', 'B': '#C77C1E', 'C': '#7A5CB0', 'D': '#1F7A4D', 'E': '#1F6FA8', 'F': '#8A6A2A'}

md = open(SRC, encoding='utf-8').read()

def tag_markers(t):
    t = t.replace('文件未载明', '<span class="tag na">文件未载明</span>')
    t = t.replace('背景（非申报文件）', '<span class="tag bg">背景（非申报文件）</span>')
    return t

def conv(t):
    t = tag_markers(t)
    h = markdown.markdown(t, extensions=['tables', 'sane_lists'])
    h = re.sub(r'<table>', '<div class="tw"><table>', h)
    h = h.replace('</table>', '</table></div>')
    # K numbers -> mono chips
    h = re.sub(r'(?<![\w/])((?:K|DEN|P)\d{6})(?![\w/])', r'<code class="k">\1</code>', h)
    return h

# split overview / groups / tail
parts = re.split(r'^# 第 ([A-F]) 组 · ([^\n]+)\n', md, flags=re.M)
overview = parts[0]
groups = []
for i in range(1, len(parts), 3):
    groups.append((parts[i], parts[i+1].strip(), parts[i+2]))
# the tail (appendix) is at the end of the last group's text: split at '\n## 附录 A'
last_letter, last_title, last_body = groups[-1]
m = re.search(r'\n---\n\n## 附录 A', last_body)
tail = ''
if m:
    tail = last_body[m.start():]
    groups[-1] = (last_letter, last_title, last_body[:m.start()])

# overview: drop the H1 and the generated TOC block (we build our own nav)
overview = re.sub(r'^# [^\n]+\n', '', overview, count=1)
overview = re.sub(r'\n## 目录（分节）.*$', '\n', overview, flags=re.S)
ov_html = conv(overview)

nav = []
body = []
target_index = []  # for search
for letter, gtitle, gbody in groups:
    color = GROUP_COLORS[letter]
    tsplit = re.split(r'^## ([^\n]+)\n', gbody, flags=re.M)
    intro = tsplit[0].strip()
    nav_items = []
    sections_html = []
    n = 0
    for j in range(1, len(tsplit), 2):
        ttitle = tsplit[j].strip().replace('*', '')
        tbody = tsplit[j+1]
        if ttitle.startswith('附'):
            xid = f'x-{letter}'
            nav_items.append(f'<li><a href="#{xid}">{html.escape(ttitle.split("（")[0])}</a></li>')
            sections_html.append(f'<section class="target extra" id="{xid}"><header class="th"><span class="chip" style="--gc:{color}">{letter} 组</span><h2>{html.escape(ttitle)}</h2></header>{conv(tbody)}</section>')
            continue
        n += 1
        tid = f't-{letter}-{n}'
        short = re.sub(r'（.*$', '', ttitle).strip()
        nav_items.append(f'<li><a href="#{tid}" data-target="{tid}">{html.escape(short)}</a></li>')
        # subsections
        sub = re.split(r'^### ([^\n]+)\n', tbody, flags=re.M)
        lead = conv(sub[0]) if sub[0].strip() else ''
        subs = []
        for k in range(1, len(sub), 2):
            stitle = sub[k].strip()
            sbody = conv(sub[k+1])
            num = re.match(r'(\d+)\.', stitle)
            is_key = num and num.group(1) in ('4', '8')
            subs.append(f'<details class="sub{" key" if is_key else ""}" open><summary>{html.escape(stitle)}</summary><div class="sb">{sbody}</div></details>')
        sections_html.append(
            f'<section class="target" id="{tid}" data-group="{letter}">'
            f'<header class="th"><span class="chip" style="--gc:{color}">{letter} 组</span>'
            f'<h2>{html.escape(ttitle)}</h2>'
            f'<div class="tools"><button type="button" data-act="open">全部展开</button><button type="button" data-act="close">全部折叠</button></div></header>'
            f'{lead}{"".join(subs)}</section>')
        target_index.append({'id': tid, 'title': ttitle, 'group': letter, 'text': re.sub(r'\s+', ' ', tbody)[:200000]})
    gid = f'g-{letter}'
    nav.append(f'<li class="grp" style="--gc:{color}"><a href="#{gid}" class="gl"><span class="dot"></span>{letter}. {html.escape(gtitle)}</a><ul>{"".join(nav_items)}</ul></li>')
    body.append(f'<div class="group" id="{gid}" style="--gc:{color}"><h1 class="gh"><span class="dot"></span>第 {letter} 组 · {html.escape(gtitle)}</h1>{conv(intro) if intro else ""}{"".join(sections_html)}</div>')

tail_html = conv(tail.replace('\n---\n', '\n', 1)) if tail else ''

n_targets = len(target_index)
page = f'''<title>510(k) 靶点 Cutoff 图谱</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:ital,wght@0,400;0,500;0,600;1,400&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
:root{{
  --ground:#F6F7F6; --paper:#FFFFFF; --ink:#1B2229; --ink-2:#4A555C; --ink-3:#7A858B;
  --line:#D9DFDF; --line-2:#EAEEEE; --accent:#146C7A; --accent-soft:#E1F0F2; --accent-ink:#0E4F5A;
  --na-bg:#FFF3D1; --na-ink:#7A5200; --bg-bg:#ECE8F8; --bg-ink:#4B3B8E; --key:#F3F7F7;
  --shadow:0 1px 2px rgba(20,40,45,.06);
  --sans:"IBM Plex Sans","PingFang SC","Hiragino Sans GB","Microsoft YaHei","Noto Sans CJK SC","Source Han Sans SC",system-ui,sans-serif;
  --mono:"IBM Plex Mono",ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;
}}
@media (prefers-color-scheme: dark){{ :root:not([data-theme="light"]){{
  color-scheme:dark; --ground:#0F1416; --paper:#151B1E; --ink:#E4E9E9; --ink-2:#AEB8BC; --ink-3:#7E8A8F;
  --line:#2A3438; --line-2:#1F2729; --accent:#5EC0CD; --accent-soft:#12333A; --accent-ink:#9EDDE5;
  --na-bg:#3A2E10; --na-ink:#F0C86A; --bg-bg:#2A2450; --bg-ink:#C9BDF5; --key:#131A1C; --shadow:none; }} }}
:root[data-theme="dark"]{{
  color-scheme:dark; --ground:#0F1416; --paper:#151B1E; --ink:#E4E9E9; --ink-2:#AEB8BC; --ink-3:#7E8A8F;
  --line:#2A3438; --line-2:#1F2729; --accent:#5EC0CD; --accent-soft:#12333A; --accent-ink:#9EDDE5;
  --na-bg:#3A2E10; --na-ink:#F0C86A; --bg-bg:#2A2450; --bg-ink:#C9BDF5; --key:#131A1C; --shadow:none; }}
*{{box-sizing:border-box}}
body{{margin:0;background:var(--ground);color:var(--ink);font-family:var(--sans);font-size:15px;line-height:1.65;-webkit-font-smoothing:antialiased}}
a{{color:var(--accent);text-decoration:none}} a:hover{{text-decoration:underline}}
.top{{position:sticky;top:env(safe-area-inset-top,0px);z-index:20;background:var(--paper);border-bottom:1px solid var(--line);display:flex;align-items:center;gap:14px;padding:10px 16px;flex-wrap:wrap}}
.top .brand{{font-weight:600;letter-spacing:.01em}} .top .brand small{{display:block;font-weight:400;color:var(--ink-3);font-size:12px}}
.top .search{{margin-left:auto;display:flex;gap:8px;align-items:center;flex:1 1 240px;max-width:460px}}
.top input{{flex:1;min-width:0;font:inherit;padding:7px 10px;border:1px solid var(--line);border-radius:6px;background:var(--ground);color:var(--ink)}}
.top input:focus{{outline:2px solid var(--accent);outline-offset:1px}}
.top .cnt{{font-size:12px;color:var(--ink-3);white-space:nowrap;font-variant-numeric:tabular-nums}}
.top button, .tools button{{font:inherit;font-size:12px;padding:5px 9px;border:1px solid var(--line);border-radius:6px;background:var(--paper);color:var(--ink-2);cursor:pointer}}
.top button:hover, .tools button:hover{{border-color:var(--accent);color:var(--accent)}}
.wrap{{display:grid;grid-template-columns:270px minmax(0,1fr);gap:0;max-width:1440px;margin:0 auto}}
nav.side{{position:sticky;top:calc(env(safe-area-inset-top,0px) + 58px);align-self:start;height:calc(100vh - 58px);overflow:auto;padding:14px 12px 40px 16px;border-right:1px solid var(--line-2);font-size:13.5px}}
nav.side ul{{list-style:none;margin:0;padding:0}}
nav.side > ul > li{{margin-bottom:10px}}
nav.side a{{display:block;padding:3px 8px;border-radius:5px;color:var(--ink-2)}}
nav.side a:hover{{background:var(--accent-soft);color:var(--accent-ink);text-decoration:none}}
nav.side a.on{{background:var(--accent-soft);color:var(--accent-ink);font-weight:500}}
nav.side .gl{{font-weight:600;color:var(--ink);display:flex;align-items:center;gap:7px}}
nav.side li.grp ul{{margin-left:6px;border-left:2px solid var(--gc,var(--line));padding-left:6px;margin-top:2px}}
nav.side li.hide{{display:none}}
.dot{{display:inline-block;width:9px;height:9px;border-radius:50%;background:var(--gc,var(--accent));flex:none}}
main{{padding:22px 28px 80px;min-width:0}}
.ov{{max-width:980px}}
.ov h2, .group h2{{font-size:20px;margin:34px 0 10px;text-wrap:balance}}
.ov h3{{font-size:16px;margin:24px 0 8px}}
.ov > p:first-child{{color:var(--ink-2)}}
blockquote{{margin:0 0 18px;padding:10px 16px;border-left:3px solid var(--accent);background:var(--paper);color:var(--ink-2);font-size:14px}}
blockquote p{{margin:6px 0}}
.gh{{font-size:24px;margin:56px 0 12px;padding-top:22px;border-top:2px solid var(--gc);display:flex;align-items:center;gap:10px;text-wrap:balance}}
.gh .dot{{width:12px;height:12px}}
.target{{background:var(--paper);border:1px solid var(--line);border-radius:8px;box-shadow:var(--shadow);margin:18px 0 26px;padding:6px 22px 14px;scroll-margin-top:70px}}
.th{{display:flex;align-items:center;gap:12px;flex-wrap:wrap;padding:12px 0 8px;border-bottom:1px solid var(--line-2);margin-bottom:8px}}
.th h2{{margin:0;font-size:19px;flex:1 1 300px}}
.chip{{font-size:11.5px;font-weight:600;letter-spacing:.04em;color:#fff;background:var(--gc);padding:2px 8px;border-radius:999px;white-space:nowrap}}
.tools{{display:flex;gap:6px}}
details.sub{{border-top:1px solid var(--line-2);padding:4px 0}}
details.sub summary{{cursor:pointer;font-weight:600;font-size:15px;padding:8px 0;color:var(--ink);list-style:none;display:flex;gap:8px;align-items:center}}
details.sub summary::-webkit-details-marker{{display:none}}
details.sub summary::before{{content:"";width:7px;height:7px;border-right:1.5px solid var(--ink-3);border-bottom:1.5px solid var(--ink-3);transform:rotate(-45deg);transition:transform .15s;flex:none}}
details.sub[open] summary::before{{transform:rotate(45deg)}}
details.sub.key summary{{color:var(--accent-ink)}}
details.sub.key .sb{{background:var(--key);border-radius:6px;padding:2px 14px 8px;margin-bottom:8px}}
.sb{{font-size:14.5px}} .sb > :first-child{{margin-top:4px}}
.sb h4{{font-size:14px;margin:16px 0 6px}}
p{{margin:8px 0}} ul,ol{{padding-left:22px;margin:6px 0}} li{{margin:3px 0}}
.tw{{overflow-x:auto;margin:10px 0 14px;border:1px solid var(--line);border-radius:6px}}
table{{border-collapse:collapse;font-size:13.2px;min-width:100%;font-variant-numeric:tabular-nums}}
th,td{{padding:6px 10px;border-bottom:1px solid var(--line-2);vertical-align:top;text-align:left}}
th{{background:var(--ground);font-weight:600;white-space:nowrap;position:sticky;top:0}}
tr:last-child td{{border-bottom:none}}
code{{font-family:var(--mono);font-size:.92em;background:var(--ground);padding:1px 4px;border-radius:4px}}
code.k{{color:var(--accent-ink);background:var(--accent-soft);white-space:nowrap}}
pre{{overflow-x:auto;background:var(--ground);border:1px solid var(--line);border-radius:6px;padding:10px 12px;font-size:13px}}
pre code{{background:none;padding:0}}
.tag{{display:inline-block;font-size:12px;line-height:1.3;padding:1px 6px;border-radius:4px;font-weight:500;white-space:nowrap}}
.tag.na{{background:var(--na-bg);color:var(--na-ink)}} .tag.bg{{background:var(--bg-bg);color:var(--bg-ink)}}
mark.hit{{background:var(--na-bg);color:inherit;padding:0 1px}}
.legend{{display:flex;gap:14px;flex-wrap:wrap;font-size:13px;color:var(--ink-2);margin:10px 0 0}}
.nores{{display:none;color:var(--ink-3);padding:6px 8px;font-size:13px}}
.menu-toggle{{display:none}}
:focus-visible{{outline:2px solid var(--accent);outline-offset:2px}}
@media (max-width:900px){{
  .wrap{{grid-template-columns:1fr}}
  nav.side{{position:static;height:auto;max-height:none;border-right:none;border-bottom:1px solid var(--line);padding:8px 16px 12px}}
  nav.side:not(.open){{display:none}}
  .menu-toggle{{display:inline-block}}
  main{{padding:16px 16px 60px}}
  .target{{padding:4px 14px 10px}}
}}
@media (prefers-reduced-motion:reduce){{ details.sub summary::before{{transition:none}} }}
</style>

<div class="top">
  <div class="brand">510(k) 靶点 Cutoff 图谱<small>FDA 决策摘要 · 25 个靶点 · 6 组 · 2026-09-22</small></div>
  <div class="search">
    <button type="button" class="menu-toggle" id="menu">目录</button>
    <input id="q" type="search" placeholder="搜靶点、K 号、厂家、cutoff…" aria-label="搜索">
    <span class="cnt" id="cnt">{n_targets} 个靶点</span>
  </div>
</div>
<div class="wrap">
<nav class="side" id="side">
  <ul>
    <li><a href="#overview" class="gl">总览与方法</a>
      <ul>
        <li><a href="#ov-0">0. 这份整理回答什么问题</a></li>
        <li><a href="#ov-1">1. 510(k) 与 PMA 边界</a></li>
        <li><a href="#ov-2">2. cutoff 逻辑的八种原型</a></li>
        <li><a href="#ov-21">2.1 靶点 × 原型速览表</a></li>
        <li><a href="#ov-3">3. 分析性能标准对照</a></li>
        <li><a href="#ov-4">4. 临床验证设计范式</a></li>
        <li><a href="#ov-5">5. 样本类型</a></li>
      </ul></li>
    {"".join(nav)}
    <li><a href="#appendix" class="gl">附录 · 复现方法与指南</a></li>
  </ul>
  <div class="nores" id="nores">没有匹配的靶点</div>
</nav>
<main>
  <div class="ov" id="overview">
    <div class="legend"><span><span class="tag na">文件未载明</span> FDA 文件里没有这项信息</span><span><span class="tag bg">背景（非申报文件）</span> 来自指南/文献，非 FDA 原文</span><span><code class="k">K123456</code> 510(k) 号，§10 有原文链接</span></div>
    {ov_html}
  </div>
  {"".join(body)}
  <div class="ov" id="appendix">{tail_html}</div>
</main>
</div>

<script>
(function(){{
  // anchor ids for overview h2
  var map={{'0.':'ov-0','1.':'ov-1','2.':'ov-2','3.':'ov-3','4.':'ov-4','5.':'ov-5'}};
  document.querySelectorAll('#overview h2').forEach(function(h){{var k=h.textContent.trim().slice(0,2); if(map[k]) h.id=map[k];}});
  document.querySelectorAll('#overview h3').forEach(function(h){{ if(h.textContent.trim().indexOf('2.1')===0) h.id='ov-21'; }});
  // expand/collapse per target
  document.querySelectorAll('.target .tools button').forEach(function(b){{
    b.addEventListener('click',function(){{
      var open=b.dataset.act==='open';
      b.closest('.target').querySelectorAll('details.sub').forEach(function(d){{d.open=open;}});
    }});
  }});
  // menu toggle (mobile)
  var side=document.getElementById('side');
  document.getElementById('menu').addEventListener('click',function(){{side.classList.toggle('open');}});
  side.addEventListener('click',function(e){{ if(e.target.tagName==='A' && window.innerWidth<=900) side.classList.remove('open'); }});
  // search: filter targets + nav, highlight hits inside open targets
  var q=document.getElementById('q'), cnt=document.getElementById('cnt'), nores=document.getElementById('nores');
  var targets=Array.prototype.slice.call(document.querySelectorAll('.target'));
  var total=targets.length;
  var cache=targets.map(function(t){{return t.textContent.toLowerCase();}});
  function clearMarks(root){{ root.querySelectorAll('mark.hit').forEach(function(m){{ m.replaceWith(document.createTextNode(m.textContent)); }}); root.normalize(); }}
  function markText(root, needle){{
    var walker=document.createTreeWalker(root, NodeFilter.SHOW_TEXT, null);
    var nodes=[]; while(walker.nextNode()) nodes.push(walker.currentNode);
    nodes.forEach(function(n){{
      var s=n.nodeValue, i=s.toLowerCase().indexOf(needle); if(i<0) return;
      var frag=document.createDocumentFragment(), pos=0;
      while(i>=0){{ frag.appendChild(document.createTextNode(s.slice(pos,i))); var m=document.createElement('mark'); m.className='hit'; m.textContent=s.slice(i,i+needle.length); frag.appendChild(m); pos=i+needle.length; i=s.toLowerCase().indexOf(needle,pos); }}
      frag.appendChild(document.createTextNode(s.slice(pos))); n.parentNode.replaceChild(frag,n);
    }});
  }}
  var timer=null;
  q.addEventListener('input',function(){{ clearTimeout(timer); timer=setTimeout(run,180); }});
  function run(){{
    var v=q.value.trim().toLowerCase(); var shown=0;
    targets.forEach(function(t,i){{
      clearMarks(t);
      var hit=!v || cache[i].indexOf(v)>=0;
      t.hidden=!hit;
      var navA=side.querySelector('a[data-target="'+t.id+'"]'); if(navA) navA.parentElement.classList.toggle('hide',!hit);
      if(hit){{ shown++; if(v && v.length>=2) markText(t,v); }}
    }});
    document.querySelectorAll('.group').forEach(function(g){{ g.hidden = v && !g.querySelector('.target:not([hidden])'); }});
    document.getElementById('overview').hidden=!!v; document.getElementById('appendix').hidden=!!v;
    cnt.textContent = v ? (shown+' / '+total+' 个靶点') : (total+' 个靶点');
    nores.style.display = (v && shown===0) ? 'block' : 'none';
  }}
  // active nav on scroll
  var links=Array.prototype.slice.call(side.querySelectorAll('a[data-target]'));
  if('IntersectionObserver' in window){{
    var io=new IntersectionObserver(function(es){{
      es.forEach(function(e){{ if(e.isIntersecting){{ links.forEach(function(l){{l.classList.toggle('on', l.dataset.target===e.target.id);}}); }} }});
    }},{{rootMargin:'-20% 0px -70% 0px'}});
    targets.forEach(function(t){{io.observe(t);}});
  }}
  try{{ var s=localStorage.getItem('atlas-q'); if(s){{ q.value=s; run(); }} }}catch(e){{}}
  q.addEventListener('change',function(){{ try{{ localStorage.setItem('atlas-q', q.value); }}catch(e){{}} }});
}})();
</script>
'''
open(OUT, 'w', encoding='utf-8').write(page)
print('written', OUT, len(page), 'chars;', n_targets, 'targets')
