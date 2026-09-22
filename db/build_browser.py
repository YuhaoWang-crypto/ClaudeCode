"""Build docs/FDA_510k_marker_db.html: a single-file browser over db/export/bundle.json.

Usage: python3 db/build_browser.py [bundle.json] [out.html]
"""
import sys, os, json
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
BUNDLE = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'export', 'bundle.json')
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, 'docs', 'FDA_510k_marker_db.html')
b = json.loads(open(BUNDLE, encoding='utf-8').read())
# pre-render curated markdown sections to HTML so the page needs no markdown library at runtime
import markdown, re
sys.path.insert(0, os.path.join(ROOT, 'docs', 'tools'))
try:
    from mdnorm import normalise
except Exception:
    normalise = lambda s: s
for c in b['curated']:
    for k, v in list(c['sections'].items()):
        h = markdown.markdown(normalise(v), extensions=['tables', 'sane_lists'])
        h = h.replace('文件未载明', '<span class="tag na">文件未载明</span>').replace('背景（非申报文件）', '<span class="tag bg">背景（非申报文件）</span>')
        c['sections'][k] = h
data = json.dumps(b, ensure_ascii=False, separators=(',', ':')).replace('<', '\\u003c').replace('\ufffd', '')
n_m, n_s, n_e, n_c = len(b['markers']), len(b['submissions']), len(b['extraction']), len(b['curated'])
n_ds = sum(1 for s in b['submissions'].values() if s.get('ds'))

page = r'''<title>FDA IVD Marker 数据库</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
:root{
  --ground:#F5F7F7; --paper:#FFFFFF; --ink:#1B2229; --ink-2:#4A555C; --ink-3:#7A858B; --line:#D9DFDF; --line-2:#EAEEEE;
  --accent:#146C7A; --accent-soft:#E1F0F2; --accent-ink:#0E4F5A; --ok:#2E7D4F; --ok-soft:#E3F2E8; --warn:#9A6A00; --warn-soft:#FFF3D1;
  --mute:#8A9399; --mute-soft:#EEF1F1; --cur:#5B4B9E; --cur-soft:#ECE8F8; --shadow:0 1px 2px rgba(20,40,45,.06);
  --sans:"IBM Plex Sans","PingFang SC","Hiragino Sans GB","Microsoft YaHei","Noto Sans CJK SC",system-ui,sans-serif;
  --mono:"IBM Plex Mono",ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;
}
@media (prefers-color-scheme: dark){ :root:not([data-theme="light"]){ color-scheme:dark;
  --ground:#0F1416; --paper:#151B1E; --ink:#E4E9E9; --ink-2:#AEB8BC; --ink-3:#7E8A8F; --line:#2A3438; --line-2:#1F2729;
  --accent:#5EC0CD; --accent-soft:#12333A; --accent-ink:#9EDDE5; --ok:#7CCB98; --ok-soft:#173424; --warn:#F0C86A; --warn-soft:#3A2E10;
  --mute:#7E8A8F; --mute-soft:#1E2628; --cur:#C9BDF5; --cur-soft:#2A2450; --shadow:none; } }
:root[data-theme="dark"]{ color-scheme:dark;
  --ground:#0F1416; --paper:#151B1E; --ink:#E4E9E9; --ink-2:#AEB8BC; --ink-3:#7E8A8F; --line:#2A3438; --line-2:#1F2729;
  --accent:#5EC0CD; --accent-soft:#12333A; --accent-ink:#9EDDE5; --ok:#7CCB98; --ok-soft:#173424; --warn:#F0C86A; --warn-soft:#3A2E10;
  --mute:#7E8A8F; --mute-soft:#1E2628; --cur:#C9BDF5; --cur-soft:#2A2450; --shadow:none; }
*{box-sizing:border-box}
html,body{height:100%}
body{margin:0;background:var(--ground);color:var(--ink);font-family:var(--sans);font-size:14px;line-height:1.5;-webkit-font-smoothing:antialiased}
a{color:var(--accent);text-decoration:none} a:hover{text-decoration:underline}
button,select,input{font:inherit;color:var(--ink)}
.top{position:sticky;top:env(safe-area-inset-top,0px);z-index:20;background:var(--paper);border-bottom:1px solid var(--line);display:flex;align-items:center;gap:14px;padding:10px 16px;flex-wrap:wrap}
.brand{font-weight:600} .brand small{display:block;font-weight:400;color:var(--ink-3);font-size:12px}
.stats{display:flex;gap:14px;font-size:12px;color:var(--ink-3);flex-wrap:wrap;font-variant-numeric:tabular-nums}
.stats b{color:var(--ink);font-weight:600}
.top input[type=search]{margin-left:auto;flex:1 1 220px;max-width:440px;padding:7px 10px;border:1px solid var(--line);border-radius:6px;background:var(--ground)}
.top input:focus{outline:2px solid var(--accent);outline-offset:1px}
.wrap{display:grid;grid-template-columns:250px minmax(0,1fr);max-width:1600px;margin:0 auto;min-height:calc(100% - 58px)}
aside{border-right:1px solid var(--line-2);padding:14px 14px 40px;font-size:13px}
aside h4{margin:14px 0 6px;font-size:11.5px;letter-spacing:.06em;text-transform:uppercase;color:var(--ink-3)}
aside label{display:flex;gap:6px;align-items:center;padding:2px 0;color:var(--ink-2);cursor:pointer}
aside label span.n{margin-left:auto;color:var(--ink-3);font-variant-numeric:tabular-nums;font-size:12px}
aside select,aside input[type=text]{width:100%;padding:5px 8px;border:1px solid var(--line);border-radius:5px;background:var(--paper)}
.row2{display:flex;gap:6px}
.btn{padding:5px 10px;border:1px solid var(--line);border-radius:6px;background:var(--paper);color:var(--ink-2);cursor:pointer;font-size:12.5px}
.btn:hover{border-color:var(--accent);color:var(--accent)}
.btn.primary{background:var(--accent);border-color:var(--accent);color:#fff}
main{padding:16px 20px 60px;min-width:0;display:grid;grid-template-columns:minmax(0,1fr);gap:16px}
.panel{background:var(--paper);border:1px solid var(--line);border-radius:8px;box-shadow:var(--shadow)}
.ph{display:flex;align-items:center;gap:10px;padding:10px 14px;border-bottom:1px solid var(--line-2);flex-wrap:wrap}
.ph h2{margin:0;font-size:15px} .ph .sub{color:var(--ink-3);font-size:12.5px}
.tw{overflow-x:auto}
table{border-collapse:collapse;width:100%;font-size:13px;font-variant-numeric:tabular-nums}
th,td{padding:6px 10px;border-bottom:1px solid var(--line-2);text-align:left;vertical-align:top}
th{font-weight:600;color:var(--ink-2);white-space:nowrap;background:var(--ground);position:sticky;top:0;cursor:pointer;user-select:none}
th.on::after{content:" ▾";color:var(--accent)} th.on.asc::after{content:" ▴"}
tr.m{cursor:pointer} tr.m:hover td{background:var(--accent-soft)} tr.m.sel td{background:var(--accent-soft);font-weight:500}
.chip{display:inline-block;font-size:11.5px;padding:1px 7px;border-radius:999px;background:var(--mute-soft);color:var(--ink-2);white-space:nowrap;margin:1px 2px 1px 0}
.chip.pc{background:var(--accent-soft);color:var(--accent-ink);font-family:var(--mono);cursor:pointer}
.chip.cur{background:var(--cur-soft);color:var(--cur)}
.dot{display:inline-block;width:8px;height:8px;border-radius:50%;margin-right:5px;vertical-align:middle}
.dot.ds{background:var(--ok)} .dot.sm{background:var(--warn)} .dot.no{background:var(--mute)}
code,.k{font-family:var(--mono);font-size:.92em}
.k{color:var(--accent-ink);background:var(--accent-soft);padding:1px 5px;border-radius:4px;white-space:nowrap}
.pager{display:flex;gap:8px;align-items:center;padding:8px 14px;font-size:12.5px;color:var(--ink-3)}
.detail{padding:0 14px 14px}
.detail h3{font-size:14px;margin:16px 0 6px}
.kv{display:grid;grid-template-columns:150px minmax(0,1fr);gap:4px 12px;font-size:13px}
.kv div:nth-child(odd){color:var(--ink-3)}
details.sub{border-top:1px solid var(--line-2);padding:4px 0}
details.sub summary{cursor:pointer;font-weight:600;font-size:13.5px;padding:6px 0;list-style:none;display:flex;gap:8px;align-items:center;flex-wrap:wrap}
details.sub summary::-webkit-details-marker{display:none}
details.sub summary::before{content:"";width:7px;height:7px;border-right:1.5px solid var(--ink-3);border-bottom:1.5px solid var(--ink-3);transform:rotate(-45deg);flex:none}
details.sub[open] summary::before{transform:rotate(45deg)}
.fields{display:grid;grid-template-columns:170px minmax(0,1fr);gap:4px 12px;font-size:12.8px;padding:4px 0 8px 15px}
.fields div:nth-child(odd){color:var(--ink-3)} .fields div:nth-child(even){white-space:pre-wrap}
.md{font-size:13.5px;padding:0 0 8px 15px} .md table{font-size:12.5px;display:block;overflow-x:auto} .md th,.md td{border:1px solid var(--line-2)} .md th{position:static}
.md p{margin:6px 0} .md ul{padding-left:20px}
.tag{display:inline-block;font-size:11.5px;line-height:1.3;padding:1px 6px;border-radius:4px;font-weight:500;white-space:nowrap}
.tag.na{background:var(--warn-soft);color:var(--warn)} .tag.bg{background:var(--cur-soft);color:var(--cur)}
.legend{font-size:12px;color:var(--ink-3);display:flex;gap:14px;flex-wrap:wrap;padding:8px 14px}
.empty{padding:30px;text-align:center;color:var(--ink-3)}
.menu-toggle{display:none}
:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
@media (max-width:900px){ .wrap{grid-template-columns:1fr} aside{border-right:none;border-bottom:1px solid var(--line)} aside:not(.open){display:none} .menu-toggle{display:inline-block} main{padding:12px 16px 50px} .kv,.fields{grid-template-columns:1fr} }
</style>

<div class="top">
  <div class="brand">FDA IVD Marker 数据库<small>510(k) · De Novo · PMA 目录 + 决策摘要结构化字段 + 25 靶点人工整理层</small></div>
  <div class="stats" id="stats"></div>
  <button class="btn menu-toggle" id="menu">筛选</button>
  <input id="q" type="search" placeholder="搜 marker、K 号、厂家、产品代码…" aria-label="搜索">
</div>
<div class="wrap">
<aside id="side">
  <h4>FDA 专业分类</h4><div id="f-cat"></div>
  <h4>标签层级</h4><div id="f-tier"></div>
  <h4>申报路径</h4><div id="f-path"></div>
  <h4>证据层级</h4>
  <label><input type="checkbox" id="f-ds"> 至少一份决策摘要已结构化</label>
  <label><input type="checkbox" id="f-cur"> 有人工整理层（25 靶点）</label>
  <h4>清关年份（任一关联申报）</h4>
  <div class="row2"><input type="text" id="f-y1" placeholder="起 如 2015" inputmode="numeric"><input type="text" id="f-y2" placeholder="止" inputmode="numeric"></div>
  <h4>产品代码</h4><input type="text" id="f-pc" placeholder="如 MMI">
  <div style="margin-top:12px" class="row2"><button class="btn" id="reset">清除筛选</button></div>
  <p style="color:var(--ink-3);font-size:12px;margin-top:16px">数据：openFDA 快照 2026-09-14；FDA 决策摘要抓取 2026-09-22；目录层来自 FDA_Marker_Catalog_CN v1（CLIA 文档号关联）。目录中的 510(k) 关联到 2014 年为止，2015 年后由产品代码补齐（in_catalog=0）。</p>
</aside>
<main>
  <div class="panel">
    <div class="ph"><h2>Marker 列表</h2><span class="sub" id="mcount"></span></div>
    <div class="tw"><table id="mt"><thead><tr>
      <th data-k="label">Marker</th><th data-k="category_cn">分类</th><th data-k="tier">层级</th><th data-k="n_sub">关联申报</th><th data-k="n_ds">有决策摘要</th><th data-k="latest">最新清关</th><th data-k="codes">产品代码</th><th data-k="cur">人工层</th>
    </tr></thead><tbody id="mb"></tbody></table></div>
    <div class="pager"><button class="btn" id="prev">上一页</button><span id="pg"></span><button class="btn" id="next">下一页</button></div>
    <div class="legend"><span><span class="dot ds"></span>决策摘要已抓取并解析</span><span><span class="dot sm"></span>仅 510(k) summary / 元数据</span><span><span class="dot no"></span>仅目录记录</span><span><span class="chip cur">人工层</span> 25 靶点图谱中有完整 10 节整理</span></div>
  </div>
  <div class="panel" id="dp" hidden>
    <div class="ph"><h2 id="dtitle"></h2><span class="sub" id="dsub"></span><button class="btn" id="dclose" style="margin-left:auto">关闭</button></div>
    <div class="detail" id="dbody"></div>
  </div>
</main>
</div>
<script id="data" type="application/json">__DATA__</script>
<script>
(function(){
  var D=JSON.parse(document.getElementById('data').textContent);
  var SUB=D.submissions, EXT=D.extraction, PCS=D.product_codes, LINKS=D.links;
  var CUR={}; D.curated.forEach(function(c){ (c.marker_ids||'').split('; ').forEach(function(m){ if(m) (CUR[m]=CUR[m]||[]).push(c); }); });
  var pcIndex={}; Object.keys(SUB).forEach(function(k){ var pc=SUB[k].pc; if(pc) (pcIndex[pc]=pcIndex[pc]||[]).push(k); });
  var esc=function(s){ return String(s==null?'':s).replace(/[&<>"]/g,function(c){return {'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c];}); };
  var md=function(s){ return s; };  // sections are pre-rendered HTML
  // per-marker derived stats; submissions linked = catalog links ∪ (by product code of catalog-linked submissions, post-2014 fill)
  D.markers.forEach(function(m){
    var subs={}; (LINKS[m.marker_id]||[]).forEach(function(l){ subs[l[1].split('/')[0]]=l[0]; });
    var codes={}, ncoded=0; Object.keys(subs).forEach(function(k){ var s=SUB[k]; if(s&&s.pc){ codes[s.pc]=(codes[s.pc]||0)+1; ncoded++; } });
    // post-2014 fill by product code, only for codes that carry the marker (>=3 catalog links or >=25% of its coded links),
    // plus the product codes named in the curated atlas for this marker (covers newer codes such as PMT/PRI for procalcitonin)
    var fill={}; Object.keys(codes).forEach(function(pc){ if(codes[pc]>=3 || codes[pc]>=0.25*ncoded) fill[pc]=1; });
    (CUR[m.marker_id]||[]).forEach(function(c){ (c.product_codes||'').split('; ').forEach(function(pc){ if(pc){ fill[pc]=1; if(!codes[pc]) codes[pc]=0; } }); });
    Object.keys(fill).forEach(function(pc){ (pcIndex[pc]||[]).forEach(function(k){ if(!subs[k]) subs[k]='510(k)*'; }); });
    m.subs=subs; m.codes=Object.keys(codes).sort(function(a,b){return codes[b]-codes[a];}); m.codeCounts=codes; m.n_sub=Object.keys(subs).length;
    var ds=0, latest=0, paths={};
    Object.keys(subs).forEach(function(k){ var s=SUB[k]; paths[subs[k].replace('*','')]=1; if(s){ if(s.ds&&EXT[k]) ds++; if(s.y&&s.y>latest) latest=s.y; } else { paths[subs[k].replace('*','')]=1; } });
    m.n_ds=ds; m.latest=latest||''; m.paths=Object.keys(paths); m.cur=CUR[m.marker_id]?1:0;
    m.hay=(m.label+' '+m.category_cn+' '+m.codes.join(' ')+' '+Object.keys(subs).join(' ')).toLowerCase();
  });
  document.getElementById('stats').innerHTML='<span><b>'+D.markers.length+'</b> markers</span><span><b>'+Object.keys(SUB).length+'</b> 申报号</span><span><b>'+Object.keys(EXT).length+'</b> 份决策摘要已结构化</span><span><b>'+D.curated.length+'</b> 靶点人工层</span><span>生成 '+D.generated+'</span>';
  // filters
  function counts(key){ var c={}; D.markers.forEach(function(m){ var v=m[key]||'—'; c[v]=(c[v]||0)+1; }); return c; }
  function boxes(id,key,vals){ var el=document.getElementById(id); el.innerHTML=vals.map(function(v){ return '<label><input type="checkbox" data-f="'+key+'" value="'+esc(v[0])+'"> '+esc(v[0])+'<span class="n">'+v[1]+'</span></label>'; }).join(''); }
  var cc=counts('category_cn'); boxes('f-cat','category_cn',Object.keys(cc).sort(function(a,b){return cc[b]-cc[a];}).map(function(k){return [k,cc[k]];}));
  var tc=counts('tier'); boxes('f-tier','tier',Object.keys(tc).sort(function(a,b){return tc[b]-tc[a];}).map(function(k){return [k,tc[k]];}));
  var pc2={}; D.markers.forEach(function(m){ m.paths.forEach(function(p){ pc2[p]=(pc2[p]||0)+1; }); }); boxes('f-path','paths',Object.keys(pc2).map(function(k){return [k,pc2[k]];}));
  var state={q:'',cat:{},tier:{},paths:{},ds:false,cur:false,y1:null,y2:null,pc:'',sort:'label',asc:true,page:0,sel:null};
  var PAGE=60;
  function apply(){
    var q=state.q.toLowerCase(), pc=state.pc.toUpperCase();
    var rows=D.markers.filter(function(m){
      if(q && m.hay.indexOf(q)<0) return false;
      if(Object.keys(state.cat).length && !state.cat[m.category_cn]) return false;
      if(Object.keys(state.tier).length && !state.tier[m.tier]) return false;
      if(Object.keys(state.paths).length && !m.paths.some(function(p){return state.paths[p];})) return false;
      if(state.ds && !m.n_ds) return false;
      if(state.cur && !m.cur) return false;
      if(pc && m.codes.indexOf(pc)<0) return false;
      if(state.y1||state.y2){ var ok=Object.keys(m.subs).some(function(k){ var s=SUB[k]; return s&&s.y&&(!state.y1||s.y>=state.y1)&&(!state.y2||s.y<=state.y2); }); if(!ok) return false; }
      return true;
    });
    var k=state.sort, asc=state.asc?1:-1;
    rows.sort(function(a,b){ var x=a[k],y=b[k]; if(k==='codes'){x=a.codes.join(','),y=b.codes.join(',');} if(typeof x==='number'||typeof y==='number'){ return ((x||0)-(y||0))*asc; } return String(x||'').localeCompare(String(y||''),'zh')*asc; });
    var n=rows.length, pages=Math.max(1,Math.ceil(n/PAGE)); if(state.page>=pages) state.page=pages-1;
    var slice=rows.slice(state.page*PAGE,(state.page+1)*PAGE);
    document.getElementById('mcount').textContent=n+' / '+D.markers.length;
    document.getElementById('pg').textContent='第 '+(state.page+1)+' / '+pages+' 页';
    document.getElementById('mb').innerHTML=slice.map(function(m){
      var dot=m.n_ds?'ds':(m.n_sub?'sm':'no');
      return '<tr class="m'+(state.sel===m.marker_id?' sel':'')+'" data-id="'+m.marker_id+'"><td><span class="dot '+dot+'"></span>'+esc(m.label)+'</td><td>'+esc(m.category_cn)+'</td><td>'+esc((m.tier||'').split('／')[0])+'</td><td>'+m.n_sub+'</td><td>'+m.n_ds+'</td><td>'+(m.latest||'')+'</td><td>'+m.codes.slice(0,6).map(function(c){return '<span class="chip pc" data-pc="'+c+'" title="'+esc((PCS[c]||{}).device_name||'')+'">'+c+'</span>';}).join('')+(m.codes.length>6?'<span class="chip">+'+(m.codes.length-6)+'</span>':'')+'</td><td>'+(m.cur?'<span class="chip cur">人工层</span>':'')+'</td></tr>';
    }).join('') || '<tr><td colspan="8" class="empty">没有匹配的 marker</td></tr>';
    document.querySelectorAll('#mt th').forEach(function(th){ th.classList.toggle('on',th.dataset.k===state.sort); th.classList.toggle('asc',th.dataset.k===state.sort&&state.asc); });
  }
  document.getElementById('side').addEventListener('change',function(e){
    var t=e.target; if(t.dataset.f){ var set=state[t.dataset.f==='category_cn'?'cat':(t.dataset.f==='tier'?'tier':'paths')]; if(t.checked) set[t.value]=1; else delete set[t.value]; }
    state.ds=document.getElementById('f-ds').checked; state.cur=document.getElementById('f-cur').checked;
    state.y1=parseInt(document.getElementById('f-y1').value)||null; state.y2=parseInt(document.getElementById('f-y2').value)||null; state.pc=document.getElementById('f-pc').value.trim();
    state.page=0; apply();
  });
  ['f-y1','f-y2','f-pc'].forEach(function(id){ document.getElementById(id).addEventListener('input',function(){ state.y1=parseInt(document.getElementById('f-y1').value)||null; state.y2=parseInt(document.getElementById('f-y2').value)||null; state.pc=document.getElementById('f-pc').value.trim(); state.page=0; apply(); }); });
  var qt=null; document.getElementById('q').addEventListener('input',function(e){ clearTimeout(qt); qt=setTimeout(function(){ state.q=e.target.value.trim(); state.page=0; apply(); },150); });
  document.getElementById('reset').addEventListener('click',function(){ state.cat={};state.tier={};state.paths={};state.ds=false;state.cur=false;state.y1=state.y2=null;state.pc='';state.q='';state.page=0; document.querySelectorAll('#side input').forEach(function(i){ if(i.type==='checkbox') i.checked=false; else i.value=''; }); document.getElementById('q').value=''; apply(); });
  document.getElementById('prev').addEventListener('click',function(){ if(state.page>0){state.page--;apply();} });
  document.getElementById('next').addEventListener('click',function(){ state.page++; apply(); });
  document.querySelector('#mt thead').addEventListener('click',function(e){ var th=e.target.closest('th'); if(!th) return; if(state.sort===th.dataset.k) state.asc=!state.asc; else { state.sort=th.dataset.k; state.asc=true; } apply(); });
  document.getElementById('menu').addEventListener('click',function(){ document.getElementById('side').classList.toggle('open'); });
  document.getElementById('mb').addEventListener('click',function(e){
    var pc=e.target.closest('.chip.pc'); if(pc){ document.getElementById('f-pc').value=pc.dataset.pc; state.pc=pc.dataset.pc; state.page=0; apply(); e.stopPropagation(); return; }
    var tr=e.target.closest('tr.m'); if(tr){ show(tr.dataset.id); }
  });
  document.getElementById('dclose').addEventListener('click',function(){ document.getElementById('dp').hidden=true; state.sel=null; apply(); });
  // detail
  var FIELD_LABEL={measurand:'Measurand',type_of_test:'检测类型',intended_use:'Intended use',indications:'Indications for use',specimen_types:'样本类型（推导）',assay_cutoff:'Assay cut-off',clinical_cutoff:'Clinical cut-off',cutoff_numbers:'cutoff 数值（推导）',reference_range:'Expected values / reference range',standards:'Standards / guidance referenced',clsi_codes:'CLSI/ISO 编号（推导）',precision:'Precision',detection_limit:'Detection limit',traceability:'Traceability / stability',method_comparison:'Method comparison',mc_slope:'斜率（推导）',mc_r:'r（推导）',clinical_studies:'Clinical studies',clinical_sensitivity:'Clinical sensitivity',clinical_specificity:'Clinical specificity',sens_pct:'灵敏度 %（推导）',spec_pct:'特异度 %（推导）',sample_n:'n（推导）',predicate:'Predicate',instrument:'Instrument',conclusion:'Conclusion'};
  var ORDER=['measurand','type_of_test','intended_use','specimen_types','assay_cutoff','clinical_cutoff','cutoff_numbers','reference_range','clsi_codes','detection_limit','mc_slope','mc_r','clinical_studies','sens_pct','spec_pct','sample_n','predicate'];
  function subRow(k,pathTag){
    var s=SUB[k]||{}, e=EXT[k]; var dot=e?'ds':(s.sm||s.ds?'sm':'no');
    var url=k.indexOf('DEN')===0?'https://www.accessdata.fda.gov/scripts/cdrh/cfdocs/cfpmn/denovo.cfm?id='+k:(k.indexOf('P')===0?'https://www.accessdata.fda.gov/scripts/cdrh/cfdocs/cfpma/pma.cfm?id='+k:'https://www.accessdata.fda.gov/scripts/cdrh/cfdocs/cfpmn/pmn.cfm?ID='+k);
    var dsurl='https://www.accessdata.fda.gov/cdrh_docs/reviews/'+k+'.pdf';
    var head='<summary><span class="dot '+dot+'"></span><span class="k">'+k+'</span> <span>'+esc(s.d||s.y||'')+'</span> <span>'+esc(s.a||'')+'</span> <span style="color:var(--ink-2)">'+esc(s.n||'')+'</span>'+(s.pc?'<span class="chip pc">'+s.pc+'</span>':'')+(pathTag==='510(k)*'?'<span class="chip" title="目录未收录，按产品代码从 openFDA 补入">补入</span>':'')+(s.cat===0?'':'')+'</summary>';
    var body='<div class="fields"><div>FDA 记录</div><div><a href="'+url+'" target="_blank" rel="noopener">'+url.replace('https://www.accessdata.fda.gov','…')+'</a>'+(e||s.ds?' · <a href="'+dsurl+'" target="_blank" rel="noopener">决策摘要 PDF</a>':'')+'</div>';
    if(s.reg) body+='<div>21 CFR</div><div>'+esc(s.reg)+'</div>';
    if(e){ body+='<div>模板</div><div>'+esc(e.tpl)+'（页面只显示截断的关键字段；完整段落在 SQLite/CSV 的 extraction 表）</div>'; ORDER.forEach(function(f){ if(e[f]) body+='<div>'+FIELD_LABEL[f]+'</div><div>'+esc(e[f])+'</div>'; }); }
    else body+='<div>结构化字段</div><div style="color:var(--ink-3)">'+(s.ds?'决策摘要已抓取，未识别模板':'无决策摘要文本（'+(s.sm?'仅 510(k) summary':'目录/元数据')+'）')+'</div>';
    body+='</div>';
    return '<details class="sub">'+head+body+'</details>';
  }
  function show(id){
    var m=D.markers.filter(function(x){return x.marker_id===id;})[0]; if(!m) return; state.sel=id; apply();
    document.getElementById('dtitle').textContent=m.label; document.getElementById('dsub').textContent=m.marker_id+' · '+m.category_cn+' · '+m.tier;
    var h='<div class="kv"><div>关联申报</div><div>'+m.n_sub+'（目录 '+Object.keys(m.subs).filter(function(k){return m.subs[k]!=='510(k)*';}).length+'，按产品代码补入 '+Object.keys(m.subs).filter(function(k){return m.subs[k]==='510(k)*';}).length+'）</div><div>产品代码</div><div>'+(m.codes.map(function(c){ var p=PCS[c]||{}; return '<span class="chip pc">'+c+'</span> '+esc(p.device_name||'')+(p.regulation_number?' · 21 CFR '+p.regulation_number:'')+(p.device_class?' · Class '+p.device_class:'')+'<br>'; }).join('')||'—')+'</div>'+(m.clia_analyte_id?'<div>CLIA analyte ID</div><div>'+esc(m.clia_analyte_id)+'</div>':'')+'</div>';
    (CUR[id]||[]).forEach(function(c){
      h+='<h3><span class="chip cur">人工整理层</span> '+esc(c.title)+'</h3><div class="kv"><div>cutoff 原型</div><div>'+esc(c.archetype||'—')+'</div><div>代表性 cutoff</div><div>'+esc(c.representative_cutoff||'—')+'</div><div>产品代码</div><div>'+esc(c.product_codes)+'</div></div>';
      [4,8,6,7,3,9].forEach(function(n){ var sec=c.sections[n]||c.sections[String(n)]; if(!sec) return; var title={3:'3. 预期用途与声明类型',4:'4. 阳性/阴性判定与 cutoff 逻辑',6:'6. 样本类型',7:'7. 分析性能标准',8:'8. 临床验证设计与结果',9:'9. 厂家间差异'}[n]; h+='<details class="sub"'+(n===4?' open':'')+'><summary>'+title+'</summary><div class="md">'+md(sec)+'</div></details>'; });
    });
    var ks=Object.keys(m.subs).sort(function(a,b){ var x=(SUB[a]||{}).d||'', y=(SUB[b]||{}).d||''; return y.localeCompare(x)||b.localeCompare(a); });
    var withE=ks.filter(function(k){return EXT[k];}), rest=ks.filter(function(k){return !EXT[k];});
    h+='<h3>申报记录（'+ks.length+'；'+withE.length+' 份有结构化决策摘要）</h3>'+withE.map(function(k){return subRow(k,m.subs[k]);}).join('');
    if(rest.length) h+='<details class="sub"><summary>其余 '+rest.length+' 条（无决策摘要文本）</summary>'+rest.slice(0,400).map(function(k){return subRow(k,m.subs[k]);}).join('')+(rest.length>400?'<div class="empty">仅显示前 400 条</div>':'')+'</details>';
    document.getElementById('dbody').innerHTML=h; document.getElementById('dp').hidden=false;
    document.getElementById('dp').scrollIntoView({behavior:'smooth',block:'start'});
    try{ history.replaceState(null,'','#'+id); }catch(e){}
  }
  apply();
  var hash=(location.hash||'').slice(1); if(/^M\d{4}$/.test(hash)) show(hash);
})();
</script>
'''
page = page.replace('__DATA__', data)
open(OUT, 'w', encoding='utf-8').write(page)
print('written', OUT, round(len(page.encode('utf-8')) / 1e6, 1), 'MB', f'| markers {n_m} submissions {n_s} (decision summaries {n_ds}) extracted {n_e} curated {n_c}')
