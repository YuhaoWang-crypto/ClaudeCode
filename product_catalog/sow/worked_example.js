
/* ═══════════ 参考数据（内置样本，演示用） ═══════════ */
const REF_DRUGS=[
 ["Imatinib","Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1","BCR-ABL / 酪氨酸激酶抑制"],
 ["Gefitinib","COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1","EGFR 抑制"],
 ["Erlotinib","C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1","EGFR 抑制"],
 ["Lapatinib","CS(=O)(=O)CCNCc1ccc(-c2ccc3ncnc(Nc4ccc(OCc5cccc(F)c5)c(Cl)c4)c3c2)o1","EGFR/HER2 抑制"],
 ["Sorafenib","CNC(=O)c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(C(F)(F)F)c3)cc2)ccn1","多激酶抑制"],
 ["Sunitinib","CCN(CC)CCNC(=O)c1c(C)[nH]c(C=C2C(=O)Nc3ccc(F)cc32)c1C","多激酶抑制"],
 ["Dasatinib","Cc1nc(Nc2ncc(C(=O)Nc3c(C)cccc3Cl)s2)cc(N2CCN(CCO)CC2)n1","SRC/ABL 抑制"],
 ["Colchicine","COc1cc2c(c(OC)c1OC)-c1ccc(OC)c(=O)cc1C(NC(C)=O)CC2","微管蛋白结合"],
 ["Nocodazole","COC(=O)Nc1nc2ccc(C(=O)c3cccs3)cc2[nH]1","微管蛋白结合"],
 ["Vinblastine-like","CCC1(O)CC2CN(CCc3c([nH]c4ccccc34)C(C(=O)OC)(c3cc4c(cc3OC)N(C)C3C(O)(C(=O)OC)C(OC(C)=O)C5(CC)C=CCN6CCC43C65)C2)C1","微管蛋白结合"],
 ["Vorinostat","O=C(NO)CCCCCCC(=O)Nc1ccccc1","HDAC 抑制"],
 ["Trichostatin A","CC(C=CC=CC(=O)NO)C(C)C(=O)c1ccc(N(C)C)cc1","HDAC 抑制"],
 ["Entinostat","NCc1cccc(NC(=O)c2ccc(CNC(=O)OCc3cccnc3)cc2)c1N","HDAC 抑制"],
 ["Bortezomib","CC(C)CC(NC(=O)C(Cc1ccccc1)NC(=O)c1cnccn1)B(O)O","蛋白酶体抑制"],
 ["MG-132","CC(C)CC(NC(=O)C(CC(C)C)NC(=O)C(CC(C)C)NC(=O)OCc1ccccc1)C=O","蛋白酶体抑制"],
 ["Methotrexate","CN(Cc1cnc2nc(N)nc(N)c2n1)c1ccc(C(=O)NC(CCC(=O)O)C(=O)O)cc1","抗叶酸 / DHFR"],
 ["5-Fluorouracil","O=c1[nH]cc(F)c(=O)[nH]1","抗代谢 / TS 抑制"],
 ["Gemcitabine","NC1=NC(=O)N(C2OC(CO)C(O)C2(F)F)C=C1","抗代谢 / 核苷类"],
 ["Simvastatin","CCC(C)(C)C(=O)OC1CC(C)C=C2C=CC(C)C(CCC3CC(O)CC(=O)O3)C12","HMG-CoA 还原酶抑制"],
 ["Lovastatin","CCC(C)C(=O)OC1CC(C)C=C2C=CC(C)C(CCC3CC(O)CC(=O)O3)C12","HMG-CoA 还原酶抑制"],
 ["Dexamethasone","CC1CC2C3CCC4=CC(=O)C=CC4(C)C3(F)C(O)CC2(C)C1(O)C(=O)CO","糖皮质激素受体激动"],
 ["Tamoxifen","CCC(=C(c1ccccc1)c1ccc(OCCN(C)C)cc1)c1ccccc1","雌激素受体调节"],
 ["Raloxifene","OCCN1CCCCC1CCOc1ccc(C(=O)c2c(-c3ccc(O)cc3)sc3cc(O)ccc23)cc1","雌激素受体调节"],
 ["Thalidomide","O=C1CCC(N2C(=O)c3ccccc3C2=O)C(=O)N1","CRBN 分子胶"],
 ["Lenalidomide","NC1=CC=CC2=C1C(=O)N(C1CCC(=O)NC1=O)C2","CRBN 分子胶"],
 ["Chloroquine","CCN(CC)CCCC(C)Nc1ccnc2cc(Cl)ccc12","溶酶体 / 自噬阻断"],
 ["Metformin","CN(C)C(=N)N=C(N)N","线粒体复合物 I / 代谢"],
 ["Rotenone","COc1cc2c(cc1OC)C1COc3cc(OC)c(OC)cc3C1OC2=O","线粒体复合物 I"],
 ["Curcumin","COc1cc(C=CC(=O)CC(=O)C=Cc2ccc(O)c(OC)c2)ccc1O","频繁命中 / PAINS 风险"],
 ["Quercetin","O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12","频繁命中 / PAINS 风险"]
];
const AW={C:12.011,N:14.007,O:15.999,S:32.06,P:30.974,F:18.998,Cl:35.45,Br:79.904,I:126.90,B:10.81,Si:28.09,Se:78.97,H:1.008,Pt:195.08,Na:22.99,K:39.098,Zn:65.38,Fe:55.845,Mg:24.305,Ca:40.078,Mn:54.938,Cu:63.546,Li:6.94};
const VAL={C:4,N:3,O:2,S:2,P:3,F:1,Cl:1,Br:1,I:1,B:3,Si:4,Se:2};
const LOGP={C:0.40,N:-0.85,O:-0.55,S:0.30,P:-0.20,F:0.30,Cl:0.75,Br:0.95,I:1.10,B:0.0,Si:0.5,Se:0.4};
const ALERTS=[
 [/c1ccc\(O\)c\(O\)c1|c\(O\)c\(O\)c/,"儿茶酚（PAINS 频繁命中）"],
 [/C\(=O\)C=CC\(=O\)|C=CC\(=O\)C/,"Michael 受体 / 迈克尔加成"],
 [/\[N\+\]\(=O\)\[O-\]|N\(=O\)=O/,"硝基芳香"],
 [/N=N/,"偶氮基团"],
 [/C\(=S\)|S\(=O\)\(=O\)N.*N.*S\(=O\)/,"硫代羰基 / 多磺酰胺"],
 [/O=C1CSC\(=N1\)|C1SC\(=O\)N.*1/,"罗丹宁类支架"],
 [/c1ccc2c\(c1\)C\(=O\)c1ccccc1C2=O/,"醌类"],
 [/B\(O\)O/,"硼酸（共价弹头，非警报但需标注）"]
];
const GENES=(()=>{ // symbol, pathway, essentiality tier, in-corpus, base |effect|
 const raw=[
 ["KRAS","MAPK","selective",1,.82],["NRAS","MAPK","selective",1,.61],["HRAS","MAPK","non",1,.24],
 ["BRAF","MAPK","selective",1,.66],["MAP2K1","MAPK","selective",1,.58],["MAPK1","MAPK","selective",1,.63],
 ["MAPK3","MAPK","non",1,.31],["EGFR","MAPK",  "selective",1,.55],["SHP2","MAPK","selective",0,.44],
 ["PTPN11","MAPK","selective",1,.57],["SOS1","MAPK","selective",1,.49],["RAF1","MAPK","selective",1,.52],
 ["PIK3CA","PI3K-AKT","selective",1,.60],["AKT1","PI3K-AKT","selective",1,.47],["MTOR","PI3K-AKT","essential",1,.78],
 ["PTEN","PI3K-AKT","non",1,.28],["RPTOR","PI3K-AKT","essential",1,.71],["RICTOR","PI3K-AKT","selective",1,.40],
 ["TP53","细胞周期/凋亡","non",1,.35],["MDM2","细胞周期/凋亡","selective",1,.54],["RB1","细胞周期/凋亡","non",1,.30],
 ["CDK1","细胞周期/凋亡","essential",1,.93],["CDK2","细胞周期/凋亡","selective",1,.42],["CDK4","细胞周期/凋亡","selective",1,.51],
 ["CDK6","细胞周期/凋亡","selective",1,.45],["CCND1","细胞周期/凋亡","selective",1,.48],["CCNE1","细胞周期/凋亡","selective",1,.44],
 ["PLK1","细胞周期/凋亡","essential",1,.95],["AURKA","细胞周期/凋亡","essential",1,.80],["AURKB","细胞周期/凋亡","essential",1,.86],
 ["WEE1","细胞周期/凋亡","selective",1,.53],["KIF11","细胞周期/凋亡","essential",1,.91],["TOP2A","细胞周期/凋亡","essential",1,.74],
 ["BCL2","细胞周期/凋亡","selective",1,.39],["MCL1","细胞周期/凋亡","essential",1,.76],["BAX","细胞周期/凋亡","non",1,.22],
 ["CASP3","细胞周期/凋亡","non",1,.19],["CASP9","细胞周期/凋亡","non",0,.18],
 ["BRCA1","DNA 损伤修复","selective",1,.58],["BRCA2","DNA 损伤修复","selective",1,.55],["PARP1","DNA 损伤修复","selective",1,.47],
 ["ATM","DNA 损伤修复","selective",1,.43],["ATR","DNA 损伤修复","essential",1,.79],["CHEK1","DNA 损伤修复","essential",1,.83],
 ["RAD51","DNA 损伤修复","essential",1,.81],["WRN","DNA 损伤修复","selective",1,.50],["MLH1","DNA 损伤修复","non",1,.26],
 ["MSH2","DNA 损伤修复","non",1,.27],["POLQ","DNA 损伤修复","selective",0,.41],["FANCD2","DNA 损伤修复","selective",1,.46],
 ["RPL5","核心essential","essential",1,.97],["RPS6","核心essential","essential",1,.98],["RPL11","核心essential","essential",1,.96],
 ["POLR2A","核心essential","essential",1,.99],["EEF2","核心essential","essential",1,.97],["PSMA1","核心essential","essential",1,.95],
 ["PSMD1","核心essential","essential",1,.94],["SNRPD1","核心essential","essential",1,.93],["RAN","核心essential","essential",1,.92],
 ["ESR1","核受体/表观","selective",1,.52],["AR","核受体/表观","selective",1,.38],["NR3C1","核受体/表观","selective",1,.36],
 ["HDAC1","核受体/表观","selective",1,.49],["HDAC2","核受体/表观","non",1,.29],["EZH2","核受体/表观","selective",1,.45],
 ["KAT2A","核受体/表观","selective",0,.34],["SMARCA4","核受体/表观","selective",1,.56],["ARID1A","核受体/表观","non",1,.25],
 ["MYC","转录/代谢","essential",1,.85],["STK11","转录/代谢","non",1,.23],["KEAP1","转录/代谢","non",1,.21],
 ["NFE2L2","转录/代谢","selective",1,.37],["HMGCR","转录/代谢","selective",1,.44],["SLC7A11","转录/代谢","selective",0,.33],
 ["OR2L13","阴性对照","non",1,.03],["OR5A1","阴性对照","non",0,.02],["GFP-ctrl","阴性对照","non",1,.01]];
 const m={}; raw.forEach(r=>m[r[0].toUpperCase()]={sym:r[0],path:r[1],ess:r[2],seen:!!r[3],eff:r[4]}); return m;
})();
const PATH_SET=[...new Set(Object.values(GENES).map(g=>g.path))];

/* ═══════════ 真实计算：SMILES 解析 ═══════════ */
function parseSMILES(s){
  s=s.trim(); if(!s) return null;
  const atoms=[]; let i=0, ringOpen={}, ringPairs=0, branchStack=[], prev=-1, pendingBond=1, ok=true;
  const two=["Cl","Br","Si","Se","Na"];
  while(i<s.length){
    const c=s[i];
    if(c==='('){branchStack.push(prev); i++; continue;}
    if(c===')'){prev=branchStack.pop(); if(prev===undefined){ok=false;prev=-1;} i++; continue;}
    if(c==='-'){pendingBond=1;i++;continue;} if(c==='='){pendingBond=2;i++;continue;}
    if(c==='#'){pendingBond=3;i++;continue;} if(c===':'){pendingBond=1.5;i++;continue;}
    if(c==='/'||c==='\\'||c==='.'||c==='@'||c==='+'||c==='%'){ if(c==='.')prev=-1; i++; continue;}
    if(/[0-9]/.test(c)){
      const lbl=c;
      if(ringOpen[lbl]!==undefined){ const a=ringOpen[lbl]; delete ringOpen[lbl];
        if(a>=0&&prev>=0){atoms[a].bonds+=pendingBond;atoms[prev].bonds+=pendingBond;atoms[a].ring=1;atoms[prev].ring=1;ringPairs++;}
      } else ringOpen[lbl]=prev;
      pendingBond=1; i++; continue;
    }
    let el=null, arom=false, hExp=0, chg=0;
    if(c==='['){
      const j=s.indexOf(']',i); if(j<0){ok=false;break;}
      const inner=s.slice(i+1,j);
      const em=inner.match(/^[0-9]*([A-Za-z][a-z]?)/); if(!em){i=j+1;continue;}
      el=em[1]; arom=/^[a-z]/.test(el); el=el[0].toUpperCase()+el.slice(1);
      const hm=inner.match(/H([0-9]*)/); if(hm) hExp=hm[1]===''?1:+hm[1];
      const cm=inner.match(/([+-])([0-9]*)/); if(cm) chg=(cm[1]==='+'?1:-1)*(cm[2]===''?1:+cm[2]);
      i=j+1;
    } else {
      const t2=s.slice(i,i+2);
      if(two.includes(t2)){el=t2;i+=2;}
      else if(/[A-Za-z]/.test(c)){ arom=/[a-z]/.test(c); el=c.toUpperCase(); i++; }
      else {i++;continue;}
    }
    const idx=atoms.length;
    atoms.push({el,arom,hExp,chg,bonds:0,ring:0});
    if(prev>=0){atoms[prev].bonds+=pendingBond;atoms[idx].bonds+=pendingBond;}
    pendingBond=1; prev=idx;
  }
  if(!atoms.length) return null;
  let mw=0,hc=0,hbd=0,hba=0,tpsa=0,logp=0,arC=0,rot=0;
  atoms.forEach(a=>{
    const w=AW[a.el]; if(w===undefined) {ok=false; return;}
    const v=VAL[a.el]||0;
    let h=a.hExp;
    if(!h && v){ const used=a.arom? a.bonds+1 : a.bonds; h=Math.max(0,Math.round(v-used+a.chg*(a.el==='N'?1:0))); }
    mw+=w+h*AW.H; hc+=h;
    logp+=(a.arom&&a.el==='N'?-0.45:(LOGP[a.el]!==undefined?LOGP[a.el]:0)+(a.arom?-0.10:0))-h*0.02;
    if(a.el==='N'||a.el==='O'){ hba++; if(h>0) hbd+=h;
      if(a.el==='O') tpsa+= h>0?20.2:(a.arom?13.1:9.2);
      else tpsa+= h>=2?26.0:(h===1?12.0:(a.arom?12.9:3.2));
    }
    if(a.el==='S') tpsa+=25.3; if(a.el==='P') tpsa+=13.6;
    if(a.arom&&a.el==='C') arC++;
    if(!a.ring&&a.bonds>=2&&a.bonds<=4) rot++;
  });
  const rings=ringPairs;
  return {ok, n:atoms.length, mw:+mw.toFixed(1), hbd, hba, tpsa:+tpsa.toFixed(1),
    logp:+(logp).toFixed(2), rings, arRings:+(arC/6).toFixed(1), rot:Math.max(0,rot-rings*2),
    lip:(mw<=500?1:0)+(logp<=5?1:0)+(hbd<=5?1:0)+(hba<=10?1:0)};
}
function fp(s){ // 字符 n-gram 哈希指纹（演示用）
  const bits=new Uint8Array(512); const t=s.replace(/\s/g,'');
  for(let n=2;n<=4;n++) for(let i=0;i+n<=t.length;i++){
    let h=2166136261; const g=t.slice(i,i+n);
    for(let k=0;k<g.length;k++){h^=g.charCodeAt(k);h=Math.imul(h,16777619);}
    bits[(h>>>0)%512]=1;
  }
  return bits;
}
function tanimoto(a,b){let i=0,u=0;for(let k=0;k<512;k++){const x=a[k],y=b[k];if(x&&y)i++;if(x||y)u++;}return u?i/u:0;}
function alertsOf(s){const out=[];ALERTS.forEach(([re,name])=>{if(re.test(s))out.push(name);});return out;}

/* 种子随机（可复现） */
function rng(seed){let s=seed>>>0||1;return()=>{s^=s<<13;s>>>=0;s^=s>>17;s^=s<<5;s>>>=0;return s/4294967296;};}
const hashStr=s=>{let h=2166136261;for(let i=0;i<s.length;i++){h^=s.charCodeAt(i);h=Math.imul(h,16777619);}return h>>>0;};

/* ═══════════ 板图生成（真实布局算法） ═══════════ */
function makePlate(items,fmt,seed){
  const R=fmt===384?16:8, C=fmt===384?24:12;
  const rows="ABCDEFGHIJKLMNOP".slice(0,R).split("");
  const rand=rng(seed);
  const wells=[];
  for(let r=0;r<R;r++)for(let c=0;c<C;c++){
    const edge=(r===0||r===R-1||c===0||c===C-1);
    wells.push({r,c,id:rows[r]+String(c+1).padStart(2,"0"),edge,role:edge?"empty":null});
  }
  const inner=wells.filter(w=>!w.edge);
  // 对照列：内区第一列 = DMSO，最后一列 = 阳性/阴性交替
  const cols=[...new Set(inner.map(w=>w.c))].sort((a,b)=>a-b);
  const dmsoCol=cols[0], ctrlCol=cols[cols.length-1];
  inner.forEach(w=>{
    if(w.c===dmsoCol) w.role="dmso";
    else if(w.c===ctrlCol) w.role=(w.r%2===0)?"pos":"neg";
  });
  const slots=inner.filter(w=>w.role===null);
  // 种子化随机布位，避免同类相邻聚集
  const order=slots.map((w,i)=>({w,k:rand()})).sort((a,b)=>a.k-b.k).map(o=>o.w);
  const placed=items.slice(0,order.length);
  placed.forEach((it,i)=>{order[i].role="sample";order[i].item=it;order[i].rank=i+1;});
  return {R,C,rows,wells,capacity:order.length,placed:placed.length,
    nDmso:inner.filter(w=>w.role==="dmso").length,
    nPos:inner.filter(w=>w.role==="pos").length,
    nNeg:inner.filter(w=>w.role==="neg").length,
    nEdge:wells.filter(w=>w.edge).length};
}

/* ═══════════ 图表 ═══════════ */
function enrichmentChart(scores,baseScores,hitFlag){
  const n=scores.length; if(n<4) return "";
  const W=760,H=250,mL=52,mR=18,mT=14,mB=36;
  const iw=W-mL-mR, ih=H-mT-mB;
  const curve=(sc)=>{
    const idx=sc.map((s,i)=>[s,i]).sort((a,b)=>b[0]-a[0]).map(p=>p[1]);
    const tot=hitFlag.reduce((a,b)=>a+b,0)||1; let c=0;
    return idx.map((i,k)=>{c+=hitFlag[i];return [(k+1)/n,c/tot];});
  };
  const m=curve(scores), b=curve(baseScores);
  const X=v=>mL+v*iw, Y=v=>mT+(1-v)*ih;
  const path=(pts)=>"M"+X(0)+","+Y(0)+"L"+pts.map(p=>X(p[0])+","+Y(p[1])).join("L");
  const ticks=[0,.25,.5,.75,1];
  const at=(pts,f)=>{let best=pts[0];for(const p of pts){if(p[0]<=f)best=p;}return best[1];}
  const m20=at(m,.2), b20=at(b,.2);
  return `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="富集曲线：模型 vs 结构指纹基线">
    ${ticks.map(t=>`<line x1="${mL}" x2="${W-mR}" y1="${Y(t)}" y2="${Y(t)}" stroke="#0d172014"/>
      <text x="${mL-8}" y="${Y(t)+4}" text-anchor="end" font-size="10" font-family="ui-monospace,monospace" fill="#667078">${Math.round(t*100)}%</text>`).join("")}
    ${ticks.map(t=>`<text x="${X(t)}" y="${H-12}" text-anchor="middle" font-size="10" font-family="ui-monospace,monospace" fill="#667078">${Math.round(t*100)}%</text>`).join("")}
    <line x1="${mL}" x2="${W-mR}" y1="${Y(0)}" y2="${Y(0)}" stroke="#0d172040"/>
    <line x1="${X(0)}" y1="${Y(0)}" x2="${X(1)}" y2="${Y(1)}" stroke="#0d172033" stroke-dasharray="4 4"/>
    <text x="${X(.62)}" y="${Y(.55)}" font-size="10.5" fill="#667078" transform="rotate(-24 ${X(.62)} ${Y(.55)})">随机排序</text>
    <path d="${path(b)}" fill="none" stroke="var(--c2)" stroke-width="2" stroke-linejoin="round"/>
    <path d="${path(m)}" fill="none" stroke="var(--c1)" stroke-width="2" stroke-linejoin="round"/>
    <line x1="${X(.2)}" y1="${mT}" x2="${X(.2)}" y2="${Y(0)}" stroke="#0d172033" stroke-dasharray="3 3"/>
    <circle cx="${X(.2)}" cy="${Y(m20)}" r="4.5" fill="var(--c1)" stroke="#fffdf8" stroke-width="2"/>
    <circle cx="${X(.2)}" cy="${Y(b20)}" r="4.5" fill="var(--c2)" stroke="#fffdf8" stroke-width="2"/>
    <text x="${X(.2)+10}" y="${Y(m20)-6}" font-size="11" font-weight="700" fill="var(--c1)">模型 ${Math.round(m20*100)}%</text>
    <text x="${X(.2)+10}" y="${Y(b20)+16}" font-size="11" font-weight="700" fill="var(--c2)">基线 ${Math.round(b20*100)}%</text>
    <text x="${mL-40}" y="${mT-2}" font-size="10" font-family="ui-monospace,monospace" fill="#667078">累计召回</text>
    <text x="${W-mR}" y="${H-12}" text-anchor="end" font-size="10" font-family="ui-monospace,monospace" fill="#667078">筛选比例 →</text>
  </svg>`;
}
function rankBar(rows,key,label){
  const n=rows.length, bh=17, gap=3, W=760, mL=132, mR=64, H=n*(bh+gap)+22;
  const max=Math.max(...rows.map(r=>r[key]))||1;
  return `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="${label}">
   ${rows.map((r,i)=>{const y=i*(bh+gap)+8, w=(r[key]/max)*(W-mL-mR);
     return `<text x="${mL-10}" y="${y+bh-4}" text-anchor="end" font-size="11.5" font-weight="600" fill="#0d1720">${r.sym}</text>
     <rect x="${mL}" y="${y}" width="${Math.max(2,w)}" height="${bh}" rx="3" fill="${r.seen?'var(--c1)':'var(--c3)'}"/>
     <text x="${mL+Math.max(2,w)+8}" y="${y+bh-4}" font-size="11" font-family="ui-monospace,monospace" fill="#667078">${r[key].toFixed(2)}</text>`;}).join("")}
  </svg>`;
}


// 追加 20 个化合物：含刻意的近似物（用于演示冗余聚类）与频繁命中子（用于演示干扰剔除）
const EXTRA=[
 ["Paracetamol","CC(=O)Nc1ccc(O)cc1","COX / 解热镇痛"],
 ["Naproxen","COc1ccc2cc(C(C)C(=O)O)ccc2c1","COX 抑制"],
 ["Diclofenac","OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl","COX 抑制"],
 ["Celecoxib","Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1","COX-2 抑制"],
 ["Ciprofloxacin","O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O","细菌拓扑异构酶"],
 ["Vemurafenib","CCCS(=O)(=O)Nc1ccc(F)c(C(=O)c2c[nH]c3ncc(-c4ccc(Cl)cc4)cc23)c1F","BRAF 抑制"],
 ["Crizotinib","CC(Oc1cc(-c2cnn(C3CCNCC3)c2)cnc1N)c1c(Cl)ccc(F)c1Cl","ALK/MET 抑制"],
 ["Ruxolitinib","N#CCC(C1CCCC1)n1cc(-c2ncnc3[nH]ccc23)cn1","JAK 抑制"],
 ["Panobinostat","Cc1[nH]c2ccccc2c1CCNCc1ccc(C=CC(=O)NO)cc1","HDAC 抑制"],
 ["Belinostat","O=C(NO)C=Cc1cccc(S(=O)(=O)Nc2ccccc2)c1","HDAC 抑制"],
 ["Olaparib","O=C(c1ccc(CC2=NNC(=O)c3ccccc32)cc1F)N1CCN(C(=O)C2CC2)CC1","PARP 抑制"],
 ["Resveratrol","Oc1ccc(C=Cc2cc(O)cc(O)c2)cc1","频繁命中 / 多酚"],
 ["Genistein","O=c1c(-c2ccc(O)cc2)coc2cc(O)cc(O)c12","频繁命中 / 多酚"],
 // 刻意加入的近似物（内部历史库的典型情况）
 ["SAHA-analog-1","O=C(NO)CCCCCCC(=O)Nc1ccc(Cl)cc1","HDAC 类似物"],
 ["SAHA-analog-2","O=C(NO)CCCCCCC(=O)Nc1ccc(F)cc1","HDAC 类似物"],
 ["SAHA-analog-3","O=C(NO)CCCCCCC(=O)Nc1ccc(C)cc1","HDAC 类似物"],
 ["Gefitinib-analog","COc1cc2ncnc(Nc3ccc(F)c(Br)c3)c2cc1OCCCN1CCOCC1","EGFR 类似物"],
 ["Erlotinib-analog","C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCC)cc23)c1","EGFR 类似物"],
 ["Curcumin-analog","COc1cc(C=CC(=O)CC(=O)C=Cc2ccc(O)c(F)c2)ccc1O","多酚类似物"],
 ["Nitro-frag","O=[N+]([O-])c1ccc(C(=O)Nc2ccccc2)cc1","硝基芳香片段"]
];
const LIB=[...REF_DRUGS,...EXTRA];
const parsed=LIB.map(([n,s,m])=>({n,s,m,d:parseSMILES(s),f:fp(s),al:alertsOf(s)})).filter(o=>o.d&&o.d.ok);
console.log("库容量:",LIB.length,"可解析:",parsed.length);

// Butina 聚类（真实算法）：邻居数排序 → 取中心 → 吸收阈值内成员 → 移除 → 重复
function butina(items,cut){
  const N=items.length, sim=Array.from({length:N},()=>new Float32Array(N));
  for(let i=0;i<N;i++)for(let j=i+1;j<N;j++){const t=tanimoto(items[i].f,items[j].f);sim[i][j]=t;sim[j][i]=t;}
  const nb=items.map((_,i)=>{const a=[];for(let j=0;j<N;j++)if(j!==i&&sim[i][j]>=cut)a.push(j);return a;});
  const order=items.map((_,i)=>i).sort((a,b)=>nb[b].length-nb[a].length);
  const assigned=new Set(), clusters=[];
  for(const c of order){ if(assigned.has(c))continue;
    const mem=[c,...nb[c].filter(j=>!assigned.has(j))];
    mem.forEach(j=>assigned.add(j)); clusters.push(mem); }
  return {clusters,sim};
}
for(const cut of [0.55,0.65,0.75]){
  const {clusters}=butina(parsed,cut);
  const single=clusters.filter(c=>c.length===1).length;
  console.log(` Tc≥${cut}: ${clusters.length} 簇（${single} 单例, 最大簇 ${Math.max(...clusters.map(c=>c.length))}） → 压缩 ${(parsed.length/clusters.length).toFixed(2)}×`);
}
const {clusters,sim}=butina(parsed,0.65);
console.log("\n多成员簇（冗余候选）:");
clusters.filter(c=>c.length>1).forEach(c=>console.log("  ["+c.map(i=>parsed[i].n).join(" | ")+"]"));

// 干扰子
const flagged=parsed.filter(o=>o.al.length);
console.log("\n结构警报命中:",flagged.length+"/"+parsed.length,`(${(flagged.length/parsed.length*100).toFixed(0)}%)`);
flagged.forEach(o=>console.log("  "+o.n.padEnd(18)+o.al.join("; ")));

// MoA 归类：留一法（leave-drug-out）—— 每个化合物用其余全部做参考集
const refIdx=parsed.map((_,i)=>i);
let hit1=0,hit3=0,ood=0,scored=0;
const norm=s=>s.replace(/类似物|-analog.*/,"").trim();
parsed.forEach((q,qi)=>{
  const cand=refIdx.filter(i=>i!==qi).map(i=>({i,t:tanimoto(q.f,parsed[i].f)})).sort((a,b)=>b.t-a.t);
  if(cand[0].t<0.28){ood++;return;}
  scored++;
  const truth=norm(q.m);
  const top1=norm(parsed[cand[0].i].m), top3=cand.slice(0,3).map(c=>norm(parsed[c.i].m));
  if(top1===truth)hit1++; if(top3.includes(truth))hit3++;
});
console.log("\nleave-drug-out MoA 归类（"+scored+" 条可评分, "+ood+" 条 OOD）:");
console.log("  top-1 准确率 "+(hit1/scored*100).toFixed(0)+"%   top-3 "+(hit3/scored*100).toFixed(0)+"%");
const nClass=new Set(parsed.map(o=>norm(o.m))).size;
console.log("  MoA 类别数 "+nClass+" → 随机 top-1 期望 "+(100/nClass).toFixed(0)+"%  ⇒ 富集 "+((hit1/scored)/(1/nClass)).toFixed(1)+"×");
// 板容量下的覆盖
const P=makePlate(parsed.map((o,i)=>({...o,_q:0})),384,1);
console.log("\n384 板: 样品孔 "+P.capacity+" → 本库 "+parsed.length+" 条全部可入板；覆盖 "+nClass+" 个 MoA 类别");
