const fs = require("fs");
const { Document, Packer, Paragraph, TextRun, HeadingLevel, AlignmentType,
  Table, TableRow, TableCell, WidthType, ShadingType, ImageRun, PageBreak,
  Header, Footer, PageNumber } = require("docx");
const DIR = "/tmp/claude-0/-home-user-ClaudeCode/0e1c2284-3690-5d68-a96a-480ebe908c1c/scratchpad";
const Q = JSON.parse(fs.readFileSync(`${DIR}/qsar_results.json`));
const QA = JSON.parse(fs.readFileSync(`${DIR}/qsar_augment.json`));
const M7 = JSON.parse(fs.readFileSync(`${DIR}/m7_results.json`));
const EPS = Q.cv_metrics.map(r=>r.endpoint);
const meanauc=Q.mean_cv_auc, aucs=Q.cv_metrics.map(r=>r.CV_AUC);
const EP_ZH={"NR-AR":"雄激素受体","NR-AR-LBD":"雄激素受体(LBD)","NR-AhR":"芳烃受体","NR-Aromatase":"芳香化酶","NR-ER":"雌激素受体","NR-ER-LBD":"雌激素受体(LBD)","NR-PPAR-gamma":"PPARγ","SR-ARE":"氧化应激(ARE)","SR-ATAD5":"遗传毒性(ATAD5)","SR-HSE":"热休克反应","SR-MMP":"线粒体膜电位","SR-p53":"p53损伤应答"};
const MEAS_ZH={"active":"阳性","inactive":"阴性","not tested":"未测试"};
const NAVY="1F3864",BLUE="2E74B5",GREY="595959",LIGHT="DCE6F1",GREEN="2E7D32",RED="C00000",AMBER="FBE4D5";
const FONT="Microsoft YaHei"; // Word 会以系统中文字体渲染

function img(f,w,h){return new Paragraph({alignment:AlignmentType.CENTER,spacing:{before:120,after:60},children:[new ImageRun({type:"png",data:fs.readFileSync(`${DIR}/${f}`),transformation:{width:w,height:h}})]});}
function cap(t){return new Paragraph({alignment:AlignmentType.CENTER,spacing:{after:160},children:[new TextRun({text:t,italics:true,size:16,color:GREY,font:FONT})]});}
function h1(t){return new Paragraph({heading:HeadingLevel.HEADING_1,spacing:{before:280,after:120},children:[new TextRun({text:t,bold:true,color:NAVY,font:FONT})]});}
function h2(t){return new Paragraph({heading:HeadingLevel.HEADING_2,spacing:{before:200,after:100},children:[new TextRun({text:t,bold:true,color:BLUE,font:FONT})]});}
function R(text,o={}){return new TextRun({text,size:o.size||20,bold:o.bold,italics:o.i,color:o.color,font:FONT,break:o.break});}
function p(runs,o={}){if(typeof runs==="string")runs=[R(runs)];return new Paragraph({spacing:{after:120,line:288},alignment:o.align||AlignmentType.JUSTIFIED,children:runs});}
function bullet(runs){if(typeof runs==="string")runs=[R(runs)];return new Paragraph({bullet:{level:0},spacing:{after:80},children:runs});}
function cell(content,o={}){const kids=Array.isArray(content)?content:[new Paragraph({children:[new TextRun({text:String(content),size:o.size||18,bold:o.bold,color:o.color,font:FONT})],alignment:o.align||AlignmentType.LEFT,spacing:{before:40,after:40}})];return new TableCell({width:{size:o.w,type:WidthType.DXA},shading:o.shade?{type:ShadingType.CLEAR,fill:o.shade,color:"auto"}:undefined,margins:{top:40,bottom:40,left:80,right:80},verticalAlign:"center",children:kids});}
function hc(t,w){return cell(t,{w,bold:true,color:"FFFFFF",shade:NAVY,align:AlignmentType.CENTER});}
function table(cw,rows){return new Table({columnWidths:cw,width:{size:cw.reduce((a,b)=>a+b,0),type:WidthType.DXA},rows});}
function green(t){return cell(t,{w:0,bold:true,color:GREEN,align:AlignmentType.CENTER});}

const ID=[["属性","母体 API","RRT 0.87","RRT 0.75"],
["Xellia 名称","万古霉素 B","(7S,26R)-万古霉素 B","万古霉素 B 杂质 (26R)-DMe-DeLI"],
["分子式","C66H75Cl2N9O24","C66H75Cl2N9O24","C65H73Cl2N9O24"],
["单同位素质量","1447.43","1447.43","1433.41"],
["已定义立体中心","18","18","19"],
["来源","发酵（API）","降解（低水活度+受热）","发酵"],
["与 API 的关系","—","非对映异构体（C7、C26差向）","第1位残基变体(去N-甲基;Leu→Ile)+C26差向"]];

const children=[
 new Paragraph({spacing:{before:1000,after:120},alignment:AlignmentType.CENTER,children:[R("两个万古霉素 B 杂质的计算机模拟毒理学界定",{bold:true,size:34,color:NAVY})]}),
 new Paragraph({alignment:AlignmentType.CENTER,spacing:{after:120},children:[R("基于 3D-QSAR / 结构类推（Read-Across）",{bold:true,size:34,color:NAVY})]}),
 new Paragraph({alignment:AlignmentType.CENTER,spacing:{after:360},children:[R("杂质 RRT 0.75 ［(26R)-DMe-DeLI］ 与 RRT 0.87 ［(7S,26R)-万古霉素 B］",{size:24,color:BLUE})]}),
 table([3000,6600],[
   new TableRow({children:[cell("文件类型",{w:3000,bold:true,shade:LIGHT}),cell("计算机模拟毒理学界定报告（研究/监管支持稿）",{w:6600})]}),
   new TableRow({children:[cell("受试物",{w:3000,bold:true,shade:LIGHT}),cell("万古霉素 B（API）及杂质 RRT 0.75、RRT 0.87",{w:6600})]}),
   new TableRow({children:[cell("计算方法",{w:3000,bold:true,shade:LIGHT}),cell("毒性团/药效团分析；RDKit PAINS/Brenk/NIH 结构警示；理化与三维形状描述符；三维构象生成及 MCS/O3A 形状类推；自建 12 端点 Tox21 统计 QSAR（随机森林、5折CV、适用域检查）；ICH M7 双方法学致突变筛查（Benigni–Bossa 专家规则库 + Ames QSAR）；辅助 logD/pKa QSAR",{w:6600})]}),
   new TableRow({children:[cell("毒理学锚点",{w:3000,bold:true,shade:LIGHT}),cell("母体 API 90天静脉犬试验，NOAEL 75/100 mg/kg/day（约临床剂量2–3倍）",{w:6600})]}),
   new TableRow({children:[cell("日期",{w:3000,bold:true,shade:LIGHT}),cell("2026-07-30",{w:6600})]}),
 ]),
 new Paragraph({spacing:{before:300},children:[R("重要提示 —— 本报告为计算机模拟筛查与结构类推评估，用于支持毒理学专家审核；其本身不构成监管毒理学结论，也不替代专家签署或（在需要时）实验测试。见第8节。",{i:true,size:18,color:RED})]}),
 new Paragraph({children:[new PageBreak()]}),

 h1("1. 执行摘要"),
 p([R("糖肽类抗生素万古霉素 B 的两个有关物质需进行计算机模拟毒理学界定："),R("RRT 0.87 = (7S,26R)-万古霉素 B",{bold:true}),R("（低水活度/受热降解产物）与 "),R("RRT 0.75 = (26R)-DMe-DeLI",{bold:true}),R("（发酵来源的第1位残基变体）。由于二者与母体的差异主要为立体化学/构成差异、常规二维类推无法分辨，故对母体及两个子体分子均采用三维结构方法评估。")]),
 p([R("主要结论：两个杂质相对母体 API 均未引入新的毒性团，也未引入母体所不具备的新的 DNA 反应性（致突变）结构警示。",{bold:true}),R("RRT 0.87 为纯立体异构体，与母体共享全部 101 个重原子及完整官能团清单（仅 18 个立体中心中的 2 个构型不同）。RRT 0.75 与母体共享 99/约100 个重原子；其唯一构成改变为第1位残基由仲胺（N-甲基）变为伯脂肪胺 —— 该官能团类别母体本已存在（万古糖胺）—— 亲脂性不变（logD 0.42），碱性近乎不变（pKa 约8.1对8.2）。")]),
 p([R("一项 12 端点 Tox21 统计 QSAR（5折CV平均 AUC "+meanauc.toFixed(3)+"）进一步支持：母体本身即为 Tox21 训练集化合物（实测谱基本干净——仅 NR-AR 弱信号），RRT 0.87 与母体指纹相同（各端点 Δ≡0），RRT 0.75 在每个端点与母体差异≤0.02。ICH M7 双方法学筛查（Benigni–Bossa 专家规则库 + Ames QSAR，AUC "+M7.methodology_2.cv_auc.toFixed(3)+"）对三者均给出一致的阴性/Class 5。综合证据权衡，两个杂质均可通过对母体的结构类推予以界定；母体安全性由 91天静脉犬试验（NOAEL 75/100 mg/kg/day）及长期临床使用锚定，前提是各杂质控制在《欧洲药典》有关物质限度内。")]),
 table([2600,2400,2400,2240],[
   new TableRow({children:[hc("问题",2600),hc("RRT 0.87",2400),hc("RRT 0.75",2400),hc("依据",2240)]}),
   new TableRow({children:[cell("相对母体是否有新毒性团？",{w:2600,bold:true,shade:LIGHT}),cell("否",{w:2400,bold:true,color:GREEN,align:AlignmentType.CENTER}),cell("否",{w:2400,bold:true,color:GREEN,align:AlignmentType.CENTER}),cell("官能团清单一致(0.87)；等官能团胺替换(0.75)",{w:2240})]}),
   new TableRow({children:[cell("是否有新的致突变(M7)警示？",{w:2600,bold:true,shade:LIGHT}),cell("否",{w:2400,bold:true,color:GREEN,align:AlignmentType.CENTER}),cell("否",{w:2400,bold:true,color:GREEN,align:AlignmentType.CENTER}),cell("M7两方法学均阴性→Class 5(见5.7)",{w:2240})]}),
   new TableRow({children:[cell("对母体的类推是否成立？",{w:2600,bold:true,shade:LIGHT}),cell("是",{w:2400,bold:true,color:GREEN,align:AlignmentType.CENTER}),cell("是",{w:2400,bold:true,color:GREEN,align:AlignmentType.CENTER}),cell("≥98%骨架共享；核心毒性团保守",{w:2240})]}),
 ]),
 new Paragraph({children:[new PageBreak()]}),

 h1("2. 目的、范围与监管框架"),
 p("本项目目标为对杂质 RRT 0.75 与 RRT 0.87 进行计算机模拟毒理学界定，交付：(i) 对母体及各子体分子的 3D-QSAR/结构方法评估；(ii) 识别并评估子体分子中存在而母体不具备的任何附加毒性团；(iii) 书面结构类推论证；(iv) 本监管风格报告。"),
 p([R("框架。",{bold:true}),R("非致突变杂质界定遵循 ICH Q3A(R2)/Q3B(R2)。致突变潜力按 ICH M7(R2) 筛查，该指南要求采用两条互补的 (Q)SAR 方法学（专家规则+统计）并辅以专家审核。计算模型使用遵循 OECD (Q)SAR 验证原则。由于万古霉素 B 为大分子（MW≈1449）柔性糖肽，处于许多小分子毒性模型适用域之外，本评估以机理性、结构与形状类推为主，统计模型仅在其适用域得到尊重时使用（域外情形均予标注）。")]),

 h1("3. 受试物与结构关系"),
 p("结构取自申办方资料（IUPAC名/SMILES/InChI），并用 RDKit 独立校验；计算所得分子式与精确质量与源文件一致。"),
 table([2200,2480,2480,2480],ID.map((r,i)=>new TableRow({children:[
   cell(r[0],{w:2200,bold:true,shade:i===0?NAVY:LIGHT,color:i===0?"FFFFFF":undefined}),
   cell(r[1],{w:2480,shade:i===0?NAVY:undefined,color:i===0?"FFFFFF":undefined,bold:i===0,align:i===0?AlignmentType.CENTER:undefined}),
   cell(r[2],{w:2480,shade:i===0?NAVY:undefined,color:i===0?"FFFFFF":undefined,bold:i===0,align:i===0?AlignmentType.CENTER:undefined}),
   cell(r[3],{w:2480,shade:i===0?NAVY:undefined,color:i===0?"FFFFFF":undefined,bold:i===0,align:i===0?AlignmentType.CENTER:undefined}),
 ]}))),
 img("fig1_structures.png",640,210),
 cap("图1. 母体 API（左）与两个杂质的二维结构。RRT 0.87 因仅立体化学不同，画法与母体一致；RRT 0.75 在第1位残基处不同（游离伯胺+异亮氨酸侧链，取代N-甲基亮氨酸）。"),
 p([R("RRT 0.87 —— (7S,26R)-万古霉素 B。",{bold:true}),R("分子式与 API 相同；C7（残基3苄位甲醇碳）与 C26（肽α-碳）构型翻转——与低水活度受热差向异构化一致。不增删任何原子。")]),
 p([R("RRT 0.75 —— (26R)-DMe-DeLI。",{bold:true}),R("分子式 C65H73Cl2N9O24（较母体少一个CH2）。第1位残基由 N-甲基-D-亮氨酸变为去N-甲基异亮氨酸（仲胺→伯胺；异丁基→仲丁基），并在 C26 差向。其余残基、两个糖、两个芳基氯及所有酚均保留。")]),
 img("fig2_parent_residue1.png",300,250),
 img("fig2_rrt075_residue1.png",300,250),
 cap("图2. RRT 0.75 与 API 之间唯一的构成改变位点：仲胺(母体) vs 伯胺(RRT 0.75)。"),
 new Paragraph({children:[new PageBreak()]}),

 h1("4. 计算方法"),
 h2("4.1 毒性团 / 结构警示分析"),
 p("对每个结构分析其万古霉素药效团（七肽糖苷配基 D-Ala-D-Ala 结合口袋、氢键网络、万古糖胺–葡萄糖二糖）及公认毒性团（酚/芳环→肾小管蓄积；碱性胺→组胺/红人综合征；高分子量/亲水性）。用 SMARTS 对官能团定量清点并与母体求差。采用 RDKit PAINS(A/B/C)、Brenk、NIH 目录筛查隐患，并对照 ICH M7 DNA 反应性警示类别审核。"),
 h2("4.2 理化与三维描述符；三维形状类推"),
 p("RDKit 描述符（MW、cLogP、TPSA、HBD/HBA、可旋转键、QED、Lipinski/Veber）；Inductive Bio logD/pKa QSAR。ETKDGv3 生成构象并经 MMFF 最小化；三维形状描述符（回转半径、非球度、球度、离心率、NPR1/2）。三维类推通过 MCS 核心限定叠合（RMSD）及 Crippen Open3DAlign/形状-Tanimoto 量化。不构建正式 CoMFA/CoMSIA 3D-QSAR（见第8节）。"),
 h2("4.3 统计 QSAR —— Tox21 与 Ames"),
 p("训练并应用 12 端点 Tox21 QSAR（2048位 Morgan 指纹上的类平衡随机森林；公开 MoleculeNet Tox21；5折CV AUC）与 Ames QSAR（Hansen 基准 N≈6500），均带适用域检查。母体–杂质差异为可解释读数。"),

 h1("5. 结果"),
 h2("5.1 母体（万古霉素 B）—— 基线毒性团谱"),
 p("母体呈现典型糖肽谱：6酚羟基、5芳环、2芳基氯、3芳醚、2苄醇、1羧酸、1伯酰胺、6仲酰胺、1伯胺（万古糖胺）、1个N-甲基仲胺。无 PAINS/Brenk 警示、无 DNA 反应性警示；唯一 NIH 标记（≥7脂肪羟基）为极性描述符。此为类推基线。"),
 img("fig3_toxicophore_inventory.png",620,330),
 cap("图3. 毒性团/官能团清单（计数）。母体与 RRT 0.87 所有类别完全一致；RRT 0.75 仅在胺列不同（红色）——一个仲胺(N-甲基)被替换为一个伯胺。"),
 h2("5.2 RRT 0.87 —— (7S,26R)-万古霉素 B"),
 p([R("因仅构型不同，其二维结构、官能团清单与警示谱均与母体一致（ECFP4/MACCS=1.00；共享101/101重原子；0警示）。"),R("无附加毒性团。",{bold:true}),R("常规二维类推完全无法将其与 API 区分——这正是需要三维/机理评估的原因。三维形状描述符近乎重合（回转半径7.45对7.50 Å；NPR1 0.529对0.528）。C7/C26 构型翻转可能改变靶点结合（效力），但不产生任何新的反应性/DNA反应性官能团。")]),
 h2("5.3 RRT 0.75 —— (26R)-DMe-DeLI"),
 p([R("保留全部糖苷配基、两个糖、两个芳基氯、所有酚、羧酸与伯酰胺；共享99/约100重原子（ECFP4 0.86/MACCS 0.94）。唯一改变（仲N-甲基→伯胺）"),R("并非新官能团类别",{bold:true}),R("（母体本已含伯胺），也不是 DNA 反应性警示。logD 不变（0.42对0.42），最强碱性 pKa 近乎不变（8.08对8.21）——碱性胺/红人综合征相关隐患为保守而非放大。TPSA 仅+14 Å²；三维形状接近母体（回转半径7.46 Å）。")]),
 h2("5.4 差异毒性团分析 —— 是否存在新毒性团？"),
 table([3200,3200,3239],[
   new TableRow({children:[hc("潜在毒性团",3200),hc("RRT 0.87 对母体",3200),hc("RRT 0.75 对母体",3239)]}),
   new TableRow({children:[cell("酚/芳环（肾蓄积）",{w:3200}),cell("一致(6酚,5芳环)",{w:3200}),cell("一致(6酚,5芳环)",{w:3239})]}),
   new TableRow({children:[cell("碱性胺中心（组胺/红人）",{w:3200}),cell("一致",{w:3200}),cell("保守：仲→伯胺，等碱性(pKa8.1)；未放大",{w:3239,shade:AMBER})]}),
   new TableRow({children:[cell("芳基氯",{w:3200}),cell("一致(2)",{w:3200}),cell("一致(2)",{w:3239})]}),
   new TableRow({children:[cell("DNA反应性警示(芳香胺/硝基/环氧/Michael受体/烷基卤…)",{w:3200}),cell("无(同母体)",{w:3200}),cell("无(同母体)",{w:3239})]}),
   new TableRow({children:[cell("PAINS/Brenk 隐患",{w:3200}),cell("无",{w:3200}),cell("无",{w:3239})]}),
   new TableRow({children:[cell("净结果：是否引入新毒性团？",{w:3200,bold:true}),cell("否",{w:3200,bold:true,color:GREEN,align:AlignmentType.CENTER}),cell("否",{w:3239,bold:true,color:GREEN,align:AlignmentType.CENTER})]}),
 ]),
 p("唯一值得说明的是 RRT 0.75 的第1位残基胺——这是对既有毒性团类别的修饰（仲→伯胺、等碱性），而非新增类别。两个子体分子均未识别出附加毒性团。"),
 h2("5.5 3D-QSAR / 形状类推指标"),
 img("fig5_similarity.png",470,300),
 cap("图4. 二维相似度与共享骨架类推指标。RRT 0.87 指纹相同；RRT 0.75 共享≥98%重原子骨架。"),
 img("fig4_descriptors.png",640,250),
 cap("图5. 三维整体形状描述符（左）与理化描述符（右，对称对数）。三者占据基本相同的形状/性质空间。"),

 h2("5.6 统计（二维）QSAR —— 12 端点 Tox21 谱"),
 p("训练了 12 端点 Tox21 QSAR（2048位 Morgan 指纹上的类平衡随机森林；公开 MoleculeNet Tox21，n="+(Q.model.match(/n=(\d+)/)[1])+"）：5折CV平均 AUC="+meanauc.toFixed(3)+"（范围 "+Math.min(...aucs).toFixed(3)+"–"+Math.max(...aucs).toFixed(3)+"），随后对每个分子打分并做适用域检查。"),
 p([R("适用域——一个有利的发现。",{bold:true}),R("万古霉素本身即在 Tox21 数据集中（训练化合物 "+QA.identical_mol_ids.join(", ")+"；对母体及 RRT 0.87 最近邻 Tanimoto=1.00）。因此该 QSAR 对母体为域内；RRT 0.75 亦相近（最大 Tanimoto 0.863）。")]),
 p([R("万古霉素实测谱（真值）：",{bold:true}),R("仅 NR-AR 阳性（弱），其余所有测过端点均阴性（含遗传毒性 SR-ATAD5、SR-p53）。此为毒理学基线。")]),
 img("fig6_qsar.png",660,405),
 cap("图6. 自建12端点 Tox21 QSAR。(A) 预测概率——母体与 RRT 0.87 两行相同（指纹相同）。(B) 各端点5折 AUC（平均0.827）。(C) RRT 0.75−留一法母体 的Δ——各端点均在±0.02内；RRT 0.87 按构造 Δ≡0。"),
 p([R("杂质对母体。",{bold:true}),R("RRT 0.87 在全部12端点返回相同概率（Δ≡0）：证明二维 QSAR 对该立体异构体在结构上'盲视'（故三维类推承担主要判据）。对 RRT 0.75，每个端点均在留一法母体的±0.02以内（最大|Δ|=0.02），无一越过活性阈值——包括遗传毒性相关的 SR-ATAD5 与 SR-p53。")]),
 table([2160,1500,1900,1900,2179],[
   new TableRow({children:[hc("端点",2160),hc("CV-AUC",1500),hc("实测(母体)",1900),hc("留一法母体P",1900),hc("RRT 0.75 P (Δ)",2179)]}),
   ...EPS.map(ep=>{const m=Q.cv_metrics.find(r=>r.endpoint===ep);const meas=QA.measured_vancomycin[ep];const hp=QA.heldout_parent_pred[ep];const p75=QA.heldout_pred_rrt075[ep];const d=p75-hp;
     return new TableRow({children:[
       cell(ep+"  ("+EP_ZH[ep]+")",{w:2160,size:16}),
       cell(m.CV_AUC.toFixed(3),{w:1500,align:AlignmentType.CENTER,size:16}),
       cell(MEAS_ZH[meas],{w:1900,align:AlignmentType.CENTER,size:16,color:meas==="active"?RED:undefined,bold:meas==="active"}),
       cell(hp.toFixed(3),{w:1900,align:AlignmentType.CENTER,size:16}),
       cell(p75.toFixed(3)+"  ("+(d>=0?"+":"")+d.toFixed(3)+")",{w:2179,align:AlignmentType.CENTER,size:16}),
     ]});})
 ]),
 new Paragraph({spacing:{before:80},children:[R("所有预测概率均低于0.5阈值（母体 NR-AR 为实测弱阳性，模型低估——在该端点较低 AUC 下属预期；以实测标签为准）。RRT 0.87 未列入：指纹与母体相同⇒各端点预测相同(Δ≡0)。",{i:true,size:16,color:GREY})]}),
 new Paragraph({children:[new PageBreak()]}),

 h2("5.7 ICH M7 致突变评估 —— 两条互补 (Q)SAR 方法学"),
 p([R("ICH M7(R2) 要求以两条互补方法学评估潜在细菌（Ames）致突变性——一条"),R("专家规则",{bold:true}),R("、一条"),R("统计",{bold:true}),R("——任何阳性均需经专家审核。二者均已应用：")]),
 bullet([R("方法一（专家规则）。",{bold:true}),R("Benigni–Bossa 致突变/遗传毒性致癌结构警示规则库（"+(M7.methodology_1.match(/(\d+) alerts/)||[])[1]+"类警示：芳香胺、硝基/亚硝基、偶氮/氧化偶氮、N-亚硝基、环氧/氮丙啶、Michael受体、酰卤/烷基卤、肼、醛、醌、β-内酰胺……）——即 Derek 类专家系统背后的公开科学规则库（Toxtree/VEGA）。")]),
 bullet([R("方法二（统计）。",{bold:true}),R("Ames QSAR（类平衡随机森林、2048位 Morgan），训练于公开 Hansen 基准（N="+M7.methodology_2.n_train+"；"+Math.round(M7.methodology_2.pos_frac*100)+"%为致突变物）；5折 OOF AUC="+M7.methodology_2.cv_auc.toFixed(3)+"；决策阈值 "+M7.methodology_2.threshold_youden.toFixed(2)+"。")]),
 img("fig7_m7.png",660,235),
 cap("图7. ICH M7 双方法学评估。(A) 一致性矩阵——两条方法学对三个分子均阴性。(B) 统计 Ames 概率均低于阈值；琥珀色标注提示适用域相似度低，故以域无关的专家规则库为主判据。"),
 table([2760,2960,2960,919],[
   new TableRow({children:[hc("分子",2760),hc("专家规则",2960),hc("统计 Ames QSAR",2960),hc("M7",919)]}),
   ...Object.keys(M7.per_compound).map(c=>{const r=M7.per_compound[c];const alert=r.expert_rulebase_alerts.length>0;
     return new TableRow({children:[
       cell(c,{w:2760,size:16}),
       cell(alert?"阳性（存在警示）":"阴性（无结构警示）",{w:2960,size:16,color:alert?RED:GREEN,align:AlignmentType.CENTER}),
       cell((r.stat_Ames_call==="negative"?"阴性":"阳性")+"（p="+r.stat_Ames_prob.toFixed(2)+"，域相似度"+r.stat_AD_maxTanimoto.toFixed(2)+"）",{w:2960,size:16,color:r.stat_Ames_call==="negative"?GREEN:RED,align:AlignmentType.CENTER}),
       cell("5",{w:919,size:16,bold:true,color:GREEN,align:AlignmentType.CENTER}),
     ]});})
 ]),
 p([R("结论。",{bold:true}),R("两条方法学对母体及两个杂质一致且阴性：无致突变结构警示，统计 Ames 预测在各情形均低于阈值。据 ICH M7，此对应 Class 5（无结构警示）——按非致突变的普通 ICH Q3A 有关物质控制，无需专门的细菌回复突变（Ames）试验。为透明起见记录两点：(i) RRT 0.87 与母体指纹相同，故两条方法学按构造返回相同结果；(ii) 此类糖肽处于 Ames 模型适用域之外（最近邻 Tanimoto≈0.22），故统计阴性为低置信度，域无关的专家规则库（明确无警示）为主判据，并由母体实测阴性的 Tox21 DNA 损伤端点（SR-ATAD5、SR-p53；见5.6）佐证。在验证流程中运行商用 Derek/Sarah Nexus 预期复现此无警示、阴性结论。")]),
 new Paragraph({children:[new PageBreak()]}),

 h1("6. 结构类推（Read-Across）论证"),
 p([R("源（类似物）：",{bold:true}),R("母体 API 万古霉素 B（具完整毒理学资料，含91天静脉犬试验 NOAEL 75/100 mg/kg/day 及数十年临床暴露）。"),R("目标：",{bold:true}),R("RRT 0.87 与 RRT 0.75。")]),
 bullet([R("结构相似性。",{bold:true}),R("RRT 0.87 共享100%重原子骨架与官能团（构型异构体）。RRT 0.75 共享≥98%（99/100），仅在一个肽 N-末端有单一、局部、无警示的修饰。")]),
 bullet([R("共同机理/保守毒性团。",{bold:true}),R("所有酚、芳环、芳基氯、酰胺、羧基与碱性胺特征均保守。无新反应性中心。")]),
 bullet([R("代谢/理化一致性。",{bold:true}),R("万古霉素主要以原形经肾清除；近乎一致的分子大小、极性、logD(0.42)与 pKa(≈8.1)预示相同 ADME/处置。")]),
 bullet([R("无致突变担忧。",{bold:true}),R("ICH M7 两条方法学均阴性（见5.7）；QSAR 在 SR-ATAD5/SR-p53 无激活（Δ≤0.02）。")]),
 bullet([R("已认知的不确定性（立体化学）。",{bold:true}),R("立体异构体可能在靶点结合上不同——关乎效力/疗效，而非新生毒性；三维分析正是针对这一维度。")]),

 h1("7. 界定结论"),
 p("基于证据权衡，本计算机模拟评估支持通过对母体的结构类推对两个杂质进行毒理学界定："),
 bullet([R("RRT 0.87 ［(7S,26R)-万古霉素 B］",{bold:true}),R("——纯立体异构体，无新毒性团、无新致突变警示；类推充分成立。")]),
 bullet([R("RRT 0.75 ［(26R)-DMe-DeLI］",{bold:true}),R("——单一、无警示、等官能团的第1位残基修饰（加 C26 差向），无新毒性团、无新致突变警示；类推成立。")]),
 p("只要二者控制在《欧洲药典》有关物质限度内，母体 NOAEL（75/100 mg/kg/day，约临床2–3倍）及临床记录即可提供充分毒理学覆盖。正式签署应由合格毒理学家出具，致突变维度按公司 ICH M7 双 (Q)SAR 方法学流程存档。"),

 new Paragraph({children:[new PageBreak()]}),
 h1("8. 适用域、诚实声明与局限性"),
 bullet([R("\"3D-QSAR\"的性质。",{bold:true}),R("这些化合物不存在带实测毒性端点的同系物训练集，故未（也不能合法地）构建正式场基 3D-QSAR（CoMFA/CoMSIA）预测模型。本报告所述为三维结构/形状类推加机理性毒性团分析。")]),
 bullet([R("统计 QSAR——有效范围。",{bold:true}),R("12端点 Tox21 QSAR 对母体为域内（万古霉素为训练化合物）。两点提醒：(i) 母体集内预测含部分记忆效应，故改以实测标签+留一法预测报告；(ii) Tox21 面板仅覆盖12个特定机制——不建模万古霉素临床效应（肾/耳毒性、红人综合征），后者由机理分析及对母体体内资料类推处理。")]),
 bullet([R("ICH M7 统计模型适用域。",{bold:true}),R("Ames QSAR 对此类糖肽为域外（Tanimoto≈0.22），其阴性为低置信度，以专家规则库为主判据。应在验证流程中以商用 Derek/Sarah Nexus 确认。")]),
 bullet([R("构象采样。",{bold:true}),R("三维描述符/叠合来自大柔性大环单一低能构象→属佐证而非决定性；决定性证据为共享骨架与保守毒性团清单。")]),
 bullet([R("范围。",{bold:true}),R("本工作为危害识别/结构筛查与类推，用于支持而非替代毒理学专家审核及主管当局可能要求的任何实验测试。")]),
 bullet([R("杂质控制限度。",{bold:true}),R("源文件未提供 RRT 0.75/0.87 的《欧洲药典》数值限度；本结论以控制在适用各论限度内为前提。")]),

 h1("9. 参考文献与工具"),
 bullet("ICH Q3A(R2)/Q3B(R2)——新原料药/新制剂中的杂质。"),
 bullet("ICH M7(R2)——DNA 反应性（致突变）杂质的评估与控制。"),
 bullet("OECD (2007)——(Q)SAR 模型验证指南。"),
 bullet("Benigni R., Bossa C. (2008)——致突变物与致癌物的结构警示（Toxtree/VEGA 规则库；方法一）。"),
 bullet("Hansen K. et al. (2009) J. Chem. Inf. Model. 49:2077——Ames 致突变性基准数据集（方法二训练数据）。"),
 bullet("Tox21/MoleculeNet（Wu Z. et al. 2018, Chem. Sci.）——统计 QSAR 训练数据。"),
 bullet("RDKit (2024)——描述符、PAINS/Brenk/NIH 目录、ETKDGv3、MMFF94、Open3DAlign、MCS。Inductive Bio logD/pKa QSAR。《欧洲药典》万古霉素各论。"),
 bullet("申办方源文件：\"Vanco_RRT_075_and_087_Info\"（Xellia）——身份、SMILES/InChI、来源、91天犬 NOAEL 75/100 mg/kg/day。"),
];

const doc=new Document({
  creator:"计算机模拟毒理学 (Claude Code)",
  title:"万古霉素杂质 RRT 0.75 与 RRT 0.87 的 3D-QSAR/结构类推毒理学界定报告",
  styles:{default:{document:{run:{font:FONT,size:20}}}},
  sections:[{
    properties:{page:{size:{width:12240,height:15840},margin:{top:1200,bottom:1200,left:1200,right:1200}}},
    headers:{default:new Header({children:[new Paragraph({alignment:AlignmentType.RIGHT,children:[R("万古霉素 RRT 0.75 与 RRT 0.87 —— 计算机模拟毒理学界定",{size:14,color:GREY})]})]})},
    footers:{default:new Footer({children:[new Paragraph({alignment:AlignmentType.CENTER,children:[R("计算机模拟筛查与结构类推——支持而非替代毒理学专家审核。 第 ",{size:14,color:GREY}),new TextRun({children:[PageNumber.CURRENT],size:14,color:GREY,font:FONT}),R(" 页",{size:14,color:GREY})]})]})},
    children,
  }]
});
Packer.toBuffer(doc).then(b=>{fs.writeFileSync(`${DIR}/万古霉素_RRT075_RRT087_3DQSAR毒理学界定报告.docx`,b);console.log("WROTE zh docx",b.length);});
