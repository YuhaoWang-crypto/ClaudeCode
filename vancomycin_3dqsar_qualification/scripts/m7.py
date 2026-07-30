"""
ICH M7 mutagenicity assessment by TWO complementary (Q)SAR methodologies:
  (1) EXPERT RULE-BASED  - Benigni-Bossa structural-alert rulebase for
      mutagenicity/genotoxic carcinogenicity (the public scientific basis of
      Derek-type systems; as implemented in Toxtree/VEGA). RDKit SMARTS.
  (2) STATISTICAL QSAR   - RandomForest Ames model trained on the public
      Hansen benchmark (N=6512), Morgan r2 2048-bit, 5-fold OOF AUC + AD.
Applied to parent vancomycin B and impurities RRT 0.87 and RRT 0.75.
"""
import json, numpy as np, pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
OUT="/tmp/claude-0/-home-user-ClaudeCode/0e1c2284-3690-5d68-a96a-480ebe908c1c/scratchpad"

MOLS={
 "Parent (Vancomycin B)":"OC1=CC([C@@H](C(O)=O)NC2=O)=C(C3=C(O)C=CC([C@@H](NC([C@H](NC([C@H](CC(N)=O)NC([C@@H]([C@H](O)C(C=C4Cl)=CC=C4O5)NC([C@H](NC)CC(C)C)=O)=O)=O)C6=CC5=C(O[C@H]7[C@H](O[C@]8([H])O[C@@H](C)[C@@H](O)[C@](N)(C)C8)[C@@H](O)[C@H](O)[C@@H](CO)O7)C9=C6)=O)C(N[C@H]2[C@H](O)C%10=CC(Cl)=C(O9)C=C%10)=O)=C3)C(O)=C1",
 "RRT 0.87 (7S,26R)-Vanco B":"OC1=CC([C@H](C(O)=O)NC2=O)=C(C3=C(O)C=CC([C@@H](NC([C@H](NC([C@H](CC(N)=O)NC([C@@H]([C@@H](O)C(C=C4Cl)=CC=C4O5)NC([C@H](NC)CC(C)C)=O)=O)=O)C6=CC5=C(O[C@H]7[C@H](O[C@]8([H])O[C@@H](C)[C@@H](O)[C@](N)(C)C8)[C@@H](O)[C@H](O)[C@@H](CO)O7)C9=C6)=O)C(N[C@H]2[C@H](O)C%10=CC(Cl)=C(O9)C=C%10)=O)=C3)C(O)=C1",
 "RRT 0.75 (26R)-DMe-DeLI":"OC1=CC([C@H](C(O)=O)NC2=O)=C(C3=C(O)C=CC([C@@H](NC([C@H](NC([C@H](CC(N)=O)NC([C@@H]([C@H](O)C(C=C4Cl)=CC=C4O5)NC([C@H](N)[C@H](CC)C)=O)=O)=O)C6=CC5=C(O[C@H]7[C@H](O[C@]8([H])O[C@@H](C)[C@@H](O)[C@](N)(C)C8)[C@@H](O)[C@H](O)[C@@H](CO)O7)C9=C6)=O)C(N[C@H]2[C@H](O)C%10=CC(Cl)=C(O9)C=C%10)=O)=C3)C(O)=C1",
}

# ---------- (1) Benigni-Bossa mutagenicity/genotoxic-carcinogenicity rulebase ----------
# Representative implementation of the Benigni-Bossa structural alerts for
# mutagenicity & genotoxic carcinogenicity (Toxtree SA1-SA30 family).
BB_ALERTS = {
 "SA1 Acyl halide":"[CX3](=[OX1])[F,Cl,Br,I]",
 "SA2 Alkyl/benzyl ester of sulfonic/phosphonic acid":"[$([SX4](=O)(=O)[OX2][CX4]),$([PX4](=O)([OX2][CX4])[OX2][CX4])]",
 "SA3 N-methylol":"[NX3][CX4;H2][OX2H]",
 "SA4 Monohaloalkene":"[CX3]=[CX3][F,Cl,Br,I]",
 "SA5 S/N mustard":"[Cl,Br,I][CX4][CX4][NX3,SX2]",
 "SA6 propiolactone/propiosultone":"[#6]1[#6][OX2][CX3]1=O",
 "SA7 epoxide":"[OX2r3]1[#6r3][#6r3]1",
 "SA7 aziridine":"[NX3r3]1[#6r3][#6r3]1",
 "SA8 aliphatic halide (reactive)":"[CX4;!$(C(F)F);!$(C[!#1;!#6])][Cl,Br,I]",
 "SA9 alkyl nitrite/nitrate":"[CX4][OX2][NX2]=O",
 "SA10 alpha,beta-unsat carbonyl (Michael acceptor)":"[CX3]=[CX3][CX3]=[OX1]",
 "SA11 aliphatic aldehyde":"[CX3H1](=O)[#6X4]",
 "SA12 aromatic nitroso":"[c][NX2]=O",
 "SA13 aromatic ring N-oxide":"[n+][O-]",
 "SA14 aromatic mono/dialkyl-nitro":"[c][NX3+](=O)[O-]",
 "SA15 nitro-aromatic":"[$([c][N+](=O)[O-]),$([c][N](=O)=O)]",
 "SA16 aromatic azo":"[c][NX2]=[NX2][c]",
 "SA17 aromatic azoxy":"[c][NX2]=[NX2+]([O-])[c]",
 "SA18 aromatic hydroxylamine":"[c][NX3;H1,H2][OX2H]",
 "SA19 aromatic amine (primary/secondary/tertiary)":"[c][NX3;!$([NX3][CX3]=[OX1]);!$([NX3+])]([#1,#6])[#1,#6]",
 "SA21 aromatic diazonium":"[c][NX2+]#[NX1]",
 "SA22 hydrazine":"[NX3;!$(N=O)][NX3;!$(N=O)]",
 "SA23 aliphatic N-nitroso":"[NX3][NX2]=O",
 "SA24 aromatic N-nitroso":"[c][NX3][NX2]=O",
 "SA25 epoxide/aziridine on N":"[NX3][OX2r3]1[#6][#6]1",
 "SA27 polycyclic aromatic (>=3 fused benzene)":"c1ccc2ccc3ccccc3c2c1",
 "SA28 quinone":"O=C1C=CC(=O)C=C1",
 "SA_alkyl_sulfonate":"[CX4][OX2][SX4](=O)(=O)[#6]",
 "SA_beta_lactam":"O=C1CCN1",
}
def rulebase(smi):
    m=Chem.MolFromSmiles(smi); hits=[]
    for name,sma in BB_ALERTS.items():
        q=Chem.MolFromSmarts(sma)
        if q is None:
            hits.append(name+" [BAD_SMARTS]"); continue
        if m.HasSubstructMatch(q): hits.append(name)
    return hits

# ---------- (2) Statistical Ames QSAR (Hansen N6512) ----------
gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
def fp(s):
    m=Chem.MolFromSmiles(s); return gen.GetFingerprint(m) if m else None
def arr(f):
    a=np.zeros((2048,),dtype=np.int8); DataStructs.ConvertToNumpyArray(f,a); return a

df=pd.read_csv(f"{OUT}/ames_data.csv")
smi_col="Canonical_Smiles"; y_col="Activity"
df=df.dropna(subset=[smi_col,y_col])
fps=[fp(s) for s in df[smi_col]]
keep=[i for i,f in enumerate(fps) if f is not None]
df=df.iloc[keep].reset_index(drop=True); fps=[fps[i] for i in keep]
X=np.vstack([arr(f) for f in fps]); y=df[y_col].astype(int).values
print(f"Ames training: {X.shape}, positives={y.sum()} ({y.mean()*100:.1f}%)")

skf=StratifiedKFold(n_splits=5,shuffle=True,random_state=42)
oof=np.zeros(len(y))
for tr,te in skf.split(X,y):
    clf=RandomForestClassifier(n_estimators=400,class_weight="balanced",n_jobs=-1,random_state=42)
    clf.fit(X[tr],y[tr]); oof[te]=clf.predict_proba(X[te])[:,1]
auc=roc_auc_score(y,oof)
# operating threshold: Youden J on OOF
from sklearn.metrics import roc_curve
fpr,tpr,thr=roc_curve(y,oof); J=tpr-fpr; thr_opt=float(thr[np.argmax(J)])
print(f"Ames QSAR 5-fold OOF AUC={auc:.3f}  Youden thr={thr_opt:.3f}")
final=RandomForestClassifier(n_estimators=400,class_weight="balanced",n_jobs=-1,random_state=42).fit(X,y)

results={}
for name,smi in MOLS.items():
    qf=fp(smi); qx=arr(qf).reshape(1,-1)
    prob=float(final.predict_proba(qx)[0,1])
    sims=np.array(DataStructs.BulkTanimotoSimilarity(qf,fps))
    top=np.sort(sims)[::-1][:5]
    alerts=rulebase(smi)
    call_stat = "positive" if prob>=thr_opt else "negative"
    results[name]={
      "expert_rulebase_alerts":alerts,
      "expert_call":"positive (alert present)" if alerts else "negative (no structural alert)",
      "stat_Ames_prob":round(prob,3),
      "stat_Ames_call":call_stat,
      "stat_AD_maxTanimoto":round(float(top[0]),3),
      "stat_AD_meanTop5":round(float(top.mean()),3),
    }
    print(f"\n{name}\n  expert: {results[name]['expert_call']}  alerts={alerts}"
          f"\n  statistical Ames: p={prob:.3f} -> {call_stat}  (AD maxTanimoto={top[0]:.3f})")

# ICH M7 class per compound (both negative + no alert => Class 5)
for name in MOLS:
    r=results[name]
    both_neg = (not r["expert_rulebase_alerts"]) and (r["stat_Ames_call"]=="negative")
    r["ICH_M7_class"]= "Class 5 (no structural alert; treat as non-mutagenic)" if both_neg else \
                       "Class 3/2 - expert review required (alert or positive prediction)"

out={"methodology_1":"Benigni-Bossa expert structural-alert rulebase (Toxtree/VEGA basis; %d alerts)"%len(BB_ALERTS),
     "methodology_2":{"model":"RandomForest(400,balanced) Morgan r2 2048 on Hansen Ames N6512",
                      "n_train":int(X.shape[0]),"pos_frac":round(float(y.mean()),3),
                      "cv_auc":round(float(auc),3),"threshold_youden":round(thr_opt,3)},
     "per_compound":results}
json.dump(out,open(f"{OUT}/m7_results.json","w"),indent=2)
print("\nsaved m7_results.json")
