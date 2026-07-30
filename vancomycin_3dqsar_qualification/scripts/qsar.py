"""
Genuine 12-endpoint Tox21 QSAR (2D-QSAR) for vancomycin B + impurities.
- RandomForest (300 trees, class-balanced) on 2048-bit radius-2 Morgan fingerprints
- 5-fold stratified out-of-fold ROC-AUC per endpoint (honest performance)
- Predictions for parent, RRT 0.87, RRT 0.75
- Applicability-domain (AD): nearest-neighbour Tanimoto similarity to training set
Reproduces the molecular-toxicology skill's model spec, trained from public
MoleculeNet Tox21 (7,831 compounds).
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
RNG=42
gen=rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

def fp(smi):
    m=Chem.MolFromSmiles(smi)
    if m is None: return None
    return gen.GetFingerprint(m)
def arr(f):
    a=np.zeros((2048,),dtype=np.int8); DataStructs.ConvertToNumpyArray(f,a); return a

df=pd.read_csv(f"{OUT}/tox21.csv.gz")
ENDPOINTS=['NR-AR','NR-AR-LBD','NR-AhR','NR-Aromatase','NR-ER','NR-ER-LBD',
           'NR-PPAR-gamma','SR-ARE','SR-ATAD5','SR-HSE','SR-MMP','SR-p53']
DESC={
 "NR-AR":"Androgen receptor (full)","NR-AR-LBD":"Androgen receptor (LBD)",
 "NR-AhR":"Aryl hydrocarbon receptor","NR-Aromatase":"Aromatase (CYP19A1)",
 "NR-ER":"Estrogen receptor (full)","NR-ER-LBD":"Estrogen receptor (LBD)",
 "NR-PPAR-gamma":"PPAR-gamma","SR-ARE":"Oxidative-stress response (ARE)",
 "SR-ATAD5":"Genotoxicity / DNA damage (ATAD5)","SR-HSE":"Heat-shock response",
 "SR-MMP":"Mitochondrial membrane potential","SR-p53":"p53 DNA-damage response"}

# featurize training set once
df=df.dropna(subset=["smiles"]).reset_index(drop=True)
fps=[fp(s) for s in df["smiles"]]
ok=[i for i,f in enumerate(fps) if f is not None]
df=df.iloc[ok].reset_index(drop=True); fps=[fps[i] for i in ok]
X=np.vstack([arr(f) for f in fps])
train_fps=fps
print(f"Training matrix: {X.shape}")

# queries
Q={
 "Parent (Vancomycin B)":"OC1=CC([C@@H](C(O)=O)NC2=O)=C(C3=C(O)C=CC([C@@H](NC([C@H](NC([C@H](CC(N)=O)NC([C@@H]([C@H](O)C(C=C4Cl)=CC=C4O5)NC([C@H](NC)CC(C)C)=O)=O)=O)C6=CC5=C(O[C@H]7[C@H](O[C@]8([H])O[C@@H](C)[C@@H](O)[C@](N)(C)C8)[C@@H](O)[C@H](O)[C@@H](CO)O7)C9=C6)=O)C(N[C@H]2[C@H](O)C%10=CC(Cl)=C(O9)C=C%10)=O)=C3)C(O)=C1",
 "RRT 0.87 (7S,26R)-Vanco B":"OC1=CC([C@H](C(O)=O)NC2=O)=C(C3=C(O)C=CC([C@@H](NC([C@H](NC([C@H](CC(N)=O)NC([C@@H]([C@@H](O)C(C=C4Cl)=CC=C4O5)NC([C@H](NC)CC(C)C)=O)=O)=O)C6=CC5=C(O[C@H]7[C@H](O[C@]8([H])O[C@@H](C)[C@@H](O)[C@](N)(C)C8)[C@@H](O)[C@H](O)[C@@H](CO)O7)C9=C6)=O)C(N[C@H]2[C@H](O)C%10=CC(Cl)=C(O9)C=C%10)=O)=C3)C(O)=C1",
 "RRT 0.75 (26R)-DMe-DeLI":"OC1=CC([C@H](C(O)=O)NC2=O)=C(C3=C(O)C=CC([C@@H](NC([C@H](NC([C@H](CC(N)=O)NC([C@@H]([C@H](O)C(C=C4Cl)=CC=C4O5)NC([C@H](N)[C@H](CC)C)=O)=O)=O)C6=CC5=C(O[C@H]7[C@H](O[C@]8([H])O[C@@H](C)[C@@H](O)[C@](N)(C)C8)[C@@H](O)[C@H](O)[C@@H](CO)O7)C9=C6)=O)C(N[C@H]2[C@H](O)C%10=CC(Cl)=C(O9)C=C%10)=O)=C3)C(O)=C1",
}
qfps={k:fp(v) for k,v in Q.items()}
qX={k:arr(f).reshape(1,-1) for k,f in qfps.items()}

# ---- Applicability domain: NN Tanimoto to training set ----
ad={}
for k,qf in qfps.items():
    sims=np.array(DataStructs.BulkTanimotoSimilarity(qf, train_fps))
    top=np.sort(sims)[::-1][:5]
    ad[k]={"max_Tanimoto_to_train":round(float(top[0]),3),
           "mean_top5_Tanimoto":round(float(top.mean()),3)}
print("AD:",json.dumps(ad,indent=2))

# ---- Train + 5-fold OOF CV + predict ----
cv_rows=[]; preds={k:{} for k in Q}
for ep in ENDPOINTS:
    y=df[ep].values
    mask=~pd.isna(y)
    Xe=X[mask]; ye=y[mask].astype(int)
    # 5-fold OOF AUC
    skf=StratifiedKFold(n_splits=5,shuffle=True,random_state=RNG)
    oof=np.zeros(len(ye))
    for tr,te in skf.split(Xe,ye):
        clf=RandomForestClassifier(n_estimators=300,class_weight="balanced",
                                   n_jobs=-1,random_state=RNG)
        clf.fit(Xe[tr],ye[tr])
        oof[te]=clf.predict_proba(Xe[te])[:,1]
    auc=roc_auc_score(ye,oof)
    # final model on all labelled data
    final=RandomForestClassifier(n_estimators=300,class_weight="balanced",
                                 n_jobs=-1,random_state=RNG).fit(Xe,ye)
    for k in Q:
        preds[k][ep]=round(float(final.predict_proba(qX[k])[0,1]),3)
    cv_rows.append({"endpoint":ep,"description":DESC[ep],"n":int(mask.sum()),
                    "n_active":int(ye.sum()),"CV_AUC":round(auc,3)})
    print(f"{ep:14s} AUC={auc:.3f}  n={mask.sum():5d}  pred(parent/087/075)="
          f"{preds['Parent (Vancomycin B)'][ep]}/{preds['RRT 0.87 (7S,26R)-Vanco B'][ep]}/{preds['RRT 0.75 (26R)-DMe-DeLI'][ep]}")

cv=pd.DataFrame(cv_rows)
print(f"\nMean CV-AUC = {cv['CV_AUC'].mean():.3f}  (range {cv['CV_AUC'].min():.3f}-{cv['CV_AUC'].max():.3f})")

# deltas vs parent
P="Parent (Vancomycin B)"
delta={k:{ep:round(preds[k][ep]-preds[P][ep],3) for ep in ENDPOINTS} for k in Q}

out={"cv_metrics":cv_rows,"predictions":preds,"delta_vs_parent":delta,
     "applicability_domain":ad,
     "mean_cv_auc":round(float(cv['CV_AUC'].mean()),3),
     "model":"RandomForest(300,balanced) on Morgan r2 2048-bit; 5-fold stratified OOF AUC; "
             "trained on public MoleculeNet Tox21 (n=%d)"%X.shape[0]}
json.dump(out,open(f"{OUT}/qsar_results.json","w"),indent=2)
cv.to_csv(f"{OUT}/qsar_cv_metrics.csv",index=False)
pd.DataFrame(preds).to_csv(f"{OUT}/qsar_predictions.csv")
print("saved qsar_results.json / qsar_cv_metrics.csv / qsar_predictions.csv")
