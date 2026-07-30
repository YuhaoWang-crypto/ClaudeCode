"""Augment QSAR: (1) identify the training compound identical (Tanimoto=1.0) to the
parent and report its MEASURED Tox21 labels (ground truth for vancomycin);
(2) honest held-out prediction for the parent by refitting each endpoint with all
Tanimoto=1.0 training rows removed."""
import json, numpy as np, pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import RandomForestClassifier
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
OUT="/tmp/claude-0/-home-user-ClaudeCode/0e1c2284-3690-5d68-a96a-480ebe908c1c/scratchpad"
gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
def fp(s):
    m=Chem.MolFromSmiles(s); return gen.GetFingerprint(m) if m else None
def arr(f):
    a=np.zeros((2048,),dtype=np.int8); DataStructs.ConvertToNumpyArray(f,a); return a
ENDPOINTS=['NR-AR','NR-AR-LBD','NR-AhR','NR-Aromatase','NR-ER','NR-ER-LBD',
           'NR-PPAR-gamma','SR-ARE','SR-ATAD5','SR-HSE','SR-MMP','SR-p53']
PARENT="OC1=CC([C@@H](C(O)=O)NC2=O)=C(C3=C(O)C=CC([C@@H](NC([C@H](NC([C@H](CC(N)=O)NC([C@@H]([C@H](O)C(C=C4Cl)=CC=C4O5)NC([C@H](NC)CC(C)C)=O)=O)=O)C6=CC5=C(O[C@H]7[C@H](O[C@]8([H])O[C@@H](C)[C@@H](O)[C@](N)(C)C8)[C@@H](O)[C@H](O)[C@@H](CO)O7)C9=C6)=O)C(N[C@H]2[C@H](O)C%10=CC(Cl)=C(O9)C=C%10)=O)=C3)C(O)=C1"
RRT075="OC1=CC([C@H](C(O)=O)NC2=O)=C(C3=C(O)C=CC([C@@H](NC([C@H](NC([C@H](CC(N)=O)NC([C@@H]([C@H](O)C(C=C4Cl)=CC=C4O5)NC([C@H](N)[C@H](CC)C)=O)=O)=O)C6=CC5=C(O[C@H]7[C@H](O[C@]8([H])O[C@@H](C)[C@@H](O)[C@](N)(C)C8)[C@@H](O)[C@H](O)[C@@H](CO)O7)C9=C6)=O)C(N[C@H]2[C@H](O)C%10=CC(Cl)=C(O9)C=C%10)=O)=C3)C(O)=C1"

df=pd.read_csv(f"{OUT}/tox21.csv.gz").dropna(subset=["smiles"]).reset_index(drop=True)
fps=[fp(s) for s in df["smiles"]]
keep=[i for i,f in enumerate(fps) if f is not None]
df=df.iloc[keep].reset_index(drop=True); fps=[fps[i] for i in keep]
pf=fp(PARENT)
sims=np.array(DataStructs.BulkTanimotoSimilarity(pf,fps))
identical=np.where(sims>=0.999)[0]
print("Training rows identical to parent (Tanimoto>=0.999):")
measured={}
for idx in identical:
    row=df.iloc[idx]
    print("  mol_id=",row["mol_id"],"smiles_len",len(row["smiles"]))
    for ep in ENDPOINTS:
        v=row[ep]
        measured.setdefault(ep,[]).append(None if pd.isna(v) else int(v))
# consolidate measured label per endpoint (take max coverage / any active)
measured_final={}
for ep in ENDPOINTS:
    vals=[v for v in measured.get(ep,[]) if v is not None]
    measured_final[ep]= ("active" if any(v==1 for v in vals) else "inactive") if vals else "not tested"
print("\nMEASURED Tox21 profile for vancomycin (parent):")
for ep in ENDPOINTS: print(f"  {ep:14s}: {measured_final[ep]}")

# Honest held-out parent prediction: remove identical rows, refit, predict parent & 075
X=np.vstack([arr(f) for f in fps])
q_parent=arr(pf).reshape(1,-1); q075=arr(fp(RRT075)).reshape(1,-1)
mask_out=np.ones(len(df),bool); mask_out[identical]=False
heldout={}; pred075_clean={}
for ep in ENDPOINTS:
    y=df[ep].values; m=~pd.isna(y)
    m_train=m & mask_out
    clf=RandomForestClassifier(n_estimators=300,class_weight="balanced",n_jobs=-1,random_state=42)
    clf.fit(X[m_train], y[m_train].astype(int))
    heldout[ep]=round(float(clf.predict_proba(q_parent)[0,1]),3)
    pred075_clean[ep]=round(float(clf.predict_proba(q075)[0,1]),3)
print("\nParent held-out (vancomycin removed from training) predictions:")
for ep in ENDPOINTS: print(f"  {ep:14s}: heldout_parent={heldout[ep]}  pred_075={pred075_clean[ep]}  delta={round(pred075_clean[ep]-heldout[ep],3)}")

json.dump({"measured_vancomycin":measured_final,
           "heldout_parent_pred":heldout,
           "heldout_pred_rrt075":pred075_clean,
           "n_identical_train_rows":int(len(identical)),
           "identical_mol_ids":[str(df.iloc[i]["mol_id"]) for i in identical]},
          open(f"{OUT}/qsar_augment.json","w"),indent=2)
print("\nsaved qsar_augment.json")
