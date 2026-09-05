"""Faza 2b: porzadny klasyfikator zamiast argmax ze score_genes."""
import scanpy as sc, numpy as np, pandas as pd, json
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import balanced_accuracy_score, confusion_matrix, classification_report
P="/home/mikolajrurad/omics-data-highway/ichb-purkinje-analysis/pc-zebrin-sorting"
A=sc.read_h5ad(f"{P}/processed/purkinje_cells_v2_processed.h5ad")
y=A.obs["subcluster"].astype(str).values
SIG=json.load(open(f"{P}/processed/subtype_signatures.json"))
genes=sorted({g for v in SIG.values() for g in v})
genes=[g for g in genes if g in A.var_names]
print(f"genow sygnaturowych (unia, obecnych): {len(genes)}")
X=np.asarray(A[:,genes].X.todense()) if hasattr(A[:,genes].X,'todense') else np.asarray(A[:,genes].X)
print("macierz cech:", X.shape)
tr,te=train_test_split(np.arange(len(y)),test_size=0.3,random_state=0,stratify=y)

def run(tag,Xf):
    sc_=StandardScaler().fit(Xf[tr])
    clf=LogisticRegression(max_iter=3000,C=1.0,class_weight="balanced",n_jobs=-1)
    clf.fit(sc_.transform(Xf[tr]),y[tr])
    p=clf.predict(sc_.transform(Xf[te]))
    acc=(p==y[te]).mean(); bal=balanced_accuracy_score(y[te],p)
    print(f"\n=== {tag} ===")
    print(f"dokladnosc {acc:.3f}  zbalansowana {bal:.3f}")
    labs=sorted(set(y)); short=[l.replace("Purkinje_","") for l in labs]
    print(classification_report(y[te],p,labels=labs,target_names=short,zero_division=0))
    ap=np.where(pd.Series(p).str.contains("Anti_Aldoc"),"Aldoc-","Aldoc+")
    at=np.where(pd.Series(y[te]).str.contains("Anti_Aldoc"),"Aldoc-","Aldoc+")
    print(f"po redukcji do osi Aldoc: dokladnosc {(ap==at).mean():.3f}  zbalansowana {balanced_accuracy_score(at,ap):.3f}")
    return clf,sc_,acc,bal

run("regresja logistyczna na genach sygnaturowych",X)
Xp=A.obsm["X_pca"][:,:30]
run("regresja logistyczna na 30 PC",Xp)
