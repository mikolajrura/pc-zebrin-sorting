"""Faza 2: sygnatury 9 podtypow + walidacja na trzymanym zbiorze testowym."""
import scanpy as sc, numpy as np, pandas as pd, json, psutil
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix, classification_report, balanced_accuracy_score
GiB=1024**3
P="/home/mikolajrurad/omics-data-highway/ichb-purkinje-analysis/pc-zebrin-sorting"
A=sc.read_h5ad(f"{P}/processed/purkinje_cells_v2_processed.h5ad")
print(A.shape, f"RAM {psutil.Process().memory_info().rss/GiB:.2f} GiB")
A.obs["sub"]=A.obs["subcluster"].astype(str)

idx=np.arange(A.n_obs)
tr,te=train_test_split(idx,test_size=0.3,random_state=0,stratify=A.obs["sub"].values)
print(f"trening {len(tr)}  test {len(te)}")
Atr=A[tr].copy()

sc.tl.rank_genes_groups(Atr,"sub",method="wilcoxon",n_genes=200)
res=Atr.uns["rank_genes_groups"]
groups=list(res["names"].dtype.names)
SIG={}
print("\n=== sygnatury (top geny po filtrze logFC>0.5 i pct.1>0.25) ===")
for g in groups:
    d=sc.get.rank_genes_groups_df(Atr,group=g)
    d=d[(d.logfoldchanges>0.5)&(d.pvals_adj<0.01)]
    genes=d.head(40)["names"].tolist()
    SIG[g]=genes
    print(f"  {g:<24} n={len(genes):>3}  {', '.join(genes[:8])}")
json.dump(SIG,open(f"{P}/processed/subtype_signatures.json","w"),indent=1)

# scoring na CALYM zbiorze, klasyfikacja = argmax score
for g,genes in SIG.items():
    if genes: sc.tl.score_genes(A,genes,score_name=f"sc_{g}")
cols=[f"sc_{g}" for g in groups if SIG[g]]
S=A.obs[cols].to_numpy()
pred=np.array([groups[i] for i in S.argmax(1)])
truth=A.obs["sub"].values

print("\n=== WALIDACJA na trzymanym zbiorze testowym (30%) ===")
yt,yp=truth[te],pred[te]
print(f"dokladnosc: {(yt==yp).mean():.3f}   zbalansowana: {balanced_accuracy_score(yt,yp):.3f}")
print("\nmacierz pomylek (wiersz=prawda, kolumna=przewidziane):")
labs=sorted(set(truth))
cm=confusion_matrix(yt,yp,labels=labs)
short=[l.replace("Purkinje_","") for l in labs]
print("      "+"".join(f"{s[:9]:>10}" for s in short))
for r,l in zip(cm,short):
    print(f"{l[:9]:>9} "+"".join(f"{v:>10}" for v in r))
print()
print(classification_report(yt,yp,labels=labs,target_names=short,zero_division=0))

# uproszczenie do osi Aldoc
a_t=np.where(pd.Series(yt).str.contains("Anti_Aldoc"),"Aldoc-","Aldoc+")
a_p=np.where(pd.Series(yp).str.contains("Anti_Aldoc"),"Aldoc-","Aldoc+")
print(f"=== ta sama walidacja zredukowana do osi Aldoc+/- ===")
print(f"dokladnosc: {(a_t==a_p).mean():.3f}   zbalansowana: {balanced_accuracy_score(a_t,a_p):.3f}")
print(confusion_matrix(a_t,a_p,labels=["Aldoc+","Aldoc-"]))
A.write_h5ad(f"{P}/processed/purkinje_cells_v2_scored.h5ad",compression="gzip")
print(f"\nzapisano purkinje_cells_v2_scored.h5ad")
