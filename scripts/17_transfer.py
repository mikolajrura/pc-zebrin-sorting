"""Transfer 9 podtypow Kozareva -> Hao, z dopasowaniem glebokosci.

Uczciwosc polega na tym, ze klasyfikator walidujemy na danych Kozarevy
SZTUCZNIE SPLYCONYCH do glebokosci Stereo-seq. Dopiero jesli tam dziala,
wolno go puscic na Hao.
"""
import numpy as np, scanpy as sc, h5py, json, glob, psutil
import scipy.sparse as sp
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import balanced_accuracy_score, classification_report
GiB=1024**3
rng=np.random.default_rng(0)
P="/home/mikolajrurad/omics-data-highway/ichb-purkinje-analysis/pc-zebrin-sorting"
def mem(t): print(f"  [RAM {psutil.Process().memory_info().rss/GiB:.2f} GiB] {t}",flush=True)

A=sc.read_h5ad(f"{P}/processed/purkinje_cells_v2.h5ad"); mem("Kozareva wczytana")
X=A.layers["counts"].tocsr()
koz_genes=np.array(A.var_names)
y=A.obs["subcluster"].astype(str).values
koz_umi=np.asarray(X.sum(1)).ravel()
print(f"Kozareva: {X.shape}, mediana UMI {np.median(koz_umi):.0f}")

f0="/mnt/data1t/hao_stereoseq/Mouse1_T175.h5ad"
with h5py.File(f0,'r') as h:
    hao_genes=np.array([x.decode() if isinstance(x,bytes) else str(x) for x in h['var']['features'][:]])
    hao_umi=h['obs']['nCount_RNA'][:]
print(f"Hao: {len(hao_genes)} genow, mediana UMI {np.median(hao_umi):.0f}")
shared=np.intersect1d(koz_genes,hao_genes)
print(f"genow wspolnych: {len(shared):,}")

SIG=json.load(open(f"{P}/processed/subtype_signatures.json"))
sig=sorted({g for v in SIG.values() for g in v})
sig=[g for g in sig if g in set(shared)]
print(f"genow sygnaturowych obecnych u Hao: {len(sig)} (z {len({g for v in SIG.values() for g in v})})")

gi={g:i for i,g in enumerate(koz_genes)}
cols=np.array([gi[g] for g in sig])
Xs=X[:,cols].toarray().astype(np.float32)
print("macierz sygnaturowa:",Xs.shape)

TARGET=float(np.median(hao_umi))
print(f"\n=== splycenie Kozarevy do mediany Hao ({TARGET:.0f} UMI) ===")
def downsample(Xfull_counts, total_umi, target):
    """binomialne splycenie: kazda komorka zachowuje ulamek p=target/umi"""
    out=np.zeros_like(Xfull_counts)
    p=np.clip(target/np.maximum(total_umi,1),0,1)
    for i in range(Xfull_counts.shape[0]):
        out[i]=rng.binomial(Xfull_counts[i].astype(np.int64),p[i])
    return out.astype(np.float32)
Xd=downsample(Xs,koz_umi,TARGET)
print(f"  po splyceniu: mediana sumy po genach sygnaturowych "
      f"{np.median(Xd.sum(1)):.1f} (przed: {np.median(Xs.sum(1)):.1f})")
print(f"  mediana genow niezerowych: {np.median((Xd>0).sum(1)):.0f} (przed: {np.median((Xs>0).sum(1)):.0f})")

def norm(M,umi):
    M=M/np.maximum(umi,1)[:,None]*1e4
    return np.log1p(M)
tr,te=train_test_split(np.arange(len(y)),test_size=0.3,random_state=0,stratify=y)
def fit_eval(M,umi,tag):
    Z=norm(M,umi); s=StandardScaler().fit(Z[tr])
    clf=LogisticRegression(max_iter=4000,class_weight="balanced")
    clf.fit(s.transform(Z[tr]),y[tr])
    p=clf.predict(s.transform(Z[te]))
    acc=(p==y[te]).mean(); bal=balanced_accuracy_score(y[te],p)
    import pandas as pd
    at=np.where(pd.Series(y[te]).str.contains("Anti_Aldoc"),"Aldoc-","Aldoc+")
    ap=np.where(pd.Series(p).str.contains("Anti_Aldoc"),"Aldoc-","Aldoc+")
    print(f"  {tag:<34} 9 podtypow: acc {acc:.3f} bal {bal:.3f} | os Aldoc: acc {(at==ap).mean():.3f}")
    return clf,s,bal
print("\n=== WALIDACJA ===")
c_full,s_full,b_full=fit_eval(Xs,koz_umi,"pelna glebokosc")
c_dn,  s_dn,  b_dn  =fit_eval(Xd,np.full(len(y),TARGET),"splycona do glebokosci Hao")
np.save("/mnt/data1t/pc_rebuild/sig_cols.npy",cols)
json.dump({"sig":sig,"target_umi":TARGET,"bal_full":float(b_full),"bal_down":float(b_dn)},
          open("/mnt/data1t/pc_rebuild/transfer_meta.json","w"),indent=1)
import pickle
pickle.dump({"clf":c_dn,"scaler":s_dn,"genes":sig},open("/mnt/data1t/pc_rebuild/clf_down.pkl","wb"))
print("\nzapisano klasyfikator splycony -> clf_down.pkl")
