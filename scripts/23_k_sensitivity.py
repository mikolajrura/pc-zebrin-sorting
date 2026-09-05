"""Czy rozmiar domen to biologia, czy artefakt wygladzania?

Jesli domeny rosna proporcjonalnie do K (liczby laczonych sasiadow),
to je sam tworze. Jesli od pewnego K sie stabilizuja - sa prawdziwe.
"""
import numpy as np, h5py, glob, json
from sklearn.neighbors import NearestNeighbors
from sklearn.mixture import GaussianMixture
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
def gene(name,feats,ip,ix,dt,n):
    if name not in feats: return np.zeros(n)
    g=feats.index(name); pos=np.flatnonzero(ix==g)
    r=np.searchsorted(ip,pos,side='right')-1; v=np.zeros(n); v[r]=dt[pos]; return v

def run(f,K):
    with h5py.File(f,'r') as h:
        feats=[x.decode() if isinstance(x,bytes) else str(x) for x in h['var']['features'][:]]
        ann=h['obs']['annotation']; ac=[c.decode() for c in ann['categories'][:]]; cd=ann['codes'][:]
        sp=h['obsm']['spatial'][:]; umi=h['obs']['nCount_RNA'][:]; Pk=h['obs']['Purkinje'][:]
        ip=h['X']['indptr'][:]; ix=h['X']['indices'][:]; dt=h['X']['data'][:]; n=len(ip)-1
        A0=gene('Aldoc',feats,ip,ix,dt,n); B0=gene('Plcb4',feats,ip,ix,dt,n)
    pl=cd==ac.index('purkinje layer'); gr=cd==ac.index('granular layer'); mo=cd==ac.index('molecular layer')
    sel=pl&(Pk>np.percentile(Pk[pl],75))
    X=sp[sel]; u=umi[sel]; A=A0[sel]; B=B0[sel]; N=len(X)
    grX=sp[gr]; moX=sp[mo]
    _,gi=NearestNeighbors(n_neighbors=min(12,len(grX))).fit(grX).kneighbors(X)
    _,mi=NearestNeighbors(n_neighbors=min(12,len(moX))).fit(moX).kneighbors(X)
    nor=moX[mi].mean(1)-grX[gi].mean(1); nor/=np.maximum(np.linalg.norm(nor,axis=1,keepdims=True),1e-9)
    KG=max(K,7)
    dd,ii=NearestNeighbors(n_neighbors=KG+1).fit(X).kneighbors(X)
    r=np.repeat(np.arange(N),KG); c=ii[:,1:].ravel(); w=dd[:,1:].ravel()
    keep=((nor[r]*nor[c]).sum(1)>0.5)&(w<np.percentile(w,90))
    ncomp,cl=connected_components(coo_matrix((np.ones(keep.sum()),(r[keep],c[keep])),shape=(N,N)),directed=False)
    Ap=np.zeros(N);Bp=np.zeros(N);Up=np.zeros(N)
    for i in range(N):
        cand=ii[i][:K+1]; cand=cand[cl[cand]==cl[i]]
        Ap[i]=A[cand].sum();Bp[i]=B[cand].sum();Up[i]=u[cand].sum()
    idx2=np.log2((Ap/np.maximum(Up,1)*1e4+1)/(Bp/np.maximum(Up,1)*1e4+1))
    gm=GaussianMixture(2,random_state=0).fit(idx2.reshape(-1,1))
    lab=gm.predict(idx2.reshape(-1,1))
    if gm.means_.ravel()[0]>gm.means_.ravel()[1]: lab=1-lab
    runs=[]
    for k in range(ncomp):
        m=np.flatnonzero(cl==k)
        if len(m)<20: continue
        P0=X[m]-X[m].mean(0); _,_,vt=np.linalg.svd(P0,full_matrices=False)
        ls=lab[m][np.argsort(P0@vt[0])]; cur=1
        for t in range(1,len(ls)):
            if ls[t]==ls[t-1]: cur+=1
            else: runs.append(cur); cur=1
        runs.append(cur)
    runs=np.array(runs); big=runs[runs>=3]
    return dict(K=K,frac=float(lab.mean()),n_dom=len(runs),n_big=len(big),
                med_big=float(np.median(big)) if len(big) else 0,
                max_run=int(runs.max()) if len(runs) else 0)

print(f"{'przekroj':<16}{'K':>4}{'domen>=3':>10}{'mediana':>9}{'najdluzsza':>12}{'%wysokich':>11}")
agg={}
for f in sorted(glob.glob("/mnt/data1t/hao_stereoseq/Mouse*.h5ad")):
    nm=f.split('/')[-1].replace('.h5ad','')
    for K in [0,2,4,7,12,20]:
        r=run(f,K)
        agg.setdefault(K,[]).append(r['med_big'])
        print(f"{nm if K==0 else '':<16}{K:>4}{r['n_big']:>10}{r['med_big']:>9.0f}{r['max_run']:>12}{100*r['frac']:>10.1f}%")
    print()
print("=== mediana rozmiaru domeny usredniona po 4 przekrojach ===")
print(f"{'K (laczonych sasiadow)':<26}{'mediana domeny':>16}{'stosunek do K':>16}")
for K in [0,2,4,7,12,20]:
    m=np.mean(agg[K])
    print(f"{K:<26}{m:>16.1f}{m/max(K,1):>16.2f}")
print("\nJesli 'stosunek do K' spada wraz z K -> domeny NIE sa prostym artefaktem wygladzania.")
