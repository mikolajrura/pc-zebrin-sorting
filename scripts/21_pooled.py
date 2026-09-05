"""Domkniecie: wskaznik dwugenowy liczony na POLACZONYCH sasiadach WZDLUZ faldu.

Laczenie sasiadow rozwiazuje oba problemy naraz:
 - Plcb4 przestaje byc rzadki (5 sasiadow = 104% glebokosci Kozarevy)
 - etykiety przestaja migotac, wiec dlugosci domen staja sie mierzalne
Laczymy TYLKO w obrebie tego samego faldu, zeby nie mieszac przez szczeline.
"""
import numpy as np, h5py, glob, json
from sklearn.neighbors import NearestNeighbors
from sklearn.mixture import GaussianMixture
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
rng=np.random.default_rng(0)
def gene(name,feats,ip,ix,dt,n):
    if name not in feats: return np.zeros(n)
    g=feats.index(name); pos=np.flatnonzero(ix==g)
    r=np.searchsorted(ip,pos,side='right')-1
    v=np.zeros(n); v[r]=dt[pos]; return v

K=7   # sasiadow do polaczenia (z obliczen: 5 = 104% glebokosci Kozarevy)
RES={}
for f in sorted(glob.glob("/mnt/data1t/hao_stereoseq/Mouse*.h5ad")):
    nm=f.split('/')[-1].replace('.h5ad','')
    with h5py.File(f,'r') as h:
        feats=[x.decode() if isinstance(x,bytes) else str(x) for x in h['var']['features'][:]]
        ann=h['obs']['annotation']; ac=[c.decode() for c in ann['categories'][:]]; cd=ann['codes'][:]
        if 'purkinje layer' not in ac: continue
        sp=h['obsm']['spatial'][:]; umi=h['obs']['nCount_RNA'][:]; Pk=h['obs']['Purkinje'][:]
        ip=h['X']['indptr'][:]; ix=h['X']['indices'][:]; dt=h['X']['data'][:]; n=len(ip)-1
        ald=gene('Aldoc',feats,ip,ix,dt,n); plc=gene('Plcb4',feats,ip,ix,dt,n)
    pl=cd==ac.index('purkinje layer'); gr=cd==ac.index('granular layer'); mo=cd==ac.index('molecular layer')
    sel=pl&(Pk>np.percentile(Pk[pl],75))
    X=sp[sel]; u=umi[sel]; A=ald[sel]; B=plc[sel]; N=len(X)

    grX=sp[gr]; moX=sp[mo]
    _,gi=NearestNeighbors(n_neighbors=min(12,len(grX))).fit(grX).kneighbors(X)
    _,mi=NearestNeighbors(n_neighbors=min(12,len(moX))).fit(moX).kneighbors(X)
    nor=moX[mi].mean(1)-grX[gi].mean(1); nor/=np.maximum(np.linalg.norm(nor,axis=1,keepdims=True),1e-9)
    dd,ii=NearestNeighbors(n_neighbors=K+1).fit(X).kneighbors(X)
    r=np.repeat(np.arange(N),K); c=ii[:,1:].ravel(); w=dd[:,1:].ravel()
    dot=(nor[r]*nor[c]).sum(1); thr=np.percentile(w,90)
    keep=(dot>0.5)&(w<thr)
    G=coo_matrix((np.ones(keep.sum()),(r[keep],c[keep])),shape=(N,N))
    ncomp,cl=connected_components(G,directed=False)

    # LACZENIE tylko w obrebie tego samego faldu
    Ap=np.zeros(N); Bp=np.zeros(N); Up=np.zeros(N); npool=np.zeros(N,int)
    for i in range(N):
        cand=ii[i]; cand=cand[cl[cand]==cl[i]]
        Ap[i]=A[cand].sum(); Bp[i]=B[cand].sum(); Up[i]=u[cand].sum(); npool[i]=len(cand)
    Acp=Ap/np.maximum(Up,1)*1e4; Bcp=Bp/np.maximum(Up,1)*1e4
    idx2=np.log2((Acp+1)/(Bcp+1))
    det_b_single=100*np.mean(B>0); det_b_pool=100*np.mean(Bp>0)
    r_dep=abs(np.corrcoef(np.log(np.maximum(Up,1)),idx2)[0,1])

    gm=GaussianMixture(2,random_state=0).fit(idx2.reshape(-1,1))
    lab=gm.predict(idx2.reshape(-1,1)); mu=gm.means_.ravel()
    if mu[0]>mu[1]: lab=1-lab; mu=mu[::-1]
    nb=ii[:,1:]
    same=(lab[nb]==lab[:,None]).mean()
    perm=np.array([((lp:=rng.permutation(lab))[nb]==lp[:,None]).mean() for _ in range(200)])
    z=(same-perm.mean())/perm.std()

    runs=[];segs=0;ordered=0
    for k in range(ncomp):
        m=np.flatnonzero(cl==k)
        if len(m)<20: continue
        segs+=1; P0=X[m]-X[m].mean(0)
        _,_,vt=np.linalg.svd(P0,full_matrices=False)
        o=np.argsort(P0@vt[0]); ls=lab[m][o]; ordered+=len(ls)
        cur=1
        for t in range(1,len(ls)):
            if ls[t]==ls[t-1]: cur+=1
            else: runs.append(cur); cur=1
        runs.append(cur)
    runs=np.array(runs) if len(runs) else np.array([0])
    big=runs[runs>=3]
    RES[nm]=dict(n=int(N),K=K,med_pool=float(np.median(npool)),
        det_plcb4_single=det_b_single,det_plcb4_pool=det_b_pool,r_depth=float(r_dep),
        frac_high=float(lab.mean()),conc=float(same),z=float(z),segments=int(segs),
        n_domains=int(len(runs)),run_med=float(np.median(runs)),
        n_dom_ge3=int(len(big)),run_med_ge3=float(np.median(big)) if len(big) else 0.0)
    print(f"=== {nm} (n={N}, laczonych srednio {np.median(npool):.0f} komorek) ===")
    print(f"  Plcb4 wykryty: pojedyncza kom. {det_b_single:.1f}%  ->  po polaczeniu {det_b_pool:.1f}%")
    print(f"  zaleznosc od glebokosci |r| = {r_dep:.3f}")
    print(f"  zgodnosc sasiadow {100*same:.1f}% vs {100*perm.mean():.1f}% losowo, z={z:.1f}")
    print(f"  faldow: {segs} | domen: {len(runs)}, mediana {np.median(runs):.0f} kom.")
    print(f"  domeny >=3 komorki: {len(big)}, mediana {np.median(big) if len(big) else 0:.0f} kom.,"
          f" najdluzsza {runs.max()}")
    print(f"  udzial Aldoc-wysokich: {100*lab.mean():.1f}%\n")
json.dump(RES,open("/mnt/data1t/pc_rebuild/pooled.json","w"),indent=1)
print("zapisano pooled.json")
