"""Punkt 2 + 3: wskaznik dwugenowy Aldoc/Plcb4 oraz rozciecie wstegi na faldy.

Punkt 2 — kontrast miedzy dwoma genami zamiast progu na jednym, zeby wynik
          nie zalezal od jakosci wychwytu komorki.
Punkt 3 — sasiednie faldy sa BLISKO w przestrzeni, ale DALEKO wzdluz wstegi.
          Rozdzielamy je po orientacji: normalna warstwy Purkinjego wskazuje
          od warstwy ziarnistej ku drobinowej i po dwoch stronach szczeliny
          jest przeciwna.
"""
import numpy as np, h5py, glob, json
from sklearn.neighbors import NearestNeighbors
from sklearn.mixture import GaussianMixture
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
rng=np.random.default_rng(0)

def gene(h,name,feats,ip,ix,dt,n):
    if name not in feats: return None
    g=feats.index(name); pos=np.flatnonzero(ix==g)
    rows=np.searchsorted(ip,pos,side='right')-1
    v=np.zeros(n); v[rows]=dt[pos]; return v

RES={}
for f in sorted(glob.glob("/mnt/data1t/hao_stereoseq/Mouse*.h5ad")):
    nm=f.split('/')[-1].replace('.h5ad','')
    with h5py.File(f,'r') as h:
        feats=[x.decode() if isinstance(x,bytes) else str(x) for x in h['var']['features'][:]]
        ann=h['obs']['annotation']; ac=[c.decode() for c in ann['categories'][:]]; cd=ann['codes'][:]
        if 'purkinje layer' not in ac: continue
        sp=h['obsm']['spatial'][:]; umi=h['obs']['nCount_RNA'][:]; Pk=h['obs']['Purkinje'][:]
        ip=h['X']['indptr'][:]; ix=h['X']['indices'][:]; dt=h['X']['data'][:]
        n=len(ip)-1
        ald=gene(h,'Aldoc',feats,ip,ix,dt,n); plc=gene(h,'Plcb4',feats,ip,ix,dt,n)
    pl=cd==ac.index('purkinje layer')
    gr=cd==ac.index('granular layer'); mo=cd==ac.index('molecular layer')
    sel=pl&(Pk>np.percentile(Pk[pl],75))
    X=sp[sel]; u=umi[sel]
    A=ald[sel]/np.maximum(u,1)*1e4; B=plc[sel]/np.maximum(u,1)*1e4
    N=len(X)

    # ---- PUNKT 2: wskaznik dwugenowy ----
    idx2=np.log2((A+1)/(B+1))
    single=np.log1p(A)
    r_single=abs(np.corrcoef(np.log(u),single)[0,1])
    r_two   =abs(np.corrcoef(np.log(u),idx2)[0,1])

    # ---- PUNKT 3: normalne i rozciecie na faldy ----
    grX=sp[gr]; moX=sp[mo]
    kg=NearestNeighbors(n_neighbors=min(12,len(grX))).fit(grX)
    km=NearestNeighbors(n_neighbors=min(12,len(moX))).fit(moX)
    _,gi=kg.kneighbors(X); _,mi=km.kneighbors(X)
    nor=moX[mi].mean(1)-grX[gi].mean(1)
    nor/=np.maximum(np.linalg.norm(nor,axis=1,keepdims=True),1e-9)

    nn=NearestNeighbors(n_neighbors=7).fit(X); dd,ii=nn.kneighbors(X)
    r=np.repeat(np.arange(N),6); c=ii[:,1:].ravel(); w=dd[:,1:].ravel()
    dot=(nor[r]*nor[c]).sum(1)
    thr=np.percentile(w,90)
    keep=(dot>0.5)&(w<thr)
    G=coo_matrix((np.ones(keep.sum()),(r[keep],c[keep])),shape=(N,N))
    ncomp,cl=connected_components(G,directed=False)

    # ---- klasyfikacja na wskazniku dwugenowym ----
    gm=GaussianMixture(2,random_state=0).fit(idx2.reshape(-1,1))
    lab=gm.predict(idx2.reshape(-1,1)); mu=gm.means_.ravel()
    if mu[0]>mu[1]: lab=1-lab; mu=mu[::-1]

    # zgodnosc przestrzenna
    nb=ii[:,1:]
    same=(lab[nb]==lab[:,None]).mean()
    perm=np.array([((lp:=rng.permutation(lab))[nb]==lp[:,None]).mean() for _ in range(200)])
    z=(same-perm.mean())/perm.std()

    # ---- domeny wzdluz kazdego faldu ----
    runs=[]; segs=0; ordered=0
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
    RES[nm]=dict(n=int(N),ncomp=int(ncomp),segments=int(segs),ordered=int(ordered),
        r_single=float(r_single),r_two=float(r_two),
        frac_high=float(lab.mean()),conc=float(same),z=float(z),
        n_domains=int(len(runs)),run_med=float(np.median(runs)),
        run_p25=float(np.percentile(runs,25)),run_p75=float(np.percentile(runs,75)))
    print(f"=== {nm}  (n={N} komorek Purkinjego) ===")
    print(f"  PUNKT 2  korelacja |r| z glebokoscia (log UMI):")
    print(f"           sam Aldoc      : {r_single:.3f}")
    print(f"           Aldoc/Plcb4    : {r_two:.3f}   {'LEPIEJ' if r_two<r_single else 'GORZEJ'}"
          f"  (spadek o {100*(1-r_two/max(r_single,1e-9)):.0f}%)")
    print(f"  PUNKT 3  faldow (spojnych segmentow >=20 kom.): {segs} z {ncomp} skladowych")
    print(f"           uporzadkowano {ordered}/{N} komorek ({100*ordered/N:.0f}%)")
    print(f"  DOMENY   {len(runs)} domen, mediana {np.median(runs):.0f} komorek "
          f"(IQR {np.percentile(runs,25):.0f}-{np.percentile(runs,75):.0f})")
    print(f"           udzial Aldoc-wysokich {100*lab.mean():.1f}% | zgodnosc sasiadow "
          f"{100*same:.1f}% vs {100*perm.mean():.1f}% losowo, z={z:.1f}\n")
json.dump(RES,open("/mnt/data1t/pc_rebuild/ribbon.json","w"),indent=1)
print("zapisano ribbon.json")
