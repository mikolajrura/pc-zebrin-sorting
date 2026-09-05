"""Domeny Aldoc w warstwie Purkinjego — analiza NIEZALEZNA OD SKALI.

Wszystko wyrazone w liczbie komorek albo we frakcji dlugosci wstegi,
bo fizycznej skali jednostek Hao nie udalo sie ustalic.
"""
import numpy as np, h5py, glob, json
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import NearestNeighbors
from scipy.sparse.csgraph import connected_components, minimum_spanning_tree
from scipy.sparse import coo_matrix
rng=np.random.default_rng(0)

def gene_col(h, name):
    feats=[x.decode() if isinstance(x,bytes) else str(x) for x in h['var']['features'][:]]
    if name not in feats: return None
    g=feats.index(name)
    X=h['X']; indptr=X['indptr'][:]; idx=X['indices'][:]; dat=X['data'][:]
    pos=np.flatnonzero(idx==g); rows=np.searchsorted(indptr,pos,side='right')-1
    v=np.zeros(len(indptr)-1); v[rows]=dat[pos]; return v

out={}
for f in sorted(glob.glob("/mnt/data1t/hao_stereoseq/Mouse*.h5ad")):
    nm=f.split('/')[-1].replace('.h5ad','')
    with h5py.File(f,'r') as h:
        ann=h['obs']['annotation']; ac=[c.decode() for c in ann['categories'][:]]; cd=ann['codes'][:]
        if 'purkinje layer' not in ac: continue
        pl=cd==ac.index('purkinje layer')
        sp=h['obsm']['spatial'][:][pl]
        umi=h['obs']['nCount_RNA'][:][pl]
        P=h['obs']['Purkinje'][:][pl]
        ald=gene_col(h,'Aldoc')[pl]
    sel=P>np.percentile(P,75)
    x,y=sp[sel,0],sp[sel,1]; a=ald[sel]/np.maximum(umi[sel],1)*1e4
    n=len(x)
    g=GaussianMixture(2,random_state=0).fit(np.log1p(a).reshape(-1,1))
    lab=g.predict(np.log1p(a).reshape(-1,1)); mu=np.expm1(g.means_.ravel())
    if mu[0]>mu[1]: lab=1-lab; mu=mu[::-1]
    bic1=GaussianMixture(1,random_state=0).fit(np.log1p(a).reshape(-1,1)).bic(np.log1p(a).reshape(-1,1))
    bic2=g.bic(np.log1p(a).reshape(-1,1))
    XY=np.c_[x,y]
    nn=NearestNeighbors(n_neighbors=9).fit(XY); _,ix=nn.kneighbors(XY); nb=ix[:,1:]
    same=(lab[nb]==lab[:,None]).mean()
    perm=np.array([( (lp:=rng.permutation(lab))[nb]==lp[:,None]).mean() for _ in range(200)])
    z=(same-perm.mean())/perm.std()
    # uporzadkowanie WZDLUZ wstegi: MST -> najdluzsza sciezka w kazdej skladowej
    d,i2=NearestNeighbors(n_neighbors=6).fit(XY).kneighbors(XY)
    r=np.repeat(np.arange(n),5); c=i2[:,1:].ravel(); w=d[:,1:].ravel()
    G=coo_matrix((w,(r,c)),shape=(n,n)); T=minimum_spanning_tree(G)
    ncomp,cl=connected_components(T,directed=False)
    runs=[]; switches=0; tot=0
    for k in range(ncomp):
        m=np.flatnonzero(cl==k)
        if len(m)<25: continue
        sub=XY[m]; li=lab[m]
        pc0=sub-sub.mean(0)
        u,s,vt=np.linalg.svd(pc0,full_matrices=False)
        o=np.argsort(pc0@vt[0])
        ls=li[o]; tot+=len(ls)
        sw=int((ls[1:]!=ls[:-1]).sum()); switches+=sw
        cur=1
        for t in range(1,len(ls)):
            if ls[t]==ls[t-1]: cur+=1
            else: runs.append(cur); cur=1
        runs.append(cur)
    runs=np.array(runs) if runs else np.array([0])
    out[nm]=dict(n=int(n), low=int((lab==0).sum()), high=int((lab==1).sum()),
        mu_low=float(mu[0]), mu_high=float(mu[1]), bic1=float(bic1), bic2=float(bic2),
        conc=float(same), z=float(z), ncomp=int(ncomp), switches=int(switches),
        run_med=float(np.median(runs)), run_mean=float(runs.mean()), n_ordered=int(tot),
        frac_high=float((lab==1).mean()))
    print(f"{nm:<16} n={n:>6} | Aldoc niskie {mu[0]:>6.1f} / wysokie {mu[1]:>6.1f} CP10K, "
          f"udzial wysokich {100*(lab==1).mean():>5.1f}%")
    print(f"{'':16} BIC 1skl={bic1:>8.0f} 2skl={bic2:>8.0f} {'(2 lepsze)' if bic2<bic1 else '(1 lepsze)'}"
          f" | zgodnosc sasiadow {100*same:>5.1f}% vs {100*perm.mean():.1f}% losowo, z={z:>6.1f}")
    print(f"{'':16} przelaczen wzdluz wstegi: {switches}, mediana dlugosci domeny {np.median(runs):.0f} komorek "
          f"({100*np.median(runs)/max(tot,1)*len(runs):.2f}% wstegi na domene)")
json.dump(out,open("/mnt/data1t/pc_rebuild/domains.json","w"),indent=1)
print("\nzapisano domains.json")
