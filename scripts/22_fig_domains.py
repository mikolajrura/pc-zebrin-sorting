import numpy as np, h5py, glob, json, matplotlib
PROJ="/home/mikolajrurad/omics-data-highway/ichb-purkinje-analysis/pc-zebrin-sorting"
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
from sklearn.mixture import GaussianMixture
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
rng=np.random.default_rng(0)
SURF="#fcfcfb"; INK="#0b0b0b"; INK2="#52514e"; HI="#2a78d6"; LO="#eb6834"
def gene(name,feats,ip,ix,dt,n):
    if name not in feats: return np.zeros(n)
    g=feats.index(name); pos=np.flatnonzero(ix==g)
    r=np.searchsorted(ip,pos,side='right')-1; v=np.zeros(n); v[r]=dt[pos]; return v
K=7; allruns=[]; panels=[]
for f in sorted(glob.glob("/mnt/data1t/hao_stereoseq/Mouse*.h5ad")):
    nm=f.split('/')[-1].replace('.h5ad','')
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
    dd,ii=NearestNeighbors(n_neighbors=K+1).fit(X).kneighbors(X)
    r=np.repeat(np.arange(N),K); c=ii[:,1:].ravel(); w=dd[:,1:].ravel()
    keep=((nor[r]*nor[c]).sum(1)>0.5)&(w<np.percentile(w,90))
    ncomp,cl=connected_components(coo_matrix((np.ones(keep.sum()),(r[keep],c[keep])),shape=(N,N)),directed=False)
    Ap=np.zeros(N);Bp=np.zeros(N);Up=np.zeros(N)
    for i in range(N):
        cand=ii[i]; cand=cand[cl[cand]==cl[i]]
        Ap[i]=A[cand].sum();Bp[i]=B[cand].sum();Up[i]=u[cand].sum()
    idx2=np.log2((Ap/np.maximum(Up,1)*1e4+1)/(Bp/np.maximum(Up,1)*1e4+1))
    gm=GaussianMixture(2,random_state=0).fit(idx2.reshape(-1,1))
    lab=gm.predict(idx2.reshape(-1,1)); mu=gm.means_.ravel()
    if mu[0]>mu[1]: lab=1-lab
    for k in range(ncomp):
        m=np.flatnonzero(cl==k)
        if len(m)<20: continue
        P0=X[m]-X[m].mean(0); _,_,vt=np.linalg.svd(P0,full_matrices=False)
        ls=lab[m][np.argsort(P0@vt[0])]; cur=1
        for t in range(1,len(ls)):
            if ls[t]==ls[t-1]: cur+=1
            else: allruns.append(cur); cur=1
        allruns.append(cur)
    panels.append((nm,sp,X,lab))
allruns=np.array(allruns); big=allruns[allruns>=3]
fig=plt.figure(figsize=(20,5.6),facecolor=SURF)
gs=fig.add_gridspec(1,5,width_ratios=[1,1,1,1,1.05],wspace=0.06)
for i,(nm,sp,X,lab) in enumerate(panels):
    ax=fig.add_subplot(gs[0,i]); ax.set_facecolor(SURF)
    ax.scatter(sp[:,0],sp[:,1],s=0.7,c="#ecebe6",linewidths=0,rasterized=True)
    for k,(cc,t) in enumerate([(LO,"Aldoc-niskie"),(HI,"Aldoc-wysokie")]):
        m=lab==k
        ax.scatter(X[m,0],X[m,1],s=8,c=cc,linewidths=0,rasterized=True,
                   label=f"{t} ({100*m.mean():.0f}%)" if i==0 else None)
    ax.set_title(nm.replace('_',' '),fontsize=11.5,color=INK,loc='left',pad=6)
    ax.set_aspect('equal'); ax.invert_yaxis(); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
    if i==0: ax.legend(loc='lower left',frameon=False,fontsize=9,markerscale=1.6,labelcolor=INK2)
ax=fig.add_subplot(gs[0,4]); ax.set_facecolor(SURF)
ax.hist(big,bins=np.arange(3,60,2),color=HI,edgecolor=SURF,linewidth=0.5)
ax.axvline(np.median(big),color=LO,lw=1.6)
ax.text(np.median(big)+1.5,ax.get_ylim()[1]*0.9,f"mediana {np.median(big):.0f}",color=LO,fontsize=10)
ax.set_title("Rozmiary domen ≥3 komórek",fontsize=11.5,color=INK,loc='left',pad=6)
ax.set_xlabel("komórek Purkinjego w domenie",color=INK2,fontsize=10)
ax.set_ylabel("liczba domen",color=INK2,fontsize=10)
ax.tick_params(colors=INK2,labelsize=9)
for s in ['top','right']: ax.spines[s].set_visible(False)
for s in ['left','bottom']: ax.spines[s].set_color("#d9d7d0")
ax.grid(axis='y',color="#eceae5",lw=0.8); ax.set_axisbelow(True)
fig.suptitle("Domeny Aldoc w warstwie Purkinjego — wskaźnik dwugenowy Aldoc/Plcb4, "
             f"łączony po ~8 sąsiadach w obrębie fałdu  ·  {len(big)} domen ≥3 kom. w 4 przekrojach",
             fontsize=11.5,color=INK2,x=0.008,ha='left',y=1.0)
fig.savefig(f"{PROJ}/figures/domains_pooled.png",dpi=145,facecolor=SURF,bbox_inches='tight')
print(f"domen >=3: {len(big)}, mediana {np.median(big):.0f}, p75 {np.percentile(big,75):.0f}, max {big.max()}")
print("zapisano figures/domains_pooled.png")
