"""Cux2 i Ebf2 w SCA7 wobec WT — czy tracą wzajemnie wykluczający się charakter?"""
import rdata, numpy as np, pandas as pd, scipy.sparse as sp, psutil, sys
from rdata.conversion import _conversion as C
GiB=1024**3
_o=C.SimpleConverter._convert_next
def _p(self,obj):
    i=getattr(obj,'info',None)
    if i is not None and i.type==rdata.parser.RObjectType.CLO: return None
    return _o(self,obj)
C.SimpleConverter._convert_next=_p

def load(path):
    o=rdata.conversion.convert(rdata.parser.parse_file(path))
    md=o.__dict__["meta.data"]
    rna=o.__dict__["assays"]["RNA"]
    c=rna.__dict__["counts"] if hasattr(rna,'__dict__') else rna["counts"]
    d=c.__dict__ if hasattr(c,'__dict__') else c
    i=np.asarray(d["i"]); p=np.asarray(d["p"]); x=np.asarray(d["x"]); dim=np.asarray(d["Dim"])
    dn=d["Dimnames"]
    genes=np.array([str(g) for g in dn[0]]); cells=np.array([str(g) for g in dn[1]])
    M=sp.csc_matrix((x,i,p),shape=(int(dim[0]),int(dim[1])))   # geny x komorki
    return md,M,genes,cells

def analyse(tag,path):
    md,M,genes,cells=load(path)
    print(f"\n########## {tag} ##########")
    print(f"macierz {M.shape} (geny x komorki), nnz {M.nnz:,}, RAM {psutil.Process().memory_info().rss/GiB:.2f} GiB")
    gi={g:k for k,g in enumerate(genes)}
    umi=np.asarray(M.sum(0)).ravel()
    def vec(g):
        if g not in gi: return None
        return np.asarray(M[gi[g],:].todense()).ravel()
    # --- identyfikacja komorek Purkinjego po markerach ---
    PC=["Calb1","Car8","Pcp2","Ppp1r17","Itpr1","Slc1a6"]
    have=[g for g in PC if g in gi]
    S=np.zeros(M.shape[1])
    for g in have: S+=vec(g)/np.maximum(umi,1)*1e4
    cl=md["seurat_clusters"].astype(str).values
    sc_by=pd.Series(S).groupby(cl).median().sort_values(ascending=False)
    print(f"\nmediana sygnatury Purkinjego ({'+'.join(have)}) wg klastra:")
    for k,v in sc_by.items(): print(f"   klaster {k:>3}: {v:8.1f}   n={int((cl==k).sum())}")
    thr=sc_by.iloc[0]*0.35
    pk=[k for k,v in sc_by.items() if v>=thr]
    isPC=np.isin(cl,pk)
    print(f"\nklastry uznane za Purkinje: {pk}  ->  {isPC.sum()} komorek z {M.shape[1]}")
    typ=md["Type"].astype(str).values
    for g in ["Cux2","Ebf2","Aldoc","Plcb4","Calb1"]:
        v=vec(g)
        if v is None: print(f"  {g}: BRAK W MACIERZY"); continue
        cp=v/np.maximum(umi,1)*1e4
        print(f"\n=== {g} w komorkach Purkinjego ===")
        for t in ["WT","SCA7"]:
            m=isPC&(typ==t)
            print(f"   {t:<5} n={int(m.sum()):>5}  %dodatnich {100*np.mean(v[m]>0):>5.1f}%  "
                  f"mediana {np.median(cp[m]):>6.2f}  srednia {cp[m].mean():>6.2f} CP10K")
    # wzajemne wykluczanie
    cx=vec("Cux2"); eb=vec("Ebf2")
    if cx is not None and eb is not None:
        print(f"\n=== WZAJEMNE WYKLUCZANIE Cux2/Ebf2 (komorki Purkinjego) ===")
        for t in ["WT","SCA7"]:
            m=isPC&(typ==t)
            a=cx[m]>0; b=eb[m]>0
            both=np.mean(a&b); nei=np.mean(~a&~b)
            print(f"   {t:<5} tylko Cux2 {100*np.mean(a&~b):>5.1f}%  tylko Ebf2 {100*np.mean(b&~a):>5.1f}%  "
                  f"OBA {100*both:>5.1f}%  zadne {100*nei:>5.1f}%")
            from scipy.stats import fisher_exact
            tab=[[int(np.sum(a&b)),int(np.sum(a&~b))],[int(np.sum(~a&b)),int(np.sum(~a&~b))]]
            orr,pv=fisher_exact(tab)
            print(f"         iloraz szans (OR) = {orr:.3f}  p = {pv:.2e}   "
                  f"{'(wykluczanie)' if orr<1 else '(wspolwystepowanie)'}")
    return dict(tag=tag,n=int(isPC.sum()))

for tag,f in [("5 TYGODNI","GSM8315551_SCA7.5wk.seurat.rds.gz"),
              ("8 TYGODNI","GSM8315552_SCA7.8wk.seurat.rds.gz")]:
    try: analyse(tag,f)
    except Exception as e:
        import traceback; print(f"{tag} BLAD:"); traceback.print_exc()
