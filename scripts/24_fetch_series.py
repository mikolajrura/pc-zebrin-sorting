"""Pobranie serii 61 sekcji myszy + ekstrakcja kompaktowych tabel.

Tabela A (wszystkie komorki): geometria -> do rejestracji szeregowej
Tabela B (tylko warstwa Purkinjego): panel genow -> do analizy domen
Duze .h5ad zostaja na dysku; do paczki ida tylko tabele.
"""
import numpy as np, pandas as pd, h5py, os, json, subprocess, time, sys
S="/tmp/claude-1000/-home-mikolajrurad/1c9f8792-8a0d-4e75-8a05-a9d370bd7ee3/scratchpad"
D="/mnt/data1t/hao_stereoseq"; OUT="/mnt/data1t/hao_pack"
BASE="https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000244"
os.makedirs(OUT,exist_ok=True); os.makedirs(f"{OUT}/sections",exist_ok=True)
P="/home/mikolajrurad/omics-data-highway/ichb-purkinje-analysis/pc-zebrin-sorting"
SIG=json.load(open(f"{P}/processed/subtype_signatures.json"))
PANEL=sorted({g for v in SIG.values() for g in v} |
             {'Aldoc','Plcb4','Car8','Calb1','Pcp2','Grid2','Gpr176','Tox2','Slc1a6','Itpr1'})
print(f"panel genow: {len(PANEL)}")
CT=['Astrocyte','Bergmann','Choroid','Endothelial_mural','Endothelial_stalk','Ependymal',
    'Fibroblast','Golgi','Granule','MLI1','MLI2','Macrophage','Microglia','ODC','OPC','PLI',
    'Purkinje','UBC']
files=[l.strip() for l in open(f"{S}/mouse_sections.txt") if l.strip()]
print(f"sekcji do przetworzenia: {len(files)}")
man=[]
t0=time.time()
for k,fn in enumerate(files,1):
    path=f"{D}/{fn}"; nm=fn.replace('.h5ad','')
    outA=f"{OUT}/sections/{nm}_geom.parquet"; outB=f"{OUT}/sections/{nm}_pc.parquet"
    if os.path.exists(outA) and os.path.exists(outB):
        print(f"[{k}/{len(files)}] {nm} juz gotowe"); continue
    if not os.path.exists(path):
        r=subprocess.run(["curl","-sSL","--retry","3","-o",path,f"{BASE}/{fn}"],
                         capture_output=True)
        if r.returncode!=0 or os.path.getsize(path)<1000:
            print(f"[{k}/{len(files)}] {nm} BLAD POBIERANIA"); continue
    try:
        with h5py.File(path,'r') as h:
            feats=[x.decode() if isinstance(x,bytes) else str(x) for x in h['var']['features'][:]]
            fi={g:i for i,g in enumerate(feats)}
            o=h['obs']; sp=h['obsm']['spatial'][:]
            n=sp.shape[0]
            ann=o['annotation']; ac=[c.decode() for c in ann['categories'][:]]; acode=ann['codes'][:]
            reg=o['region']; rc=[c.decode() for c in reg['categories'][:]]; rcode=reg['codes'][:]
            umi=o['nCount_RNA'][:]; nf=o['nFeature_RNA'][:]
            scores={c:o[c][:] for c in CT if c in o}
            ip=h['X']['indptr'][:]; ix=h['X']['indices'][:]; dt=h['X']['data'][:]
        A=pd.DataFrame({"section":nm,"x":sp[:,0].astype(np.float32),"y":sp[:,1].astype(np.float32),
                        "layer":pd.Categorical([ac[c] if c>=0 else "NA" for c in acode]),
                        "region":pd.Categorical([rc[c] if c>=0 else "NA" for c in rcode]),
                        "nCount":umi.astype(np.int32),"nFeature":nf.astype(np.int32),
                        "Purkinje":scores.get('Purkinje',np.zeros(n)).astype(np.float32)})
        A.to_parquet(outA,compression="zstd",index=False)
        pl=np.array([ac[c]=="purkinje layer" if c>=0 else False for c in acode])
        cols=[g for g in PANEL if g in fi]
        gidx=np.array([fi[g] for g in cols])
        M=np.zeros((n,len(cols)),dtype=np.uint16)
        m=np.isin(ix,gidx)
        pos=np.flatnonzero(m); rows=np.searchsorted(ip,pos,side='right')-1
        remap={int(gi):j for j,gi in enumerate(gidx)}
        M[rows,[remap[int(v)] for v in ix[pos]]]=np.minimum(dt[pos],65535).astype(np.uint16)
        B=pd.DataFrame(M[pl],columns=cols)
        B.insert(0,"section",nm); B.insert(1,"x",sp[pl,0].astype(np.float32))
        B.insert(2,"y",sp[pl,1].astype(np.float32))
        B.insert(3,"nCount",umi[pl].astype(np.int32))
        for c,v in scores.items(): B[f"sc_{c}"]=v[pl].astype(np.float32)
        B.to_parquet(outB,compression="zstd",index=False)
        man.append(dict(section=nm,series=nm.split('_')[0],
                        T=int(nm.split('_T')[1]),n_cells=int(n),n_pc=int(pl.sum()),
                        genes=len(feats),bytes=os.path.getsize(path)))
        print(f"[{k}/{len(files)}] {nm}  komorek {n:>6}  PC {int(pl.sum()):>6}  "
              f"({time.time()-t0:.0f}s)",flush=True)
    except Exception as e:
        print(f"[{k}/{len(files)}] {nm} BLAD: {type(e).__name__} {e}",flush=True)
json.dump(man,open(f"{OUT}/manifest.json","w"),indent=1)
print(f"\nGOTOWE: {len(man)} sekcji w {time.time()-t0:.0f} s")
