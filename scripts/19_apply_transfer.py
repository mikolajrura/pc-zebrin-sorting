"""Nalozenie klasyfikatora (walidowanego na splyconych danych) na sekcje Hao."""
import numpy as np, h5py, glob, pickle, json
P="/home/mikolajrurad/omics-data-highway/ichb-purkinje-analysis/pc-zebrin-sorting"
M=pickle.load(open("/mnt/data1t/pc_rebuild/clf_down.pkl","rb"))
clf,scaler,genes=M["clf"],M["scaler"],M["genes"]
meta=json.load(open("/mnt/data1t/pc_rebuild/transfer_meta.json"))
print(f"klasyfikator: {len(genes)} genow, zbalansowana trafnosc na splyconych "
      f"danych Kozarevy = {meta['bal_down']:.3f}")
res={}
for f in sorted(glob.glob("/mnt/data1t/hao_stereoseq/Mouse*.h5ad")):
    nm=f.split('/')[-1].replace('.h5ad','')
    with h5py.File(f,'r') as h:
        feats=[x.decode() if isinstance(x,bytes) else str(x) for x in h['var']['features'][:]]
        fi={g:i for i,g in enumerate(feats)}
        cols=[fi.get(g,-1) for g in genes]
        ann=h['obs']['annotation']; ac=[c.decode() for c in ann['categories'][:]]; cd=ann['codes'][:]
        if 'purkinje layer' not in ac: continue
        pl=cd==ac.index('purkinje layer')
        Pk=h['obs']['Purkinje'][:]
        umi=h['obs']['nCount_RNA'][:]
        sp=h['obsm']['spatial'][:]
        X=h['X']; indptr=X['indptr'][:]; idx=X['indices'][:]; dat=X['data'][:]
        ncell=len(indptr)-1
        want={c:k for k,c in enumerate(cols) if c>=0}
        Mx=np.zeros((ncell,len(genes)),dtype=np.float32)
        mask=np.isin(idx,np.array(sorted(want)))
        pos=np.flatnonzero(mask)
        rows=np.searchsorted(indptr,pos,side='right')-1
        for p,r in zip(pos,rows):
            Mx[r,want[idx[p]]]=dat[p]
    sel=pl & (Pk>np.percentile(Pk[pl],75))
    Z=np.log1p(Mx[sel]/np.maximum(umi[sel],1)[:,None]*1e4)
    pred=clf.predict(scaler.transform(Z))
    proba=clf.predict_proba(scaler.transform(Z)).max(1)
    u,c=np.unique(pred,return_counts=True)
    frac={a:int(b) for a,b in zip(u,c)}
    ald=sum(v for k,v in frac.items() if 'Anti' not in k)
    print(f"\n=== {nm}  (n={int(sel.sum())} komorek Purkinjego) ===")
    print(f"  mediana pewnosci klasyfikacji: {np.median(proba):.3f}")
    print(f"  Aldoc+ {100*ald/sel.sum():.1f}%  |  Anti-Aldoc {100*(1-ald/sel.sum()):.1f}%")
    for k,v in sorted(frac.items(),key=lambda t:-t[1]):
        print(f"    {k.replace('Purkinje_',''):<16} {v:>6}  {100*v/sel.sum():>5.1f}%")
    np.savez_compressed(f"/mnt/data1t/pc_rebuild/transfer_{nm}.npz",
        x=sp[sel,0],y=sp[sel,1],pred=pred,proba=proba)
    res[nm]={"n":int(sel.sum()),"frac":frac,"med_proba":float(np.median(proba))}
json.dump(res,open("/mnt/data1t/pc_rebuild/transfer_results.json","w"),indent=1)
print("\nzapisano transfer_*.npz i transfer_results.json")
