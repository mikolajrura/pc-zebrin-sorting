"""Warstwa Purkinjego @10um jako granica gr|mo.

Atlas nie zawiera etykiet *pu (0 wokseli), ale zawiera *gr i *mo.
Ciala komorek Purkinjego leza na granicy tych warstw - liczymy ja morfologicznie.
"""
import gzip, numpy as np, time, json, psutil
from scipy import ndimage as ndi
GiB=1024**3
def mem(): return psutil.Process().memory_info().rss/GiB
F="/mnt/data1t/atlas/annotv3a_bbp_10.nrrd"; OUT="/mnt/data1t/pc_rebuild"
off=381; NX,NY,NZ=1415,800,1140; total=NX*NY*NZ

LOB={'LING':(10705,10707),'CENT2':(10708,10710),'CENT3':(10711,10713),
     'CUL4, 5':(10720,10722),'DEC':(10723,10725),'FOTU':(10726,10728),
     'PYR':(10729,10731),'UVU':(10732,10734),'NOD':(10735,10737),
     'SIM':(10672,10674),'ANcr1':(10675,10677),'ANcr2':(10678,10680),
     'PRM':(10681,10683),'COPY':(10684,10686),'PFL':(10687,10689),'FL':(10690,10692)}
names=sorted(LOB)
gr2lob={LOB[n][0]:i+1 for i,n in enumerate(names)}
mo2lob={LOB[n][1]:i+1 for i,n in enumerate(names)}
targets=np.array(sorted(list(gr2lob)+list(mo2lob)),dtype=np.uint32)
print(f"placikow: {len(names)}  szukanych etykiet: {len(targets)}")

t0=time.time(); done=0; CH=40_000_000
lin_p=[]; lab_p=[]
with open(F,'rb') as fh:
    fh.seek(off); dec=gzip.GzipFile(fileobj=fh,mode='rb')
    while done<total:
        n=min(CH,total-done); buf=dec.read(n*4)
        if not buf: break
        a=np.frombuffer(buf,dtype='<u4')
        m=np.isin(a,targets)
        if m.any():
            lin_p.append(np.flatnonzero(m).astype(np.int64)+done)
            lab_p.append(a[m].copy())
        done+=len(a)
    dec.close()
lin=np.concatenate(lin_p); lab=np.concatenate(lab_p); del lin_p,lab_p
print(f"przelot {done:,} w {time.time()-t0:.0f} s | wokseli kory mozdzku: {len(lin):,} | RAM {mem():.2f} GiB")

i=(lin%NX).astype(np.int32); j=((lin//NX)%NY).astype(np.int32); k=(lin//(NX*NY)).astype(np.int32)
del lin
bb=[(int(v.min()),int(v.max())) for v in (i,j,k)]
shape=tuple(h-l+1 for l,h in bb)
print(f"bbox (i,j,k): {bb}  ksztalt {shape} = {np.prod(shape):,} wokseli")

lobv=np.zeros(shape,dtype=np.uint8); layv=np.zeros(shape,dtype=np.uint8)
ii=i-bb[0][0]; jj=j-bb[1][0]; kk=k-bb[2][0]
isgr=np.isin(lab,np.array(sorted(gr2lob),dtype=np.uint32))
lobcode=np.zeros(len(lab),dtype=np.uint8)
for lid,c in gr2lob.items(): lobcode[lab==lid]=c
for lid,c in mo2lob.items(): lobcode[lab==lid]=c
lobv[ii,jj,kk]=lobcode
layv[ii,jj,kk]=np.where(isgr,1,2).astype(np.uint8)
print(f"gesta objetosc zbudowana, RAM {mem():.2f} GiB "
      f"({(lobv.nbytes+layv.nbytes)/GiB:.2f} GiB same tablice)")
del i,j,k,ii,jj,kk,lab,lobcode,isgr

st=ndi.generate_binary_structure(3,1)   # 6-sasiedztwo
mo=(layv==2); gr=(layv==1)
mo_d=ndi.binary_dilation(mo,structure=st)
pc = gr & mo_d
print(f"\nwokseli granicy gr|mo (warstwa Purkinjego): {int(pc.sum()):,}  RAM {mem():.2f} GiB")
del mo,mo_d

zz,yy,xx=np.nonzero(pc)
lobs=lobv[zz,yy,xx]
print(f"\n=== wokseli warstwy Purkinjego wg placika ===")
print(f"{'placik':<10}{'wokseli':>10}{'gr wokseli':>12}{'udzial %':>10}")
res={}
for idx,nm in enumerate(names,1):
    n=int((lobs==idx).sum()); g=int(((layv==1)&(lobv==idx)).sum())
    res[nm]=n
    print(f"{nm:<10}{n:>10,}{g:>12,}{100*n/max(g,1):>9.1f}%")
tot=sum(res.values())
print(f"\nRAZEM {tot:,} wokseli = {tot*1e-9*1000:.4f} mm3")
print(f"grubosc srednia: warstwa ma {tot/ (sum(1 for _ in names)):,.0f} wokseli/placik")

np.savez_compressed(f"{OUT}/pc_layer_10um.npz",
    i=(zz+bb[0][0]).astype(np.int32), j=(yy+bb[1][0]).astype(np.int32),
    k=(xx+bb[2][0]).astype(np.int32), lobule=lobs,
    names=np.array(names), bbox=np.array(bb), sizes=np.array([NX,NY,NZ]))
print(f"\nzapisano -> {OUT}/pc_layer_10um.npz")
json.dump(res, open(f"{OUT}/pc_layer_counts.json","w"), indent=1)
