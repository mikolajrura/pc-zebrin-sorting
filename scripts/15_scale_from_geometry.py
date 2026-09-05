"""Skala jednostki Hao wyznaczona z grubosci warstwy drobinowej.

Ta sama wielkosc anatomiczna mierzona po obu stronach:
 - atlas CCFv3a 10 um  -> grubosc w mikrometrach (skala znana)
 - sekcja Hao          -> grubosc w jednostkach Hao
Iloraz = um na jednostke Hao.
"""
import numpy as np, h5py
from scipy.spatial import cKDTree
rng=np.random.default_rng(0)
OUT="/mnt/data1t/pc_rebuild"

# ---------- ATLAS ----------
print("=== ATLAS CCFv3a 10 um ===")
d=np.load(f"{OUT}/pc_layer_10um.npz",allow_pickle=True)
names=[str(x) for x in d["names"]]
pci,pcj,pck,pclob=d["i"],d["j"],d["k"],d["lobule"]
# woksele warstwy drobinowej: odtwarzamy z NRRD strumieniowo dla wybranych placikow
import gzip
F="/mnt/data1t/atlas/annotv3a_bbp_10.nrrd"; off=381; NX,NY,NZ=1415,800,1140; total=NX*NY*NZ
LOB={'DEC':(10723,10725),'UVU':(10732,10734),'CENT3':(10711,10713),
     'ANcr1':(10675,10677),'PYR':(10729,10731)}
mo_ids=np.array([v[1] for v in LOB.values()],dtype=np.uint32)
lin=[];lab=[]
done=0;CH=40_000_000
with open(F,'rb') as fh:
    fh.seek(off); dec=gzip.GzipFile(fileobj=fh,mode='rb')
    while done<total:
        n=min(CH,total-done); buf=dec.read(n*4)
        if not buf: break
        a=np.frombuffer(buf,dtype='<u4'); m=np.isin(a,mo_ids)
        if m.any(): lin.append(np.flatnonzero(m).astype(np.int64)+done); lab.append(a[m].copy())
        done+=len(a)
    dec.close()
lin=np.concatenate(lin); lab=np.concatenate(lab)
mi=(lin%NX).astype(np.int32); mj=((lin//NX)%NY).astype(np.int32); mk=(lin//(NX*NY)).astype(np.int32)
print(f"wokseli warstwy drobinowej (5 placikow): {len(mi):,}")

atlas_th={}
for nm,(gid,moid) in LOB.items():
    li=names.index(nm)+1
    pc=np.c_[pci[pclob==li],pcj[pclob==li],pck[pclob==li]]
    sel=lab==moid
    mo=np.c_[mi[sel],mj[sel],mk[sel]]
    if len(pc)<100 or len(mo)<100: continue
    q=mo[rng.choice(len(mo),min(120000,len(mo)),replace=False)]
    dd,_=cKDTree(pc).query(q,k=1)
    th=np.percentile(dd,95)*10.0     # woksel = 10 um
    atlas_th[nm]=th
    print(f"  {nm:<7} n_mo={len(mo):>9,}  grubosc(p95 odl. do warstwy PC) = {th:>7.1f} um")

# ---------- HAO ----------
print("\n=== SEKCJE HAO ===")
import glob
hao_th={}
for f in sorted(glob.glob("/mnt/data1t/hao_stereoseq/Mouse*.h5ad")):
    nm=f.split('/')[-1].replace('.h5ad','')
    with h5py.File(f,'r') as h:
        sp=h['obsm']['spatial'][:]
        ann=h['obs']['annotation']; ac=[c.decode() for c in ann['categories'][:]]; cd=ann['codes'][:]
    if 'purkinje layer' not in ac or 'molecular layer' not in ac: continue
    pc=sp[cd==ac.index('purkinje layer')]; mo=sp[cd==ac.index('molecular layer')]
    dd,_=cKDTree(pc).query(mo,k=1)
    th=np.percentile(dd,95)
    hao_th[nm]=th
    print(f"  {nm:<16} n_pc={len(pc):>6,} n_mo={len(mo):>7,}  grubosc(p95) = {th:>6.2f} jednostek")

print("\n=== PRZELICZNIK: um na jednostke Hao ===")
a=np.mean(list(atlas_th.values())); print(f"  atlas, srednia z {len(atlas_th)} placikow: {a:.1f} um")
for nm,th in hao_th.items():
    print(f"  {nm:<16} {a/th:>6.2f} um/jednostke   (przy grubosci {th:.2f})")
print(f"\n  UWAGA: p95 to arbitralny wybor progu grubosci; ten sam prog po obu")
print(f"  stronach, wiec iloraz jest odporny na wybor, ale wartosc bezwzgledna nie.")
