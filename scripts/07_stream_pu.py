"""Warstwa Purkinjego @10um wyciagnieta STRUMIENIOWO z NRRD.

Nigdy nie materializujemy pelnych 4.81 GiB - dekompresujemy porcjami
i zatrzymujemy tylko woksele o docelowych etykietach.
"""
import gzip, numpy as np, time, os, json, psutil
GiB=1024**3
F="/mnt/data1t/atlas/annotv3a_bbp_10.nrrd"
OUT="/mnt/data1t/pc_rebuild"
def mem(): return psutil.Process().memory_info().rss/GiB

# --- naglowek ---
hdr={}; off=0
with open(F,'rb') as fh:
    while True:
        line=fh.readline(); off+=len(line)
        s=line.decode('ascii','replace').strip()
        if s=='' : break
        if ':' in s and not s.startswith('#'):
            k,v=s.split(':',1); hdr[k.strip()]=v.strip()
print("naglowek:", {k:hdr[k] for k in ['type','sizes','encoding','endian'] if k in hdr})
print("offset danych:", off)
sizes=[int(x) for x in hdr['sizes'].split()]
NX,NY,NZ=sizes                      # kolejnosc NRRD: pierwsza os najszybsza
total=NX*NY*NZ
print(f"wymiary {sizes}, wokseli {total:,}, pelny rozmiar {total*4/GiB:.2f} GiB")

# ID warstwy Purkinjego (z hierarchii CCFv3a) + kontrola
PU={1145:'CBXpu',10673:'SIMpu',10676:'ANcr1pu',10679:'ANcr2pu',10682:'PRMpu',
    10685:'COPYpu',10688:'PFLpu',10691:'FLpu',10706:'LINGpu',10709:'CENT2pu',
    10712:'CENT3pu',10715:'CUL4pu',10718:'CUL5pu',10721:'CUL45pu',10724:'DECpu',
    10727:'FOTUpu',10730:'PYRpu',10733:'UVUpu',10736:'NODpu'}
CTRL={10723:'DECgr',10725:'DECmo'}          # kontrola: znane liczby przy 25um
targets=np.array(sorted(set(PU)|set(CTRL)),dtype=np.uint32)
print(f"szukam {len(targets)} etykiet ({len(PU)} PU + {len(CTRL)} kontrolnych)")

CH=40_000_000        # wokseli na porcje -> 160 MB uint32
idx_parts=[]; lab_parts=[]; counts={}
t0=time.time(); done=0
with gzip.open(F,'rb') as gz:
    gz.read(0)
    with open(F,'rb') as raw: pass
    # gzip.open na pliku NRRD: dane zaczynaja sie po naglowku -> otwieramy recznie
with open(F,'rb') as fh:
    fh.seek(off)
    dec=gzip.GzipFile(fileobj=fh, mode='rb')
    while done<total:
        n=min(CH,total-done)
        buf=dec.read(n*4)
        if not buf: break
        arr=np.frombuffer(buf,dtype='<u4')
        m=np.isin(arr,targets)
        if m.any():
            li=np.flatnonzero(m).astype(np.int64)+done
            idx_parts.append(li); lab_parts.append(arr[m].copy())
        for t in np.unique(arr[m]) if m.any() else []:
            counts[int(t)]=counts.get(int(t),0)+int((arr[m]==t).sum())
        done+=len(arr)
        if done % (CH*4) < CH:
            print(f"  {100*done/total:5.1f}%  RAM {mem():.2f} GiB  t={time.time()-t0:.0f}s", flush=True)
    dec.close()
print(f"przeczytano {done:,} z {total:,} wokseli w {time.time()-t0:.0f} s, RAM {mem():.2f} GiB")

lin=np.concatenate(idx_parts) if idx_parts else np.array([],np.int64)
lab=np.concatenate(lab_parts) if lab_parts else np.array([],np.uint32)
del idx_parts, lab_parts
print(f"trafionych wokseli: {len(lin):,}")

print("\n=== KONTROLA: znane struktury (10um vs 25um x 2.5^3=15.625) ===")
ref25={10723:122208,10725:91558}
for t,nm in CTRL.items():
    n=counts.get(t,0); exp=ref25[t]*15.625
    print(f"  {nm:>8} @10um={n:>9,}  @25um={ref25[t]:>8,}  oczekiwane~{exp:,.0f}  stosunek={n/max(ref25[t],1):.2f}")

print("\n=== WARSTWA PURKINJEGO @10um ===")
tot=0
for t in sorted(PU):
    n=counts.get(t,0); tot+=n
    print(f"  {'OK ' if n else 'BRAK'} id={t:>6} {PU[t]:<9} wokseli={n:>9,}")
print(f"\nlacznie (bez rodzica CBXpu): {tot-counts.get(1145,0):,} wokseli"
      f" = {(tot-counts.get(1145,0))*1e-9*1000:.4f} mm3")

# rozklad na wspolrzedne
i=(lin%NX).astype(np.int32); j=((lin//NX)%NY).astype(np.int32); k=(lin//(NX*NY)).astype(np.int32)
keep=np.isin(lab,np.array(sorted(PU),dtype=np.uint32))&(lab!=1145)
np.savez_compressed(f"{OUT}/pu_voxels_10um.npz",
    i=i[keep],j=j[keep],k=k[keep],label=lab[keep],sizes=np.array(sizes))
print(f"zapisano {int(keep.sum()):,} wokseli PU -> {OUT}/pu_voxels_10um.npz")
json.dump({str(k):v for k,v in counts.items()}, open(f"{OUT}/pu_counts.json","w"), indent=1)
