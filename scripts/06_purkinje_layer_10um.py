"""Faza 5b: warstwa Purkinjego z adnotacji 10 um (OME-Zarr, czytana leniwie)."""
import numpy as np, zarr, os, time, json
from brainglobe_atlasapi import BrainGlobeAtlas
GiB=1024**3; OUT="/mnt/data1t/pc_rebuild"
Z=os.path.expanduser("~/.brainglobe/brainglobe-atlasapi/annotation-sets/"
                     "ccfv3augmented_mouse-annotation/3_0/annotations_compressed.ome.zarr")
g=zarr.open(Z, mode='r'); s0=g['s0']; s1=g['s1']
print(f"s0 (10um): {s0.shape}  s1 (25um): {s1.shape}")

atlas=BrainGlobeAtlas("ccfv3augmented_mouse_25um", check_latest=False)
# wszystkie struktury bedace potomkami Cerebellum (512)
cb_ids=set()
for sid,s in atlas.structures.items():
    if 512 in s['structure_id_path']: cb_ids.add(s['id'])
print(f"struktur w drzewie Cerebellum: {len(cb_ids)}")

print("\nliczę bounding box mozdzku na 25 um ...", flush=True)
a25=np.asarray(s1[:])
mask=np.isin(a25, list(cb_ids))
print(f"wokseli mozdzku @25um: {mask.sum():,}")
idx=np.array(np.nonzero(mask))
bb25=[(int(idx[d].min()), int(idx[d].max())) for d in range(3)]
print("bbox @25um (z,y,x):", bb25)
del a25, mask, idx

# przeskalowanie na 10 um, z marginesem
sc=2.5
bb10=[]
for d,(lo,hi) in enumerate(bb25):
    L=max(0,int(np.floor(lo*sc))-3); H=min(s0.shape[d], int(np.ceil((hi+1)*sc))+3)
    bb10.append((L,H))
print("bbox @10um (z,y,x):", bb10)
shape=[h-l for l,h in bb10]
print(f"podobjetosc: {shape} = {np.prod(shape):,} wokseli = "
      f"{np.prod(shape)*4/GiB:.2f} GiB jako uint32")

t=time.time()
sub=np.asarray(s0[bb10[0][0]:bb10[0][1], bb10[1][0]:bb10[1][1], bb10[2][0]:bb10[2][1]])
print(f"wczytane w {time.time()-t:.0f} s, faktycznie {sub.nbytes/GiB:.2f} GiB")

uniq=set(np.unique(sub).tolist())
print(f"\nunikalnych etykiet w podobjetosci: {len(uniq)}")

# szukamy ID warstwy Purkinjego: wzorzec gr/pu/mo -> pu = gr_id+1
pu_expect={}
for sid,s in atlas.structures.items():
    if s['acronym'].endswith('gr') and 512 in s.get('structure_id_path',[]):
        pu_expect[s['id']+1]=s['acronym'][:-2]+'pu ('+s['name'].replace('granular layer','Purkinje layer')+')'
print(f"\n=== szukam {len(pu_expect)} spodziewanych ID warstwy Purkinjego (gr_id+1) ===")
found={}
for pid,nm in sorted(pu_expect.items()):
    n=int((sub==pid).sum())
    found[pid]=n
    flag="OK " if n>0 else "BRAK"
    print(f"  {flag} id={pid:>6} wokseli={n:>8}  {nm}")
tot=sum(found.values())
print(f"\nlacznie wokseli warstwy Purkinjego @10um: {tot:,}  ({tot*1e-3:.1f} x 10^3)")
print(f"objetosc: {tot*(10**3)/1e9:.4f} mm3")
if tot>0:
    np.savez_compressed(f"{OUT}/pu_10um.npz",
        bbox=np.array(bb10), ids=np.array(sorted([k for k,v in found.items() if v>0])),
        counts=np.array([found[k] for k in sorted([k for k,v in found.items() if v>0])]))
    # zapis samych wspolrzednych wokseli PU (rzadkie -> tanie)
    pu_ids=np.array([k for k,v in found.items() if v>0], dtype=np.uint32)
    m=np.isin(sub, pu_ids)
    zz,yy,xx=np.nonzero(m)
    lab=sub[zz,yy,xx]
    zz=zz.astype(np.int32)+bb10[0][0]; yy=yy.astype(np.int32)+bb10[1][0]; xx=xx.astype(np.int32)+bb10[2][0]
    np.savez_compressed(f"{OUT}/pu_voxels_10um.npz", z=zz,y=yy,x=xx,label=lab)
    print(f"zapisano {len(zz):,} wokseli warstwy Purkinjego -> pu_voxels_10um.npz")
