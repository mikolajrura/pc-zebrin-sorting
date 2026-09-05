"""Dopasowanie sekcji Hao do atlasu: plaszczyzna + pozycja + skala.

Deskryptory NIEZALEZNE OD SKALI (stosunek osi glownych, solidity, wypelnienie)
sluza do znalezienia plaszczyzny; dopiero potem stosunek pol daje skale.
"""
import numpy as np, h5py, glob, json
from scipy import ndimage as ndi
from brainglobe_atlasapi import BrainGlobeAtlas
rng=np.random.default_rng(0)

def descr(mask):
    """deskryptory ksztaltu niezalezne od skali"""
    if mask.sum()<50: return None
    yy,xx=np.nonzero(mask)
    pts=np.c_[yy,xx].astype(float); pts-=pts.mean(0)
    C=np.cov(pts.T); ev=np.sort(np.linalg.eigvalsh(C))[::-1]
    aspect=np.sqrt(ev[0]/max(ev[1],1e-9))
    area=mask.sum()
    fill=area/mask[yy.min():yy.max()+1, xx.min():xx.max()+1].size
    per=np.abs(np.diff(mask.astype(np.int8),axis=0)).sum()+np.abs(np.diff(mask.astype(np.int8),axis=1)).sum()
    circ=4*np.pi*area/max(per**2,1)
    lab,n=ndi.label(mask)
    return dict(aspect=float(aspect),fill=float(fill),circ=float(circ),
                ncomp=int(n),area=int(area))

print("=== ATLAS: maska mozdzku @25 um ===")
a=BrainGlobeAtlas("ccfv3augmented_mouse_25um",check_latest=False)
cb=set(s['id'] for s in a.structures.values() if 512 in s['structure_id_path'])
ann=a.annotation
M=np.isin(ann,list(cb))
print("wokseli mozdzku @25um:",M.sum())
nz=np.array(np.nonzero(M)); bb=[(int(nz[d].min()),int(nz[d].max())) for d in range(3)]
M=M[bb[0][0]:bb[0][1]+1, bb[1][0]:bb[1][1]+1, bb[2][0]:bb[2][1]+1]
print("bbox:",bb,"ksztalt:",M.shape,"=", tuple(round(s*0.025,2) for s in M.shape),"mm")
del ann

AX={0:"strzalkowa? (os 0 = przod-tyl)",1:"poprzeczna? (os 1 = gora-dol)",2:"czolowa? (os 2 = lewo-prawo)"}
atlas_slices={}
for ax in range(3):
    for idx in range(0,M.shape[ax],2):
        sl=np.take(M,idx,axis=ax)
        d=descr(sl)
        if d and d['area']>200: atlas_slices[(ax,idx)]=d
print("przekrojow atlasu do porownania:",len(atlas_slices))

print("\n=== SEKCJE HAO: deskryptory ===")
hao={}
for f in sorted(glob.glob("/mnt/data1t/hao_stereoseq/Mouse*.h5ad")):
    nm=f.split('/')[-1].replace('.h5ad','')
    with h5py.File(f,'r') as h: sp=h['obsm']['spatial'][:]
    x=sp[:,0]-sp[:,0].min(); y=sp[:,1]-sp[:,1].min()
    G=1.0
    W=int(x.max()/G)+2; H=int(y.max()/G)+2
    m=np.zeros((H,W),bool); m[(y/G).astype(int),(x/G).astype(int)]=True
    m=ndi.binary_closing(m,np.ones((3,3)))
    m=ndi.binary_fill_holes(m)
    d=descr(m); hao[nm]=(d,m)
    print(f"  {nm:<16} aspect={d['aspect']:.2f} fill={d['fill']:.2f} circ={d['circ']:.3f} "
          f"ncomp={d['ncomp']} area={d['area']}")

print("\n=== NAJLEPSZE DOPASOWANIE PLASZCZYZNY (deskryptory bezwymiarowe) ===")
res={}
for nm,(d,m) in hao.items():
    best=[]
    for (ax,idx),ad in atlas_slices.items():
        cost=(abs(np.log(d['aspect']/ad['aspect']))*1.0
              + abs(d['fill']-ad['fill'])*1.5
              + abs(np.log(max(d['circ'],1e-6)/max(ad['circ'],1e-6)))*0.7)
        best.append((cost,ax,idx,ad))
    best.sort(key=lambda t:t[0])
    c,ax,idx,ad=best[0]
    scale_um = np.sqrt(ad['area']/d['area'])*25.0
    top3=[(f"os{b[1]}@{b[2]}",round(b[0],3)) for b in best[:3]]
    res[nm]=dict(ax=ax,idx=idx,cost=float(c),scale_um=float(scale_um))
    print(f"  {nm:<16} os={ax} ({AX[ax].split('?')[0]}) plaster={idx}  koszt={c:.3f}")
    print(f"      pole atlas={ad['area']} vs Hao={d['area']}  ->  skala = {scale_um:.2f} um/jednostke")
    print(f"      3 najlepsze: {top3}")
sc=[v['scale_um'] for v in res.values()]
print(f"\nskale: {[round(s,2) for s in sc]}  mediana {np.median(sc):.2f}  rozrzut {max(sc)-min(sc):.2f} um")
json.dump(res,open("/mnt/data1t/pc_rebuild/plane_fit.json","w"),indent=1)
