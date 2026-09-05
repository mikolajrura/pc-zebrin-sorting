"""Faza 3 (BRAMKA): czy sekcje Hao dziela uklad wspolrzednych?"""
import h5py, numpy as np, glob, os
D="/mnt/data1t/hao_stereoseq"
files=sorted(glob.glob(f"{D}/Mouse*.h5ad"))
print(f"sekcji do sprawdzenia: {len(files)}\n")
rows=[]
for f in files:
    n=os.path.basename(f)
    with h5py.File(f,'r') as h:
        sp=h['obsm']['spatial'][:]
        ann=h['obs']['annotation']
        acat=[c.decode() for c in ann['categories'][:]]; acode=ann['codes'][:]
        pl=acode==acat.index('purkinje layer') if 'purkinje layer' in acat else np.zeros(len(acode),bool)
        reg=h['obs']['region']
        rcat=[c.decode() for c in reg['categories'][:]]
        X=h['X']; shp=X.attrs.get('shape')
    x,y=sp[:,0],sp[:,1]
    rows.append(dict(name=n, ncell=sp.shape[0], npu=int(pl.sum()),
        xmin=x.min(), xmax=x.max(), ymin=y.min(), ymax=y.max(),
        xc=x.mean(), yc=y.mean(), regs=len(rcat),
        genes=int(shp[1]) if shp is not None else -1))
print(f"{'sekcja':<22}{'komorek':>9}{'PC-layer':>9}{'genow':>7}{'regionow':>9}"
      f"{'x_min':>8}{'x_max':>8}{'y_min':>8}{'y_max':>8}{'x_srodek':>10}{'y_srodek':>10}")
for r in rows:
    print(f"{r['name']:<22}{r['ncell']:>9}{r['npu']:>9}{r['genes']:>7}{r['regs']:>9}"
          f"{r['xmin']:>8.0f}{r['xmax']:>8.0f}{r['ymin']:>8.0f}{r['ymax']:>8.0f}"
          f"{r['xc']:>10.1f}{r['yc']:>10.1f}")
print()
xs=np.array([r['xmax']-r['xmin'] for r in rows]); ys=np.array([r['ymax']-r['ymin'] for r in rows])
xc=np.array([r['xc'] for r in rows]); yc=np.array([r['yc'] for r in rows])
print("=== OCENA BRAMKI ===")
print(f"rozpietosc X: min {xs.min():.0f} max {xs.max():.0f} (rozrzut {xs.max()-xs.min():.0f})")
print(f"rozpietosc Y: min {ys.min():.0f} max {ys.max():.0f} (rozrzut {ys.max()-ys.min():.0f})")
print(f"srodek X: {xc.min():.0f}..{xc.max():.0f} (rozrzut {xc.max()-xc.min():.0f})")
print(f"srodek Y: {yc.min():.0f}..{yc.max():.0f} (rozrzut {yc.max()-yc.min():.0f})")
print(f"czy geny identyczne we wszystkich: {len(set(r['genes'] for r in rows))==1}")
