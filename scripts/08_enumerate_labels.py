"""Jakie etykiety FAKTYCZNIE wystepuja w annotv3a_bbp_10.nrrd - pelny przelot."""
import gzip, numpy as np, time, json, psutil
GiB=1024**3
F="/mnt/data1t/atlas/annotv3a_bbp_10.nrrd"
off=381; NX,NY,NZ=1415,800,1140; total=NX*NY*NZ
acc=None
t0=time.time(); done=0; CH=40_000_000
with open(F,'rb') as fh:
    fh.seek(off); dec=gzip.GzipFile(fileobj=fh,mode='rb')
    while done<total:
        n=min(CH,total-done); buf=dec.read(n*4)
        if not buf: break
        a=np.frombuffer(buf,dtype='<u4')
        u,c=np.unique(a,return_counts=True)
        if acc is None: acc=dict(zip(u.tolist(),c.tolist()))
        else:
            for k,v in zip(u.tolist(),c.tolist()): acc[k]=acc.get(k,0)+v
        done+=len(a)
    dec.close()
print(f"przelot {done:,} wokseli w {time.time()-t0:.0f} s, RAM {psutil.Process().memory_info().rss/GiB:.2f} GiB")
print(f"unikalnych etykiet: {len(acc)}")
json.dump({str(k):v for k,v in acc.items()}, open("/mnt/data1t/pc_rebuild/labels_10um.json","w"))

# krzyzowo z hierarchia z Zenodo
import urllib.request, os
H="/mnt/data1t/atlas/hierarchy.json"
if not os.path.exists(H):
    urllib.request.urlretrieve("https://zenodo.org/records/15176439/files/hierarchy_bbp_atlas_pipeline.json?download=1",H)
d=json.load(open(H)); root=d['msg'][0] if isinstance(d,dict) and 'msg' in d else (d[0] if isinstance(d,list) else d)
nodes=[]
def walk(n):
    nodes.append(n)
    for c in n.get('children',[]): walk(c)
walk(root)
byid={n['id']:n for n in nodes}
print(f"struktur w hierarchii Zenodo: {len(nodes)}")

cb=[n for n in nodes if 'cerebell' in n.get('name','').lower() or
    any(n.get('acronym','').startswith(p) for p in
        ['SIM','ANcr','PRM','COPY','PFL','FL','LING','CENT','CUL','DEC','FOTU','PYR','UVU','NOD'])]
print(f"\n=== struktury mozdzku w hierarchii vs obecnosc w wolumenie 10um ===")
print(f"{'id':>8} {'akronim':<12} {'wokseli@10um':>13}  nazwa")
for n in sorted(cb,key=lambda x:x['id']):
    v=acc.get(n['id'],0)
    if 'Purkinje' in n.get('name','') or v>0:
        print(f"{n['id']:>8} {n.get('acronym',''):<12} {v:>13,}  {n.get('name','')}")
pu=[n for n in nodes if 'Purkinje layer' in n.get('name','')]
print(f"\nstruktur 'Purkinje layer' w hierarchii: {len(pu)}, z tego obecnych w wolumenie: "
      f"{sum(1 for n in pu if acc.get(n['id'],0)>0)}")
