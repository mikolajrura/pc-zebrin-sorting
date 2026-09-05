"""Eksport chmury punktow do strony 3D (kwantyzacja int16 + base64)."""
import numpy as np, json, base64, pandas as pd
OUT="/mnt/data1t/pc_rebuild"; P="/home/mikolajrurad/omics-data-highway/ichb-purkinje-analysis/pc-zebrin-sorting"
d=np.load(f"{OUT}/scene_pc_layer.npz",allow_pickle=True)
i,j,k=d["i"],d["j"],d["k"]; lob=d["lobule"]; names=[str(x) for x in d["names"]]
comp=pd.read_csv(f"{P}/processed/lobule_composition.csv",index_col=0)
print("wokseli:",len(i))
rng=np.random.default_rng(1); N=130000
s=rng.choice(len(i),N,replace=False)
X=i[s].astype(np.float64); Y=j[s].astype(np.float64); Z=k[s].astype(np.float64)
# wysrodkuj i przeskaluj na int16
P0=np.c_[X,Y,Z]; ctr=P0.mean(0); P0-=ctr
sc=32000.0/np.abs(P0).max()
Q=np.round(P0*sc).astype(np.int16)
print("skala:",sc,"zakres:",Q.min(),Q.max())
lb=lob[s].astype(np.uint8)
b64=lambda a: base64.b64encode(a.tobytes()).decode()
meta={}
for idx,nm in enumerate(names,1):
    if nm in comp.index:
        r=comp.loc[nm]
        subs={c:float(r[c]) for c in comp.columns if c not in ("n_kom","Aldoc+_frac")}
        meta[idx]={"name":nm,"n":int(r["n_kom"]),"aldoc":float(r["Aldoc+_frac"]),
                   "dominant":max(subs,key=subs.get),"subs":subs}
KOZ={'LING':'I','CENT2':'II','CENT3':'III','CUL4, 5':'CUL (IV–V)','DEC':'VI','FOTU':'VII',
     'PYR':'VIII','UVU':'IX','NOD':'X','SIM':'SIM','ANcr1':'AN1 (Crus 1)','ANcr2':'AN2 (Crus 2)',
     'PRM':'PRM','COPY':'COP','PFL':'PF','FL':'F'}
for kk,v in meta.items(): v["koz"]=KOZ.get(v["name"],v["name"])
payload={"n":int(N),"scale":float(sc),"pos":b64(Q),"lob":b64(lb),"meta":meta,
         "subtypes":[c for c in comp.columns if c not in ("n_kom","Aldoc+_frac")]}
json.dump(payload,open(f"{OUT}/web_cloud.json","w"))
import os; print("json:",os.path.getsize(f'{OUT}/web_cloud.json'),"B")
