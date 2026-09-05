"""Faza 5/6: polaczenie geometrii warstwy Purkinjego (atlas 10um)
z kompozycja podtypow zmierzona u Kozarevy (per placik)."""
import numpy as np, pandas as pd, scanpy as sc, json
OUT="/mnt/data1t/pc_rebuild"
P="/home/mikolajrurad/omics-data-highway/ichb-purkinje-analysis/pc-zebrin-sorting"
# mapowanie zweryfikowane 16/16
KOZ2ATLAS={'I':'LING','II':'CENT2','III':'CENT3','CUL':'CUL4, 5','VI':'DEC','VII':'FOTU',
           'VIII':'PYR','IX':'UVU','X':'NOD','SIM':'SIM','AN1':'ANcr1','AN2':'ANcr2',
           'PRM':'PRM','COP':'COPY','PF':'PFL','F':'FL'}
A=sc.read_h5ad(f"{P}/processed/purkinje_cells_v2.h5ad")
obs=A.obs[["regions","subcluster","aldoc_group"]].copy()
del A
obs["atlas"]=obs["regions"].astype(str).map(KOZ2ATLAS)
assert obs["atlas"].notna().all(), "brak mapowania dla jakiegos regionu"
subs=sorted(obs["subcluster"].astype(str).unique())
print(f"podtypow: {len(subs)}  placikow: {obs['atlas'].nunique()}")

comp=(obs.groupby(["atlas","subcluster"],observed=True).size()
        .unstack(fill_value=0))
comp=comp.div(comp.sum(1),axis=0)
ald=(obs.assign(pos=obs["aldoc_group"].astype(str).eq("Aldoc+"))
       .groupby("atlas",observed=True)["pos"].mean())
n_cells=obs.groupby("atlas",observed=True).size()
print("\n=== kompozycja per placik (frakcje, z danych Kozarevy) ===")
tab=comp.copy(); tab.columns=[c.replace("Purkinje_","") for c in tab.columns]
tab.insert(0,"n_kom",n_cells); tab.insert(1,"Aldoc+_frac",ald.round(3))
print(tab.round(3).to_string())
tab.to_csv(f"{P}/processed/lobule_composition.csv")

d=np.load(f"{OUT}/pc_layer_10um.npz",allow_pickle=True)
names=[str(x) for x in d["names"]]; lob=d["lobule"]
i,j,k=d["i"],d["j"],d["k"]
print(f"\nwokseli warstwy Purkinjego: {len(i):,}")
name_of={idx+1:nm for idx,nm in enumerate(names)}
lob_name=np.array([name_of[x] for x in lob])
missing=set(lob_name)-set(comp.index)
print("placiki bez danych Kozarevy:", missing if missing else "brak")

aldv=np.array([ald.get(n,np.nan) for n in lob_name],dtype=np.float32)
rng=np.random.default_rng(0)
sub_idx=np.full(len(i),-1,dtype=np.int8)
for nm in comp.index:
    m=lob_name==nm
    if m.sum()==0: continue
    p=comp.loc[nm,subs].to_numpy(dtype=float); p=p/p.sum()
    sub_idx[m]=rng.choice(len(subs),size=int(m.sum()),p=p).astype(np.int8)
np.savez_compressed(f"{OUT}/scene_pc_layer.npz",
    i=i,j=j,k=k,lobule=lob,lobule_name=lob_name,aldoc_frac=aldv,
    subtype_idx=sub_idx,subtypes=np.array(subs),names=np.array(names))
print(f"zapisano scene_pc_layer.npz")
print("\nUWAGA metodologiczna: subtype_idx to LOSOWANIE z proporcji zmierzonych")
print("dla danego placika. Proporcje sa prawdziwe, przypisanie pojedynczego")
print("woksela - nie. Kolor ciagly aldoc_frac jest w pelni pomiarowy.")
