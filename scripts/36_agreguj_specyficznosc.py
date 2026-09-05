"""Ekspresja 7 genow markerowych we WSZYSTKICH typach komorek atlasu (611 034 jader).

Odpowiada na pytanie: czy marker wyrozni Purkinje sposrod jader calego mozdzku,
czy tylko rozdziela Aldoc+/- WEWNATRZ juz wyizolowanych PC.
Mapowanie barcodow: ta sama logika co scripts/00_prep_map.py (12 przemianowanych probek).
"""
import csv, numpy as np, pandas as pd
P="/home/mikolajrurad/omics-data-highway/ichb-purkinje-analysis/pc-zebrin-sorting"
GEN={2302:"Plcb4",4548:"Calb1",9615:"Syne1",10307:"Nav3",13976:"Aldoc",15708:"Syne2",18983:"Gnaq"}
REPO2GEO={'DECa_F002':'VI2a_F002','DECb_F002':'VI2b_F002','DECc_F002':'VI2c_F002',
          'DECa_M006':'VI2a_M006','DECb_M006':'VI2b_M006','DECc_M006':'VI2c_M006',
          'Va_M002':'CULa_M002','Vb_M002':'CULb_M002','Vc_M002':'CULc_M002',
          'Vd_M002':'CULd_M002','IV_M006':'CUL_M006','Vl_M006':'VI_M006'}
def to_geo(bc):
    a=bc.split('_'); pre='_'.join(a[:2]); suf='_'.join(a[2:])
    return f"{REPO2GEO[pre]}_{suf}" if pre in REPO2GEO else bc

col={}
with open(f"{P}/raw/cb_adult_mouse_barcodes.txt") as fh:
    for i,l in enumerate(fh,1): col[l.strip()]=i
print("barcodow w mtx:",len(col))

NC=len(col)+1
typ=np.zeros(NC,dtype=np.int16); umi=np.zeros(NC,dtype=np.float32)
names=["<brak>"]; nidx={}
nreg=nmiss=0
with open(f"{P}/cerebellum-atlas-analysis/data/full_cb_metadata.csv",newline='') as fh:
    r=csv.reader(fh); h=next(r); ix={n:i for i,n in enumerate(h[1:],1)}; ix['barcode']=0
    for row in r:
        cl=row[ix['final_annotation_cluster']]
        if cl=="REMOVED": nreg+=1; continue
        c=col.get(to_geo(row[0]))
        if not c: nmiss+=1; continue
        if cl not in nidx: nidx[cl]=len(names); names.append(cl)
        typ[c]=nidx[cl]; umi[c]=float(row[ix['nUMI']])
print(f"REMOVED pominietych: {nreg}   NIETRAFIONYCH barcodow: {nmiss}")
print(f"komorek z przypisanym typem: {int((typ>0).sum())}   typow: {len(names)-1}")
assert nmiss==0, "mapowanie barcodow niekompletne"

ncell=np.bincount(typ,minlength=len(names))
S={g:np.zeros(len(names)) for g in GEN.values()}   # suma CP10K
N={g:np.zeros(len(names)) for g in GEN.values()}   # liczba komorek dodatnich
nline=0
with open("/mnt/data1t/pc_rebuild/geny_atlas.tsv") as fh:
    for line in fh:
        a,b,v=line.split()
        c=int(b); t=typ[c]
        if t==0: continue
        g=GEN[int(a)]
        S[g][t]+= int(v)/umi[c]*1e4
        N[g][t]+= 1
        nline+=1
print(f"wpisow przetworzonych: {nline}")

rows=[]
for i in range(1,len(names)):
    n=ncell[i]
    d=dict(typ=names[i], n_komorek=int(n))
    for g in GEN.values():
        d[f"{g}_cp10k"]=round(S[g][i]/n,3)
        d[f"{g}_pct"]=round(N[g][i]/n*100,1)
    rows.append(d)
T=pd.DataFrame(rows).sort_values("n_komorek",ascending=False)
T.to_csv(f"{P}/processed/specyficznosc_typy_komorek.csv",index=False)
pd.set_option("display.width",260)
print(f"\n=== EKSPRESJA W {len(T)} TYPACH KOMOREK ATLASU (suma {T.n_komorek.sum()} jader) ===")
print(T[["typ","n_komorek","Nav3_cp10k","Nav3_pct","Syne1_cp10k","Syne1_pct",
         "Aldoc_cp10k","Aldoc_pct","Calb1_cp10k","Calb1_pct"]].to_string(index=False))
print("\n=== KROTNOSC Purkinje vs najwyzszy INNY typ ===")
pk=T[T.typ=="Purkinje"].iloc[0]
oth=T[T.typ!="Purkinje"]
for g in GEN.values():
    mx=oth[f"{g}_cp10k"].max(); nm=oth.loc[oth[f"{g}_cp10k"].idxmax(),"typ"]
    kr=pk[f"{g}_cp10k"]/mx if mx>0 else float("inf")
    print(f"  {g:7s} Purkinje {pk[f'{g}_cp10k']:8.3f}  |  max inny: {nm} {mx:.3f}  |  krotnosc {kr:.1f}x")
print(f"\nzapisano processed/specyficznosc_typy_komorek.csv")
