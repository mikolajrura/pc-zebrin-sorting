"""Przeglad 9 fenotypow Purkinje: markery per subklaster + reprezentatywne komorki (medoidy PCA).
Wejscie: processed/purkinje_cells_v2_processed.h5ad (lognorm CP10K) + _v2.h5ad (surowe zliczenia).
Wyjscie: processed/fenotypy_markery.csv, processed/fenotypy_medoidy.csv,
         processed/fenotypy_medoidy_ekspresja.csv
"""
import scanpy as sc, numpy as np, pandas as pd, psutil, h5py, scipy.sparse as sp
GiB=1024**3
P="/home/mikolajrurad/omics-data-highway/ichb-purkinje-analysis/pc-zebrin-sorting"
def mem(t): print(f"  [RAM {psutil.Process().memory_info().rss/GiB:.2f} GiB] {t}", flush=True)

A=sc.read_h5ad(f"{P}/processed/purkinje_cells_v2_processed.h5ad")
mem("wczytane")
print("ksztalt:",A.shape)
A.obs["sub"]=A.obs["subcluster"].astype(str)
subs=sorted(A.obs["sub"].unique())
print("liczności:"); print(A.obs["sub"].value_counts().sort_index().to_string())

# ---- 1. markery per subklaster (Wilcoxon, jeden-vs-reszta, PELNY zbior) ----
print("\n=== rank_genes_groups Wilcoxon, 9 grup vs reszta ===", flush=True)
sc.tl.rank_genes_groups(A,"sub",method="wilcoxon",n_genes=100,pts=True)
rows=[]
for g in subs:
    d=sc.get.rank_genes_groups_df(A,group=g)
    d["subcluster"]=g
    rows.append(d)
M=pd.concat(rows)
M.to_csv(f"{P}/processed/fenotypy_markery.csv",index=False)
print(f"zapisano fenotypy_markery.csv  ({len(M)} wierszy)")
mem("po markerach")

# ---- 2. medoidy: 5 komorek najblizszych centroidowi w 30 PC ----
NPC=30; K=5
X=A.obsm["X_pca"][:,:NPC]
sel_idx=[]; sel_sub=[]; sel_dist=[]
for g in subs:
    m=(A.obs["sub"]==g).values
    ii=np.where(m)[0]
    c=X[ii].mean(0)
    d=np.linalg.norm(X[ii]-c,axis=1)
    o=np.argsort(d)[:K]
    sel_idx.extend(ii[o]); sel_sub.extend([g]*len(o)); sel_dist.extend(d[o])
sel_idx=np.array(sel_idx)
med=A.obs.iloc[sel_idx][["geo_barcode","orig_ident","regions","subcluster","aldoc_group",
                          "nGene","nUMI","pct_counts_mt","leiden_0.2"]].copy()
med["dist_do_centroidu_30PC"]=sel_dist
med["idx_processed"]=sel_idx
med.to_csv(f"{P}/processed/fenotypy_medoidy.csv")
print(f"\nzapisano fenotypy_medoidy.csv ({len(med)} komorek)")

# ---- 3. transkryptom kazdej wybranej komorki: top 25 genow po CP10K ----
V=A.var_names.to_numpy()
sub_X=A.X[sel_idx].toarray()          # 45 x 22325 lognorm
cp10k=np.expm1(sub_X)                 # z powrotem na CP10K
recs=[]
for r,(i,g) in enumerate(zip(sel_idx,sel_sub)):
    v=cp10k[r]
    top=np.argsort(-v)[:25]
    for rank,j in enumerate(top,1):
        recs.append(dict(subcluster=g, barcode=A.obs["geo_barcode"].iloc[i],
                         ranga=rank, gen=V[j], cp10k=round(float(v[j]),2)))
E=pd.DataFrame(recs)
E.to_csv(f"{P}/processed/fenotypy_medoidy_ekspresja.csv",index=False)
print(f"zapisano fenotypy_medoidy_ekspresja.csv ({len(E)} wierszy)")

# ---- 4. panel markerow kanonicznych: srednia CP10K + % komorek dodatnich per subklaster ----
PANEL=["Aldoc","Plcb4","Plcb3","Calb1","Car8","Pcp2","Pcp4","Itpr1","Grid2","Slc1a6",
       "Atxn1","Atxn7","Nrgn","Cck","Gpr176","Tox2","Ebf2","Cux2","Casq2","Eaat4"]
have=[g for g in PANEL if g in set(V)]
miss=[g for g in PANEL if g not in set(V)]
print(f"\npanel: obecne {len(have)}, brak {miss}")
cols={g:A.var_names.get_loc(g) for g in have}
sub_codes=A.obs["sub"].values
tabM=pd.DataFrame(index=subs,columns=have,dtype=float)
tabP=pd.DataFrame(index=subs,columns=have,dtype=float)
for g in subs:
    m=(sub_codes==g)
    blk=A.X[m][:,[cols[x] for x in have]].toarray()
    lin=np.expm1(blk)
    tabM.loc[g]=lin.mean(0).round(2)
    tabP.loc[g]=(blk>0).mean(0).round(3)*100
tabM.to_csv(f"{P}/processed/fenotypy_panel_srednia_cp10k.csv")
tabP.to_csv(f"{P}/processed/fenotypy_panel_pct_dodatnich.csv")
print("\n=== srednia CP10K, panel kanoniczny ===")
print(tabM.to_string())
print("\n=== % komorek z niezerowym odczytem ===")
print(tabP.to_string())
mem("koniec")
