"""Faza 0, krok 3: zlozenie h5ad + WALIDACJA per komorka.

Walidacja jest istotą tego kroku: jesli dla kazdej z 16 634 komorek liczba
wyekstrahowanych wpisow rowna sie nGene z metadanych, to nic nie zostalo
pominiete - niezaleznie od tego, czy plik mtx byl posortowany.
"""
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, sys, os
OUT="/mnt/data1t/pc_rebuild"
P="/home/mikolajrurad/omics-data-highway/ichb-purkinje-analysis/pc-zebrin-sorting"
GiB=1024**3

obs=pd.read_csv(f"{OUT}/obs.tsv",sep="\t")
n_cells=len(obs); print("komorek w obs:", n_cells)
lo,hi,ncol,expected=open(f"{OUT}/block.txt").read().split()
expected=int(expected)

print("wczytuje triplety ...", flush=True)
t=pd.read_csv(f"{OUT}/triplets.tsv",sep="\t",header=None,
              names=["cell","gene","val"],
              dtype={"cell":np.int32,"gene":np.int32,"val":np.int32})
print(f"wpisow: {len(t):,}   pamiec: {t.memory_usage(deep=True).sum()/GiB:.2f} GiB")

# ---- WALIDACJA ----
print("\n=== WALIDACJA per komorka ===")
got=np.bincount(t["cell"].to_numpy(), minlength=n_cells)
want=obs["nGene"].to_numpy()
ok=(got==want)
print(f"komorek zgodnych nnz==nGene: {ok.sum()} / {n_cells}")
if not ok.all():
    bad=np.flatnonzero(~ok)
    print(f"NIEZGODNYCH: {len(bad)}")
    d=pd.DataFrame({"row":bad,"barcode":obs['barcode'].to_numpy()[bad],
                    "nGene":want[bad],"wyekstrahowano":got[bad],
                    "roznica":got[bad]-want[bad]})
    print(d.head(20).to_string(index=False))
    d.to_csv(f"{OUT}/walidacja_niezgodne.csv",index=False)
    print(f"pelna lista: {OUT}/walidacja_niezgodne.csv")
    print(f"suma roznic: {int((got-want).sum())}")
else:
    print("WSZYSTKIE ZGODNE")
print(f"suma wpisow: {len(t):,}   oczekiwana (suma nGene): {expected:,}   roznica: {len(t)-expected:,}")

genes=[l.strip() for l in open(f"{P}/raw/cb_adult_mouse_genes.txt")]
n_genes=len(genes); print(f"\ngenow: {n_genes}")

X=sp.csr_matrix((t["val"].to_numpy(np.int32),
                 (t["cell"].to_numpy(np.int32), t["gene"].to_numpy(np.int32))),
                shape=(n_cells,n_genes), dtype=np.int32)
X.sort_indices()
print(f"macierz: {X.shape}, nnz={X.nnz:,}, {(X.data.nbytes+X.indices.nbytes+X.indptr.nbytes)/GiB:.2f} GiB")
del t

o=obs.set_index("barcode")
o.index.name=None
for c in ["orig_ident","regions","subcluster","geo_barcode"]:
    o[c]=o[c].astype("category")
o["aldoc_group"]=np.where(o["subcluster"].astype(str).str.contains("Anti_Aldoc"),
                          "Aldoc-","Aldoc+")
o["aldoc_group"]=o["aldoc_group"].astype("category")
var=pd.DataFrame(index=pd.Index(genes,name=None)); var["gene"]=genes

A=ad.AnnData(X=X, obs=o, var=var)
A.layers["counts"]=X.copy()
A.uns["provenance"]={
 "source":"GSE165371 / Kozareva et al. 2021 Nature 598:214-219",
 "mtx":"cb_adult_mouse.mtx.gz (24409 x 611034, nnz 1026400102)",
 "block_columns":f"{lo}-{hi}",
 "rebuild_reason":"oryginalny purkinje_cells.h5ad miał 15467 komórek; "
                  "naiwny join gubił 1167 (CUL -456, VI -711)",
 "sample_rename_applied":"DEC*->VI2*, Va-Vd_M002->CULa-d_M002, IV_M006->CUL_M006, Vl_M006->VI_M006",
 "validation":f"nnz==nGene dla {int(ok.sum())}/{n_cells} komórek",
}
dest=f"{P}/processed/purkinje_cells_v2.h5ad"
A.write_h5ad(dest, compression="gzip")
print(f"\nzapisano {dest}")
print(f"rozmiar: {os.path.getsize(dest):,} B")
print(A)
print("\n=== rozklad po regionach (kontrola CUL i VI) ===")
print(A.obs["regions"].value_counts().sort_index().to_string())
print("\n=== subklastry ===")
print(A.obs["subcluster"].value_counts().to_string())
