"""Faza 1: QC + normalizacja + PCA/UMAP/Leiden na przebudowanym v2, z walidacja ARI."""
import scanpy as sc, numpy as np, pandas as pd, anndata as ad, os, psutil
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
GiB=1024**3
def mem(tag): print(f"  [RAM {psutil.Process().memory_info().rss/GiB:.2f} GiB] {tag}", flush=True)
P="/home/mikolajrurad/omics-data-highway/ichb-purkinje-analysis/pc-zebrin-sorting"
sc.settings.verbosity=1
A=sc.read_h5ad(f"{P}/processed/purkinje_cells_v2.h5ad"); mem("wczytane")
print(A)
print("\n=== QC: rozklady (progi pokazane, nie narzucone po cichu) ===")
for c in ["nGene","nUMI","percent_mito"]:
    v=A.obs[c].astype(float)
    print(f"  {c:14s} min={v.min():.4g} p1={np.percentile(v,1):.4g} mediana={v.median():.4g} "
          f"p99={np.percentile(v,99):.4g} max={v.max():.4g}")
A.var["mt"]=A.var_names.str.startswith(("mt-","Mt-","MT-"))
print(f"  genow mitochondrialnych w macierzy: {int(A.var['mt'].sum())}")
sc.pp.calculate_qc_metrics(A, qc_vars=["mt"], inplace=True, log1p=False, percent_top=None)
print(f"  policzone pct_counts_mt: mediana={A.obs['pct_counts_mt'].median():.3f} "
      f"max={A.obs['pct_counts_mt'].max():.3f}")

n0=A.n_obs
keep=(A.obs["n_genes_by_counts"]>=500)&(A.obs["pct_counts_mt"]<5)
print(f"\n=== filtr: nGene>=500 i pct_mt<5 ===")
print(f"  zostaje {int(keep.sum())} z {n0}  (odrzucone {n0-int(keep.sum())})")
A=A[keep].copy(); mem("po filtrze")
sc.pp.filter_genes(A, min_cells=3)
print(f"  genow po filtrze min_cells=3: {A.n_vars}")

A.layers["counts"]=A.X.copy()
sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
A.raw=A
sc.pp.highly_variable_genes(A, n_top_genes=2000, batch_key=None)
print(f"  HVG: {int(A.var.highly_variable.sum())}")
A=A[:,A.var.highly_variable].copy(); mem("po HVG")
sc.pp.scale(A, max_value=10)
sc.tl.pca(A, n_comps=50, svd_solver="arpack")
sc.pp.neighbors(A, n_neighbors=15, n_pcs=30)
sc.tl.umap(A); mem("po UMAP")

truth=A.obs["subcluster"].astype(str).values
print("\n=== Leiden vs anotacja Kozarevy (walidacja, nie ozdoba) ===")
print(f"{'rozdzielczosc':>14}{'klastrow':>10}{'ARI':>8}{'NMI':>8}")
best=None
for r in [0.2,0.3,0.4,0.5,0.6,0.8,1.0]:
    key=f"leiden_{r}"
    sc.tl.leiden(A, resolution=r, key_added=key, flavor="igraph", n_iterations=2, directed=False)
    lab=A.obs[key].astype(str).values
    ari=adjusted_rand_score(truth,lab); nmi=normalized_mutual_info_score(truth,lab)
    k=A.obs[key].nunique()
    print(f"{r:>14}{k:>10}{ari:>8.3f}{nmi:>8.3f}")
    if best is None or ari>best[1]: best=(r,ari,nmi,k)
print(f"\nnajlepsze ARI: rozdzielczosc {best[0]}, ARI={best[1]:.3f}, NMI={best[2]:.3f}, klastrow={best[3]}")
A.uns["leiden_validation"]={"best_resolution":best[0],"ARI":float(best[1]),"NMI":float(best[2]),"n_clusters":int(best[3])}

# przywrocenie pelnego zestawu genow (raw) do zapisu
full=A.raw.to_adata()
full.obs=A.obs.copy()
for k in ["X_pca","X_umap"]: full.obsm[k]=A.obsm[k]
full.uns["leiden_validation"]=A.uns["leiden_validation"]
full.uns["provenance"]=A.uns.get("provenance",{})
dest=f"{P}/processed/purkinje_cells_v2_processed.h5ad"
full.write_h5ad(dest, compression="gzip")
print(f"\nzapisano {dest} ({os.path.getsize(dest):,} B)")
print(full)
mem("koniec")
