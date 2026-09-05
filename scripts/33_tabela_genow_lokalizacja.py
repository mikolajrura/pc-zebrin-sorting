"""Tabela 24 409 genow x (lokalizacja UniProt + profil w 9 fenotypach Purkinje).

Wejscie:  processed/purkinje_cells_v2.h5ad          (16 634 x 24 409, surowe zliczenia)
          references/uniprot_mouse_subcell.tsv      (Swiss-Prot mysz, reviewed)
Wyjscie:  processed/geny_purkinje_lokalizacja.csv

Wszystkie kolumny liczone na PELNYM zbiorze 16 634 komorek Purkinje (bez filtra QC,
ktory odrzucal 15 komorek) i PELNYM zestawie 24 409 genow.
Kategoryzacja lokalizacji: ta sama logika co scripts/29_lokalizacja.py.
"""
import scanpy as sc, numpy as np, pandas as pd, re, psutil, scipy.sparse as sp
GiB=1024**3
P="/home/mikolajrurad/omics-data-highway/ichb-purkinje-analysis/pc-zebrin-sorting"
def mem(t): print(f"  [RAM {psutil.Process().memory_info().rss/GiB:.2f} GiB] {t}", flush=True)

# ---------- 1. dane ----------
A=sc.read_h5ad(f"{P}/processed/purkinje_cells_v2.h5ad"); mem("wczytane")
print("ksztalt:",A.shape,"nnz:",A.X.nnz)
A.obs["sub"]=A.obs["subcluster"].astype(str).str.replace("Purkinje_","")
subs=sorted(A.obs["sub"].unique())
genes=A.var_names.to_numpy()

X=A.X.astype(np.float32)
tot=np.asarray(X.sum(1)).ravel()
print(f"suma zliczen/komorke: min={tot.min():.0f} mediana={np.median(tot):.0f} max={tot.max():.0f}")
sc.pp.normalize_total(A,target_sum=1e4); A.X=A.X.astype(np.float32); mem("CP10K")
CP=A.X.tocsc()                       # CSC = szybki dostep po genach
mem("CSC")

# ---------- 2. profil per fenotyp ----------
out=pd.DataFrame({"gen":genes})
masks={s:(A.obs["sub"]==s).values for s in subs}
ncell=A.n_obs
nz_glob=np.diff(CP.indptr)                        # w ilu komorkach niezerowy
out["n_komorek_wykryty"]=nz_glob
out["pct_komorek"]= (nz_glob/ncell*100).round(2)
out["srednia_cp10k_all"]=np.asarray(CP.sum(0)).ravel()/ncell
for s in subs:
    m=masks[s]; n=int(m.sum())
    B=CP[m]                                       # podzbior komorek
    out[f"mean_{s}"]=(np.asarray(B.sum(0)).ravel()/n)
    out[f"pct_{s}"]=(np.diff(B.tocsc().indptr)/n*100)
for c in out.columns:
    if c.startswith(("mean_","srednia")): out[c]=out[c].round(3)
    if c.startswith("pct_"): out[c]=out[c].round(2)
mem("profile per fenotyp")

# ---------- 3. markery (Wilcoxon, pelny ranking) ----------
sc.pp.log1p(A)
print("\nrank_genes_groups Wilcoxon, 9 grup, pelny ranking...", flush=True)
sc.tl.rank_genes_groups(A,"sub",method="wilcoxon",n_genes=A.n_vars)
best={}
for s in subs:
    d=sc.get.rank_genes_groups_df(A,group=s).set_index("names")
    d=d.reindex(genes)
    ok=(d.pvals_adj<0.01)&(d.logfoldchanges>0.5)
    best[s]=(d["scores"].where(ok,-np.inf).to_numpy(),
             d["logfoldchanges"].to_numpy(), d["pvals_adj"].to_numpy())
S=np.vstack([best[s][0] for s in subs])
arg=S.argmax(0); mx=S.max(0)
out["marker_dla"]=[subs[i] if np.isfinite(mx[j]) else "" for j,i in enumerate(arg)]
out["marker_logFC"]=[round(float(best[subs[i]][1][j]),3) if np.isfinite(mx[j]) else np.nan
                     for j,i in enumerate(arg)]
out["marker_padj"]=[float(best[subs[i]][2][j]) if np.isfinite(mx[j]) else np.nan
                    for j,i in enumerate(arg)]
out["marker_score"]=[round(float(mx[j]),2) if np.isfinite(mx[j]) else np.nan
                     for j in range(len(genes))]
print("genow bedacych markerem >=1 fenotypu:", int((out.marker_dla!="").sum()))
mem("markery")

# ---------- 4. UniProt ----------
u=pd.read_csv(f"{P}/references/uniprot_mouse_subcell.tsv",sep="\t")
u.columns=["gen","syn","loc"]; u["loc"]=u["loc"].fillna("")
EVID=re.compile(r"\{[^}]*\}"); NOTE=re.compile(r"Note=.*?(?=(?:SUBCELLULAR LOCATION:)|$)",re.S)
def terms(t):
    t=EVID.sub("",t); t=NOTE.sub("",t)
    t=t.replace("SUBCELLULAR LOCATION:"," "); t=re.sub(r"\[[^\]]*\]:"," ",t)
    return t.lower()
OTOCZKA=re.compile(r"nucleus (inner |outer )?membrane|nucleus envelope|nucleus lamina|nuclear pore|nucleus, nuclear pore")
JADRO  =re.compile(r"\bnucleus\b|\bchromosome\b|\bnucleolus\b|\bnucleoplasm\b")
INNE   =re.compile(r"\bcytoplasm\b|cell membrane|\bsecreted\b|mitochondri|endoplasmic reticulum|"
                   r"golgi|cell projection|cell surface|lysosome|endosome|peroxisome|"
                   r"cytoskeleton|cell junction|synapse|\bvesicle\b|apical|basolateral")
idx={}
for _,r in u.iterrows():
    t=terms(r["loc"])
    oto=bool(OTOCZKA.search(t)); jad=bool(JADRO.search(t)); inn=bool(INNE.search(t))
    kat="otoczka" if oto else ("jadro" if jad else ("reszta" if t.strip() else "brak_opisu"))
    rec=(kat,jad,oto,inn,r["loc"][:400],r["gen"])
    if isinstance(r["gen"],str): idx.setdefault(r["gen"].lower(),rec)
    if isinstance(r["syn"],str):
        for s2 in r["syn"].split(): idx.setdefault(s2.lower(),rec)
print(f"\nUniProt: {len(u)} rekordow -> {len(idx)} kluczy (nazwy + synonimy)")

hit=[idx.get(g.lower()) for g in genes]
out["lokalizacja"]=[h[0] if h else "nieznane" for h in hit]
out["w_jadrze"]  =[h[1] if h else False for h in hit]
out["otoczka"]   =[h[2] if h else False for h in hit]
out["tez_poza"]  =[h[3] if h else False for h in hit]
out["tylko_jadro"]=out["w_jadrze"]&~out["tez_poza"]
out["uniprot_gen"]=[h[5] if h else "" for h in hit]
out["dopasowanie"]=["nazwa" if h and h[5].lower()==g.lower() else ("synonim" if h else "brak")
                    for g,h in zip(genes,hit)]
out["loc_txt"]=[h[4] if h else "" for h in hit]

# ---------- 5. os Aldoc ----------
Ap=[f"mean_Aldoc_{i}" for i in range(1,8)]; An=["mean_Anti_Aldoc_1","mean_Anti_Aldoc_2"]
out["mean_Aldoc_plus"]=out[Ap].mean(1).round(3)
out["mean_Aldoc_minus"]=out[An].mean(1).round(3)
out["log2FC_plus_vs_minus"]=np.log2((out["mean_Aldoc_plus"]+0.01)/(out["mean_Aldoc_minus"]+0.01)).round(3)

cols=(["gen","lokalizacja","w_jadrze","otoczka","tylko_jadro","tez_poza","dopasowanie","uniprot_gen",
       "n_komorek_wykryty","pct_komorek","srednia_cp10k_all",
       "mean_Aldoc_plus","mean_Aldoc_minus","log2FC_plus_vs_minus",
       "marker_dla","marker_logFC","marker_padj","marker_score"]
      +[f"mean_{s}" for s in subs]+[f"pct_{s}" for s in subs]+["loc_txt"])
out=out[cols]
dest=f"{P}/processed/geny_purkinje_lokalizacja.csv"
out.to_csv(dest,index=False)
import os; print(f"\nzapisano {dest}  ({os.path.getsize(dest):,} B, {len(out)} wierszy, {len(cols)} kolumn)")

print("\n=== POKRYCIE ADNOTACJA ===")
print(out["lokalizacja"].value_counts().to_string())
print("\n=== sposob dopasowania ===")
print(out["dopasowanie"].value_counts().to_string())
print("\n=== pokrycie wsrod genow FAKTYCZNIE wykrytych (>=1% komorek) ===")
w=out[out.pct_komorek>=1]
print(f"genow wykrytych w >=1% komorek: {len(w)}")
print(w["lokalizacja"].value_counts().to_string())
print("\n=== pokrycie wsrod markerow ===")
mk=out[out.marker_dla!=""]
print(f"markerow: {len(mk)}")
print(mk["lokalizacja"].value_counts().to_string())
mem("koniec")
