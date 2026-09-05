"""AUDYT claimu: EBF2 rozroznia Purkinje Aldoc+ od Aldoc-.

Kazdy test moze claim OBALIC. Testy:
 A. proweniencja etykiet
 B. czy EBF2 koreluje z samym Aldoc na poziomie POJEDYNCZYCH komorek
    (globalnie i WEWNATRZ kazdego subklastra - to usuwa artefakt klastrowania)
 C. efekt probki: test parowany po probkach zamiast po komorkach
 D. czy rozdziela WEWNATRZ jednego plaika
 E. czy wynik zalezy od glebokosci sekwencjonowania
"""
import scanpy as sc, numpy as np, pandas as pd
from scipy import stats
P="/home/mikolajrurad/omics-data-highway/ichb-purkinje-analysis/pc-zebrin-sorting"
pd.set_option("display.width",250)

A=sc.read_h5ad(f"{P}/processed/purkinje_cells_v2.h5ad")
A.obs["grp"]=A.obs["aldoc_group"].astype(str)
A.obs["sub"]=A.obs["subcluster"].astype(str).str.replace("Purkinje_","")
raw_umi=np.asarray(A.X.sum(1)).ravel()
sc.pp.normalize_total(A,target_sum=1e4)
g={n:np.asarray(A.X[:,A.var_names.get_loc(n)].todense()).ravel()
   for n in ["Ebf2","Aldoc","Cux2","Plcb4","Ebf1"]}
pos=(A.obs["grp"]=="Aldoc+").values
o=A.obs

print("="*78)
print("A. PROWENIENCJA ETYKIET")
print("="*78)
print("  aldoc_group = 'Anti_Aldoc' w nazwie subklastra Kozarevy (02_build_h5ad.py:57)")
print("  Kozareva nazwala klastry PO klastrowaniu calego transkryptomu.")
print("  RYZYKO: klastrowanie moglo byc napedzane m.in. przez sam EBF2.")
print("  Test rozstrzygajacy -> sekcja B.\n")

print("="*78)
print("B. CZY EBF2 KORELUJE Z SAMYM Aldoc, NIEZALEZNIE OD KLASTROWANIA?")
print("="*78)
r,p=stats.spearmanr(g["Aldoc"],g["Ebf2"])
print(f"  wszystkie {A.n_obs} komorek:  Spearman r = {r:+.3f}  (p = {p:.2e})")
print("\n  WEWNATRZ kazdego subklastra osobno (etykieta stala -> brak artefaktu klastrowania):")
print(f"  {'subklaster':<16}{'n':>6}{'r(Aldoc,Ebf2)':>16}{'p':>12}{'sr.Aldoc':>10}{'sr.Ebf2':>10}")
zgodne=0; istotne=0
for s in sorted(o["sub"].unique()):
    m=(o["sub"]==s).values
    if m.sum()<30: continue
    rr,pp=stats.spearmanr(g["Aldoc"][m],g["Ebf2"][m])
    zgodne+= (rr<0); istotne+= (pp<0.05 and rr<0)
    print(f"  {s:<16}{m.sum():>6}{rr:>16.3f}{pp:>12.2e}{g['Aldoc'][m].mean():>10.2f}{g['Ebf2'][m].mean():>10.2f}")
print(f"\n  ujemna korelacja w {zgodne}/9 subklastrow, istotna (p<0.05) w {istotne}/9")

print("\n"+"="*78)
print("C. EFEKT PROBKI - test PAROWANY po probkach, nie po komorkach")
print("="*78)
rows=[]
for smp,d in o.groupby("orig_ident",observed=True):
    idx=d.index; m=o.index.isin(idx)
    mp=m&pos; mn=m&~pos
    if mp.sum()>=10 and mn.sum()>=10:
        rows.append(dict(probka=smp,n_pos=int(mp.sum()),n_neg=int(mn.sum()),
                         ebf2_pos=g["Ebf2"][mp].mean(), ebf2_neg=g["Ebf2"][mn].mean()))
D=pd.DataFrame(rows)
print(f"  probek z >=10 komorkami W OBU grupach: {len(D)} (z {o['orig_ident'].nunique()} wszystkich)")
D["roznica"]=D.ebf2_neg-D.ebf2_pos
w=stats.wilcoxon(D.ebf2_neg,D.ebf2_pos)
t=stats.ttest_rel(D.ebf2_neg,D.ebf2_pos)
print(f"  w ilu probkach Aldoc- ma WYZSZY EBF2 niz Aldoc+ w TEJ SAMEJ probce: {(D.roznica>0).sum()}/{len(D)}")
print(f"  Wilcoxon parowany: W={w.statistic:.0f}, p={w.pvalue:.3e}")
print(f"  t-test parowany:   t={t.statistic:.2f}, p={t.pvalue:.3e}")
print(f"  mediana roznicy: {D.roznica.median():.2f} CP10K  (zakres {D.roznica.min():.2f} do {D.roznica.max():.2f})")
print("\n  5 probek o najmniejszej roznicy (najslabszy sygnal):")
print(D.nsmallest(5,"roznica")[["probka","n_pos","n_neg","ebf2_pos","ebf2_neg","roznica"]].to_string(index=False))

print("\n"+"="*78)
print("D. CZY ROZDZIELA WEWNATRZ JEDNEGO PLAIKA?")
print("="*78)
print(f"  {'plaik':<8}{'n Aldoc+':>10}{'n Aldoc-':>10}{'EBF2 pos':>10}{'EBF2 neg':>10}{'AUC':>8}{'p':>12}")
from sklearn.metrics import roc_auc_score
aucs=[]
for reg in sorted(o["regions"].astype(str).unique()):
    m=(o["regions"].astype(str)==reg).values
    mp=m&pos; mn=m&~pos
    if mp.sum()<20 or mn.sum()<20: continue
    v=g["Ebf2"][m]; y=(~pos[m]).astype(int)
    auc=roc_auc_score(y,v)
    u=stats.mannwhitneyu(g["Ebf2"][mn],g["Ebf2"][mp])
    aucs.append(auc)
    print(f"  {reg:<8}{mp.sum():>10}{mn.sum():>10}{g['Ebf2'][mp].mean():>10.2f}{g['Ebf2'][mn].mean():>10.2f}{auc:>8.3f}{u.pvalue:>12.2e}")
print(f"\n  plaikow z obiema grupami >=20 kom.: {len(aucs)}, AUC od {min(aucs):.3f} do {max(aucs):.3f}, mediana {np.median(aucs):.3f}")

print("\n"+"="*78)
print("E. CZY TO NIE JEST EFEKT GLEBOKOSCI SEKWENCJONOWANIA?")
print("="*78)
print(f"  mediana UMI: Aldoc+ {np.median(raw_umi[pos]):.0f}, Aldoc- {np.median(raw_umi[~pos]):.0f}")
u=stats.mannwhitneyu(raw_umi[pos],raw_umi[~pos])
print(f"  Mann-Whitney na UMI: p={u.pvalue:.2e}")
r2,p2=stats.spearmanr(raw_umi,g["Ebf2"])
print(f"  korelacja EBF2 (CP10K) z glebokoscia UMI: r={r2:+.3f} (p={p2:.2e})")
print("\n  AUC EBF2 w kwartylach glebokosci (jesli stabilne -> nie jest to efekt glebokosci):")
q=np.quantile(raw_umi,[0,.25,.5,.75,1.0])
for i in range(4):
    m=(raw_umi>=q[i])&(raw_umi<=q[i+1])
    if (pos[m].sum()<20) or ((~pos[m]).sum()<20): continue
    print(f"    Q{i+1} (UMI {q[i]:.0f}-{q[i+1]:.0f}, n={m.sum()}): AUC={roc_auc_score((~pos[m]).astype(int),g['Ebf2'][m]):.3f}")
