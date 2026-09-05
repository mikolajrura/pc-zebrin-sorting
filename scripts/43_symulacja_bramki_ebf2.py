"""Symulacja dot plotu SSC-A / anty-EBF2 dla bramki R3, na danych Kozarevy.

ZALOZENIA (jawne, bo to symulacja, nie pomiar):
 1. Sygnal przeciwciala jest PROPORCJONALNY do mRNA EBF2 (CP10K) - to zalozenie
    NAJKORZYSTNIEJSZE z mozliwych. Realnie bialko uśrednia po czasie, ma tlo
    i nie ma dropoutu, ktory tu sztucznie rozsuwa rozklady.
 2. Brak szumu technicznego, brak autofluorescencji, brak tla przeciwciala.
 3. Skala logarytmiczna cytometru: sygnal = 10^3 * (CP10K + 0.3), czyli offset
    odpowiadajacy podlozu; wartosc arbitralna, wplywa tylko na polozenie osi.
 4. SSC-A losowane z rozkladu lognormalnego dopasowanego do mediany 185 251
    (R3, cerebellum-1-test) i rozrzutu obserwowanego na wykresie.
 5. Liczba zdarzen = 3 152 (najlepsza z trzech probek uzytkownika).

Wynik jest GORNA GRANICA tego, co da sie zobaczyc - nie prognoza.
"""
import scanpy as sc, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
rng=np.random.default_rng(0)
P="/home/mikolajrurad/omics-data-highway/ichb-purkinje-analysis/pc-zebrin-sorting"
BLUE="#2a78d6"; ORANGE="#eb6834"; INK="#1a1a19"; MUTED="#6b6a63"; GRID="#e6e5e0"

A=sc.read_h5ad(f"{P}/processed/purkinje_cells_v2.h5ad")
A.obs["grp"]=A.obs["aldoc_group"].astype(str)
sc.pp.normalize_total(A,target_sum=1e4)
i=A.var_names.get_loc("Ebf2")
v=np.asarray(A.X[:,i].todense()).ravel()
pos=(A.obs["grp"]=="Aldoc+").values

N=3152
idx=rng.choice(len(v),size=N,replace=False)
vv=v[idx]; pp=pos[idx]
sig=1000.0*(vv+0.3)
ssc=rng.lognormal(np.log(185251),0.28,size=N)

fig,ax=plt.subplots(1,2,figsize=(12.5,5.2),gridspec_kw=dict(width_ratios=[1,0.85]))
fig.patch.set_facecolor("#fcfcfb")
for a in ax:
    a.set_facecolor("#fcfcfb")
    for s in ("top","right"): a.spines[s].set_visible(False)
    for s in ("left","bottom"): a.spines[s].set_color(GRID)
    a.tick_params(colors=MUTED,labelsize=9)

a=ax[0]
a.scatter(sig[pp],ssc[pp],s=5,c=BLUE,alpha=.45,lw=0,label=f"Aldoc+  (n={pp.sum()})")
a.scatter(sig[~pp],ssc[~pp],s=5,c=ORANGE,alpha=.45,lw=0,label=f"Aldoc−  (n={(~pp).sum()})")
a.set_xscale("log"); a.set_yscale("log")
a.set_xlabel("anty-EBF2  (jednostki arbitralne, skala log)",color=INK,fontsize=10)
a.set_ylabel("SSC-A",color=INK,fontsize=10)
a.axvline(1000*(2.37+0.3),color=INK,lw=1.5,ls=(0,(5,3)))
a.text(1000*(2.37+0.3)*1.06,ssc.min()*1.1,"próg 2.37\n(optymalny)",fontsize=8.5,color=INK)
a.legend(frameon=False,fontsize=9,labelcolor=INK,markerscale=2.5,loc="upper left")
a.set_title(f"Symulowany dot plot bramki R3, n={N}\n(kolor = prawda, której na cytometrze NIE widzisz)",
            color=INK,fontsize=11.5,loc="left",pad=10)

a=ax[1]
bins=np.logspace(np.log10(sig.min()*0.9),np.log10(sig.max()*1.1),60)
a.hist(sig[pp],bins=bins,color=BLUE,alpha=.6,label="Aldoc+")
a.hist(sig[~pp],bins=bins,color=ORANGE,alpha=.6,label="Aldoc−")
a.hist(sig,bins=bins,histtype="step",color=INK,lw=2,label="co zobaczysz (suma)")
a.set_xscale("log")
a.axvline(1000*(2.37+0.3),color=INK,lw=1.5,ls=(0,(5,3)))
a.set_xlabel("anty-EBF2 (skala log)",color=INK,fontsize=10)
a.set_ylabel("liczba jąder",color=INK,fontsize=10)
a.legend(frameon=False,fontsize=9,labelcolor=INK)
a.set_title("Rzut na oś EBF2: czarny obrys to\nwszystko, co widzi cytometr",
            color=INK,fontsize=11.5,loc="left",pad=10)
fig.text(0.008,0.012,"SYMULACJA na danych Kozarevy (16 634 PC), n=3152 jak w R3/cerebellum-1-test. "
        "Zalozenie: sygnal bialka ∝ mRNA, zero tla i szumu — to GORNA GRANICA, nie prognoza.",
        fontsize=8,color=MUTED)
plt.tight_layout(rect=[0,0.035,1,1])
out=f"{P}/figures/symulacja_bramka_ebf2.png"
plt.savefig(out,dpi=160,facecolor=fig.get_facecolor())
print("zapisano",out)

t=2.37
hi=vv>=t
print(f"\n=== przy progu {t} CP10K, n={N} ===")
print(f"  bramka WYSOKA: {hi.sum():4d} zdarzen, w tym Aldoc- {(~pp & hi).sum():4d} -> czystosc {(~pp & hi).sum()/hi.sum()*100:.1f}%")
print(f"  bramka NISKA : {(~hi).sum():4d} zdarzen, w tym Aldoc+ {(pp & ~hi).sum():4d} -> czystosc {(pp & ~hi).sum()/(~hi).sum()*100:.1f}%")
for n in (3152,948,580):
    k=rng.choice(len(v),size=n,replace=False)
    h=v[k]>=t
    print(f"  n={n:4d}: bramka wysoka {h.sum():4d} zdarzen, niska {(~h).sum():4d}")
