"""Klasyfikacja topologiczna 222 genow otoczkowych: czy epitop jest DOSTEPNY dla
przeciwciala na jadrach NIEUTRWALONYCH i NIEPERMEABILIZOWANYCH.

Kryterium: epitop musi byc po stronie cytoplazmatycznej ZEWNETRZNEJ blony jadrowej (ONM).
  ONM        - nucleus outer membrane / cytoplasmic side          -> DOSTEPNY
  POR_CYT    - nuclear pore, strona cytoplazmatyczna              -> czesciowo dostepny
  INM        - nucleus inner membrane / nucleoplasmic side/lamina -> NIEDOSTEPNY
  NIEOKR     - "nucleus membrane"/"envelope" bez podania strony   -> nierozstrzygniete

Zrodlo tekstu: references/uniprot_mouse_subcell.tsv (PELNY, nie uciety loc_txt).
"""
import pandas as pd, numpy as np, re
P="/home/mikolajrurad/omics-data-highway/ichb-purkinje-analysis/pc-zebrin-sorting"
pd.set_option("display.width",260); pd.set_option("display.max_colwidth",120)

u=pd.read_csv(f"{P}/references/uniprot_mouse_subcell.tsv",sep="\t")
u.columns=["gen","syn","loc"]; u["loc"]=u["loc"].fillna("")
full={}
for _,r in u.iterrows():
    if isinstance(r["gen"],str): full.setdefault(r["gen"].lower(),r["loc"])
    if isinstance(r["syn"],str):
        for s in r["syn"].split(): full.setdefault(s.lower(),r["loc"])

G=pd.read_csv(f"{P}/processed/geny_purkinje_lokalizacja.csv")
O=G[G.lokalizacja=="otoczka"].copy()
O["loc_pelny"]=[full.get(g.lower(),"") for g in O.gen]
print(f"genow otoczkowych: {len(O)}   z pelnym tekstem: {(O.loc_pelny!='').sum()}")
print(f"tekst dluzszy niz 400 znakow (uciety w tabeli): {(O.loc_pelny.str.len()>400).sum()}")

ONM=re.compile(r"nucleus outer membrane|cytoplasmic side")
INM=re.compile(r"nucleus inner membrane|nucleoplasmic side|nucleus lamina")
POR=re.compile(r"nuclear pore")
def kat(t):
    tl=t.lower()
    o,i,p=bool(ONM.search(tl)),bool(INM.search(tl)),bool(POR.search(tl))
    if o and not i: return "ONM_dostepny"
    if o and i:     return "ONM+INM_mieszany"
    if i:           return "INM_niedostepny"
    if p:           return "POR_niepewny"
    return "NIEOKRESLONY"
O["topologia"]=O.loc_pelny.apply(kat)
print("\n=== TOPOLOGIA 222 GENOW OTOCZKOWYCH ===")
print(O.topologia.value_counts().to_string())

O.to_csv(f"{P}/processed/otoczka_topologia.csv",index=False)
print(f"\nzapisano processed/otoczka_topologia.csv")

print("\n=== KANDYDACI: ONM_dostepny lub ONM+INM, posortowane po wykrywalnosci ===")
k=O[O.topologia.str.startswith("ONM")].sort_values("pct_komorek",ascending=False)
print(k[["gen","topologia","pct_komorek","srednia_cp10k_all","mean_Aldoc_plus","mean_Aldoc_minus",
         "log2FC_plus_vs_minus","marker_dla"]].to_string(index=False))

print("\n=== ci sami kandydaci: pelny opis UniProt ===")
for _,r in k.head(12).iterrows():
    print(f"\n{r.gen}  [{r.topologia}, wykryty w {r.pct_komorek}% komorek]")
    print("   "+r.loc_pelny[:400].replace("SUBCELLULAR LOCATION:","").strip())

print("\n=== dla porownania: co odpada jako INM (top 10 po wykrywalnosci) ===")
i=O[O.topologia=="INM_niedostepny"].sort_values("pct_komorek",ascending=False).head(10)
print(i[["gen","pct_komorek","log2FC_plus_vs_minus","marker_dla"]].to_string(index=False))
