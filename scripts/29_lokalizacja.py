"""Przypisanie lokalizacji subkomorkowej wszystkim genom z markers_all.csv.

Zrodlo: UniProt/Swiss-Prot, mysz (organism_id:10090, reviewed:true), pole
cc_subcellular_location. Kuratorowane, nie predykcja.

Kategorie sa podyktowane eksperymentem (sortowanie WYIZOLOWANYCH JADER):
  1 otoczka  - blona/otoczka jadrowa, por jadrowy, lamina -> dziala na jadrach
               NIEUTRWALONYCH, bez permeabilizacji
  2 jadro    - nukleoplazma, chromatyna, jaderko -> wymaga utrwalenia+permeabilizacji
  3 reszta   - cytoplazma/blona kom./wydzielnicze -> odpada przy izolacji jader
"""
import re, sys, pandas as pd

UNI = sys.argv[1]; MARK = sys.argv[2]; OUT = sys.argv[3]

u = pd.read_csv(UNI, sep="\t")
u.columns = ["gen", "syn", "loc"]
u["loc"] = u["loc"].fillna("")

EVID = re.compile(r"\{[^}]*\}")
NOTE = re.compile(r"Note=.*?(?=(?:SUBCELLULAR LOCATION:)|$)", re.S)

def terms(txt):
    """Wyciaga same terminy lokalizacji: bez kodow ECO i bez prozy w Note=."""
    t = EVID.sub("", txt)
    t = NOTE.sub("", t)
    t = t.replace("SUBCELLULAR LOCATION:", " ")
    t = re.sub(r"\[[^\]]*\]:", " ", t)          # nazwy izoform/lancuchow
    return t.lower()

OTOCZKA = re.compile(
    r"nucleus (inner |outer )?membrane|nucleus envelope|nucleus lamina|"
    r"nuclear pore|nucleus, nuclear pore")
JADRO   = re.compile(r"\bnucleus\b|\bchromosome\b|\bnucleolus\b|\bnucleoplasm\b")
INNE    = re.compile(
    r"\bcytoplasm\b|cell membrane|\bsecreted\b|mitochondri|endoplasmic reticulum|"
    r"golgi|cell projection|cell surface|lysosome|endosome|peroxisome|"
    r"cytoskeleton|cell junction|synapse|\bvesicle\b|apical|basolateral")

rows = []
for _, r in u.iterrows():
    t = terms(r["loc"])
    oto = bool(OTOCZKA.search(t)); jad = bool(JADRO.search(t)); inn = bool(INNE.search(t))
    kat = "otoczka" if oto else ("jadro" if jad else ("reszta" if t.strip() else "brak"))
    rows.append((r["gen"], r["syn"], kat, jad, oto, inn, r["loc"][:300]))
ann = pd.DataFrame(rows, columns=["gen","syn","kategoria","w_jadrze","otoczka","tez_poza","loc_txt"])

# indeks: nazwa glowna + synonimy, malymi literami
idx = {}
for _, r in ann.iterrows():
    if isinstance(r["gen"], str):
        idx.setdefault(r["gen"].lower(), r)
    if isinstance(r["syn"], str):
        for s in r["syn"].split():
            idx.setdefault(s.lower(), r)

m = pd.read_csv(MARK)
def look(g):
    r = idx.get(str(g).lower())
    return pd.Series([r["kategoria"], r["w_jadrze"], r["otoczka"], r["tez_poza"], r["loc_txt"]]
                     if r is not None else ["nieznane", False, False, False, ""])
m[["kategoria","w_jadrze","otoczka","tez_poza","loc_txt"]] = m["gen"].apply(look)
m["absd"] = m["d"].abs()
m["kierunek"] = m["d"].apply(lambda v: "Aldoc+" if v > 0 else "Aldoc-")
m["tylko_jadro"] = m["w_jadrze"] & ~m["tez_poza"]
m.to_csv(OUT, index=False)

print(f"UniProt: {len(u)} rekordow, unikalnych kluczy w indeksie: {len(idx)}")
print(f"markers_all: {len(m)} genow\n")
print("=== pokrycie adnotacja ===")
print(m["kategoria"].value_counts().to_string())
print(f"\nstara flaga 'jadrowy' (49 genow) vs nowa 'w_jadrze':")
print(pd.crosstab(m["jadrowy"], m["w_jadrze"]).to_string())

for kat, tytul in [("otoczka","KATEGORIA 1: OTOCZKA/BLONA JADROWA (jadra nieutrwalone)"),
                   ("jadro",  "KATEGORIA 2: JADRO (wymaga utrwalenia + permeabilizacji)")]:
    sub = m[m["kategoria"]==kat].sort_values("absd", ascending=False).head(20)
    print(f"\n=== {tytul} — top 20 wg |d| (Kozareva) ===")
    print(sub[["gen","kierunek","f_pos","f_neg","d","m_pos","m_neg","tylko_jadro"]].to_string(index=False))

print("\n=== ile genow jadrowych ma |d| powyzej progu ===")
for th in (0.7, 0.5, 0.4, 0.3, 0.2):
    n1 = ((m.kategoria=="otoczka") & (m.absd>=th)).sum()
    n2 = ((m.kategoria=="jadro")   & (m.absd>=th)).sum()
    print(f"  |d| >= {th}:  otoczka {n1}, jadro {n2}")
