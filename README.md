# pc-zebrin-sorting

Analiza danych snRNA-seq móżdżku myszy pod kątem znalezienia markera jądrowego,
który pozwoli rozdzielić komórki Purkinjego na Aldoc dodatnie i Aldoc ujemne
na sorterze przepływowym.

**Wynik główny:** gen `Ebf2` różnicuje obie grupy (AUC 0,932), jest czynnikiem
transkrypcyjnym wiążącym DNA, więc przetrwa izolację jąder w 1 % Tritonie bez
utrwalania, i jest praktycznie nieobecny we wszystkich pozostałych typach komórek
móżdżku.

Pełny opis w [`raport/main.pdf`](raport/main.pdf).

---

## Warstwa Purkinjego w 3D

![Warstwa Purkinjego myszy — skład fenotypowy płacików](figures/film_mozdzek.gif)

Warstwa komórek Purkinjego całego móżdżku myszy, obracana wokół osi pionowej.
Każdy płacik zapala się po kolei, pokolorowany podtypami Purkinjego, które
zawiera; panel po prawej podaje jego skład.

| | |
|---|---|
| geometria | **2 279 886 wokseli** przy 10 µm |
| fenotypy | **16 634 komórki**, 9 podtypów, 16 płacików |
| render | 1920×1080, 60 kl./s, 40 s |

Warstwy Purkinjego nie ma w wydanym wolumenie adnotacji CCFv3a — wszystkie 19
struktur `*pu` mają zero wokseli, podczas gdy `*gr` i `*mo` mają miliony.
Jest więc wyliczona geometrycznie, jako granica warstwy ziarnistej i drobinowej,
bo tam z definicji leżą ciała tych komórek.

- pełny film, 40 s: [`figures/film_mozdzek_web.mp4`](figures/film_mozdzek_web.mp4)
- model interaktywny (Three.js): [`figures/purkinje3d_v1.html`](figures/purkinje3d_v1.html)
- atrybucja i DOI: [`figures/film_mozdzek_CREDITS.txt`](figures/film_mozdzek_CREDITS.txt)

Zrobione w Pythonie (NumPy, Pillow) i `ffmpeg` — własny rasteryzer programowy,
bez GPU i bez silnika graficznego. Skrypt: [`scripts/58_film_mozdzek.py`](scripts/58_film_mozdzek.py).

**Źródła danych**
Kozareva V. i wsp. *Nature* 2021; 598:214–219 · [10.1038/s41586-021-03220-z](https://doi.org/10.1038/s41586-021-03220-z) · GEO [GSE165371](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165371)
Piluso S. i wsp. *Imaging Neuroscience* 2025; 3:imag_a_00565 · [10.1162/imag_a_00565](https://doi.org/10.1162/imag_a_00565) · Zenodo [10.5281/zenodo.15176439](https://doi.org/10.5281/zenodo.15176439) · CC BY 4.0
Oparte na Allen Mouse Brain Common Coordinate Framework (CCFv3), Allen Institute for Brain Science.

---

## Dane wejściowe

| co | skąd |
|---|---|
| Atlas móżdżku myszy, 611 034 jąder | [GEO GSE165371](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165371) (Kozareva i wsp., *Nature* 2021) |
| Metadane z przypisaniem podtypów | [MacoskoLab/cerebellum-atlas-analysis](https://github.com/MacoskoLab/cerebellum-atlas-analysis) |
| Adnotacja lokalizacji białek | [UniProt](https://rest.uniprot.org), Swiss-Prot mysz, 17 283 rekordy |
| Modele przestrzenne białek | [AlphaFold DB](https://alphafold.ebi.ac.uk) |

Surowych danych GEO (3,1 GB) i plików `.h5ad` (825 MB) nie ma w repozytorium.
Odtwarza je potok `00` do `02` w około dwie minuty.

## Jak odtworzyć

```bash
# 1. Pobierz z GEO trzy pliki do raw/
#    cb_adult_mouse.mtx.gz, cb_adult_mouse_barcodes.txt, cb_adult_mouse_genes.txt

# 2. Sklonuj repozytorium autorów atlasu (metadane)
git clone https://github.com/MacoskoLab/cerebellum-atlas-analysis.git

# 3. Zbuduj podzbiór Purkinjego
python scripts/00_prep_map.py      # mapowanie kodów kreskowych, 12 przemianowanych próbek
./scripts/01_extract.sh            # strumieniowa ekstrakcja kolumn z .mtx.gz
python scripts/02_build_h5ad.py    # złożenie + walidacja per komórka

# 4. Główne wyniki
python scripts/33_tabela_genow_lokalizacja.py   # tabela 24 409 genów + lokalizacja UniProt
python scripts/44_audyt_ebf2.py                 # weryfikacja twierdzenia o EBF2
python scripts/43_symulacja_bramki_ebf2.py      # symulacja bramki cytometrycznej
```

Środowisko: Python 3.12, scanpy 1.12, anndata 0.13, numpy 2.5, scipy 1.18,
scikit-learn 1.9. PyMOL potrzebny tylko do rycin strukturalnych.

## Struktura

```
scripts/     potok analityczny, numerowany w kolejności uruchamiania
raport/      dokument LaTeX (szablon sleek) + PDF
figures/     ryciny
references/struktury/   modele AlphaFold użyte w analizie strukturalnej
processed/   tabele wynikowe (CSV, JSON)
```

## Najważniejsze skrypty

| skrypt | co robi |
|---|---|
| `00_prep_map.py` | łączy kody kreskowe GEO z metadanymi; obchodzi pułapkę 12 przemianowanych próbek, która gubi 7 % komórek Purkinjego |
| `01_extract.sh` | wycina kolumny Purkinjego z macierzy 1,03 mld wpisów, strumieniowo |
| `02_build_h5ad.py` | składa `purkinje_cells_v2.h5ad` (16 634 × 24 409) i waliduje per komórkę |
| `33_tabela_genow_lokalizacja.py` | profil wszystkich 24 409 genów w 9 podtypach + lokalizacja subkomórkowa z UniProt |
| `34_topologia_otoczki.py` | rozstrzyga, po której stronie błony jądrowej leży epitop |
| `35`, `37`, `36`, `38` | ekstrakcja wybranych genów z pełnego atlasu 611 034 jąder i profil w 18 typach komórek |
| `39_ebf2_struktura_fig.py` | pLDDT i podobieństwo EBF2 do EBF1 wzdłuż sekwencji |
| `43_symulacja_bramki_ebf2.py` | symulacja wykresu punktowego dla bramki R3 |
| `44_audyt_ebf2.py` | cztery niezależne próby obalenia twierdzenia o EBF2 |
| `45_fig_epitop_czarne.pml` | rycina strukturalna do raportu |

## Uwaga o etykietach

Podział na Aldoc dodatnie i ujemne pochodzi z nazw podtypów nadanych przez
autorów atlasu, nie z progu ustawionego w tej analizie. Ograniczenia takiego
podejścia opisuje sekcja 2.5 raportu.
