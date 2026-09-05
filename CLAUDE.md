# pc-zebrin-sorting — kontekst projektu

Wszystkie liczby w tym pliku pochodzą z komend uruchomionych na tej maszynie
2026-09-03. Przy każdej jest źródło. Jeśli coś się zmieni — przelicz, nie zgaduj.

---

## 1. Cel naukowy

Charakterystyka fenotypów komórek Purkinjego (PC) w móżdżku myszy w osi
**Aldoc+ / Aldoc−** (zebrina II), na bazie snRNAseq + spatial transcriptomics,
z docelowym przełożeniem na eksperyment mokry (ICHB PAN).

---

## 2. Sprzęt — stacja robocza `mikolajrurad`

| Element | Wartość | Źródło |
|---|---|---|
| CPU | Intel i5-10400F, 6 rdzeni / 12 wątków, max 4.30 GHz | `lscpu` |
| RAM | **15.54 GiB** (`MemTotal 16290104 kB`) | `/proc/meminfo` |
| Swap | 4 GiB, **zram** (`/dev/zram0`) — kompresowany w RAM, **nie ratuje przed OOM** | `swapon --show` |
| GPU | NVIDIA RTX 3060 Ti, **8192 MiB VRAM**, sterownik 610.57.04, CUDA UMD 13.3 | `nvidia-smi` |
| Dysk `/` i `/home` | btrfs, 464 G, **74 G wolne (84 % zajęte)** | `df -hT` |
| `/mnt/data1t` | ext4, 916 G, 488 G wolne | `df -hT` |
| `/mnt/data500g` | ntfs3, 466 G, 373 G wolne | `df -hT` |
| `/mnt/macnas` | CIFS 192.168.1.19, 3.6 T, 1.4 T wolne | `df -hT` |
| OS | Arch Linux, kernel 7.1.9-arch1-2, Hyprland | `hostnamectl` |

**Konsekwencja:** to nie jest maszyna do trzymania w RAM całego atlasu.
Reguły w sekcji 5 nie są sugestią.

---

## 3. Środowisko obliczeniowe

**Zbudowane 2026-09-03.** Wcześniej na maszynie nie było NIC — conda/mamba
usunięte, `~/miniforge_envs/` nie istnieje, systemowy Python 3.14.7 miał tylko
`numpy` i `matplotlib`.

### venv

```
/home/mikolajrurad/omics-data-highway/.venv     Python 3.12.14, 1.7 GB
```

Zbudowany **`uv` 0.12.9** (zainstalowany bez roota do `~/.local/bin`;
`sudo` na tej maszynie wymaga hasła, więc pacman odpadł). `uv` ściągnął
własnego CPythona 3.12.14 — systemowy 3.14 jest za nowy dla części scverse.

```bash
export PATH="$HOME/.local/bin:$PATH"
cd ~/omics-data-highway
source .venv/bin/activate          # albo wprost: .venv/bin/python
uv pip install --python .venv/bin/python <pakiet>
```

### Co jest w środku (zweryfikowane `importlib.metadata`)

| | | | |
|---|---|---|---|
| scanpy 1.12.4 | anndata 0.13.3 | h5py 3.16.0 | zarr 3.3.0 |
| numpy 2.5.2 | scipy 1.18.1 | pandas 3.0.5 | statsmodels 0.15.0 |
| leidenalg 0.12.0 | igraph 1.0.0 | umap-learn 0.5.12 | scikit-learn 1.9.0 |
| numba 0.67.0 | pyarrow 25.0.1 | polars 1.44.1 | pynrrd 1.1.3 |
| matplotlib 3.11.1 | seaborn 0.13.2 | harmonypy 2.0.0 | decoupler 2.2.0 |
| **napari 0.9.0** | **brainglobe-atlasapi 3.0.1** | **brainrender-napari 0.2.0** | vispy 0.16.2 |
| jupyterlab 4.6.3 | ipywidgets | | |

Import `scanpy`, `anndata`, `h5py`, `nrrd`, `BrainGlobeAtlas` — sprawdzony,
działa.

### Czego świadomie NIE MA

**`torch`, `scvi-tools`, `cell2location`, `squidpy`, `celltypist`, `scrublet`.**
Stack GPU to kolejne kilka GB i wiązanie z konkretną wersją CUDA
(sterownik: 610.57.04, CUDA UMD 13.3, RTX 3060 Ti **8 GiB VRAM**, z czego
~960 MiB zjada Hyprland). Dokładać dopiero, gdy będzie potrzebny — i wtedy
pilnować, żeby wersja CUDA torcha pasowała do sterownika.

**R nadal ma tylko 29 pakietów bazowych — bez Seurata.** Pliki `.rds`
(np. `mouse.Purkinje.rds` z CBMSTA, 5.64 GiB) są tym samym nieczytelne.

### Pliki, którym NIE WOLNO wierzyć

| Plik | Status |
|---|---|
| `TOOLS.md` | **NIEAKTUALNY W CAŁOŚCI.** Opisuje envy `pharma`/`genometools` pod `~/miniforge_envs/` — ścieżka nie istnieje. Każda wersja pakietu w tym pliku to fikcja. Zastąpiony przez sekcję powyżej. |
| `pharma.sh` | symlink do modułu bash; sam przyznaje w komentarzu, że conda usunięta. Zostały `drugview`/`pdbview` wołające `pymol` z PATH. Zero związku ze scRNA-seq. |
| `genometools.sh` | czysty bash: `reverse`/`switch`/`revcomp`. Żadnego samtools/blast/seqkit. |

## 4. `ssh eagle` — co się dzieje

```
Host eagle
    HostName eagle.man.poznan.pl
    User     mateuszmac
    IdentityFile ~/.ssh/pcss_ed25519
    ControlMaster auto / ControlPersist 10m   → drugie połączenie idzie tunelem, bez uwierzytelniania
```

Lądujesz na **węźle logowania klastra Eagle w PCSS (Poznań)**:

| | |
|---|---|
| OS | AlmaLinux 9.8 |
| użytkownik | `mateuszmac`, uid 67751201 |
| grupy grantowe | **`pl0796-01`**, **`pl0845-02`** |
| scheduler | **Slurm** (`sbatch`, `squeue`, `sinfo` w `/usr/bin`) |
| konta Slurm | `pl0796-01` i `pl0845-02`, QOS: `normal`, `tesla` |
| moduły | Environment Modules 5.3.0 — **`module` to funkcja powłoki, widoczna dopiero po `bash -l`.** Przez `ssh eagle '<cmd>'` (bez logowania) `module` NIE ISTNIEJE. |

### Partycje (`sinfo`)

| Partycja | Limit czasu | Uwagi |
|---|---|---|
| `interactive` | 10:00:00 | węzły 192 GB i 384 GB RAM |
| `standard` (domyślna) | 7-00:00:00 | węzły 192/384 GB, część z `local_ssd` |

Węzły: Intel Cascade Lake, EDR InfiniBand, **192 GB lub 384 GB RAM**.
GRES pokazuje `(null)` na tych partycjach — **partycja GPU nie została
zweryfikowana**, choć QOS `tesla` sugeruje, że istnieje. NIEPOLICZONE:
sprawdzić `sinfo -o "%P %G"` na wszystkich partycjach, zanim zaplanujesz
cokolwiek na GPU.

### Storage na eagle

| Ścieżka | FS | Uwaga |
|---|---|---|
| `$HOME` = `/mnt/storage_3/home/mateuszmac` | NFS (`pacyfik-nfs:/home`) | **`df` na `$HOME` pokazuje 1.0 G / 0 użyte** — to limit per-user, nie miejsce na dane. Katalog jest praktycznie pusty (`.ssh`, `.bash_history`). |
| `~/pl0796-01`, `~/pl0845-02` | dowiązania do katalogów grantowych, `dr-xr-x---` (tylko odczyt dla usera) | |
| `/mnt/storage_5/scratch/pl0796-01` | **Lustre**, 6.5 P (1.3 P użyte) | katalog grupowy, `drwxrws---`, zużycie ~141 MB / 278 826 plików |
| `/mnt/storage_5/scratch/pl0845-02` | Lustre | praktycznie pusty (4 k, 1 plik) |
| `/mnt/storage_2/scratch`, `/mnt/storage_2/project_data` | Lustre 6.4 P | |
| `/mnt/storage_6/project_data` | NFS 14 P | |
| `/mnt/storage_4/archive` | Ceph 8 P | archiwum |

**Dane liczbowe idą na `/mnt/storage_5/scratch/pl0796-01/`, NIGDY do `$HOME`.**

### Software na eagle

`module -t avail` → **84 moduły**. Istotne:

```
python/3.13.0-gcc-14.2.0
r/4.4.1-gcc-14.2.0
samtools/1.2-gcc-14.2.0          (wersja 1.2 — antyk)
cuda/12.4.1  12.6.0  12.8.0  13.2.1
amber/2025
```

**Brak `conda`, `mamba`, `apptainer`, `singularity` w PATH** (sprawdzone
`command -v`). Jest tylko `/usr/bin/python3`. Czyli: żeby uruchomić scanpy
na eagle, trzeba zbudować venv na module `python/3.13.0` w katalogu scratch,
albo załatwić kontener. **To jest otwarta kwestia — nie zakładaj, że da się
`pip install` z węzła obliczeniowego (dostęp do sieci na węzłach niesprawdzony).**

### Kiedy używać eagle

Zasada: **wszystko, co dotyka pełnej macierzy 611 034 komórek, idzie na eagle.**
Węzeł 384 GB robi w RAM to, czego ta stacja nie zrobi nigdy (patrz sekcja 5).
Lokalnie zostają: podzbiór Purkinje, wykresy, pisanie kodu.

---

## 5. BUDŻET PAMIĘCI — reguły twarde

### Rozmiar danych (zmierzony, nie szacowany)

Nagłówek `raw/cb_adult_mouse.mtx.gz` (`zcat | head -2`):

```
24409 611034 1026400102          ← genów × komórek × niezerowych
```

Gęstość **6.88 %**. Weryfikacja krzyżowa: suma kolumny `nGene` po komórkach
non-REMOVED w `full_cb_metadata.csv` = **1 026 400 606**, czyli o **504** więcej
niż nnz w mtx (0.000049 %). Zgodność potwierdza, że macierz GEO = 611 034
komórek po QC.

### Ile to zajmuje w RAM

| Obiekt | Reprezentacja | Rozmiar |
|---|---|---|
| **Cała macierz** | CSR float32 + indices int32 | **7.65 GiB** |
| Cała macierz | CSR float64 + int64 (domyślne, jeśli nie wymusisz!) | **15.30 GiB** |
| Cała macierz | gęsta float32 | **55.56 GiB** |
| Same Granule (477 176 kom.) | CSR f32/i32 | **5.20 GiB** |
| **Same Purkinje (16 634 kom., nnz 78 152 971)** | CSR f32/i32 | **0.58 GiB** |
| Same Purkinje | gęsta float32 (24409 × 16634) | **1.51 GiB** |

Dostępny RAM w chwili pomiaru: **9.3 GiB** (`free -h`).

### Reguły

1. **NIGDY `sc.read_10x_mtx()` / `scipy.io.mmread()` na `raw/cb_adult_mouse.mtx.gz`
   na tej maszynie.** `mmread` buduje COO w float64 + dwa int64 → 15+ GiB tylko
   na macierz, plus kopia przy konwersji do CSR. To gwarantowany OOM,
   a zram tego nie uratuje.
2. **Wczytuj podzbiór, nie całość.** Purkinje to 2.72 % komórek i 7.61 % nnz.
   Filtrowanie robi się **strumieniowo po pliku `.mtx.gz`** (awk/parser
   liniowy po indeksach kolumn), a nie przez wczytanie i `adata[mask]`.
3. **Zawsze `dtype=np.float32`, indeksy `int32`.** float64+int64 podwaja koszt.
4. **`anndata.read_h5ad(..., backed='r')`** dla oglądania `.obs`/`.var`.
   Do wczytania w pamięć — dopiero po sprawdzeniu `X.shape` i `nnz`.
5. **Twardy limit: żaden obiekt w RAM > 2.33 GiB** (25 % dostępnego RAM),
   co odpowiada **max ~312 mln nnz** przy f32/i32. Powyżej tego — eagle.
6. **Zanim wczytasz cokolwiek, policz.** `nnz × (4+4) + (n_cells+1) × 4`.
   Jeśli wynik > 2.33 GiB → nie wczytujesz, tylko mówisz i proponujesz eagle.
7. **GPU ma 8 GiB VRAM i dzieli je z Hyprlandem** (~960 MiB zajęte w spoczynku).
   Realny budżet dla scvi-tools ≈ 6.5 GiB. Batch size dobierany, nie domyślny.
8. **Nie ma miejsca na dysku.** `/home` ma 74 G wolnego przy 84 % zajętości.
   Duże pliki pośrednie → `/mnt/data1t` (488 G wolne), nie do `processed/`.

---

## 6. Dane w projekcie

### Proweniencja danych — co to właściwie jest

**Kozareva V, Martin C, Osorno T, Rudolph S, Guo C, Vanderburg C, Nadaf N,
Regev A, Regehr WG, Macosko E. „A transcriptomic atlas of mouse cerebellar
cortex comprehensively defines cell types". *Nature* 2021; 598(7879): 214–219.**
DOI `10.1038/s41586-021-03220-z` · PMID `34616064` · PMCID `PMC8494635`
· GEO `GSE165371` · Single Cell Portal `SCP795`
(dostęp do źródeł: 2026-09-03)

Errata (*Author Correction*, `10.1038/s41586-021-04373-7`, PMID `35022615`)
dotyczy **wyłącznie literówki w numerze akcesyjnym GEO** w sekcji Data
availability — poprawiono na GSE165371. Nie zmienia danych ani analiz.

Metoda: **snRNA-seq, 10x Chromium V3**, jądra z kory móżdżku dorosłej myszy,
**16 płacików**: I, II, III, CUL, VI, VII, VIII, IX, X, AN1, AN2, PRM, SIM,
COP, F, PF.

#### Zgodność publikacji z plikami na dysku (weryfikacja krzyżowa)

| Wielkość | Publikacja | Zmierzone lokalnie | Zgoda |
|---|---|---|---|
| jądra pozyskane | 780 553 | 792 365 wierszy `full_cb_metadata.csv` | ✗ patrz niżej |
| jądra po QC | **611 034** | **611 034** (`cluster_metadata.csv`, `barcodes.txt`, kolumny mtx) | ✓ |
| Purkinje | **16 634** | **16 634** | ✓ |
| Granule | **477 176** | **477 176** | ✓ |
| klastry PC | **9** (7 Aldoc+, 2 Aldoc−) | **9** (`Aldoc_1..7` + `Anti_Aldoc_1..2`) | ✓ |
| mediana UMI/profil | 2 862 | **2 750** (mediana `nUMI`, n=611 034) | ✗ patrz niżej |
| płeć | 530 063 M / 250 490 F | 431 221 M / 179 813 F (po QC) | spójne |

**Dwie rozbieżności, obie wyjaśnialne, żadna nie zweryfikowana do końca:**

1. **792 365 vs 780 553.** README repo mówi: „780,553 + 11,812 FACS sorted
   putative Purkinje profiles". 780 553 + 11 812 = **792 365** — dokładnie tyle,
   ile wierszy ma `full_cb_metadata.csv`. Czyli plik zawiera atlas **plus**
   jądra Purkinje sortowane FACS-em, a liczba z abstraktu dotyczy samego atlasu.
   Publikacja (fragment, który udało się pobrać) **nie wspomina o FACS** —
   informacja pochodzi wyłącznie z README repo.
2. **2 750 vs 2 862 mediany UMI.** Nasza mediana liczona po 611 034 komórkach
   PO QC; liczba z pracy prawdopodobnie po pełnym zestawie 780 553.
   **NIEPOLICZONE** — `full_cb_metadata.csv` zawiera `nUMI` także dla komórek
   `REMOVED`, więc da się to sprawdzić; nie zrobiono tego.

Liczba unikalnych próbek (`orig.ident` w `full_cb_metadata.csv`): **77**.

#### Co praca mówi o Purkinje — i co się z tym zgadza w naszych danych

- 9 klastrów PC, **7 Aldoc-dodatnich, 2 Aldoc-ujemne** → u nas dokładnie
  `Purkinje_Aldoc_1..7` (10 427 kom.) i `Purkinje_Anti_Aldoc_1..2` (6 207 kom.).
- Cytat: „Combinatorial expression of *Aldoc* and at least one subtype-specific
  marker fully identified the Purkinje clusters". Markery przykładowe podane
  w pracy: **`Gpr176`**, **`Tox2`**.
- Cytat: „Most of this PC diversity was concentrated in the posterior
  cerebellum, particularly the uvula and nodulus" (uvula = IX, nodulus = X).
  Zgadza się z naszym rozkładem: `Aldoc_5` dominuje w IX (699) i VIII (497).
- Dane przestrzenne w pracy: **smFISH** (walidacja markerów MLI: `Sorcs3`,
  `Nxph1`, `Gjd2`, `Cdh22`, `Lgi2`, `Sst`), **Slide-seq**, **HCR**, Allen Brain
  Atlas. **MERFISH nie występuje** w pobranym fragmencie.
  Uwaga: surowych danych Slide-seq/smFISH **nie ma w tym katalogu** — w repo są
  tylko pochodne (`data/SST_GJD_LGI_ROI_intensities.csv`, `data/coords/`).

### `raw/` — 5.9 GB

| Plik | Rozmiar | Co to |
|---|---|---|
| `cb_adult_mouse.mtx.gz` | 3 119 066 566 B | macierz 24409 × 611034, nnz 1 026 400 102 |
| `cb_adult_mouse_barcodes.txt` | 611 034 linii | kolumny macierzy |
| `cb_adult_mouse_genes.txt` | 24 409 linii | wiersze macierzy, symbole MGI |

**`GSE165371_cb_adult_mouse.tar.gz` — PRZENIESIONY DO ARCHIWUM 2026-09-03.**
Był duplikatem trzech plików powyżej (`tar -tzvf` to potwierdził).
Leży teraz w `hostbrr:backup/snrnaseq/GSE165371_cb_adult_mouse.tar.gz`
(3 123 821 257 B, mtime 2021-01-22 zachowany).
Weryfikacja: md5 round-trip **`c602b2312ab48b1ebd96bcb068aa9301`** — zgodne
przed wysłaniem i po pobraniu z powrotem. Dopiero potem skasowano lokalny.
Odzyskane miejsce: 73.7 → **76.61 GiB** wolnego na `/home` (`df -B1`).

Indeksy wierszy genów markerowych (nr linii w `genes.txt` = indeks 1-based
w mtx): `Plcb4` 2302, `Car8` 4508, `Calb1` 4548, `Grid2` 7385, `Itpr1` 7670,
`Slc1a6` 10032, `Pcp2` 10502, `Atxn7` 11511, `Nrgn` 12357, `Cck` 13157,
**`Aldoc` 13976**, `Atxn1` 14945, `Pcp4` 17317, `Plcb3` 18863.

### `processed/purkinje_cells.h5ad` — 152 624 292 B — ⚠ WADLIWY, NIE UŻYWAĆ

Otwarty h5py 2026-09-03. **Zawiera błąd, który go dyskwalifikuje.**

| | |
|---|---|
| kształt | **15 467 × 24 409** (CSC, `float32` + `int32`) |
| nnz | 72 331 479 → **0.54 GiB w RAM**, gęsta float32 1.41 GiB |
| `X` | surowe zliczenia (próbka 2 mln z 72.3 mln: same całkowite, min 1, max 756) |
| `obsm`, `varm`, `obsp`, `varp`, `layers`, `uns` | **wszystkie puste** — brak PCA/UMAP/sąsiedztw |
| `var` | tylko kolumna `gene`, 24 409 pozycji |

**Powinno być 16 634 komórek. Jest 15 467. Brakuje 1 167 (7.0 %).**

Strata nie jest rozłożona losowo — siedzi wyłącznie w dwóch regionach,
dokładnie tych, których próbki są przemianowane między GEO a metadanymi
(sekcja o pułapce barcodów):

| region | metadane | h5ad | strata |
|---|---:|---:|---:|
| **CUL** | 592 | 136 | **−456 (77.0 %)** |
| **VI** | 2 285 | 1 574 | **−711 (31.1 %)** |
| pozostałe 14 regionów | — | — | 0 (0.0 %) |
| **razem** | **16 634** | **15 467** | **−1 167 (7.0 %)** |

To jest podpis naiwnego joinu `barcodes.txt` ↔ metadane, bez mapowania
`VI2*→DEC*` i `CUL*→V*/IV*`. **Culmen stracił trzy czwarte komórek Purkinjego.**
Każdy wniosek o przednim móżdżku z tego pliku jest nieważny.

Dodatkowo dwie kolumny `obs` są martwe: **`lobule` i `mouse` mają jedną
kategorię `NA` dla wszystkich 15 467 komórek.**

Kolumny `obs`: `nGene`, `nUMI`, `regions` (16 kat.), `percent_mito`,
`first_step_annotation` (1 kat.), `final_annotation_cluster` (1 kat.),
`final_annotation_subcluster` (9 kat.), `aldoc_group` (2 kat.:
Aldoc+ 10 037 / Aldoc− 5 430 — powinno być 10 427 / 6 207), `lobule`, `mouse`.

**PRZEBUDOWANY 2026-09-04 → `processed/purkinje_cells_v2.h5ad`.**
Ten plik (`purkinje_cells.h5ad`) zostaje jako historyczny; **nie używać**.

### `cerebellum-atlas-analysis/` — repo Macosko Lab

`git remote`: `https://github.com/MacoskoLab/cerebellum-atlas-analysis.git`,
HEAD `142415c`. Kod do Kozareva et al. 2021, *Nature*
(`s41586-021-03220-z`). Skrypty w **R** (Seurat/LIGER) — a Seurata tu nie ma.

Wartościowe pliki danych:

| Plik | Wiersze | Kolumny |
|---|---|---|
| `data/full_cb_metadata.csv` | 792 365 | `nGene, nUMI, orig.ident, res.0.1, regions, percent_mito, first_step_annotation, final_annotation_cluster, final_annotation_subcluster` |
| `data/cluster_metadata.csv` | 611 034 | `nGene, nUMI, sex, region, subcluster, cluster` |
| `data/coords/` | 16 plików CSV | współrzędne przestrzenne per płacik (I–X, AN1/2, COP, CUL, F, PF, PRM, SIM) |
| `data/var_genes_union.txt` | | geny zmienne |
| `data/SST_GJD_LGI_ROI_intensities.csv` | | smFISH |

`full_cb_metadata` ma **181 331 komórek `REMOVED`**; 792 365 − 181 331 = 611 034,
czyli dokładnie tyle, co w macierzy.

### ⚠ PUŁAPKA: nazwy próbek nie zgadzają się między GEO a metadanymi

Naiwny join `barcodes.txt` ↔ `cluster_metadata.csv` gubi **53 002 komórek
(8.7 %)**, w tym **1 167 Purkinje (7.0 % wszystkich PC)**.

Przyczyna: **12 próbek ma inną nazwę w GEO niż w metadanych repo.**
Zweryfikowane w całości (dla każdej pary porównano ZBIORY sufiksów barcodów,
0 różnic):

| w `barcodes.txt` (GEO) | w metadanych repo | n |
|---|---|---|
| `VI2a_F002` | `DECa_F002` | 8225 |
| `VI2b_F002` | `DECb_F002` | 6980 |
| `VI2c_F002` | `DECc_F002` | 4533 |
| `VI2a_M006` | `DECa_M006` | 3693 |
| `VI2b_M006` | `DECb_M006` | 3808 |
| `VI2c_M006` | `DECc_M006` | 3915 |
| `CULa_M002` | `Va_M002` | 4421 |
| `CULb_M002` | `Vb_M002` | 4707 |
| `CULc_M002` | `Vc_M002` | 5091 |
| `CULd_M002` | `Vd_M002` | 4432 |
| `CUL_M006` | `IV_M006` | 1824 |
| `VI_M006` | `Vl_M006` | 1373 |

**Każdy pipeline musi zmapować te nazwy przed joinem.**

**ROZSTRZYGNIĘTE 2026-09-03 — to nie jest sprzeczność, to dwie nomenklatury.**
Sprawdzone: dla każdej spornej próbki kolumna `region` w `cluster_metadata.csv`
ma wartość zgodną z anatomią i z nazewnictwem publikacji:

| prefiks w repo | `region` | prefiks w GEO | anatomia |
|---|---|---|---|
| `DECa/b/c_*` | `VI` | `VI2a/b/c_*` | declive = płacik VI |
| `Va/Vb/Vc/Vd_M002`, `IV_M006` | `CUL` | `CULa–d_M002`, `CUL_M006` | culmen = płaciki IV–V |
| `Vl_M006` | `VI` | `VI_M006` | `Vl` = literówka `VI` |

Kolumna `region` ma **dokładnie 16 unikalnych wartości**, identycznych z 16
płacikami wymienionymi w Kozareva et al. 2021. Czyli: `region` używa
nomenklatury publikacji, prefiks barcodu w repo używa klasycznych nazw
anatomicznych robaka, GEO używa nomenklatury publikacji.

**Do analizy przestrzennej używaj kolumny `region` / `regions`, nie prefiksu
barcodu.** Prefiks to tylko etykieta próbki.

### Subklastry Purkinje (`full_cb_metadata.csv`, cały plik)

| Subklaster | n | śr. nGene | śr. nUMI |
|---|---|---|---|
| `Purkinje_Anti_Aldoc_2` | 4136 | 4857 | 27854 |
| `Purkinje_Aldoc_3` | 3577 | 4877 | 28189 |
| `Purkinje_Aldoc_5` | 2238 | 4604 | 24736 |
| `Purkinje_Anti_Aldoc_1` | 2071 | 4068 | 21059 |
| `Purkinje_Aldoc_4` | 1802 | 4893 | 27723 |
| `Purkinje_Aldoc_6` | 1569 | 4890 | 27521 |
| `Purkinje_Aldoc_1` | 622 | 4365 | 21760 |
| `Purkinje_Aldoc_7` | 544 | 4177 | 20787 |
| `Purkinje_Aldoc_2` | 75 | 5524 | 33349 |
| **razem** | **16 634** | 4698 | |

Podział Aldoc / Anti-Aldoc: **10 427 vs 6 207**.
`Purkinje_Aldoc_2` (n=75) jest za mały do samodzielnych testów DE — nie
raportuj z niego statystyk bez zaznaczenia liczności.

Rozkład po płacikach jest silnie niejednorodny i zgodny z pasami zebriny,
np. `Aldoc_5` dominuje w IX (699) / VIII (497), `Anti_Aldoc`/`Aldoc_3`
w PF / AN2 / PRM. Dokładne tabele: `awk` po kolumnach 6, 9, 10.

### `references/Data_File_S1.xlsx`

**PLIK MA 0 BAJTÓW.** Pusty placeholder. Nie próbuj go czytać.

### `artykuły/` → symlink

`/home/mikolajrurad/rura/Knowledge highway/ichb-workflow/artykuły` — 9 PDF-ów
(izolacja jąder PC, zebrin-2, bartelt-et-al-2024, ataksyna-7/GCN5, nanopore).

### `drivelab/` → symlink

`/home/mikolajrurad/CloudHub/Kompotron/drivelab` — mount rclone dysku Google
laboratorium (15 G, 14 G użyte). Dokumenty Google (pptx/docx/xlsx) widnieją
jako **0 bajtów** — to normalne dla rclone, plik trzeba wyeksportować, a nie
czytać w miejscu. Zawiera m.in. `CUT&Tag Dicer`, `CUT&Tag Pol RNA II`,
`Antibodies／Taqmans.xlsx`, `Frozen samples -80.xlsx`, `Grant proposals`.
**To jest mount sieciowy — każdy `ls -R`, `grep -r`, `du` po nim jest wolny
i generuje ruch. Nie skanuj rekurencyjnie bez potrzeby.**

---

## 7. Czego w projekcie NIE MA

- **Danych spatial transcriptomics.** Deklarowany cel obejmuje ST, ale w
  katalogu nie ma ani jednego pliku Visium / MERFISH / Slide-seq. Jedyna
  informacja przestrzenna to `regions` (płacik) w metadanych i
  `data/coords/*.csv` z repo Macosko.
- **Środowiska Pythona.** Patrz sekcja 3.
- **Seurata / LIGER** — a repo `cerebellum-atlas-analysis` jest w całości w R
  i na nich oparte. Odtworzenie ich figur wymagałoby zbudowania stacku R.

---

## 8. Reguły pracy w tym katalogu

1. Zanim wczytasz plik do pamięci — **policz jego rozmiar w RAM** i powiedz
   liczbę. Powyżej 2.33 GiB: nie wczytujesz.
2. Nie cytuj `TOOLS.md`, `pharma.sh`, `genometools.sh` jako źródła.
3. Nie mów nic o zawartości `processed/purkinje_cells.h5ad` bez otwarcia go.
4. Każdy join barcodów uwzględnia mapowanie z sekcji 6.
5. Liczności subklastrów podawaj razem z `n` — cztery z dziewięciu mają
   n < 1000.
6. Nic nie instalujesz i nic nie pobierasz bez pytania.
7. `/home` ma 74 G wolnego. Duże pliki → `/mnt/data1t`.

---

## 9. Warstwy storage'u — gdzie co trzymać

| Warstwa | Ścieżka | Wolne | Rola |
|---|---|---|---|
| **robocza** | `/home/…/pc-zebrin-sorting` (btrfs `/`) | **74 G (84 % zajęte)** | tylko kod, metadane, małe wyniki |
| **duże pośrednie** | `/mnt/data1t` (ext4) | 488 G | macierze pośrednie, rozpakowane rzeczy |
| **obliczeniowa HPC** | `/mnt/storage_5/scratch/pl0796-01` na eagle (Lustre) | 4.9 P wolne na FS | wszystko, co dotyka pełnych 611 034 komórek |
| **ARCHIWUM** | `/home/mikolajrurad/CloudHub/Brrzeszczot` (HostBrr) | patrz niżej | dane archiwalne, zimne |
| NAS | `/mnt/macnas` (CIFS) | 1.4 T | |

### `Brrzeszczot` = HostBrr — archiwum

Stan zmierzony 2026-09-03:

```
rclone mount hostbrr:backup /home/mikolajrurad/CloudHub/Brrzeszczot \
  --daemon --vfs-cache-mode writes --ftp-concurrency 8 \
  --dir-cache-time 30s --timeout 60s
```

Remote `[hostbrr]`: **type = ftp**, host `ftpde2.hostypanel.com`, port 21,
explicit TLS. Zawartość: `obsidian/` i `snrnaseq/` — **`snrnaseq/` jest pusty**
(`rclone size hostbrr:backup/snrnaseq` → `Total objects: 0, Total size: 0 B`).

**Cztery rzeczy, o które się rozbijesz, jeśli ich nie wiesz:**

1. **`df` kłamie.** Pokazuje `1.0P size / 0 used / 0%`. To wartość zastępcza —
   backend FTP nie umie raportować zajętości:
   `rclone about hostbrr:backup` → `NOTICE: Failed to about: ... doesn't support about`.
   **Nigdy nie cytuj `df` dla tego mounta jako dostępnego miejsca.** Realny limit
   to quota konta hostingowego, nie 1 PB. Żeby poznać zużycie:
   `rclone size hostbrr:backup` (skanuje, trwa).

2. **To FTP, nie obiektowy storage — brak random access.** Nie da się otworzyć
   `.h5ad` w trybie `backed='r'` ani czytać `.mtx.gz` strumieniowo z tego mounta
   z sensowną wydajnością. **Każdy plik danych trzeba najpierw ściągnąć**
   (`rclone copy hostbrr:backup/… /mnt/data1t/…`), policzyć lokalnie, i ewentualnie
   odesłać wynik. Traktuj to jak taśmę, nie jak dysk.

3. **`--vfs-cache-mode writes` = zapis kończy się pozornie.** Plik ląduje
   najpierw w lokalnym cache, a wysyłka na FTP idzie **w tle**. `cp` który
   „wrócił" nie oznacza, że dane są na serwerze. Weryfikacja wyłącznie przez
   `rclone size` / `rclone lsl` na remote, nie przez `ls` na mouncie.

4. **Cache zjada `/home`, który ma 74 G wolnego.** `~/.cache/rclone` waży już
   **633 MB** i leży na tym samym btrfs co projekt. Wysłanie 3 GB archiwum
   najpierw zapisze 3 GB na dysku, który jest w 84 % pełny. Przy dużych
   transferach: `rclone copy` bezpośrednio źródło→remote (**omija mount i VFS
   cache**), zamiast `cp` przez katalog mounta.

### ⚠ PUŁAPKA: multi-thread download z tego FTP nie działa

Zmierzone 2026-09-03. `rclone copy` **z** remote'a domyślnie używa
multi-thread dla dużych plików i na tym serwerze **pada**:

```
ERROR : multi-thread copy: failed to open source: open: write tcp
        192.168.1.60:44704->5.175.233.4:21: write: broken pipe
Attempt 1/3, 2/3, 3/3 — wszystkie nieudane
```

Trwało 12 minut (21:07:31 → 21:19:28) i skończyło się zerem.
**Poprawnie** — jednym strumieniem:

```bash
rclone copy hostbrr:backup/<plik> <cel>/   --multi-thread-streams 0 --transfers 1 --checkers 1 --low-level-retries 20
```

Ta sama zawartość zeszła wtedy w **83 s** (21:19:50 → 21:21:13), 2.909 GiB.
Upload multi-thread nie dotyczy — idzie jednym strumieniem i działa
(zmierzone: 2.909 GiB w 8 min 12 s, średnio 5.96 MiB/s).

### Procedura archiwizacji (sprawdzona)

```bash
md5sum <plik>                                    # 1. suma przed
rclone copy <plik> hostbrr:backup/snrnaseq/      # 2. upload (omija mount!)
rclone lsl hostbrr:backup/snrnaseq/              # 3. rozmiar na remote
rclone copy hostbrr:backup/snrnaseq/<plik> /mnt/data1t/tmp/   --multi-thread-streams 0 --transfers 1         # 4. round-trip
md5sum /mnt/data1t/tmp/<plik>                    # 5. suma po — MUSI się zgadzać
rm <plik> /mnt/data1t/tmp/<plik>                 # 6. dopiero teraz kasuj
```

Krok 4–5 jest **obowiązkowy**: FTP nie ma hashy, więc `rclone check`
porównuje tylko rozmiary, a to nie jest dowód integralności.

Uwaga: po `rm` na btrfs `df -h` przez chwilę pokazuje stan sprzed —
miejsce zwalnia się asynchronicznie. Sprawdzaj `df -B1` albo
`btrfs filesystem usage /`, nie `df -h` od razu.

### Reguła

Archiwum jest **zimne i jednokierunkowe**: idą tam rzeczy skończone i
odtwarzalne — surowe archiwa GEO, stare wersje przetworzonych obiektów,
wyniki zamkniętych analiz. **Nie jest to storage roboczy** i nic w pipeline
nie może z niego czytać w trakcie liczenia.

Pierwsza rzecz tam przeniesiona: `GSE165371_cb_adult_mouse.tar.gz` (2026-09-03).

---

## 10. Kandydat na warstwę spatial: CBMSTA / Hao et al. 2024 (NIE POBRANE)

**Hao S. et al. „Cross-species single-cell spatial transcriptomic atlases of the
cerebellar cortex". *Science* 2024; 385: eado3927.**
DOI `10.1126/science.ado3927` (dostęp 2026-09-03; **pełny tekst za paywallem,
`science.org` zwraca HTTP 403** — poniższe pochodzi z abstraktu, strony STOmics
i listingu FTP, nie z metod pracy).

Metoda: **Stereo-seq** (STOmics/BGI), rozdzielczość przestrzenna **~500 nm** —
to jest te „500 nm". Gatunki: **mysz, marmozeta, makak**. Kluczowe twierdzenie
pracy: dwa **naczelno-specyficzne podtypy Purkinje** różniące się ekspresją
**`GRID2`**, rozmieszczone w odrębnych mikrodomenach o różnej liczności
w płacikach.

`Grid2` jest w naszym `raw/cb_adult_mouse_genes.txt` w **linii 7385**.

### Baza: CBMSTA, akcesja `STDS0000244`

Listing: `https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000244/`
(bez logowania). **468 plików** w katalogu głównym:

| Format | Liczba | Co to |
|---|---|---|
| `.gem.gz` | 231 | macierze bin1 (surowe, współrzędne × gen) |
| `.h5ad` (sekcje) | 222 | pojedyncze sekcje Stereo-seq |
| `.h5ad` (snRNAseq) | 6 | `{Mouse,Macaque,Marmoset}.sn[.norm].h5ad` |
| `.rds` | 9 | obiekty Seurat per typ komórki |

Sekcje Stereo-seq wg gatunku: **Marmoset 111, Mouse 61, Macaque 50.**

### Rozmiary — zmierzone `curl -I` (Content-Length), 2026-09-03

| Plik | Bajty | GiB |
|---|---:|---:|
| `Mouse.sn.norm.h5ad` | 9 502 522 200 | **8.85** |
| `Mouse.sn.h5ad` | 8 275 712 358 | **7.71** |
| `mouse.Purkinje.rds` | 6 058 589 878 | **5.64** |
| `Marmoset.sn.h5ad` | 4 250 115 038 | 3.96 |
| `Macaque.sn.h5ad` | 3 313 709 956 | 3.09 |
| `macaque.Purkinje.rds` | 2 913 371 046 | 2.71 |
| `marmoset.Purkinje.rds` | 1 771 161 495 | 1.65 |
| `AST_subcluster.rds` | 1 546 198 922 | 1.44 |
| `macaque_granule.rds` | 1 153 162 938 | 1.07 |
| `Mouse1_T175.gem.gz` | 1 127 467 628 | 1.05 |
| `Macaque1_T80.h5ad` | 807 059 307 | 0.75 |
| `Mouse1_T175.h5ad` | 330 089 110 | 0.31 |

### Rozmiar całego zbioru — ZMIERZONY 2026-09-03

**Wszystkie 468 plików, `curl -I` na każdym, zero braków:**

## **289 717 330 095 B = 269.82 GiB = 0.263 TiB**

| kategoria | n | suma GiB |
|---|---:|---:|
| `.gem.gz` (bin1, surowe) | 231 | **163.85** |
| `.h5ad` sekcje Stereo-seq | 222 | 58.09 |
| `.h5ad` snRNAseq | 6 | 28.61 |
| `.rds` (Seurat) | 9 | 19.27 |

| gatunek | n | suma GiB |
|---|---:|---:|
| Macaque | 104 | 114.52 |
| Marmoset | 235 | 99.47 |
| **Mouse** | **125** | **50.75** |
| wspólne | 4 | 5.07 |

**Rozbicie dla myszy — to jest realna decyzja o pobieraniu:**

| | n | GiB |
|---|---:|---:|
| `.gem.gz` bin1 (surowe, raczej niepotrzebne) | 61 | 20.63 |
| `.h5ad` snRNAseq (`Mouse.sn` + `.norm`) | 2 | 16.56 |
| **`.h5ad` sekcje Stereo-seq** | **61** | **7.92** |
| `.rds` (Seurat, nieczytelny bez R) | 1 | 5.64 |

**Punkt wejścia: 61 sekcji Stereo-seq myszy = 7.92 GiB, średnio 133 MiB/sekcję.**
To mieści się wszędzie i każda sekcja z osobna wchodzi do RAM bez problemu.

Cały zbiór (269.82 GiB) zająłby **55 % wolnego miejsca na `/mnt/data1t`**
(488 G wolne) — technicznie się mieści, ale nie ma powodu go ciągnąć w całości.
Pobierać selektywnie.

**Pułapka:** `mouse_granule.rds` **nie istnieje** — URL zwraca stronę 404
nginx o długości 146 B, którą `curl -I` raportuje jako `Content-Length: 146`.
W listingu są tylko `macaque_granule.rds` i `marmoset_granule.rds`.
Nie ufaj `Content-Length` bez sprawdzenia, czy plik jest na liście.

### Konsekwencje dla tej maszyny (15.54 GiB RAM, 9.3 GiB wolne)

1. **`Mouse.sn.norm.h5ad` (8.85 GiB skompresowany) i `Mouse.sn.h5ad`
   (7.71 GiB) są nie do wczytania lokalnie w całości.** To rozmiar **na dysku**;
   po dekompresji do RAM będzie **więcej**, o ile — **NIEPOLICZONE**.
   Te pliki idą na eagle (węzeł 384 GB) albo nigdzie.
2. **`.rds` to obiekty Seurat/R — a na tej maszynie R ma 29 pakietów bazowych,
   bez Seurata.** `mouse.Purkinje.rds` (5.64 GiB) wymagałby zbudowania stacku R.
   Współczynnik rozprężenia RDS (gzip) **NIEPOLICZONY**.
3. **Do pracy lokalnej nadają się pojedyncze sekcje `.h5ad`** —
   `Mouse1_T175.h5ad` ma 0.31 GiB. To jest właściwa jednostka pracy na tym
   sprzęcie: jedna–kilka sekcji naraz, nie cały atlas.
4. **Miejsce na dysku:** `/home` ma 74 G wolnego. Cały zbiór (468 plików,
   rozmiar nieznany) **na pewno tam nie wejdzie**. Pobieranie wyłącznie na
   `/mnt/data1t` (488 G wolne), selektywnie.

### Czego NIE wiem o tej pracy

Z pełnego tekstu (paywall) **nie mam**: liczby komórek/spotów, liczby jąder
w snRNA-seq, specyfikacji chipów, dokładnej geometrii DNB, ani nazw podtypów
Purkinje. **Nie wolno tych liczb podawać, dopóki nie zostanie zdobyty PDF.**

---

## 11. Referencja 3D móżdżku: CCFv3a (BBP) — do wizualizacji fenotypów

**Piluso S., Verasztó C., Carey H., Delattre É., L'Yvonnet T., Colnot É.,
Romani A., Bjaalie J.G., Markram H., Keller D. „An extended and improved CCFv3
annotation and Nissl atlas of the entire mouse brain". *Imaging Neuroscience*
2025; 3: imag_a_00565.** DOI `10.1162/imag_a_00565`
Dane: **Zenodo `15176439`**, licencja CC-BY-4.0 (dostęp 2026-09-03).

### Dlaczego ten, a nie zwykły CCFv3

Oryginalny CCFv3 Allena ma w móżdżku **etykiety tylko na poziomie płacików,
bez podziału na warstwy**. CCFv3a dokłada dla **wszystkich 16 płacików**
warstwę ziarnistą, drobinową i — kluczowe — **warstwę Purkinjego**
(dwie warstwy wokseli przy 10 µm).

W BrainGlobe dostępny jako **`ccfv3augmented_mouse_10um`** /
**`ccfv3augmented_mouse_25um`** (z siatkami 3D).

### ✅ Płaciki Kozarevy mapują się na CCFv3a w stosunku 1:1 — 16/16

Zweryfikowane programowo na `hierarchy_bbp_atlas_pipeline.json`
(1798 struktur, 19 struktur „Purkinje layer"):

| Kozareva `region` | CCFv3a akronim | id | nazwa CCFv3a |
|---|---|---:|---|
| `I` | `LINGpu` | 10706 | Lingula (I), Purkinje layer |
| `II` | `CENT2pu` | 10709 | Lobule II, Purkinje layer |
| `III` | `CENT3pu` | 10712 | Lobule III, Purkinje layer |
| `CUL` | `CUL4, 5pu` | 10721 | Lobules IV-V, Purkinje layer |
| `VI` | `DECpu` | 10724 | Declive (VI), Purkinje layer |
| `VII` | `FOTUpu` | 10727 | Folium-tuber vermis (VII), Purkinje layer |
| `VIII` | `PYRpu` | 10730 | Pyramus (VIII), Purkinje layer |
| `IX` | `UVUpu` | 10733 | Uvula (IX), Purkinje layer |
| `X` | `NODpu` | 10736 | Nodulus (X), Purkinje layer |
| `SIM` | `SIMpu` | 10673 | Simple lobule, Purkinje layer |
| `AN1` | `ANcr1pu` | 10676 | Crus 1, Purkinje layer |
| `AN2` | `ANcr2pu` | 10679 | Crus 2, Purkinje layer |
| `PRM` | `PRMpu` | 10682 | Paramedian lobule, Purkinje layer |
| `COP` | `COPYpu` | 10685 | Copula pyramidis, Purkinje layer |
| `PF` | `PFLpu` | 10688 | Paraflocculus, Purkinje layer |
| `F` | `FLpu` | 10691 | Flocculus, Purkinje layer |

Rodzic wszystkich: `CBXpu` (id 1145), „Cerebellar cortex, Purkinje layer".

### ⚠ KOREKTA 2026-09-04 — te struktury mają ZERO wokseli

Powyższa tabela pochodzi z **hierarchii** (`hierarchy_bbp_atlas_pipeline.json`).
**Wydany wolumen adnotacji ich nie zawiera.** Sprawdzone pełnym przelotem przez
wszystkie **1 290 480 000 wokseli** `annotv3a_bbp_10.nrrd`: każda z 19 struktur
`*pu` ma **0 wokseli**, podczas gdy każda `*gr` i `*mo` ma miliony
(np. `DECgr` 1 918 320, `DECmo` 1 437 124).

Kontrola poprawności odczytu: `DECgr` przy 10 µm / przy 25 µm = **15.70**
wobec teoretycznego 2.5³ = 15.625. Odczyt jest dobry — struktur po prostu nie ma.

W atlasie brainglobe 25 µm widać to po dziurach w numeracji: 10672 `SIMgr`,
**brak 10673**, 10674 `SIMmo`. Przyczyna jest fizyczna: warstwa Purkinjego
ma w CCFv3a 2 woksele przy 10 µm (20 µm), więc przy 25 µm jest podwokselowa.

**Rozwiązanie zastosowane:** warstwa liczona geometrycznie jako **granica
`*gr` | `*mo`** (6-sąsiedztwo) — tam z definicji leżą ciała komórek Purkinjego.
Wynik: **2 279 886 wokseli = 2.28 mm³**. Kontrola spójności: granica stanowi
**7.0–11.8 %** objętości warstwy ziarnistej w każdym z 16 płacików osobno,
mimo czterdziestokrotnych różnic objętości między nimi.
Skrypt: `scripts/09_pc_layer_geom.py`, wynik: `/mnt/data1t/pc_rebuild/pc_layer_10um.npz`.

**Uwaga:** to mapowanie potwierdza niezależnie ustalenie z sekcji 6 —
prefiksy `DEC*` i `CUL*` w barcodach to **nazwy z atlasu Allena**
(Declive = VI, Culmen = IV–V), a nie inna tkanka.

### Rozmiary i RAM — zmierzone

Nagłówki NRRD odczytane HTTP range request z Zenodo:

| Plik | Na dysku | Wymiary | Typ | RAM (uint32) | RAM (uint16) |
|---|---:|---|---|---:|---:|
| `annotv3a_bbp_10.nrrd` | 28.89 MiB | 1415×800×1140 | uint32 | **4.81 GiB** | 2.40 GiB |
| `annotv3a_bbp_25.nrrd` | 3.90 MiB | 566×320×456 | uint32 | **0.31 GiB** | 0.15 GiB |
| `hierarchy_bbp_atlas_pipeline.json` | 0.82 MiB | — | JSON | — | — |
| `arav3a_bbp_nisslCOR_10.nrrd` | 2091 MiB | — | — | NIEPOLICZONE | — |

Pliki są małe **tylko na dysku** — to gzip po danych etykietowych, które
kompresują się ekstremalnie. **28.89 MiB rozpakowuje się do 4.81 GiB.**

**Reguła dla tej maszyny:** używaj **25 µm** (0.31 GiB, mieści się w limicie
2.33 GiB z sekcji 5). Wersja **10 µm to 4.81 GiB — przekracza limit ponad
dwukrotnie**; wczytywać tylko jako uint16 (2.40 GiB, nadal na granicy) albo
po wycięciu samego móżdżku. Udział objętościowy móżdżku: **NIEPOLICZONE**.

Cały rekord Zenodo waży **21.9 GB**, ale do tego zadania potrzeba **dwóch
plików o łącznej wadze poniżej 5 MiB** (`annotv3a_bbp_25.nrrd` + hierarchia).

### ⚠ Ograniczenie merytoryczne, którego nie da się obejść danymi

Kozareva daje **etykietę płacika, nie współrzędną**. Zmapowanie fenotypów
Aldoc na CCFv3a da **16 pomalowanych płacików**, a nie pasy zebriny.
Pasy zebriny to wąskie pasma parasagittalne **wewnątrz** płacików —
są **drobniejsze niż rozdzielczość anotacji Kozarevy**.

Żeby dostać prawdziwe pasy w 3D, potrzeba danych ze współrzędnymi:
**Stereo-seq z sekcji 10 (Hao et al.)** ma `x`,`y` per sekcja i rekonstrukcję
3D — to jest właściwe źródło pasów. CCFv3a jest wtedy układem odniesienia,
do którego się rejestruje, a nie źródłem sygnału.

Nie prezentuj mapy „per płacik" jako mapy pasów zebriny. To dwa różne
poziomy rozdzielczości.


---

## 12. Stan po sesji 2026-09-04

### `processed/purkinje_cells_v2.h5ad` — PLIK ROBOCZY

**16 634 × 24 409**, CSR `int32`, nnz **78 152 875**, 0.58 GiB w RAM.
Zbudowany `scripts/00_prep_map.py` → `01_extract.sh` → `02_build_h5ad.py`.

Ekstrakcja: blok kolumn **43856–60489** (ciągły, 0 przerw), wyjście wczesne
po przeskanowaniu 182 392 892 wpisów, **0 spadków numeru kolumny** (mtx jest
posortowany). Czas: 94 s.

**Walidacja per komórka:** 16 538 / 16 634 ma `nnz == nGene` dokładnie.
Pozostałe **96 to wszystkie dokładnie −1** (pełna lista:
`/mnt/data1t/pc_rebuild/walidacja_niezgodne.csv`), rozsiane po wielu próbkach.
To artefakt w tym, jak policzono `nGene` u źródła — ta sama natura co globalna
rozbieżność 504 z sekcji 5. Dotyczy 0.577 % komórek.

**Naprawa działa** (porównanie z wadliwym oryginałem):

| region | stary | v2 | metadane |
|---|---:|---:|---:|
| CUL | 136 | **592** | 592 ✓ |
| VI | 1 574 | **2 285** | 2 285 ✓ |

Wszystkie 16 regionów i 9 subklastrów zgadza się z metadanymi co do komórki.

### `processed/purkinje_cells_v2_processed.h5ad`

QC (`nGene>=500`, `pct_mt<5`) odrzuciło **15 z 16 634** — komórki Purkinjego są
bardzo czyste (mediana `nGene` 4846, `nUMI` 25 780, `pct_counts_mt` 0.021 %).
HVG 2000, PCA 50, UMAP, Leiden przy 7 rozdzielczościach.

**Walidacja klastrowania:** przy rozdzielczości 0.2 wychodzi dokładnie
**9 klastrów, ARI = 0.499, NMI = 0.577** względem `subcluster` Kozarevy.

### Klasyfikacja podtypów — UWAGA METODOLOGICZNA

| metoda | dokładność | zbalansowana |
|---|---:|---:|
| `sc.tl.score_genes` + argmax | 0.334 | 0.520 |
| regresja logistyczna, 225 genów sygnaturowych | 0.882 | 0.893 |
| **regresja logistyczna, 30 PC** | **0.917** | **0.929** |
| to samo, zredukowane do osi Aldoc+/− | **0.974** | 0.972 |

**`score_genes` + argmax jest bezużyteczny do klasyfikacji wieloklasowej** —
wyniki nie są porównywalne między zestawami genów, klasa `Aldoc_1` zagarnia
wszystko (recall 1.00 przy precyzji 0.09). Nie używać tej metody do transferu
etykiet. Każdy z 9 podtypów osiąga F1 ≥ 0.86 przy regresji logistycznej.

### ⚠ Sekcje Hao NIE są ko-rejestrowane

Zmierzone na 4 sekcjach (`scripts/03_gate_registration.py`):

| sekcja | komórek | PC-layer | genów | środek X | środek Y |
|---|---:|---:|---:|---:|---:|
| Mouse1_T167 | 18 548 | 2 933 | 21 481 | 267.9 | 76.6 |
| Mouse1_T175 | 62 227 | 5 697 | 23 192 | 200.8 | 140.6 |
| Mouse1_T190 | 78 108 | 15 659 | 23 177 | 241.7 | 181.7 |
| Mouse2_T360 | 84 660 | 12 475 | 22 517 | 303.2 | 161.0 |

Środki rozjeżdżają się o 102 (X) i 105 (Y); **listy genów też są różne**.
**Nie stackować ich bez rejestracji.** Dlatego 3D bierze się z atlasu,
a sekcje Hao zostają jako 2D w pełnej rozdzielczości.

### Wyniki wizualne

- `figures/pc_layer_3d.png` — cztery rzuty 3D warstwy Purkinjego
- `figures/purkinje3d.html` — interaktywny model 3D (Three.js), opublikowany jako artifact
- `processed/lobule_composition.csv` — skład 9 podtypów per płacik
- `processed/subtype_signatures.json` — sygnatury genowe

### Skrypty (`scripts/`)

`00_prep_map.py` mapowanie barcodów · `01_extract.sh` ekstrakcja z mtx ·
`02_build_h5ad.py` złożenie + walidacja · `03_gate_registration.py` bramka Hao ·
`04_process.py` QC/PCA/UMAP/Leiden · `08_enumerate_labels.py` etykiety atlasu ·
`09_pc_layer_geom.py` warstwa PU jako granica gr|mo · `10_signatures.py` sygnatury ·
`12_classifier.py` klasyfikatory · `11_build_scene.py` scena · `13b_render_mpl.py` render ·
`14_export_web.py` eksport do strony

---

## 13. Domeny Aldoc — wskaźnik dwugenowy + łączenie po sąsiadach (2026-09-04)

### Dlaczego transfer podtypów Kozareva→Hao ZAWIÓDŁ

**Nie używać `sc.tl.score_genes` ani klasyfikatora trenowanego na snRNA-seq
bezpośrednio na Stereo-seq.** Zmierzone:

| frakcja zliczeń | Kozareva (jądra) | Hao (cała tkanka) |
|---|---:|---:|
| **`Malat1`** | **22.67 %** | 0.04 % |
| rybosomalne `Rps*`/`Rpl*` | 0.18 % | **4.83 %** |
| mitochondrialne `mt-*` | 0.02 % | **0.64 %** |
| `Meg3` | 0.53 % | 0.02 % |

To nie jest różnica biologiczna, tylko **przedziału komórkowego**: snRNA-seq
mierzy jądro, Stereo-seq całą tkankę. Klasyfikator (0.794 zbalansowanej na
spłyconych danych Kozarevy) na Hao dał **100 % Aldoc+ i 81–88 % w klasie
`Aldoc_2`, która u Kozarevy ma n=75 (0.45 %)**. Dowód błędu: komórki zmierzone
jako Aldoc-niskie (39.3 CP10K) i Aldoc-wysokie (102.9 CP10K) dostały **tę samą
etykietę** (`Aldoc_2` 90 % vs 87 %).

Spłycenie kontroluje głębokość, **nie kontroluje platformy**. To był mój błąd
metodologiczny.

### Korekta głębokości

Mediana 1 734 UMI dotyczy **wszystkich** komórek przekroju. Same komórki
warstwy Purkinjego: **5 324 UMI, 2 432 geny** (21 % głębokości Kozarevy).
Łączenie sąsiadów: **5 sąsiadów = 104 %**, 20 = 420 % głębokości Kozarevy.

### Co DZIAŁA — procedura

`scripts/20_axis_and_ribbon.py`, `scripts/21_pooled.py`

1. **Wskaźnik dwugenowy** `log2((Aldoc_cp10k+1)/(Plcb4_cp10k+1))` zamiast progu
   na samym Aldoc. Obniża zależność od głębokości o 23–70 %.
2. **Rozcięcie wstęgi na fałdy** po orientacji: normalna warstwy Purkinjego
   (od ziarnistej ku drobinowej) jest **przeciwna po dwóch stronach szczeliny**,
   więc łączymy tylko komórki o zgodnej normalnej (iloczyn skalarny > 0.5).
   Daje 9–45 spójnych fałdów, porządkuje 72–89 % komórek.
3. **Łączenie ~8 sąsiadów w obrębie fałdu** — bez tego `Plcb4` jest za rzadki.

### Wyniki (4 przekroje)

| przekrój | n | `Plcb4` wykryty: pojed. → łączone | zależność od głębokości \|r\| | zgodność sąsiadów | z | Aldoc-wysokie |
|---|---:|---|---:|---|---:|---:|
| Mouse1_T167 | 733 | 22.2 % → **84.7 %** | 0.051 | 84.7 % vs 50.0 % | 37.6 | 48.2 % |
| Mouse1_T175 | 1 424 | 59.7 % → **99.1 %** | 0.019 | 91.3 % vs 49.9 % | 58.4 | 50.9 % |
| Mouse1_T190 | 3 915 | 19.5 % → **81.4 %** | 0.201 | 87.2 % vs 51.3 % | 90.5 | 42.0 % |
| Mouse2_T360 | 3 119 | 13.0 % → **61.2 %** | 0.156 | 86.9 % vs 50.7 % | 79.9 | 55.9 % |

**Udział Aldoc-wysokich 42–56 % w czterech niezależnych przekrojach** — przed
poprawką było 49–84 %. Ta zbieżność to najmocniejsza kontrola wewnętrzna.

### Rozmiary domen

**525 domen ≥3 komórek w 4 przekrojach: mediana 8 komórek Purkinjego,
p75 = 17, maksimum 121.**

Domeny 1–2 komórkowe (reszta rozkładu) traktować jako szum resztkowy,
nie jako pasy.

**Skali fizycznej NIE ZNAM** — dwie metody dały niespójne wyniki
(3.98–11.04 i 6.28–13.48 µm/jednostkę, `scripts/15_*`, `scripts/16_*`).
Rozmiary domen podawać **w liczbie komórek**, nigdy w mikrometrach.

### Wynik wizualny

`figures/domains_pooled.png` — 4 przekroje + rozkład rozmiarów domen.

---

## 14. Rejestracja serii Hao do atlasu (2026-09-05)

Paczka: `/mnt/data1t/hao_pack/` (~60 MB). Dane: 60 z 61 sekcji myszy;
**`Mouse2_T359` nie ma pola `annotation`** — błąd w danych źródłowych, nie u nas.
`Mouse1` T167–T198 (32, ciągłe), `Mouse2` T345–T373 (28, luka na T359).

### Skala — pierwszy powtarzalny pomiar

| metoda | skala xy | odstęp z |
|---|---|---|
| grubość warstwy drobinowej (`15_*`) | 3.98–11.04 µm/j | — |
| kształt pojedynczego skrawka (`16_*`) | 6.28–13.48 µm/j | — |
| **cała seria dopasowana do atlasu** | **19.2–20.0 µm/j** | **~100 µm** |

Pierwsze dwie metody były bezwartościowe — opierały się na pojedynczych
skrawkach. Dopiero seria daje dość więzów. **Używać 19–20 µm/jednostkę.**

### Pułapki, w które wpadłem (wszystkie: brak kary za zwyrodnienie)

1. **ICP ze swobodną skalą** zjeżdża do progu (0.500), bo kurczenie źródła
   zmniejsza każdą odległość. Skrawki tej samej serii mają skalę **1, na sztywno**.
2. **Procrustes na nieprzeskalowanych punktach** przy odpowiednikach ze
   skalowanych → niespójność między iteracjami.
3. **Jednokierunkowy Chamfer w 3D** ścisnął chmurę do punktu: skala `1.8e-17`,
   koszt „0.00 µm" wyglądał jak sukces. Konieczny **koszt symetryczny**.
4. **Wspólna skala dla `x,y` i `z`** — `z` to numer skrawka, inna jednostka.
5. **Mieszanie radianów z mikrometrami** w jednym wektorze Nelder-Mead →
   optymalizator rusza tylko obrotem, przesunięcia wychodzą dokładnie 0.
   Parametry muszą być bezwymiarowe.
6. **ICP na kształcie ~symetrycznym** znajduje odbicie ~180°. Ograniczenie
   `|obrót| < 40°` między sąsiadami.

### Dryf akumulacyjny — zmierzony, usunięty, i NIE był przyczyną

Łańcuch 31 par: suma obrotów **+126.2°**; bezpośrednia rejestracja
`T167→T198`: **+38.8°** (błąd 4.00). **Dryf 87.4°.**
Łańcuch jest przy tym **niestabilny między uruchomieniami** (+126.2°, +23.7°,
−298.6° w trzech przebiegach) — to samo w sobie uzasadnia bundle adjustment.

Odporne BA (141 więzów o odstępach 1,2,3,5,8; IRLS Tukeya, 25 odrzuconych)
sprowadziło obrót ostatniego skrawka do −17.4°, **ale dopasowanie do atlasu
się NIE poprawiło** (101 → 108 µm). Dryf nie był przyczyną tych 100 µm.

⚠ Wariant odpornego BA dał medianę 97 µm, ale **rozciągnął bryłę do 14.33 mm
przy móżdżku 9.05 mm** — odrzucony przez kontrolę gabarytów. Zawsze sprawdzać
gabaryty, nie tylko medianę.

### Co zadziałało: poprawka per skrawek

| | mediana | w 20 µm | w 50 µm |
|---|---:|---:|---:|
| po globalnym dopasowaniu | 103 µm | ~11 % | ~26 % |
| **+ poprawka per skrawek** | **76 µm** | 15.7 % | 35.7 % |
| + wygładzenie wzdłuż serii | 85 µm | 13.2 % | 31.2 % |

Wygładzenie **pogarsza** — skrawki różnią się indywidualnie.

**Walidacja na wstrzymanych komórkach** (poprawka uczona na połowie, mierzona
na drugiej): uczące 107→81 µm, **testowe 106→81 µm, różnica 0.0 p.p.**
Zero przeuczenia, poprawa w całości prawdziwa.

### ✅ WZAJEMNA WALIDACJA DWÓCH ZBIORÓW

Po rejestracji 84.6 % komórek trafia w promień 200 µm od warstwy Purkinjego
atlasu → przypisanie do płacika. **206 807 komórek, 15 płacików.**

Porównanie udziału Aldoc+ zmierzonego **przestrzennie** (Hao) z mierzonym
**z dysekcji** (Kozareva):

**Pearson r = +0.775 · Spearman r = +0.732, p = 0.0019**

Dwa niezależne zbiory — inne zwierzęta, inna platforma, inna metoda — zgadzają
się co do **kolejności** płacików. To waliduje cały pipeline.

**ALE wartości bezwzględne są ściśnięte:** Hao 29–81 %, Kozareva 7–100 %,
średni błąd **19.7 punktu procentowego**, systematycznie: płaciki przednie
zawyżone (+20 do +35), tylne i kłaczki zaniżone (−22 do −40).
**Cytować kolejność, nie wartości bezwzględne.**

### Granica, której tymi danymi nie przekroczę

81 µm to cztery grubości warstwy Purkinjego — wystarczy na przypisanie do
płacika (rozmiary mm), za mało na posadzenie komórki na konkretnym fałdzie.
Do tego **CCFv3a to średnia z 1 675 myszy**; pojedynczy móżdżek różni się od
średniej i tej różnicy nie zmierzę. Nie wiem, ile z tych 81 µm to mój błąd,
a ile osobnicza zmienność.

### Pliki

`/mnt/data1t/hao_pack/`: `registration.py`, `derift.py`, `ba_robust.py`,
`refine2.py`, `holdout.py`, `lobules.py`, `eagle_job.sh`, `README.md`,
`sections/*.parquet`, `wyniki/`.
Figury: `figures/atlas_fit.png`, `figures/cross_validation.png`.
