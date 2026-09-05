# Dostępne środowiska / narzędzia

Conda envs pod `~/miniforge_envs/`: **cb1-gromacs**, **cytoscape**, **genometools**, **kaggle-cli**, **pharma**, **tools**.

Dla tego projektu (scRNA-seq Purkinje, `purkinje_cells.h5ad`) właściwe jest **`pharma`** — nie `genometools`.

---

## `pharma` — env do analizy (PRIMARY)

Źródło: `~/.config/bash/modules/pharma.sh` → symlink `./pharma.sh`
Ścieżka env: `/home/mikolajrurad/miniforge_envs/pharma` (Python 3.10.19)

### Aktywacja / helpery

```bash
pharma                           # source env + prompt "(pharma)"
condout                          # conda deactivate (alias)
drugview [...]                   # GUI app
pdbview <pdbid>                  # PyMOL fetch + cartoon
pdbview -compare <id1> <id2> ... # PyMOL grid, wiele struktur
pdbview -compare -align <ref> <id1> ...  # z wyrównaniem
```

### Stack scverse / single-cell

- **anndata 0.11.4**, **mudata 0.3.2**, **h5py 3.15.1**
- **scanpy 1.11.5**, **scvi-tools 1.3.3**
- **scrublet 0.2.3** (doublety), **celltypist 1.7.1**, **decoupler 2.1.2**
- **leidenalg 0.11.0**, **umap-learn 0.5.9**
- scikit-learn 1.7.2, scikit-image 0.25.2, scipy 1.15.2, pandas 2.3.3, numpy, statsmodels 0.14.6
- matplotlib 3.10.8, seaborn 0.13.2, plotly 6.5, napari 0.5.5

### Deep learning (CUDA 12.9)

- pytorch 2.7.1 (CUDA), pytorch-lightning 2.6.0
- tensorflow 2.19.1 (CUDA)
- pyro-ppl 1.9.1, numpyro 0.19.0

### Chem / strukturalne

- **rdkit 2025.09.3**, pymol-open-source 3.1 (+ pymol 3.2 pip), biopython 1.86
- geopandas-base, dask 2024.11.2, jupyterlab 4.5

### Szybki start dla projektu

```bash
pharma
python -c "import anndata; a = anndata.read_h5ad('processed/purkinje_cells.h5ad'); print(a)"
jupyter lab --no-browser
```

---

## `genometools` — env do NGS CLI (SECONDARY)

Źródło: `~/.config/bash/modules/genometools.sh` → symlink `./genometools.sh`
Ścieżka: `/home/mikolajrurad/miniforge_envs/genometools` (Python 3.11.15)

### Aktywacja

```bash
genometools              # source env
gtrun <cmd> [args]       # jednorazowe uruchomienie
revcomp ATCG             # reverse complement
```

### Zawartość

| Narzędzie     | Wersja | Zastosowanie                               |
| ------------- | ------ | ------------------------------------------ |
| samtools      | 1.23   | BAM/SAM/CRAM                               |
| bcftools      | 1.23.1 | VCF/BCF                                    |
| htslib        | 1.23   | backend                                    |
| bedtools      | 2.31.1 | genome arithmetic                          |
| blast         | 2.17.0 | blastn/p/x                                 |
| minimap2      | 2.30   | aligner long-read / spliced                |
| mafft         | 7.526  | MSA                                        |
| seqkit        | 2.13.0 | FASTA/FASTQ toolkit                        |
| entrez-direct | 24.0   | NCBI E-utilities                           |
| ncbi-vdb      | 3.3.0  | SRA backend                                |
| primer3-py    | 2.3.0  | projektowanie primerów                     |
| biopython     | 1.86   | `Bio.*`                                    |

Stack Pythona: numpy 2.4.2, jupyterlab 4.5.5, kaggle 1.8.3. **Brak** scanpy/anndata/pandas/scipy/R.

---

## Inne envy (do wglądu)

- `cb1-gromacs` — GROMACS (dynamika molekularna)
- `cytoscape` — Cytoscape (sieci biologiczne)
- `kaggle-cli` — CLI Kaggle
- `tools` — ogólny

## Pełna lista paczek w dowolnym env

```bash
mamba list -p ~/miniforge_envs/<env>
```
