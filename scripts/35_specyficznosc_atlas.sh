#!/bin/bash
# Wyciaga z pelnej macierzy 611 034 komorek wiersze 7 genow markerowych.
# Wiersze mtx: Plcb4 2302, Calb1 4548, Syne1 9615, Nav3 10307, Aldoc 13976, Syne2 15708, Gnaq 18983
set -euo pipefail
P=/home/mikolajrurad/omics-data-highway/ichb-purkinje-analysis/pc-zebrin-sorting
O=/mnt/data1t/pc_rebuild/geny_atlas.tsv
zcat "$P/raw/cb_adult_mouse.mtx.gz" | awk 'NR>2 && ($1==2302||$1==4548||$1==9615||$1==10307||$1==13976||$1==15708||$1==18983){print $1"\t"$2"\t"$3}' > "$O"
wc -l "$O"
