#!/bin/bash
# Drugi przelot: kandydaci do bramki 1 (Ranbp2 - kontrola protokolu uzytkownika)
# i bramki 2 (Cux2, Ebf2, Zbtb20, Zfpm2, Nr2f2, Prkca + nukleoporyny referencyjne).
set -euo pipefail
P=/home/mikolajrurad/omics-data-highway/ichb-purkinje-analysis/pc-zebrin-sorting
O=/mnt/data1t/pc_rebuild/geny_atlas2.tsv
zcat "$P/raw/cb_adult_mouse.mtx.gz" | awk 'NR>2 && ($1==1574||$1==6597||$1==8856||$1==9855||$1==11154||$1==11963||$1==14472||$1==16199||$1==17112||$1==18537){print $1"\t"$2"\t"$3}' > "$O"
wc -l "$O"
