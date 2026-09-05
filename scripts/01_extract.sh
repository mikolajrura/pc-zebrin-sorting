#!/bin/bash
# Faza 0, krok 2: ekstrakcja bloku Purkinje z mtx.gz, z wczesnym wyjsciem.
# Poprawnosc gwarantuje walidacja per-komorka w kroku 3, nie zalozenie o sortowaniu.
set -o pipefail
P=/home/mikolajrurad/omics-data-highway/ichb-purkinje-analysis/pc-zebrin-sorting
OUT=/mnt/data1t/pc_rebuild
read LO HI N EXPECTED < $OUT/block.txt
echo "blok kolumn $LO..$HI  komorek $N  oczekiwane nnz $EXPECTED"
echo "start $(date +%T)"
zcat $P/raw/cb_adult_mouse.mtx.gz | tail -n +3 | LC_ALL=C gawk -v LO=$LO -v HI=$HI '
BEGIN{ prev=0; drops=0; seen=0; kept=0 }
{
  seen++
  c=$2+0
  if(c<prev) drops++
  prev=c
  if(c>HI){ print "EXIT_EARLY po wpisie "seen" kolumna "c > "/dev/stderr"; exit }
  if(c>=LO){ print (c-LO)"\t"($1-1)"\t"$3; kept++ }
}
END{
  print "wpisow_przeskanowanych\t"seen  > "/dev/stderr"
  print "wpisow_zachowanych\t"kept      > "/dev/stderr"
  print "spadkow_kolumny_do_tego_miejsca\t"drops > "/dev/stderr"
}' > $OUT/triplets.tsv 2> $OUT/extract_stats.txt
echo "koniec $(date +%T)"
cat $OUT/extract_stats.txt
echo "linii w wyjsciu: $(wc -l < $OUT/triplets.tsv)"
echo "oczekiwano:      $EXPECTED"
ls -la $OUT/triplets.tsv
