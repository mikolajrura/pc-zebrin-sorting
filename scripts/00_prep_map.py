"""Faza 0, krok 1: mapa barcodow Purkinje -> kolumny mtx + tabela obs."""
import csv, sys, os
P="/home/mikolajrurad/omics-data-highway/ichb-purkinje-analysis/pc-zebrin-sorting"
OUT="/mnt/data1t/pc_rebuild"
REPO2GEO={'DECa_F002':'VI2a_F002','DECb_F002':'VI2b_F002','DECc_F002':'VI2c_F002',
          'DECa_M006':'VI2a_M006','DECb_M006':'VI2b_M006','DECc_M006':'VI2c_M006',
          'Va_M002':'CULa_M002','Vb_M002':'CULb_M002','Vc_M002':'CULc_M002',
          'Vd_M002':'CULd_M002','IV_M006':'CUL_M006','Vl_M006':'VI_M006'}
def to_geo(bc):
    a=bc.split('_'); pre='_'.join(a[:2]); suf='_'.join(a[2:])
    return f"{REPO2GEO[pre]}_{suf}" if pre in REPO2GEO else bc

col={}
with open(f"{P}/raw/cb_adult_mouse_barcodes.txt") as fh:
    for i,l in enumerate(fh,1): col[l.strip()]=i
print("barcodow w mtx:", len(col))

rows=[]
with open(f"{P}/cerebellum-atlas-analysis/data/full_cb_metadata.csv", newline='') as fh:
    r=csv.reader(fh); h=next(r); h[0]='barcode'; ix={n:i for i,n in enumerate(h)}
    for row in r:
        if row[ix['final_annotation_cluster']]=='Purkinje': rows.append(row)
print("Purkinje w metadanych:", len(rows))

out=[]; miss=[]
for row in rows:
    bc=row[0]; geo=to_geo(bc); c=col.get(geo)
    (out if c else miss).append((bc,geo,c,row) if c else (bc,geo))
print(f"trafionych: {len(out)}   NIETRAFIONYCH: {len(miss)}")
if miss:
    for m in miss[:10]: print("   ",m)
    sys.exit("PRZERWANE")

out.sort(key=lambda t:t[2])
cols=[t[2] for t in out]
assert cols[-1]-cols[0]+1==len(cols), "blok kolumn NIE jest ciagly"
print(f"blok kolumn: {cols[0]}..{cols[-1]} (ciagly, n={len(cols)})")

with open(f"{OUT}/obs.tsv","w") as fh:
    fh.write("row\tbarcode\tgeo_barcode\tmtx_col\tnGene\tnUMI\torig_ident\tregions\tpercent_mito\tsubcluster\n")
    for i,(bc,geo,c,row) in enumerate(out):
        fh.write("\t".join([str(i),bc,geo,str(c),row[ix['nGene']],row[ix['nUMI']],
                            row[ix['orig.ident']],row[ix['regions']],row[ix['percent_mito']],
                            row[ix['final_annotation_subcluster']]])+"\n")
tot=sum(int(t[3][ix['nGene']]) for t in out)
print(f"suma nGene (oczekiwane nnz): {tot}")
with open(f"{OUT}/block.txt","w") as fh: fh.write(f"{cols[0]} {cols[-1]} {len(cols)} {tot}\n")
print("zapisano obs.tsv i block.txt")
