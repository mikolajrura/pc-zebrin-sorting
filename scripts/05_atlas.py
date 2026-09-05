"""Faza 5: pobranie atlasu CCFv3a 25um i wyciecie warstwy Purkinjego per placik."""
import time, os
t0=time.time()
from brainglobe_atlasapi import BrainGlobeAtlas
import numpy as np
print("pobieram ccfv3augmented_mouse_25um ...", flush=True)
atlas=BrainGlobeAtlas("ccfv3augmented_mouse_25um", check_latest=False)
print(f"pobrane w {time.time()-t0:.0f} s")
ann=atlas.annotation
print("annotation:", ann.shape, ann.dtype, f"{ann.nbytes/1024**3:.2f} GiB w RAM")
print("rozdzielczosc (um):", atlas.resolution)
print("orientacja:", atlas.orientation)
ref=atlas.reference
print("reference:", ref.shape, ref.dtype)

# struktury warstwy Purkinjego
pu=[s for s in atlas.structures.values() if 'Purkinje layer' in s['name']]
print(f"\nstruktur 'Purkinje layer': {len(pu)}")
present=[]
uniq=set(np.unique(ann).tolist())
for s in sorted(pu,key=lambda s:s['id']):
    inv = s['id'] in uniq
    n = int((ann==s['id']).sum()) if inv else 0
    present.append((s['id'],s['acronym'],s['name'],n))
    print(f"  id={s['id']:>6} {s['acronym']:<12} wokseli={n:>8}  {s['name']}")
tot=sum(p[3] for p in present)
print(f"\nlacznie wokseli warstwy Purkinjego: {tot}  ({tot*25**3/1e9:.3f} mm3)")
np.save("/mnt/data1t/pc_rebuild/atlas_pu_ids.npy", np.array([p[0] for p in present if p[3]>0]))
print("zapisano atlas_pu_ids.npy")
