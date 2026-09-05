"""Faza 6: render 3D warstwy Purkinjego pokolorowanej zmierzona frakcja Aldoc+."""
import numpy as np, napari, warnings, os
warnings.filterwarnings("ignore")
from matplotlib.colors import LinearSegmentedColormap, Normalize
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
OUT="/mnt/data1t/pc_rebuild"; P="/home/mikolajrurad/omics-data-highway/ichb-purkinje-analysis/pc-zebrin-sorting"
BLUE=["#cde2fb","#9ec5f4","#6da7ec","#3987e5","#2a78d6","#1c5cab","#104281","#0d366b"]
seq=LinearSegmentedColormap.from_list("b",BLUE)

d=np.load(f"{OUT}/scene_pc_layer.npz",allow_pickle=True)
i,j,k=d["i"],d["j"],d["k"]; ald=d["aldoc_frac"]; lobn=d["lobule_name"].astype(str)
print(f"wokseli: {len(i):,}")
rng=np.random.default_rng(0)
N=300000
sel=rng.choice(len(i),size=min(N,len(i)),replace=False)
print(f"do renderu: {len(sel):,} punktow (podprobka)")
# atlas 10um -> mikrometry, os (i,j,k) wg NRRD sizes 1415/800/1140
pts=np.c_[k[sel],j[sel],i[sel]].astype(np.float32)   # napari: (z,y,x)
a=ald[sel]
cols=seq(Normalize(np.nanmin(ald),np.nanmax(ald))(a))

v=napari.Viewer(show=False,ndisplay=3)
v.add_points(pts,size=2.2,face_color=cols,border_width=0,name="warstwa Purkinjego / Aldoc+",
             shading='none',blending='translucent_no_depth')
v.dims.ndisplay=3
v.camera.zoom=1.0
shots={}
for nm,(ang,zoom) in {"boczny":((0,0,90),1.0),"tylny":((0,0,0),1.0),
                      "gorny":((90,0,0),1.0),"skos":((30,40,20),1.0)}.items():
    v.camera.angles=ang; v.camera.zoom=zoom
    v.reset_view()
    v.camera.angles=ang
    shots[nm]=v.screenshot(canvas_only=True,size=(1100,1100))
    print("  render:",nm,shots[nm].shape)
v.close()

fig,axes=plt.subplots(1,4,figsize=(20,5.6),facecolor="#fcfcfb")
for ax,(nm,img) in zip(axes,shots.items()):
    ax.imshow(img); ax.set_title(nm,fontsize=12,color="#0b0b0b"); ax.axis("off")
sm=plt.cm.ScalarMappable(cmap=seq,norm=Normalize(np.nanmin(ald),np.nanmax(ald)))
cb=fig.colorbar(sm,ax=axes,fraction=0.012,pad=0.01)
cb.set_label("frakcja Aldoc+ w płaciku (zmierzona, Kozareva)",color="#52514e",fontsize=10)
cb.outline.set_visible(False); cb.ax.tick_params(colors="#52514e",labelsize=9)
fig.suptitle("Warstwa Purkinjego myszy, CCFv3a 10 µm — 2 279 886 wokseli, "
             "kolor = zmierzona frakcja Aldoc+ w danym płaciku",
             fontsize=12,color="#52514e",x=0.01,ha='left',y=0.97)
fig.savefig(f"{P}/figures/pc_layer_3d.png",dpi=140,facecolor="#fcfcfb",bbox_inches='tight')
print("zapisano figures/pc_layer_3d.png")
