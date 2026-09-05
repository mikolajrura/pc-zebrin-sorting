"""Statyczne rzuty 3D warstwy Purkinjego (matplotlib - niezawodny offscreen)."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
OUT="/mnt/data1t/pc_rebuild"; P="/home/mikolajrurad/omics-data-highway/ichb-purkinje-analysis/pc-zebrin-sorting"
SURF="#fcfcfb"; INK="#0b0b0b"; INK2="#52514e"
seq=LinearSegmentedColormap.from_list("b",["#cde2fb","#9ec5f4","#6da7ec","#3987e5","#2a78d6","#1c5cab","#104281","#0d366b"])
d=np.load(f"{OUT}/scene_pc_layer.npz",allow_pickle=True)
i,j,k=d["i"].astype(np.float32),d["j"].astype(np.float32),d["k"].astype(np.float32)
ald=d["aldoc_frac"]; lobn=d["lobule_name"].astype(str)
print(f"wokseli: {len(i):,}  zakres Aldoc+: {np.nanmin(ald):.3f}..{np.nanmax(ald):.3f}")
rng=np.random.default_rng(0); N=90000
s=rng.choice(len(i),N,replace=False)
# um: i=oś 1415 (L-P-S x), j=800, k=1140
X,Y,Z=i[s]*0.01, j[s]*0.01, k[s]*0.01     # mm
A=ald[s]
norm=Normalize(np.nanmin(ald),np.nanmax(ald))
views=[("tył (widok od potylicy)",(12,-90)),("bok",(8,0)),("góra",(80,-90)),("skos",(24,-52))]
fig=plt.figure(figsize=(20,5.8),facecolor=SURF)
for n,(t,(el,az)) in enumerate(views,1):
    ax=fig.add_subplot(1,4,n,projection='3d',facecolor=SURF)
    ax.scatter(X,Y,Z,c=A,cmap=seq,norm=norm,s=0.6,linewidths=0,alpha=0.85,rasterized=True)
    ax.view_init(elev=el,azim=az)
    ax.set_title(t,fontsize=12,color=INK,pad=2)
    ax.set_box_aspect((np.ptp(X),np.ptp(Y),np.ptp(Z)))
    ax.set_axis_off()
sm=plt.cm.ScalarMappable(cmap=seq,norm=norm)
cb=fig.colorbar(sm,ax=fig.axes,fraction=0.012,pad=0.01)
cb.set_label("frakcja Aldoc+ w płaciku (zmierzona: Kozareva, n=16 634)",color=INK2,fontsize=10)
cb.outline.set_visible(False); cb.ax.tick_params(colors=INK2,labelsize=9)
fig.suptitle("Warstwa Purkinjego myszy — CCFv3a 10 µm, granica warstwy ziarnistej i drobinowej\n"
             f"2 279 886 wokseli (2.28 mm³), na rysunku podpróbka {N:,}",
             fontsize=11.5,color=INK2,x=0.008,ha='left',y=1.02)
fig.savefig(f"{P}/figures/pc_layer_3d.png",dpi=145,facecolor=SURF,bbox_inches='tight')
print("zapisano figures/pc_layer_3d.png")
