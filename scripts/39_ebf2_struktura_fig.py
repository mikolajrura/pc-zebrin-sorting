"""Figura: profil pLDDT Ebf2 vs Ebf1 + tozsamosc sekwencji, z zaznaczonym epitopem HPA.
Pokazuje, dlaczego epitop 413-550 rozroznia paralogi, a jednoczesnie nie nadaje sie do dokowania.
"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
S="/tmp/claude-1000/-home-mikolajrurad-omics-data-highway-ichb-purkinje-analysis-pc-zebrin-sorting/9d42a4c4-0d4a-4a5a-ab45-d0bef76c7937/scratchpad"
P="/home/mikolajrurad/omics-data-highway/ichb-purkinje-analysis/pc-zebrin-sorting"
BLUE="#2a78d6"; ORANGE="#eb6834"; INK="#1a1a19"; MUTED="#6b6a63"; GRID="#e6e5e0"

def plddt(f):
    r={}
    for l in open(f):
        if l.startswith("ATOM") and l[12:16].strip()=="CA": r[int(l[22:26])]=float(l[60:66])
    return np.array([r[k] for k in sorted(r)])
def rd(f): return "".join(l.strip() for l in open(f) if not l.startswith(">"))

p2=plddt(f"{S}/af_ebf2.pdb"); p1=plddt(f"{S}/af_ebf1.pdb")
e2=rd(f"{S}/ebf2.fasta"); e1=rd(f"{S}/Ebf1.fasta")

# Needleman-Wunsch (ten sam kod co w analizie) -> tozsamosc w oknie 21 wzdluz Ebf2
def nw(a,b,m=2,mm=-1,g=-2):
    n,l=len(a),len(b)
    H=np.zeros((n+1,l+1),dtype=np.int32); Pt=np.zeros((n+1,l+1),dtype=np.int8)
    H[0,:]=np.arange(l+1)*g; H[:,0]=np.arange(n+1)*g; Pt[0,1:]=3; Pt[1:,0]=2
    bb=np.frombuffer(b.encode(),dtype='S1')
    for i in range(1,n+1):
        sc=np.where(bb==a[i-1].encode(),m,mm); prev,cur=H[i-1],H[i]
        for j in range(1,l+1):
            v=(prev[j-1]+sc[j-1],prev[j]+g,cur[j-1]+g)
            k=int(np.argmax(v)); cur[j]=v[k]; Pt[i,j]=[1,2,3][k]
    i,j=n,l; A=B=""
    while i>0 or j>0:
        k=Pt[i,j]
        if k==1: A=a[i-1]+A; B=b[j-1]+B; i-=1; j-=1
        elif k==2: A=a[i-1]+A; B="-"+B; i-=1
        else: A="-"+A; B=b[j-1]+B; j-=1
    return A,B
A,B=nw(e2,e1)
match=[]; pos=0
for x,y in zip(A,B):
    if x!="-":
        pos+=1; match.append(1.0 if x==y else 0.0)
match=np.array(match); W=21
ident=np.convolve(match,np.ones(W)/W,mode="same")*100
x=np.arange(1,len(p2)+1)

fig,axes=plt.subplots(2,1,figsize=(11,6.2),sharex=True,
                      gridspec_kw=dict(height_ratios=[1,0.85],hspace=0.16))
fig.patch.set_facecolor("#fcfcfb")
DOMS=[(34,244,"domena wiążąca DNA"),(253,336,"IPT/TIG")]
EPI=(413,550); EPI2=(496,528)

for ax in axes:
    ax.set_facecolor("#fcfcfb")
    for s in ("top","right"): ax.spines[s].set_visible(False)
    for s in ("left","bottom"): ax.spines[s].set_color(GRID)
    ax.tick_params(colors=MUTED,labelsize=9); ax.grid(axis="y",color=GRID,lw=0.8)
    ax.set_axisbelow(True)
    for a,b,_ in DOMS: ax.axvspan(a,b,color=BLUE,alpha=0.06,lw=0)
    ax.axvspan(*EPI,color=ORANGE,alpha=0.11,lw=0)

ax=axes[0]
ax.axhline(70,color=MUTED,lw=1,ls=(0,(4,3)))
ax.plot(x,p2,color=BLUE,lw=2,label="EBF2 (O08792)")
ax.plot(np.arange(1,len(p1)+1),p1,color=ORANGE,lw=2,alpha=.75,label="EBF1 (Q07802)")
ax.set_ylim(0,105); ax.set_ylabel("pLDDT (AlphaFold)",color=INK,fontsize=10)
ax.text(0.5,72,"próg 70 = struktura wiarygodna",color=MUTED,fontsize=8,ha="left",va="bottom")
for a,b,nm in DOMS: ax.text((a+b)/2,101,nm,ha="center",fontsize=8.5,color=BLUE)
ax.text((EPI[0]+EPI[1])/2,101,"immunogen HPA003954",ha="center",fontsize=8.5,color=ORANGE)
ax.legend(frameon=False,fontsize=9,loc="lower left",labelcolor=INK)
ax.set_title("EBF2 mysz: gdzie jest struktura, a gdzie epitop odróżniający od EBF1",
             color=INK,fontsize=12.5,pad=26,loc="left",weight="bold")

ax=axes[1]
ax.axhline(100,color=GRID,lw=1)
ax.plot(x,ident,color=INK,lw=2)
ax.fill_between(x,ident,color=INK,alpha=0.07)
ax.set_ylim(18,108); ax.set_ylabel("tożsamość z EBF1 (%)\nokno 21 aa",color=INK,fontsize=10)
ax.set_xlabel("pozycja w sekwencji EBF2 (aa)",color=INK,fontsize=10)
ax.set_xlim(1,len(p2))
seg=ident[EPI[0]-1:EPI[1]]
ax.annotate(f"{seg.mean():.0f}% w epitopie\n→ przeciwciało rozróżnia",
            xy=(500,seg.mean()),xytext=(500,29),
            ha="center",fontsize=9.5,color=ORANGE,weight="bold",
            bbox=dict(boxstyle="round,pad=0.35",fc="#fcfcfb",ec=ORANGE,lw=1),
            arrowprops=dict(arrowstyle="-",color=ORANGE,lw=1.2))
segd=ident[33:244]
ax.annotate(f"{segd.mean():.0f}% w domenie DNA\n→ tu przeciwciało by NIE rozróżniło",
            xy=(140,segd.mean()),xytext=(140,29),ha="center",fontsize=9.5,color=BLUE,weight="bold",
            bbox=dict(boxstyle="round,pad=0.35",fc="#fcfcfb",ec=BLUE,lw=1),
            arrowprops=dict(arrowstyle="-",color=BLUE,lw=1.2))
fig.text(0.008,0.012,"Dane: AlphaFold DB v6 + UniProt (Needleman-Wunsch, match+2/mismatch-1/gap-2). "
         "Epitop HPA003954: 0/138 reszt o pLDDT≥70.",fontsize=8,color=MUTED)
plt.tight_layout(rect=[0,0.03,1,1])
out=f"{P}/figures/ebf2_epitop_struktura.png"
plt.savefig(out,dpi=160,facecolor=fig.get_facecolor())
print(f"zapisano {out}")
print(f"tozsamosc z EBF1: epitop {seg.mean():.1f}%, domena DNA {segd.mean():.1f}%")
